#!/usr/bin/env python3
"""Independent brute-force check of BraWl's hcp energies and ASRO.

Reads brawl.inp (lattice ideal_hcp or hcp, with c_over_a, epi_indexing
and shell_tolerance), the interaction file and a NetCDF configuration
written by BraWl, then

  * enumerates every pair of atoms in the periodic cell with an explicit
    minimum-image search using the real-space lattice and basis vectors,
  * classifies pairs by (in-plane separation, layers apart), orders the
    classes by distance (ties: layers apart, then in-plane separation)
    and groups them into interaction blocks exactly as documented for
    epi_indexing = 'shells' (merge within shell_tolerance) or 'classes',
  * evaluates E = 1/2 sum_i sum_g sum_{j in g(i)} V^{(g)}(s_i, s_j) and
    the conditional pair probabilities r(i, j, g) per group,

and compares them with the energy and ASRO trajectory files written by
the run. None of BraWl's neighbour tables are used.

Usage:
  hcp_brute_force_check.py brawl.inp config.nc energy_trajectory.dat
      [--asro asro_trajectory.dat] [--step N]

Exit status 0 if everything agrees, 1 otherwise.

Authors: C. D. Woodgate, Claude Fable 5.1 (2026)
"""
import argparse
import itertools
import math
import re
import sys

import numpy as np
import netCDF4 as nc

IDEAL = math.sqrt(8.0 / 3.0)
BASIS_FRAC = np.array([[1/3, 2/3, 1/4], [2/3, 1/3, 3/4]])


def read_inp(path):
    d = {}
    for line in open(path):
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            d[k.strip()] = v.strip().strip("'\"")
    return d


def read_vij(path, n_species, n_blocks):
    vals = [float(x) for x in re.split(r"\s+", open(path).read().strip())]
    need = n_species * n_species * n_blocks
    assert len(vals) >= need, (len(vals), need)
    vals = vals[:need]
    # Fortran list-directed read into V(n_species, n_species, n_blocks),
    # first index fastest
    v = np.array(vals).reshape((n_blocks, n_species, n_species))
    return v.transpose(2, 1, 0)                   # v[i, j, g]


def lattice(gamma):
    a1 = np.array([1.0, 0.0, 0.0])
    a2 = np.array([-0.5, math.sqrt(3.0) / 2.0, 0.0])
    a3 = np.array([0.0, 0.0, gamma])
    return np.array([a1, a2, a3])                 # rows


def class_key(dfrac):
    """Exact class key (9 rho^2, layers apart) from a fractional pair
    vector whose in-plane parts are multiples of 1/3 and whose c part is
    a multiple of 1/4."""
    x = int(round(3 * dfrac[0]))
    y = int(round(3 * dfrac[1]))
    z4 = int(round(4 * dfrac[2]))
    assert abs(3 * dfrac[0] - x) < 1e-9 and abs(3 * dfrac[1] - y) < 1e-9
    assert abs(4 * dfrac[2] - z4) < 1e-9 and z4 % 2 == 0
    return (x * x + y * y - x * y, abs(z4) // 2)


def class_d2(key, gamma):
    k1, zeta = key
    return k1 / 9.0 + zeta * zeta * gamma * gamma / 4.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp")
    ap.add_argument("config_nc")
    ap.add_argument("traj")
    ap.add_argument("--asro", default=None)
    ap.add_argument("--step", type=int, default=None)
    args = ap.parse_args()

    p = read_inp(args.inp)
    assert p["lattice"] in ("hcp", "ideal_hcp"), p["lattice"]
    gamma = IDEAL if p["lattice"] == "ideal_hcp" else float(p["c_over_a"])
    indexing = p.get("epi_indexing", "shells")
    tol = float(p.get("shell_tolerance", "1.0e-6"))
    n1, n2, n3 = int(p["n_1"]), int(p["n_2"]), int(p["n_3"])
    ns, nr = int(p["n_species"]), int(p["interaction_range"])
    wc = int(p.get("wc_range", "0"))
    L = lattice(gamma)
    V = read_vij(p["interaction_file"], ns, nr)

    cfg = np.array(nc.Dataset(args.config_nc).variables["configuration"][:])
    assert cfg.shape == (n3, n2, n1, 2), cfg.shape
    cfg = cfg.transpose(3, 2, 1, 0)               # cfg[b, i, j, k]
    sites = [(b, i, j, k) for b in range(2) for i in range(n1)
             for j in range(n2) for k in range(n3)]
    frac = np.array([BASIS_FRAC[b] + np.array([i, j, k])
                     for (b, i, j, k) in sites], dtype=float)
    spec = np.array([cfg[b, i, j, k] for (b, i, j, k) in sites]) - 1
    N = len(sites)
    box = np.array([n1, n2, n3], dtype=float)
    counts = np.bincount(spec, minlength=ns)

    # ---- pair classes present in the periodic cell -------------------
    # Minimum-image search: for each pair keep every image within the
    # cutoff implied by the groups needed (determined below from an
    # independent enumeration), so that periodic double counting matches
    # the Hamiltonian's convention of one interaction per offset.
    n_groups = max(nr, wc - 1)

    # Independent enumeration of classes out to a generous radius
    keys = {}
    for b0 in range(2):
        for b in range(2):
            for off in itertools.product(range(-6, 7), range(-6, 7),
                                         range(-4, 5)):
                dfrac = BASIS_FRAC[b] - BASIS_FRAC[b0] + np.array(off)
                key = class_key(dfrac)
                if key == (0, 0):
                    continue
                keys[key] = class_d2(key, gamma)
    order = sorted(keys, key=lambda k: (round(math.sqrt(keys[k]), 12),
                                        k[1], k[0]))
    # tolerance-aware ordering: re-sort treating near-equal distances
    # as ties
    def before(a, b):
        da, db = math.sqrt(keys[a]), math.sqrt(keys[b])
        if abs(da - db) <= tol * max(da, db):
            return (a[1], a[0]) < (b[1], b[0])
        return da < db
    ranked = []
    for k in order:
        i = len(ranked)
        while i > 0 and before(k, ranked[i - 1]):
            i -= 1
        ranked.insert(i, k)
    group_of = {}
    g = 0
    for i, k in enumerate(ranked):
        if i == 0:
            g = 1
        elif indexing == "classes":
            g += 1
        elif math.sqrt(keys[k]) > math.sqrt(keys[ranked[i - 1]]) * (1 + tol):
            g += 1
        group_of[k] = g
    used = {k: g for k, g in group_of.items() if g <= n_groups}
    d_max = max(math.sqrt(keys[k]) for k in used) * (1 + 2 * tol)
    print("groups used:", n_groups, " cutoff distance / a: %.4f" % d_max)
    gd = {}
    for k in ranked:
        if k in used and used[k] not in gd:
            gd[used[k]] = math.sqrt(keys[k])
    print("group distances / a:", [round(gd[g], 4) for g in sorted(gd)])

    # ---- brute force energy and pair counts ---------------------------
    images = np.array(list(itertools.product(range(-2, 3), repeat=3)),
                      dtype=float)
    energy = 0.0
    r_dens = np.zeros((ns, ns, n_groups + 1))
    coord = np.zeros(n_groups + 1)
    for a in range(N):
        d = frac - frac[a]
        d = d - np.round(d / box) * box
        r_dens[spec[a], spec[a], 0] += 1          # on-site term
        for im in images:
            dd = d + im * box
            dist = np.sqrt(np.einsum("ij,ij->i", dd @ L, dd @ L))
            for bidx in np.nonzero(dist <= d_max)[0]:
                if bidx == a and np.allclose(im, 0):
                    continue
                key = class_key(dd[bidx])
                if key not in used:
                    continue
                g = used[key]
                coord[g] += 1
                if g <= nr:
                    energy += V[spec[a], spec[bidx], g - 1]
                r_dens[spec[a], spec[bidx], g] += 1
    energy *= 0.5
    print("coordination per site by group:",
          [int(round(x)) for x in (coord / N)[1:]])

    ok = True
    traj = [l.split() for l in open(args.traj)
            if not l.strip().startswith("#")]
    traj = {int(s): float(e) for s, e in traj}
    step = args.step if args.step is not None else max(traj)
    e_code = traj[step]
    diff = abs(energy - e_code)
    print(f"energy: brute force {energy:.12e} Ry, BraWl (step {step}) "
          f"{e_code:.12e} Ry, difference {diff:.2e}")
    if diff > 1e-9:
        ok = False

    if args.asro:
        for si in range(ns):
            r_dens[si] /= counts[si]
        lines = [l for l in open(args.asro) if not l.strip().startswith("#")]
        vals = np.array([float(x) for x in lines[-1].split()[1:]])
        assert vals.size == ns * ns * wc, (vals.size, ns, wc)
        code = vals.reshape((wc, ns, ns)).transpose(2, 1, 0)
        adiff = np.abs(r_dens[:, :, :wc] - code).max()
        print(f"ASRO: max |brute force - BraWl| over (i, j, group) "
              f"= {adiff:.2e}")
        if adiff > 2e-5:                           # file holds 5 decimals
            ok = False

    print("OK" if ok else "MISMATCH")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

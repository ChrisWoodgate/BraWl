@page lattices Lattice types and neighbour tables

This page describes how `BraWl` represents a lattice internally, which lattice types are available, and how the hexagonal close-packed (hcp) lattices are handled: `ideal_hcp` with hard-coded neighbour tables, and `hcp` with an arbitrary axial ratio built on runtime *pair tables*. It is intended both for users choosing a lattice and for developers adding one.

## How a configuration is stored

A configuration is an integer array `config(b, i, j, k)` holding the species on each site: `b` is the index of the atom in the basis and `(i, j, k)` the index of a cell along the three grid directions. The array extents are `n_basis` and `n_grid_1`, `n_grid_2`, `n_grid_3`, which `initialise_function_pointers` sets from the lattice type and the input values `n_1`, `n_2`, `n_3`. Two conventions are in use:

- **Cubic lattices** (`simple_cubic`, `bcc`, `fcc`) use `n_basis = 1` on a simple cubic grid of spacing \f$a/2\f$ spanning the conventional cubic cell, so `n_grid_i = 2 n_i`. Grid points that are not lattice sites hold the value 0 and are skipped everywhere. This is the original representation of the code and keeps the hard-coded neighbour offsets of those lattices simple.
- **Lattices with a multi-atom basis** (`ideal_hcp`, `hcp`) use one grid point per primitive cell, `n_grid_i = n_i`, and the basis index distinguishes the atoms in the cell. Every grid point is a site.

Each lattice provides, through procedure pointers in `run_params`, a local energy kernel (`nbr_energy`), a random site (`rdm_site`) and a random neighbour (`rdm_nbr`), so the sampling algorithms never need to know the lattice. The total energy (`total_energy`) and the pair energy used for trial swaps (`pair_energy`) are built from the local kernel. The grid vectors (`lattice_vectors`, one per column, in units of the lattice parameter) and the Cartesian basis positions (`basis_vectors`) are used by the xyz writer.

## The cubic lattices

`n_1`, `n_2`, `n_3` count conventional cubic cells (2 atoms each for `bcc`, 4 for `fcc`; for `simple_cubic` each "cell" is a 2x2x2 block of 8 sites), and `lattice_parameter` is the cubic lattice constant. Neighbour offsets for each coordination shell are hard-coded in `bw_hamiltonian.f90` (`fcc_shellN_energy` and so on) and selected according to `interaction_range`, which supports up to 2 shells for `simple_cubic`, 10 for `bcc` and 6 for `fcc`. Shell distances for the ASRO analysis are computed in units of the grid spacing by `lattice_shells`.

## The hcp lattices

### Geometry

Both hcp lattices use the hexagonal primitive cell with lattice vectors, in units of the in-plane lattice parameter \f$a\f$,
\f[
\mathbf{a}_1 = (1, 0, 0), \quad \mathbf{a}_2 = (-\tfrac{1}{2}, \tfrac{\sqrt{3}}{2}, 0), \quad \mathbf{a}_3 = (0, 0, c/a),
\f]
and a two-atom basis at fractional coordinates \f$(\tfrac{1}{3}, \tfrac{2}{3}, \tfrac{1}{4})\f$ and \f$(\tfrac{2}{3}, \tfrac{1}{3}, \tfrac{3}{4})\f$. `n_1`, `n_2`, `n_3` count primitive cells of two atoms each and `lattice_parameter` is \f$a\f$. For `ideal_hcp` the axial ratio is \f$c/a = \sqrt{8/3} \approx 1.63299\f$; for `hcp` it is given by the input keyword `c_over_a`.

### Pair classes, and why distance shells are not enough

The Bragg-Williams energy
\f[
H = \tfrac{1}{2} \sum_{i} \sum_{n} \sum_{j \in n(i)} \sum_{\gamma\gamma'} V^{(n)}_{\gamma\gamma'} \, \xi_{i\gamma} \xi_{j\gamma'}
\f]
depends only on *which* pairs of sites share an interaction block \f$n\f$, never on the distance between them. For hcp the set of neighbour offsets is fixed by the topology of the lattice; the axial ratio only decides how those offsets are ordered by distance and which of them happen to coincide.

The natural unit is therefore the **pair class**: a set of neighbour offsets equivalent under the crystal symmetry. Every operation of the hcp space group \f$P6_3/mmc\f$ maps the \f$c\f$ axis onto \f$\pm c\f$, so the magnitude of the \f$c\f$-component of a pair vector is invariant, and the in-plane part transforms under the hexagonal point group. A class is thus identified by the in-plane separation and the number of layers (in units of \f$c/2\f$) between the two atoms. In the code the key is the exact integer pair \f$(9\rho^2/a^2, 2|\Delta z|/c)\f$, where \f$\rho\f$ is the in-plane separation; the squared distance of a class is \f$d^2/a^2 = \rho^2/a^2 + (2|\Delta z|/c)^2 (c/a)^2/4\f$.

Two consequences matter in practice:

1. **At the ideal ratio some classes coincide in distance but are not equivalent.** The six in-plane and six out-of-plane nearest neighbours of an hcp atom lie at the same distance only because \f$c/a = \sqrt{8/3}\f$; no symmetry operation connects them, and in a real alloy their effective pair interactions differ. Of the ten ideal shells, shells 1 and 4 each merge two inequivalent classes (6+6 and 6+12); the other eight are single classes.
2. **Away from the ideal ratio the ordering changes.** For \f$c/a < \sqrt{8/3}\f$ (Ti, Zr, Mg) the out-of-plane nearest neighbours come first; for \f$c/a > \sqrt{8/3}\f$ (Zn, Cd) the in-plane ones do, and the two-atom class along \f$c\f$ moves outward past several others.

The class key coincides with the symmetry orbits of pair vectors for all in-plane separations below \f$7a\f$, which is the first distance at which two inequivalent in-plane vectors share a length; the table builder refuses to tabulate beyond that.

### `ideal_hcp`: hard-coded tables

For the ideal ratio the ten distance shells are hard-coded as compile-time tables in `bw_hamiltonian.f90` (`ideal_hcp_shellN_nbrs`, one list of offsets per sublattice of the central atom), evaluated by `ideal_hcp_shell_energy`, with the periodically wrapped cell indices computed once per site by `ideal_hcp_wrap`. The tables are generated, and checked against a brute-force enumeration, by `tests/ideal_hcp_neighbour_tables.py`. The shells, their radii and coordination numbers are:

| Shell | Distance | Coordination | Classes (in-plane \f$\rho/a\f$, layers) |
|---|---|---|---|
| 1 | \f$a\f$ | 12 | (1, 0) and (\f$1/\sqrt{3}\f$, 1) |
| 2 | \f$\sqrt{2}\,a\f$ | 6 | (\f$2/\sqrt{3}\f$, 1) |
| 3 | \f$\sqrt{8/3}\,a\f$ | 2 | (0, 2) |
| 4 | \f$\sqrt{3}\,a\f$ | 18 | (\f$\sqrt{3}\f$, 0) and (\f$\sqrt{7/3}\f$, 1) |
| 5 | \f$\sqrt{11/3}\,a\f$ | 12 | (1, 2) |
| 6 | \f$2a\f$ | 6 | (2, 0) |
| 7 | \f$\sqrt{5}\,a\f$ | 12 | (\f$\sqrt{13/3}\f$, 1) |
| 8 | \f$\sqrt{17/3}\,a\f$ | 12 | (\f$\sqrt{3}\f$, 2) |
| 9 | \f$\sqrt{6}\,a\f$ | 6 | (\f$4/\sqrt{3}\f$, 1) |
| 10 | \f$\sqrt{19/3}\,a\f$ | 6 | (\f$1/\sqrt{3}\f$, 3) |

Because the kernels index the interactions by distance shell, `ideal_hcp` necessarily uses one interaction for both classes of shells 1 and 4. If class-resolved interactions are available, use `hcp` with `c_over_a` set to the ideal ratio and `epi_indexing = classes` (see below); at the ideal ratio with `epi_indexing = shells` the two lattices give identical results.

### `hcp`: runtime pair tables for any axial ratio

Lattice `hcp` builds its neighbour tables at start-up in `neighbour_tables.f90` (`build_hcp_tables`) from the geometry above and the input keywords:

| Keyword | Meaning | Default |
|---|---|---|
| `c_over_a` | Axial ratio \f$c/a\f$ (required for `hcp`; must not be given for `ideal_hcp`) | none |
| `epi_indexing` | `shells`: classes whose distances agree within `shell_tolerance` are merged into one interaction block; `classes`: every class is its own block | `shells` |
| `shell_tolerance` | Relative tolerance on distance used for the merge | `1.0e-6` |

The builder enumerates all pair vectors in a box of cells that is enlarged until enough *groups* lie safely inside, identifies the classes by their exact key, ranks them by distance for the given \f$c/a\f$ (ties broken by layers apart, then in-plane separation, so in-plane classes precede out-of-plane ones at equal distance), and assigns them to groups according to `epi_indexing`. Groups are the blocks of the interaction file and the shells of the ASRO analysis; `interaction_range` is the number of blocks the file must contain, and `wc_range - 1` the number of groups beyond the on-site term reported by the ASRO analysis. The tables hold `max(interaction_range, wc_range - 1)` groups.

The energy is evaluated by a single kernel, `hcp_energy_tables`, of the same form as the ideal one: the neighbour list of each sublattice is ordered by group and each entry carries its block index, so the kernel sums group by group with one table lookup per neighbour. In a 1024-atom benchmark of 4.1 million trial swaps it runs roughly 10 to 30 percent slower than the hard-coded `ideal_hcp` kernels, whose compile-time tables the compiler can unroll, so for the ideal ratio with shell indexing `ideal_hcp` remains the fast path. The ASRO analysis (`lattice_shells`, `radial_densities`) reads the same lists, so there is one definition of which pairs belong to which group. "Nearest-neighbour" swaps (`nbr_swap = T`) use the twelve topological nearest neighbours of the hcp lattice whatever the axial ratio, since for \f$c/a > \sqrt{8/3}\f$ the first distance shell contains only the six in-plane atoms.

At start-up the full table is printed, for example for \f$c/a = 1.587\f$ with `epi_indexing = shells`:

```
 Group    d/a   Mult.   Classes (in-plane d/a, layers of c/2, mult.)
 ------------------------------------------------------------------
     1  0.9813     6   (0.5774, 1,  6)
     2  1.0000     6   (1.0000, 0,  6)
     3  1.4011     6   (1.1547, 1,  6)
     4  1.5870     2   (0.0000, 2,  2)
     5  1.7213    12   (1.5275, 1, 12)
     6  1.7321     6   (1.7321, 0,  6)
     7  1.8758    12   (1.0000, 2, 12)
     8  2.0000     6   (2.0000, 0,  6)
     9  2.2278    12   (2.0817, 1, 12)
    10  2.3492    12   (1.7321, 2, 12)
 ------------------------------------------------------------------
 Interactions use groups 1 to  10 (distances up to  2.3492 a)
```

Match the blocks of the interaction file to the groups in this order. Groups whose distances differ by less than a thousand times `shell_tolerance` are flagged, since a ratio entered as `1.633` splits the ideal shells at the fourth decimal; use `ideal_hcp`, or `epi_indexing = classes`, when that is not intended. The interaction file is checked to contain at least `interaction_range` blocks before it is read, and a surplus is reported. Configuration files written in NetCDF record `c_over_a` and the indexing mode, and a restart refuses a file written with different values.

Limits worth knowing: `interaction_range` is a number of blocks, not a cutoff radius, so ten groups of a non-ideal lattice reach a shorter distance than ten ideal shells (the covered radius is printed); the tables are limited to in-plane separations below \f$7a\f$; and the class key is specific to hcp, so another lattice with a basis would supply its own key function to the same machinery.

### Validation

- `tests/ideal_hcp_neighbour_tables.py --check src/random_site.f90 src/bw_hamiltonian.f90` verifies the hard-coded ideal tables against a brute-force enumeration with exact arithmetic (run in CI).
- `tests/hcp_brute_force_check.py` recomputes the energy and the conditional pair probabilities of a configuration written by `BraWl` from the real-space geometry alone, applying the documented grouping rules, and compares them with the trajectory files (run in CI on the hcp regression test).
- `tests/compare_energy_trajectories.py` checks that `hcp` with `c_over_a` equal to the ideal ratio reproduces the `ideal_hcp` reference trajectory of test 05 (run in CI); the two agree bit for bit.
- Tests 05 (`ideal_hcp`) and 06 (`hcp`, \f$c/a = 1.587\f$) are fixed-seed regression tests.

## Adding another lattice

Two routes exist. A lattice with a single-atom basis on a cubic grid follows the `fcc` pattern: hard-coded shell kernels in `bw_hamiltonian.f90`, random site and neighbour routines in `random_site.f90`, an occupancy pattern in `initial_setup`, and a block in `initialise_function_pointers` that sets `n_basis`, `n_atoms`, `n_grid_i`, the grid and basis vectors and the procedure pointers. A lattice with a multi-atom basis is better served by the table engine: provide the geometry and an exact class key (the analogue of \f$(9\rho^2, 2|\Delta z|)\f$) for its symmetry, and the enumeration, ranking, grouping, kernel and ASRO code in `neighbour_tables.f90`, `bw_hamiltonian.f90` and `analytics.f90` apply unchanged. In both cases add a fixed-seed regression test with reference data and, where possible, an independent brute-force check.

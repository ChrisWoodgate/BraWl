import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from matplotlib.backends.backend_pdf import PdfPages

#font_size = 22
np.set_printoptions(suppress=True)
latex_rcparams = {
    #"font.size": font_size,
    "font.family": "serif",
    "text.usetex": True,
}
plt.rcParams.update(latex_rcparams)

# ------------------------------------------------------------
# Compute Cv curve inside a folder
# ------------------------------------------------------------
def compute_cv_curve(path):
    cwd = os.getcwd()
    try:
        os.chdir(path)

        bin_edges = np.array(nc.Dataset("data/wl_dos_bins.nc")["grid data"][:], dtype=np.float64)
        wl_logdos = np.array(nc.Dataset("data/wl_dos.nc")["grid data"][:], dtype=np.float64)

        kb_ev = 8.617333262e-5
        ev_to_ry = 13.605693122
        n_atoms = 2048

        bin_edges = bin_edges * ev_to_ry
        bin_width = bin_edges[1] - bin_edges[0]

        start_temp = 0
        end_temp = 5000
        step_size = 25
        temperatures = np.arange(start_temp, end_temp + step_size, step_size).astype(np.float64)
        temperatures[0] = 1e-6

        mean_energies = np.zeros_like(temperatures)
        heat_caps = np.zeros_like(temperatures)

        for itemp, T in enumerate(temperatures):
            beta = 1.0 / (kb_ev * T)

            prob = np.zeros(len(bin_edges) - 1, dtype=np.float64)
            for ibin, edge in enumerate(bin_edges[:-1]):
                bin_energy = edge + 0.5 * bin_width
                prob[ibin] = wl_logdos[ibin] - beta * bin_energy

            prob -= np.max(prob)
            prob = np.exp(prob)
            prob = np.nan_to_num(prob)
            prob /= np.sum(prob)

            mean_energy = np.dot(bin_edges[:-1] + 0.5 * bin_width, prob)
            mean_energies[itemp] = mean_energy

            msq_dev = np.zeros(len(bin_edges) - 1, dtype=np.float64)
            for ibin, edge in enumerate(bin_edges[:-1]):
                bin_energy = edge + 0.5 * bin_width
                msq_dev[ibin] = (bin_energy - mean_energies[itemp]) ** 2

            heat_caps[itemp] = np.dot(msq_dev, prob) / (kb_ev * T**2)

        heat_caps /= (kb_ev * n_atoms)
        return temperatures, heat_caps

    finally:
        os.chdir(cwd)

# ------------------------------------------------------------
# Shoulder candidates (curvature peaks, no smoothing)
# ------------------------------------------------------------
def shoulder_indices(T, y, prominence=None, top_n=None, min_sep=None):
    T = np.asarray(T, float)
    y = np.asarray(y, float)

    d1 = np.gradient(y, T)
    d2 = np.gradient(d1, T)
    curv = -d2

    if prominence is None:
        mad = np.median(np.abs(curv - np.median(curv)))
        prominence = 3.0 * 1.4826 * mad

    curv_pos = np.maximum(curv, 0.0)
    if np.max(curv_pos) > 0:
        curv_pos = curv_pos / np.max(curv_pos)

    idx, _ = find_peaks(curv_pos, prominence=prominence)

    if min_sep is not None and len(idx) > 0:
        min_sep = int(min_sep)
        scores = curv_pos[idx]
        order = np.argsort(scores)[::-1]
        chosen = []
        for j in order:
            s = int(idx[j])
            if all(abs(s - c) >= min_sep for c in chosen):
                chosen.append(s)
        idx = np.array(sorted(chosen), dtype=int)

    if top_n is not None and len(idx) > 0:
        n = int(top_n)
        if n > 0 and len(idx) > n:
            scores = curv_pos[idx]
            keep = np.argsort(scores)[::-1][:n]
            idx = np.sort(idx[keep])

    return idx, curv_pos

def filter_shoulders_near_peaks(shoulder_idx, peak_idx, N):
    shoulder_idx = np.asarray(shoulder_idx, dtype=int)
    peak_idx = np.asarray(peak_idx, dtype=int)
    N = int(N)

    if shoulder_idx.size == 0 or peak_idx.size == 0 or N < 0:
        return shoulder_idx

    keep = np.ones(shoulder_idx.shape[0], dtype=bool)
    for i, s in enumerate(shoulder_idx):
        if np.any(np.abs(peak_idx - s) <= N):
            keep[i] = False
    return shoulder_idx[keep]

# ------------------------------------------------------------
# MAIN: one page per Cu* folder + final phase diagram page
# ------------------------------------------------------------
base = os.getcwd()
folders = sorted([f for f in os.listdir(base) if f.startswith("Cu")])

cv_dict = {}

# phase diagram point stores (3 groups)
x_pk1, T_pk1 = [], []   # highest Cv peak per composition
x_pk2, T_pk2 = [], []   # second-highest Cv peak per composition
x_pk3, T_pk3 = [], []   # third-highest Cv peak per composition
x_other, T_other = [], []  # remaining peaks + shoulders

with PdfPages("copper_all_dirs.pdf") as pdf:
    for folder in folders:
        path = os.path.join(base, folder)
        x = float(folder.replace("Cu", ""))

        T, Cv = compute_cv_curve(path)
        cv_dict[str(x)] = {"temp": T, "cv": Cv}

        # Peaks (indices into T/Cv)
        peaks, _ = find_peaks(Cv / np.max(Cv), prominence=0.005)

        # Shoulders
        shoulders, curv = shoulder_indices(T, Cv, prominence=0.01, top_n=3, min_sep=10)
        shoulders = filter_shoulders_near_peaks(shoulders, peaks, N=5)

        # Combine peaks and shoulders for phase point grouping
        peaks = np.sort(np.concatenate((peaks, shoulders)))

        # ---- group phase points by Cv height at peak ----
        if len(peaks) > 0:
            pk_heights = Cv[peaks]
            order = np.argsort(pk_heights)[::-1]  # highest Cv first
            pk_sorted = peaks[order]

            # highest peak
            x_pk1.append(x)
            T_pk1.append(T[pk_sorted[0]])

            # second-highest peak (if exists)
            if len(pk_sorted) >= 2:
                x_pk2.append(x)
                T_pk2.append(T[pk_sorted[1]])

            # third-highest peak (if exists)
            if len(pk_sorted) >= 3:
                x_pk3.append(x)
                T_pk3.append(T[pk_sorted[2]])

            # any remaining peaks -> "other"
            if len(pk_sorted) > 3:
                rem = pk_sorted[3:]
                x_other.extend([x] * len(rem))
                T_other.extend(T[rem].tolist())


        # ---- per-folder plot page ----
        #fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
        #ax.plot(T, Cv, label="Cv(T)")

        #if np.max(curv) > 0:
        #    ax.plot(T, curv * np.max(Cv), label="Curvature (scaled)", alpha=0.6)

        #for i, pk in enumerate(peaks):
        #    ax.axvline(T[pk], color="red", linestyle="--", alpha=0.7,
        #               label="Peak" if i == 0 else None)

        #for i, sh in enumerate(shoulders):
        #    ax.axvline(T[sh], color="green", linestyle=":", alpha=0.7,
        #               label="Shoulder" if i == 0 else None)

        #ax.set_title(folder)
        #ax.set_xlabel("Temperature (K)")
        #ax.set_ylabel(r"$C_V$ ($k_B$/atom)")
        #ax.grid(True, alpha=0.3)
        #ax.legend()

        #pdf.savefig(fig)
        #plt.close(fig)

    # ---- final page: phase diagram (3 groups) ----
    #fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)

    #if len(x_pk1):
    #    ax.scatter(x_pk1, T_pk1)

    #if len(x_pk2):
    #    ax.scatter(x_pk2, T_pk2)

    #ax.set_title("Phase diagram")
    #ax.set_xlabel("Cu concentration, $x$")
    #ax.set_ylabel("Temperature (K)")
    #ax.grid(True, alpha=0.3)
    #ax.legend()

    #pdf.savefig(fig)
    #plt.close(fig)

# ---- final page: phase diagram (3 groups) ----
fig_height = 3.375; fig_width=4.5

fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)

ax.set_ylim(200, 3600)
ax.set_xlim(0, 1)

if len(x_pk1):
    ax.scatter(x_pk1, T_pk1)

if len(x_pk2):
    ax.scatter(x_pk2, T_pk2)

ax.set_title(r"Predicted $\textrm{Cu}_x (\textrm{Ni}_{3/4} \textrm{Al}_{1/4})_{1-x}$ Phase Diagram")
ax.set_xlabel("Cu concentration, $x$")
ax.set_ylabel("Temperature (K)")
fig.savefig("phase_diagram.pdf")

np.save("cv_dict.npy", cv_dict, allow_pickle=True)

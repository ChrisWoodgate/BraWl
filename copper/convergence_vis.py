import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import netCDF4 as nc
import os
import sys
from matplotlib.colors import ListedColormap
font_size = 18
np.set_printoptions(suppress=True)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": font_size})

colors = {
    "steel_blue": "#1F77B4",
    "orange": "#FF7F0E",
    "forest_green": "#2CA02C",
    "firebrick_red": "#D62728",
    "lavender": "#9467BD",
    "brown": "#8C564B",
    "orchid": "#E377C2",
    "gray": "#7F7F7F",
}
custom_cmap = ListedColormap(colors.values())
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=custom_cmap.colors)

kb_ev = 8.617333262e-5
ev_to_ry = 13.605693122
base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "convergence")
compositions = ["Cu0.00", "Cu0.50"]
fss_dirs = ["2", "4", "6", "8"]
cutoff_dirs = ["2", "4", "6", "8"]

# fixed temperature grid controls
start_temp = 0
end_temp = 5000
step_size = 25
temperatures = np.arange(start_temp, end_temp + step_size, step_size).astype(np.float64)
temperatures[0] = 1e-64

def heat_capacity(path, n_atoms):
    bin_edges = np.array(nc.Dataset(path + "data/wl_dos_bins.nc")["grid data"][:], dtype=np.float64)
    wl_logdos = np.array(nc.Dataset(path + "data/wl_dos.nc")["grid data"][:], dtype=np.float64)
    bin_edges_ev = bin_edges*ev_to_ry
    bin_width = bin_edges_ev[1] - bin_edges_ev[0]
    energies = bin_edges_ev[:-1] + 0.5*bin_width

    mean_energies = np.zeros(len(temperatures))
    for itemp, T in enumerate(temperatures):
        prob = wl_logdos - energies/(kb_ev*T)
        prob = prob - np.max(prob)
        prob = np.exp(prob)
        prob = np.nan_to_num(prob/np.sum(prob))
        mean_energies[itemp] = np.dot(energies, prob)

    # C_v = dE/dT in units of k_B per atom
    return np.gradient(mean_energies, temperatures)/kb_ev/n_atoms

def read_atoms(path):
    # Total atoms from species_numbers in brawl.inp
    with open(path + "brawl.inp", "r") as file:
        for line in file:
            if line.startswith("species_numbers"):
                return sum(int(x) for x in line.split("=")[1].split())

def plot_study(study, dirs, label, xlabel_fmt, out_name):
    # One figure per composition, styled as the specific heat panels in copper_vis.py
    label_size_offset = 3
    cv_ticks = 1000
    fig_height = 4.5
    fig_width = 4 * 1.68 * 1.5
    for composition in compositions:
        fig, ax = plt.subplots(1, 1, figsize=(fig_width, fig_height*1.2), constrained_layout=True)
        for d in dirs:
            path = os.path.join(base_path, study, composition, d) + "/"
            if not os.path.exists(path + "data/wl_dos.nc"):
                print("Skipping (no data yet):", path)
                continue
            n_atoms = read_atoms(path)
            ax.plot(temperatures, heat_capacity(path, n_atoms), '-', label=xlabel_fmt.format(d, n_atoms))
        ax.tick_params(direction="in")
        ax.tick_params(axis='both', pad=6)
        ax.grid(True, axis='x')
        ax.set_xlabel(r'Temperature, $T$ (K)', fontsize=font_size+label_size_offset)
        ax.set_ylabel(r'Specific heat, $C_V$ ($k_B$/atom)', fontsize=font_size+label_size_offset)
        ax.set_xticks(np.arange(start_temp, end_temp + cv_ticks, cv_ticks))
        ax.set_xlim(start_temp, end_temp)
        ax.set_ylim(0, ax.get_ylim()[1]*1.01)
        legend = ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), title=composition.replace("Cu", "$x_{\\mathrm{Cu}}=") + "$, " + label)
        for handle in legend.legend_handles:
            handle.set_linewidth(3)
        fig.savefig(os.path.join(base_path, out_name.format(composition)), bbox_inches='tight')
        plt.close(fig)

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "both"
    if which in ("fss", "both"):
        plot_study("fss", fss_dirs, "6 shells", "${0}^3$ cells, $N={1}$", "fss_{}.pdf")
    if which in ("cutoff", "both"):
        plot_study("cutoff", cutoff_dirs, "2048 atoms", "{0} shells", "cutoff_{}.pdf")

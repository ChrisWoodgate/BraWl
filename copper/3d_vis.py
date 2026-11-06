import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    #"font.size": font_size
    })

# ------------------------------------------------------------
# Compute Cv curve inside a folder
# ------------------------------------------------------------
def compute_cv_curve(path):
    os.chdir(path)

    bin_edges = np.array(nc.Dataset("data/wl_dos_bins.nc")["grid data"][:], dtype=np.float64)
    wl_logdos = np.array(nc.Dataset("data/wl_dos.nc")["grid data"][:], dtype=np.float64)

    kb_ev = 8.167333262e-5
    ev_to_ry = 13.605693122
    kb_ry = kb_ev / ev_to_ry

    n_atoms = 2048
    ev_to_mev = 1000

    bin_edges = bin_edges * ev_to_ry
    bin_width = bin_edges[1] - bin_edges[0]
    E = bin_edges[:-1] + 0.5 * bin_width

    start_temp = 0
    end_temp = 5000
    step_size = 25
    temperatures = np.arange(start_temp, end_temp + step_size, step_size).astype(np.float64)
    temperatures[0] = 1e-64

    mean_energies = np.zeros_like(temperatures)
    heat_caps = np.zeros_like(temperatures)

    for itemp, new_temp in enumerate(temperatures):

        beta = 1.0/(kb_ev*new_temp)

        # Reweighted histogram
        prob = np.zeros(len(bin_edges)-1)
        
        # Reweight histogram
        for ibin, edge in enumerate(bin_edges[:-1]):
            bin_energy = edge + 0.5*bin_width
            prob[ibin] = wl_logdos[ibin]-beta*bin_energy
    
        prob = prob-np.max(prob)
        prob = np.exp(prob)

        # Normalise
        prob = prob/(np.sum(prob))
        prob = np.nan_to_num(prob)

        # Mean energy
        mean_energy = np.dot(bin_edges[:-1]+0.5*bin_width, prob)
        mean_energies[itemp] = mean_energy

        # Compute heat capacity using the histogram
        msq_dev = np.zeros(len(bin_edges)-1)
        for ibin, edge in enumerate(bin_edges[:-1]):
            bin_energy = edge + 0.5*bin_width
            msq_dev[ibin] = (bin_energy - mean_energies[itemp])**2
            
        heat_caps[itemp] = np.dot(msq_dev, prob)/(kb_ev*new_temp**2)

    #heat_caps = np.gradient(mean_energies, temperatures)

    heat_caps /= kb_ev*n_atoms

    return temperatures, heat_caps

# ------------------------------------------------------------
# MAIN SCRIPT
# ------------------------------------------------------------
base = os.getcwd()
folders = sorted([f for f in os.listdir(base) if f.startswith("Cu")])

# Storage for surface plot
all_temps = []
all_cvs = []
all_conc = []

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# draw the figure so the renderer knows sizes
fig.canvas.draw()

# get the tight bounding box in figure coordinates
tight_bbox = ax.get_tightbbox(fig.canvas.get_renderer())
tight_bbox = tight_bbox.transformed(fig.transFigure.inverted())

color = "#1F77B4"
cv_dict = {}

for folder in folders:
    path = os.path.join(base, folder)

    conc = float(folder.replace("Cu", ""))     # z-coordinate (composition)
    temps, cv = compute_cv_curve(path)

    cv_dict[folder.replace("Cu", "")] = {
        "temp": temps,
        "cv": cv
    }

    z = np.full_like(temps, conc)

    # 3D curves
    ax.plot(temps, z, cv, label=folder, color=color)

    # store for surface
    all_temps.append(temps)
    all_cvs.append(cv)
    all_conc.append(np.full_like(temps, conc))

os.chdir(base)

np.save("cv_dict.npy", cv_dict)

plt.gca().invert_xaxis()

ax.xaxis.pane.fill = False
ax.yaxis.pane.fill = False
ax.zaxis.pane.fill = False

# Remove grid lines
ax.grid(False)

ax.xaxis.pane.set_edgecolor('none')
ax.yaxis.pane.set_edgecolor('none')
ax.zaxis.pane.set_edgecolor('none')

# labeling
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Cu Concentration, $x$")
ax.set_zlabel("Specific Heat, $C_V$ ($k_B$/atom)")
#ax.set_title("Specific Heat vs Temperature for Cu-alloys")
#ax.legend(loc="upper left", bbox_to_anchor=(-0.15, 0.9))

# custom padding in figure fraction units
pad_left   = -3
pad_right  = 9
pad_bottom = -0.5
pad_top    = 5  # negative trims more from top

new_bbox = tight_bbox.from_extents(
    tight_bbox.x0  - pad_left,
    tight_bbox.y0  - pad_bottom,
    tight_bbox.x1  + pad_right,
    tight_bbox.y1  + pad_top
)


plt.savefig("3d.pdf", bbox_inches=new_bbox)

print(np.max(all_cvs, axis=1))

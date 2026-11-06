import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m
import netCDF4 as nc
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import copy
import itertools
from cycler import cycler
font_size = 18
np.set_printoptions(suppress=True)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": font_size})

def flip(items, ncol):
    return list(itertools.chain(*[items[i::ncol] for i in range(ncol)]))

def inner_mod(a,b):
    res = a%b
    return res if not res else res-b if a<0 else res

colors = {
    "steel_blue": "#1F77B4",
    "light_steel_blue": "#AEC7E8",
    "orange": "#FF7F0E",
    "light_orange": "#FFBB78",
    "forest_green": "#2CA02C",
    "light_green": "#98DF8A",
    "firebrick_red": "#D62728",
    "soft_red": "#FF9896",
    "lavender": "#9467BD",
    "light_lavender": "#C5B0D5",
    "brown": "#8C564B",
    "tan": "#C49C94",
    "orchid": "#E377C2",
    "light_orchid": "#F7B6D2",
    "gray": "#7F7F7F",
    "light_gray": "#C7C7C7",
    "yellow_green": "#BCBD22",
    "light_yellow_green": "#DBDB8D",
    "turquoise": "#17BECF",
    "light_turquoise": "#9EDAE5"
}

custom_cmap = ListedColormap(colors.values())
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=custom_cmap.colors)

filename = "data/wl_dos_bins.nc"
bin_edges = nc.Dataset(filename)
bin_edges = np.array(bin_edges["grid data"][:], dtype=np.float64)

filename = "data/wl_dos.nc"
wl_logdos = nc.Dataset(filename)
wl_logdos = np.array(wl_logdos["grid data"][:], dtype=np.float64)

filename = "asro/rho_of_E.nc"  
rho_of_E = nc.Dataset(filename)

rho = rho_of_E.variables['rho data']
U_data = rho_of_E.variables['U data']

filename = "load_balance/wl_lb_bins.dat"
wl_lb_bins = nc.Dataset(filename)
wl_lb_bins = np.array(wl_lb_bins["grid data"][:], dtype=np.int32).T

try:
  idx = np.where(wl_lb_bins[:,0] == -1)[0][0]
except:
  idx = np.shape(wl_lb_bins)[0]

elements = []
concentrations = None

with open("brawl.inp", "r") as file:
    for line in file:
        if "species_names" in line:
            elements = line.split("=")[1].strip().split()
            n_species = len(elements)
        elif "!species_concentrations" in line:
            concentrations = np.array([float(x) for x in line.split("=")[1].strip().split()])
        elif "species_concentrations" in line:
            concentrations = np.array([float(x) for x in line.split("=")[1].strip().split()])

n_species = len(elements)
if concentrations is None:
    concentrations = np.ones(n_species) / n_species

print(concentrations)

E = 0.5 * (bin_edges[1:] + bin_edges[:-1])

# log DOS
plt.figure()
plt.plot(E, wl_logdos)
plt.xlabel("Energy")
plt.ylabel(r"$\log g(E)$")
plt.title("Density of States")
plt.tight_layout()
plt.savefig("1_wl_logdos.pdf")

# first derivative
dS_dE = np.gradient(wl_logdos, E)
plt.figure()
plt.plot(E, dS_dE)
plt.xlabel("Energy")
plt.ylabel(r"$\beta(E) = \frac{d\log g}{dE}$")
plt.title("Microcanonical Temperature")
plt.tight_layout()
plt.savefig("2_wl_betaE.pdf")

# second derivative
d2S_dE2 = np.gradient(dS_dE, E)
plt.figure()
plt.plot(E, d2S_dE2)
plt.xlabel("Energy")
plt.ylabel(r"$\frac{d^2\log g}{dE^2}$")
plt.title("Second Derivative of log DOS")
plt.tight_layout()
plt.savefig("3_wl_d2S_dE2.pdf")


filename = "data/wl_hist.nc"
wl_hist = nc.Dataset(filename)
wl_hist = np.array(wl_hist["grid data"][:], dtype=np.float64)

kb_ev = 8.617333262e-5
ev_to_ry = 13.605693122
kb_ry = kb_ev/ev_to_ry
n_atoms = 2048
ev_to_mev = 1000
# unit conversion
#kb_ev = kb_ry
bin_edges = bin_edges*ev_to_ry # converts from Ryd to eV

orig_temp = 1000
start_temp = 0 
end_temp = 8000

step_size = 25
temperatures = np.arange(start_temp, end_temp+step_size, step_size).astype(np.float64)
# Non-uniform plot array, denser at low temperatures
temperatures_plot = np.unique(
    np.concatenate((
        np.array([50, 100, 200, 400, 1000]),
        np.arange(start_temp, end_temp + 1, 1000)
    ))
).astype(np.int32)
temperatures[0] = 1e-64
print(temperatures_plot)

# ------------------------------------------------

bin_width = bin_edges[1] - bin_edges[0]

wcs = np.zeros((len(U_data), n_species, n_species, 2))
for i in range(len(U_data)):
    for j in range(n_species):       # central species
        for k in range(n_species):   # neighbor species
            wcs[i, j, k, 0] = 1.0 - rho[i, 1, j, k] / (12.0 * concentrations[j])
            wcs[i, j, k, 1] = 1.0 - rho[i, 2, j, k] / (6.0 * concentrations[j])

pairs = np.zeros([int(n_species*(n_species+1)/2),2], dtype=np.int16)
k = 0
for i in range(n_species):
    for j in range(i, n_species):
        pairs[k] = [i, j]
        k += 1

fig_height = 4
fig_width = fig_height * 1.68
fig_height = 4.5

cv_fig, [hist_ax, cv_ax1, cv_ax2] = plt.subplots(1, 3, figsize=(fig_width*3.6,fig_height*1.2), constrained_layout=True)
cv_ax1_asro = cv_ax1.twinx()
cv_ax2_asro = cv_ax2.twinx()

# Initialise arrays 
mean_energies   = np.zeros(len(temperatures))
gibbs_energies  = np.zeros(len(temperatures))
heat_caps       = np.zeros(len(temperatures))
entropies       = np.zeros(len(temperatures))
asr_orders_1    = np.zeros((len(pairs),len(temperatures)))
asr_orders_2    = np.zeros((len(pairs),len(temperatures)))

beta = 1.0/(kb_ev*temperatures[-1])
prob = np.zeros(len(bin_edges)-1)
# Reweight histogram
for ibin, edge in enumerate(bin_edges[:-1]):
    bin_energy = edge + 0.5*bin_width
    prob[ibin] = wl_logdos[ibin]-beta*bin_energy
prob = prob-np.max(prob)
prob = np.exp(prob)
# Normalise
prob = prob/(np.sum(prob))
zero_energy = (bin_edges[:-1][np.argmax(prob)]+0.5*bin_width)/n_atoms*ev_to_mev

hist_min = 0
hist_max = 0
prob_max = 0
# Loop over temperatures of interest
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

    # Only plot every 5th histogram to avoid crowding the axes
    if np.isin(temperatures_plot, int(new_temp)).any() == True:
        strlabel = r"$T$={} K".format(int(new_temp))
        index_to_zero = np.where(prob < 1e-8)
        hist_prob = copy.deepcopy(prob)
        hist_prob[index_to_zero] = 0
        non_zero = np.nonzero(hist_prob)
        hist_ax.stairs(prob, bin_edges/n_atoms*ev_to_mev-zero_energy, fill=True, alpha=0.1)
        hist_ax.stairs(prob, bin_edges/n_atoms*ev_to_mev-zero_energy, label=strlabel, fill=False, alpha=1)
        if (bin_edges[np.max(non_zero)+1]/n_atoms*ev_to_mev-zero_energy > hist_max):
          hist_max = bin_edges[np.max(non_zero)+1]/n_atoms*ev_to_mev-zero_energy
        if (bin_edges[np.min(non_zero)]/n_atoms*ev_to_mev-zero_energy < hist_min):
          hist_min = bin_edges[np.min(non_zero)]/n_atoms*ev_to_mev-zero_energy
        if (np.max(prob) < 0.4 and np.max(prob) > prob_max):
            prob_max = np.max(prob)

    # Mean energy
    mean_energy = np.dot(bin_edges[:-1]+0.5*bin_width, prob)
    mean_energies[itemp] = mean_energy

    for ipair, pair in enumerate(pairs):
      a = 1
      asr_order_1 = np.dot(wcs[:,pair[0],pair[1],0], prob)
      asr_orders_1[ipair,itemp] = asr_order_1
      asr_order_2 = np.dot(wcs[:,pair[0],pair[1],1], prob)
      asr_orders_2[ipair,itemp] = asr_order_2
    # Compute Specific Heat using the histogram
    msq_dev = np.zeros(len(bin_edges)-1)
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        msq_dev[ibin] = (bin_energy - mean_energies[itemp])**2
        
    heat_caps[itemp] = np.dot(msq_dev, prob)/(kb_ev*new_temp**2)

heat_caps = np.gradient(mean_energies, temperatures)

gibbs_energies[0] = mean_energies[0]
gibbs_energies[1] = mean_energies[0]
for itemp in range(2,len(temperatures)):

  beta_i = 1.0/(kb_ev*temperatures[itemp])
  beta_j = 1.0/(kb_ev*temperatures[itemp-1])

  gibbs_energies[itemp] = beta_j*gibbs_energies[itemp-1]/beta_i+((mean_energies[itemp-1]+mean_energies[itemp])*(beta_i-beta_j))/(2*beta_i)

for itemp, new_temp in enumerate(temperatures):
  entropies[itemp] = (mean_energies[itemp]-gibbs_energies[itemp])/new_temp


# Complete plots using data computed above
local_max_indices = np.where((heat_caps[1:-1] > heat_caps[:-2]) & (heat_caps[1:-1] > heat_caps[2:]))[0] + 1
local_min_indices = np.where((heat_caps[1:-1] < heat_caps[:-2]) & (heat_caps[1:-1] < heat_caps[2:]))[0] + 1

# unit conversion
mean_energies = mean_energies/n_atoms*ev_to_mev
heat_caps = heat_caps/kb_ev/n_atoms
entropies = entropies/kb_ev/n_atoms

y_offset = -0.175
x_offset = 0.5

mean_energies_diff = np.max(mean_energies) - np.min(mean_energies)
heat_cap_diff = np.max(heat_caps) - np.min(heat_caps)
entropies_diff = np.max(entropies) - np.min(entropies)

y_offset = -0.135
x_offset = 0.5

hist_ax.tick_params(direction="in")
cv_ax1.tick_params(direction="in")
cv_ax2.tick_params(direction="in")

label_size_offset = 3
hist_ax.set_xlabel(r'Internal energy, $E$ (meV/atom)', fontsize=font_size+label_size_offset)
hist_ax.set_ylabel(r'Probability density, $P(E)$ (meV$^{-1}$)', fontsize=font_size+label_size_offset)
cv_ax1.set_xlabel(r'Temperature, $T$ (K)', fontsize=font_size+label_size_offset)
cv_ax2.set_xlabel(r'Temperature, $T$ (K)', fontsize=font_size+label_size_offset)
cv_ax1.set_ylabel(r'Specific heat, $C_V$ ($k_B$/atom)', fontsize=font_size+label_size_offset)
cv_ax2.set_ylabel(r'Specific heat, $C_V$ ($k_B$/atom)', fontsize=font_size+label_size_offset)
cv_ax1_asro.set_ylabel(r'$\alpha^{pq}_1$', fontsize=font_size+label_size_offset)
cv_ax2_asro.set_ylabel(r'$\alpha^{pq}_2$', fontsize=font_size+label_size_offset)
hist_ax.grid(True, axis='x')
cv_ax1.grid(True, axis='x')
cv_ax2.grid(True, axis='x')

asro_max = 0
asro_min = 0
cv_ax1.plot(temperatures, heat_caps, '-', color='black', label=r"$C_V$")
for ipair, pair in enumerate(pairs):
  cv_ax1_asro.plot(temperatures, asr_orders_1[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])
  if (np.max(asr_orders_1[ipair]) > asro_max):
    asro_max = np.max(asr_orders_1[ipair])
  if (np.min(asr_orders_1[ipair]) < asro_min):
    asro_min = np.min(asr_orders_1[ipair])
cv_ax2.plot(temperatures, heat_caps, '-', color='black', label=r"$C_V$")
for ipair, pair in enumerate(pairs):
  cv_ax2_asro.plot(temperatures, asr_orders_2[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])
  if (np.max(asr_orders_2[ipair]) > asro_max):
    asro_max = np.max(asr_orders_2[ipair])
  if (np.min(asr_orders_2[ipair]) < asro_min):
    asro_min = np.min(asr_orders_2[ipair])

cv_ticks = 1000
x_ticks_subplots = np.arange(np.around(start_temp/cv_ticks, decimals=0)*cv_ticks, np.around(end_temp/cv_ticks, decimals=0)*cv_ticks+cv_ticks, cv_ticks)
cv_ax1.set_xticks(x_ticks_subplots)
cv_ax1.tick_params(axis='both', pad=6)
cv_ax2.set_xticks(x_ticks_subplots)
cv_ax2.tick_params(axis='both', pad=6)
x_min, x_max = hist_ax.get_xlim()
rough_interval = (x_max - x_min) / 10
magnitude = 10 ** np.floor(np.log10(rough_interval))
normalized = rough_interval / magnitude
if normalized < 1.5:
    interval = magnitude
elif normalized < 3.5:
    interval = 2 * magnitude
elif normalized < 7:
    interval = 5 * magnitude
else:
    interval = 10 * magnitude
ticks = np.arange(int(np.floor(x_min / interval) * interval), 
                   int(np.ceil(x_max / interval) * interval) + interval, 
                   int(interval))
hist_ax.set_xticks(ticks)

x_ticks_subplots = np.arange(np.around((bin_edges[0]/n_atoms*ev_to_mev-zero_energy)/10, decimals=0)*10, np.around((bin_edges[-1]/n_atoms*ev_to_mev-zero_energy)/10, decimals=0)*10+10, 10)
if len(x_ticks_subplots) > 10:
    x_ticks_subplots = x_ticks_subplots[::2]
#hist_ax.set_xticks(x_ticks_subplots)

handles, labels = hist_ax.get_legend_handles_labels()
hist_ax_legend = hist_ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(-0.6, 1), ncol=1)
cv_ax2.legend(loc='upper left', bbox_to_anchor=(1.2, 1), ncol=1)
cv_ax2_asro.legend(loc='upper left', bbox_to_anchor=(1.2, 0.875), ncol=1)

# Loop through the legend handles and change their linewidth
for handle in hist_ax_legend.legend_handles:
    handle.set_linewidth(3)  # Increase the line width for each legend line

alternate_colors = list(colors.values())[::2]
cv_ax1_asro.set_prop_cycle(cycler(color=alternate_colors))

margin_pct = 0.05
for ax in [hist_ax, cv_ax1, cv_ax1_asro]:
  # Get axis limits
  x_min, x_max = ax.get_xlim()
  y_min, y_max = ax.get_ylim()

  # Calculate the margin values (5% of the range)
  x_margin = (x_max - x_min) * margin_pct
  y_margin = (y_max - y_min) * margin_pct

  # Get the current tick positions
  x_ticks = ax.get_xticks()
  y_ticks = ax.get_yticks()

  # Filter out ticks within the margin
  x_ticks_filtered = [tick for tick in x_ticks if x_min + x_margin <= tick <= x_max - x_margin or tick == 0]
  y_ticks_filtered = [tick for tick in y_ticks if y_min + y_margin <= tick <= y_max - y_margin or tick == 0]


  # Set the new ticks
  #ax.set_xticks(x_ticks_filtered)
  #ax.set_yticks(y_ticks_filtered)

asro_diff = np.abs(asro_min - asro_max)
hist_ax.set_xlim(hist_min,hist_max)
hist_ax.set_ylim(0, prob_max*1.01)
cv_ax1.set_xlim(start_temp, end_temp)
cv_ax1.set_ylim(0, cv_ax1.get_ylim()[1]*1.01)
cv_ax1_asro.set_ylim(asro_min-0.01*asro_diff, asro_max+0.01*asro_diff)
cv_ax2.set_xlim(start_temp, end_temp)
cv_ax2.set_ylim(0, cv_ax2.get_ylim()[1]*1.01)
cv_ax2_asro.set_ylim(asro_min-0.01*asro_diff, asro_max+0.01*asro_diff)

print("Hist y-limit", prob_max*1.01)
print("SHC y-limit", cv_ax1.get_ylim()[1])
print("ASRO y-limit", asro_min-0.01*asro_diff, asro_max+0.01*asro_diff)
custom_h = prob_max*1.01
custom_l = asro_min-0.01*asro_diff
custom_u = asro_max+0.01*asro_diff
custom_c = cv_ax1.get_ylim()[1]

cv_ax1.set_ylim(0, custom_c)
cv_ax1_asro.set_ylim(custom_l, custom_u)
cv_ax2.set_ylim(0, custom_c)
cv_ax2_asro.set_ylim(custom_l, custom_u)
hist_ax.set_ylim(0, custom_h)

cv_fig.savefig('{}.pdf'.format(os.path.basename(os.getcwd())), bbox_inches='tight')

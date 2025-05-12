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
font_size = 13
np.set_printoptions(suppress=True)
plt.rcParams.update({"font.size": font_size})

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

filename = "data/wl_dos_bins.dat"
bin_edges = nc.Dataset(filename)
bin_edges = np.array(bin_edges["grid data"][:], dtype=np.float64)

filename = "data/wl_dos.dat"
wl_logdos = nc.Dataset(filename)
wl_logdos = np.array(wl_logdos["grid data"][:], dtype=np.float64)
plt.plot(wl_logdos)

filename = "asro/rho_of_E.dat"  
rho_of_E = nc.Dataset(filename)

rho = rho_of_E.variables['rho data'][:]
U_data = rho_of_E.variables['U data'][:]

elements = ""  # Initialize an empty string to store species data
with open("brawl.inp", 'r') as file:
    for line in file:
        if 'species_name' in line:  # Check if the line contains the species_name variable
            parts = line.split("=")  # Split the line by '=' and get everything after it
            if len(parts) > 1:
                elements = parts[1].strip().split()  # Get the data after '=' and strip spaces
            break  # Exit the loop after finding species_name
n_species = len(elements)
concentrations = 1.0/n_species*np.ones(n_species)

wl_logdos = wl_logdos - np.max(wl_logdos)

filename = "data/wl_hist.dat"
wl_hist = nc.Dataset(filename)
wl_hist = np.array(wl_hist["grid data"][:], dtype=np.float64)

kb_ev = 8.167333262e-5
ev_to_ry = 13.605693122
kb_ry = kb_ev/ev_to_ry
n_atoms = 128
ev_to_mev = 1000
# unit conversion
#kb_ev = kb_ry
bin_edges = bin_edges*ev_to_ry # converts from Ryd to eV

orig_temp = 3000
start_temp = 0 
end_temp = 3000

step_size = 25
temperatures = np.arange(start_temp, end_temp+step_size, step_size).astype(np.float64)
temperatures_plot = np.arange(start_temp, end_temp+200, 200).astype(np.int32)
temperatures[0] = 1e-64
print(temperatures_plot)

# ------------------------------------------------

bin_width = bin_edges[1] - bin_edges[0]

wcs = np.zeros((len(U_data),n_species,n_species,2))
for i in range(len(U_data)):
  wcoft = np.zeros((n_species,n_species,2))
  for j in range(len(elements)):
      for k in range(len(elements)):
          wcoft[j,k,0] = 1.0-1.0/8.0*rho[i,1,j,k]/concentrations[j]
          wcoft[j,k,1] = 1.0-1.0/6.0*rho[i,2,j,k]/concentrations[j]
  wcs[i] = wcoft

pairs = np.zeros([int(len(elements)*(len(elements)+1)/2),2], dtype=np.int16)
k = 0
for i in range(len(elements)):
    for j in range(i, len(elements)):
        pairs[k] = [i, j]
        k += 1

fig_height = 4; fig_width=fig_height*1.68

cv_fig, [hist_ax, cv_ax1] = plt.subplots(2, 1, figsize=(fig_width*1.2,fig_height*2.4), constrained_layout=True)
cv_ax1_asro = cv_ax1.twinx()

# Initialise arrays 
mean_energies   = np.zeros(len(temperatures))
gibbs_energies  = np.zeros(len(temperatures))
heat_caps       = np.zeros(len(temperatures))
entropies       = np.zeros(len(temperatures))
asr_orders_1      = np.zeros((len(pairs),len(temperatures)))

beta = 1.0/(kb_ev*temperatures[-1])
prob = np.zeros(len(bin_edges)-1)
# Reweight histogram
for ibin, edge in enumerate(bin_edges[:-1]):
    bin_energy = edge + 0.5*bin_width
    prob[ibin] = wl_logdos[ibin]-beta*bin_energy
prob = prob-np.max(prob)/2
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
        strlabel = "T={} K".format(int(new_temp))
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
      asr_order_1 = np.dot(wcs[:,pair[0],pair[1],0], prob)
      asr_orders_1[ipair,itemp] = asr_order_1

    # Compute heat capacity using the histogram
    msq_dev = np.zeros(len(bin_edges)-1)
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        msq_dev[ibin] = (bin_energy - mean_energies[itemp])**2
        
    heat_caps[itemp] = np.dot(msq_dev, prob)*bin_width/(kb_ev*new_temp**2)

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

label_size_offset = 3
hist_ax.set_title("Wang-Landau sampling on AlTiCrMo", fontsize=font_size+label_size_offset)
hist_ax.set_xlabel(r'Internal Energy $E$ (meV/atom)', fontsize=font_size+label_size_offset)
hist_ax.set_ylabel(r'Probability Density P($E$) (meV$^{-1}$)', fontsize=font_size+label_size_offset)
cv_ax1.set_ylabel(r'$C_v$ ($k_B$/atom)', fontsize=font_size+label_size_offset)
cv_ax1.set_xlabel(r'Temperature (K)', fontsize=font_size+label_size_offset)
cv_ax1_asro.set_ylabel(r'$\alpha^{pq}_1$', fontsize=font_size+label_size_offset)
cv_ax1.grid(True, axis='x')

asro_max = 0
asro_min = 0
cv_ax1.plot(temperatures, heat_caps, '-o', markersize=4, label="Heat Capacity")
for ipair, pair in enumerate(pairs):
  cv_ax1_asro.plot(temperatures, asr_orders_1[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])
  if (np.max(asr_orders_1[ipair]) > asro_max):
    asro_max = np.max(asr_orders_1[ipair])
  if (np.min(asr_orders_1[ipair]) < asro_min):
    asro_min = np.min(asr_orders_1[ipair])

cv_ticks = 400
x_ticks_subplots = np.arange(np.around(start_temp/cv_ticks, decimals=0)*cv_ticks, np.around(end_temp/cv_ticks, decimals=0)*cv_ticks+cv_ticks, cv_ticks)
cv_ax1.set_xticks(x_ticks_subplots)
cv_ax1.tick_params(axis='both', pad=6)

x_ticks_subplots = np.arange(np.around((bin_edges[0]/n_atoms*ev_to_mev-zero_energy)/10, decimals=0)*10, np.around((bin_edges[-1]/n_atoms*ev_to_mev-zero_energy)/10, decimals=0)*10+10, 10)
hist_ax.set_xticks(x_ticks_subplots)

handles, labels = hist_ax.get_legend_handles_labels()
hist_ax_legend = hist_ax.legend(flip(handles, 4), flip(labels, 4), loc='upper center', bbox_to_anchor=(x_offset, y_offset-0.05), ncol=4)
cv_ax1.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset-0.25))
cv_ax1_asro.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset-0.0275), ncol=int(len(pairs)/2))

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

print("Hist y-limit", prob_max*1.01)
print("SHC y-limit", cv_ax1.get_ylim()[1])
print("ASRO y-limit", asro_min-0.01*asro_diff, asro_max+0.01*asro_diff)
custom_h = 0.03
custom_l = -1.04
custom_u = 1.04
custom_c = 1

cv_ax1.set_ylim(0, custom_c)
cv_ax1_asro.set_ylim(custom_l, custom_u)
hist_ax.set_ylim(0, custom_h)

cv_fig.savefig('{}.pdf'.format(''.join(elements)), bbox_inches='tight')

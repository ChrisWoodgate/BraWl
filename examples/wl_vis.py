import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m
import netCDF4 as nc
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import copy
np.set_printoptions(suppress=True)
plt.rcParams.update({"text.usetex": True,
                     "font.size": 12})

def inner_mod(a,b):
    res = a%b
    return res if not res else res-b if a<0 else res

colors = {
    "soft_blue": "#4F81BD",
    "muted_orange": "#F7A800",
    "soft_red": "#D14F5D",
    "gentle_green": "#6DAF69",
    "soft_purple": "#9B59B6",
    "dusty_pink": "#E17D85",
    "slate_blue": "#6A7F99",
    "warm_tan": "#D3B69B",
    "earthy_brown": "#8E735B",
    "moss_green": "#6A8A3B",
    "rich_gold": "#B88A2A",
    "charcoal_gray": "#4A4A48"
}

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

indexed_colors_values = {index: value for index, (key, value) in enumerate(colors.items())}
# Convert hex to RGB, then to HSV
rgb_colors = [mcolors.hex2color(color) for color in colors.values()]
hsv_colors = [mcolors.rgb_to_hsv(rgb) for rgb in rgb_colors]
# Sort colors based on the Hue value (hsv_colors[i][0])
sorted_colors_by_hue = [indexed_colors_values[i] for i in np.argsort([hsv[0] for hsv in hsv_colors])]
custom_cmap = ListedColormap(sorted_colors_by_hue)
custom_cmap = ListedColormap(colors.values())
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=custom_cmap.colors)

subfolders = [ f.name for f in os.scandir(os.getcwd()) if f.is_dir() ]
print("Available directories:")
print(subfolders)
directory = input("Input directory to pull data from: ")

filename = "{}/wl_dos_bins.dat".format(directory)
bin_edges = nc.Dataset(filename)
bin_edges = np.array(bin_edges["grid data"][:], dtype=np.float64)

filename = "{}/wl_dos.dat".format(directory)
wl_logdos = nc.Dataset(filename)
wl_logdos = np.array(wl_logdos["grid data"][:], dtype=np.float64)
#print(wl_logdos)

filename = "{}/radial_densities/rho_of_E.dat".format(directory)  
rho_of_E = nc.Dataset(filename)

rho = rho_of_E.variables['rho data'][:]
U_data = rho_of_E.variables['U data'][:]

elements = ['Al', 'Ti', 'Cr', 'Mo']
n_species = 4
concentrations = 1.0/n_species*np.ones(n_species)

wl_logdos = wl_logdos - np.max(wl_logdos)
wl_logdos = np.log(np.exp(wl_logdos)/np.sum(np.exp(wl_logdos)))

filename = "{}/wl_hist.dat".format(directory)
wl_hist = nc.Dataset(filename)
wl_hist = np.array(wl_hist["grid data"][:], dtype=np.float64)

kb_ev = 8.167333262e-5
ev_to_ry = 13.605693122
kb_ry = kb_ev/ev_to_ry
n_atoms = 432
ev_to_mev = 1000
J = 0.001
N = 6

# unit conversion
bin_edges = bin_edges*ev_to_ry # converts from Ryd to eV

orig_temp = 3000
start_temp = 400
end_temp = 3000

step_size = 25
temperatures = np.arange(start_temp, end_temp+step_size, step_size)
temperatures_plot = np.arange(start_temp, end_temp+(end_temp-start_temp)/15, (end_temp-start_temp)/15)
temperatures_plot = np.round(temperatures_plot, -2).astype(np.int32)
plt.plot(np.exp(wl_logdos))
#plt.show()
plt.close()
# ------------------------------------------------

bin_width = bin_edges[1] - bin_edges[0]

prob = np.zeros(len(bin_edges)-1)
beta = 1.0/(kb_ry*orig_temp)
for ibin, edge in enumerate(bin_edges[:-1]):
    bin_energy = edge + 0.5*bin_width
    prob[ibin] = np.exp(wl_logdos[ibin]-beta*bin_energy)
plt.plot(prob)
#plt.show()
plt.close()

wcs = np.zeros((len(U_data),n_species,n_species,2))
for i in range(len(U_data)):
  wcoft = np.zeros((n_species,n_species,2))
  for j in range(len(elements)):
      for k in range(len(elements)):
          wcoft[j,k,0] = 1.0-1.0/8.0*rho[i,1,j,k]/concentrations[j]
          wcoft[j,k,1] = 1.0-1.0/6.0*rho[i,2,j,k]/concentrations[j]
  wcs[i] = wcoft

elements = ""  # Initialize an empty string to store species data
with open("{}/input.txt".format(directory), 'r') as file:
    for line in file:
        if 'species_name' in line:  # Check if the line contains the species_name variable
            parts = line.split("=")  # Split the line by '=' and get everything after it
            if len(parts) > 1:
                elements = parts[1].strip().split()  # Get the data after '=' and strip spaces
            break  # Exit the loop after finding species_name

pairs = np.zeros([int(len(elements)*(len(elements)+1)/2),2], dtype=np.int16)
k = 0
for i in range(len(elements)):
    for j in range(i, len(elements)):
        pairs[k] = [i, j]
        k += 1

# Setup plots
#fig, [ax1, ax2, ax3, ax4] = plt.subplots(1, 4, figsize=(24, 7), constrained_layout=True)
#fig, [ax1, ax2, ax3] = plt.subplots(1, 3, figsize=(24, 9), constrained_layout=True)

fig_height = 4; fig_width=fig_height*1.68
fig1, ax1 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)
fig2, ax2 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)
fig3, ax3 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)
ax4 = ax3.twinx()
fig4, ax5 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)
ax6 = ax5.twinx()
fig5, ax7 = plt.subplots(figsize=(fig_width,fig_height), constrained_layout=True)

ax1.set_xlabel(r'Energy $U$ (meV/atom)')
ax1.set_ylabel(r'Probability P($U$)')
#ax1.set_title(r'Energy Histograms')

ax2.set_xlabel(r'Temperature (K)')
ax2.set_ylabel(r'$\langle U \rangle$ (meV/atom)')
#ax2.set_title(r'Mean Energy $\langle U \rangle$')  

ax3.set_xlabel(r'Temperature (K)')
ax3.set_ylabel(r'$C_v$ ($k_B$/atom)')
#ax3.set_title(r'Heat Capacity $C_v$') 

ax4.set_ylabel(r'$\alpha^{pq}_1$')

ax5.set_xlabel(r'Temperature (K)')
ax5.set_ylabel(r'$C_v$ ($k_B$/atom)')
#ax5.set_title(r'Heat Capacity $C_v$') 

ax6.set_ylabel(r'$\alpha^{pq}_2$')

ax7.set_xlabel(r'Temperature (K)')
ax7.set_ylabel(r'Entropy ($k_B$/atom)')
#ax7.set_title(r'Entropy $S$')

cv_fig, [hist_ax, cv_ax1, cv_ax2] = plt.subplots(3, 1, figsize=(fig_width,fig_height*3), constrained_layout=True)
cv_ax1_asro = cv_ax1.twinx()
cv_ax2_asro = cv_ax2.twinx()

# Initialise arrays 
mean_energies   = np.zeros(len(temperatures))
gibbs_energies  = np.zeros(len(temperatures))
heat_caps       = np.zeros(len(temperatures))
entropies       = np.zeros(len(temperatures))
asr_orders_1      = np.zeros((len(pairs),len(temperatures)))
asr_orders_2      = np.zeros((len(pairs),len(temperatures)))

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

# Loop over temperatures of interest
for itemp, new_temp in enumerate(temperatures):

    beta = 1.0/(kb_ev*new_temp)

    # Reweighted histogram
    prob = np.zeros(len(bin_edges)-1)
    
    # Reweight histogram
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        prob[ibin] = wl_logdos[ibin]-beta*bin_energy
   
    prob = prob-np.max(prob)/2
    prob = np.exp(prob)

    # Normalise
    prob = prob/(np.sum(prob))
    
    # Only plot every 5th histogram to avoid crowding the axes
    if np.isin(temperatures_plot, int(new_temp)).any() == True:
        strlabel = "T={} K".format(new_temp)
        index_to_zero = np.where(prob < 1e-5)
        hist_prob = copy.deepcopy(prob)
        hist_prob[index_to_zero] = 0
        non_zero = np.nonzero(hist_prob)
        ax1.stairs(hist_prob[np.min(non_zero):np.max(non_zero)], bin_edges[np.min(non_zero):np.max(non_zero)+1]/n_atoms*ev_to_mev-zero_energy, label=strlabel, fill=True)
        hist_ax.stairs(hist_prob[np.min(non_zero):np.max(non_zero)], bin_edges[np.min(non_zero):np.max(non_zero)+1]/n_atoms*ev_to_mev-zero_energy, label=strlabel, fill=True)

    # Mean energy
    mean_energy = np.dot(bin_edges[:-1]+0.5*bin_width, prob)
    mean_energies[itemp] = mean_energy

    for ipair, pair in enumerate(pairs):
      asr_order_1 = np.dot(wcs[:,pair[0],pair[1],0], prob)
      asr_orders_1[ipair,itemp] = asr_order_1
      asr_order_2 = np.dot(wcs[:,pair[0],pair[1],1], prob)
      asr_orders_2[ipair,itemp] = asr_order_2

    # Compute heat capacity using the histogram
    msq_dev = np.zeros(len(bin_edges)-1)
    for ibin, edge in enumerate(bin_edges[:-1]):
        bin_energy = edge + 0.5*bin_width
        msq_dev[ibin] = (bin_energy - mean_energies[itemp])**2
        
    heat_caps[itemp] = np.dot(msq_dev, prob)*bin_width/(kb_ev*new_temp**2)

gibbs_energies[0] = mean_energies[0]
for itemp in range(1,len(temperatures)):

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

ax2.plot(temperatures, mean_energies, '-o', markersize=4)
ax3.plot(temperatures, heat_caps, '-o', markersize=4, label="Heat Capacity")
for ipair, pair in enumerate(pairs):
  ax4.plot(temperatures, asr_orders_1[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])
ax5.plot(temperatures, heat_caps, '-o', markersize=4, label="Heat Capacity")
for ipair, pair in enumerate(pairs):
  ax6.plot(temperatures, asr_orders_2[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])
ax7.plot(temperatures, entropies, '-o', markersize=4)

y_offset = -0.175
x_offset = 0.5
ax1.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset), ncol=4)

ax3.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset-0.225))
ax4.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset), ncol=int(len(pairs)/2))

ax5.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset-0.225))
ax6.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset), ncol=int(len(pairs)/2))

ax3.tick_params(direction="in")
ax4.tick_params(direction="in")
ax5.tick_params(direction="in")
ax6.tick_params(direction="in")

mean_energies_diff = np.max(mean_energies) - np.min(mean_energies)
heat_cap_diff = np.max(heat_caps) - np.min(heat_caps)
entropies_diff = np.max(entropies) - np.min(entropies)

for x in local_max_indices:
  ax2.vlines(x=temperatures[x], ymin=np.min(mean_energies)-0.1*mean_energies_diff, ymax=np.max(mean_energies), color=colors["firebrick_red"], linestyle='--')
  ax3.vlines(x=temperatures[x], ymin=np.min(heat_caps)-0.1*heat_cap_diff, ymax=np.max(heat_caps), color=colors["firebrick_red"], linestyle='--')
  ax5.vlines(x=temperatures[x], ymin=np.min(heat_caps)-0.1*heat_cap_diff, ymax=np.max(heat_caps), color=colors["firebrick_red"], linestyle='--')
  ax7.vlines(x=temperatures[x], ymin=np.min(entropies)-0.1*entropies_diff, ymax=np.max(entropies), color=colors["firebrick_red"], linestyle='--')

ax2.set_ylim(np.min(mean_energies)-0.025*np.abs(mean_energies_diff), np.max(mean_energies)+0.025*np.abs(mean_energies_diff))
ax3.set_ylim(np.min(heat_caps)-0.025*np.abs(heat_cap_diff), np.max(heat_caps)+0.025*np.abs(heat_cap_diff))
ax5.set_ylim(np.min(heat_caps)-0.025*np.abs(heat_cap_diff), np.max(heat_caps)+0.025*np.abs(heat_cap_diff))
ax7.set_ylim(np.min(entropies)-0.025*np.abs(entropies_diff), np.max(entropies)+0.025*np.abs(entropies_diff))

fig1.savefig('figures/{}_energy_histogram.svg'.format(''.join(elements)), bbox_inches='tight')
fig2.savefig('figures/{}_mean_energy.svg'.format(''.join(elements)), bbox_inches='tight')
fig3.savefig('figures/{}_heat_capacity_1.svg'.format(''.join(elements)), bbox_inches='tight')
fig4.savefig('figures/{}_heat_capacity_2.svg'.format(''.join(elements)), bbox_inches='tight')
fig5.savefig('figures/{}_entropy.svg'.format(''.join(elements)), bbox_inches='tight')

y_offset = -0.125
x_offset = 0.5

hist_ax.tick_params(direction="in")
cv_ax1.tick_params(direction="in")
cv_ax2.tick_params(direction="in")

hist_ax.set_xlabel(r'Energy $U$ (meV/atom)')
hist_ax.set_ylabel(r'Probability P($U$)')
cv_ax2.set_xlabel(r'Temperature (K)')
cv_ax1.set_ylabel(r'$C$ ($k_B$/atom)')
cv_ax2.set_ylabel(r'$C$ ($k_B$/atom)')
cv_ax1_asro.set_ylabel(r'$\alpha^{pq}_1$')
cv_ax2_asro.set_ylabel(r'$\alpha^{pq}_2$')
cv_ax1_asro.set_xticklabels([])
cv_ax1.grid(True, axis='x')
cv_ax2.grid(True, axis='x')

x_ticks_subplots = np.arange(np.around(start_temp/250, decimals=0)*250, np.around(end_temp/250, decimals=0)*250+250, 250)
cv_ax1.set_xticks(x_ticks_subplots)
cv_ax2.set_xticks(x_ticks_subplots)

x_ticks_subplots = np.arange(np.around((bin_edges[0]/n_atoms*ev_to_mev-zero_energy)/10, decimals=0)*10, np.around((bin_edges[-1]/n_atoms*ev_to_mev-zero_energy)/10, decimals=0)*10+10, 10)

hist_ax.set_xticks(x_ticks_subplots)

cv_ax1.plot(temperatures, heat_caps, '-o', markersize=4, label="Heat Capacity")
for ipair, pair in enumerate(pairs):
  cv_ax1_asro.plot(temperatures, asr_orders_1[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])
cv_ax2.plot(temperatures, heat_caps, '-o', markersize=4, label="Heat Capacity")
for ipair, pair in enumerate(pairs):
  cv_ax2_asro.plot(temperatures, asr_orders_2[ipair], label = elements[pair[0]] + '-' + elements[pair[1]])

hist_ax.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset), ncol=4)
cv_ax2.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset-0.2))
cv_ax2_asro.legend(loc='upper center', bbox_to_anchor=(x_offset, y_offset), ncol=int(len(pairs)/2))

margin_pct = 0.05
for ax in [hist_ax, cv_ax1, cv_ax2, cv_ax1_asro, cv_ax2_asro]:
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
  ax.set_xticks(x_ticks_filtered)
  ax.set_yticks(y_ticks_filtered)

hist_ax.set_ylim(0, hist_ax.get_ylim()[1])
cv_ax1.set_ylim(0, cv_ax1.get_ylim()[1])
cv_ax2.set_ylim(0, cv_ax2.get_ylim()[1])

cv_fig.savefig('figures/{}.svg'.format(''.join(elements)), bbox_inches='tight')

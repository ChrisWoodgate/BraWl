import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m
import netCDF4 as nc
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import itertools

font_size = 32
np.set_printoptions(suppress=True)
#plt.rcParams.update({"text.usetex": True,
#                     "font.size": font_size})
plt.rcParams.update({"font.size": font_size})
plt.rc('font', family='serif')
#plt.rc('text', usetex=True)
plt.rcParams["figure.figsize"] = (12*1.75,8*1.75)

def flip(items, ncol):
    return list(itertools.chain(*[items[i::ncol] for i in range(ncol)]))

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

#subfolders = [ f.name for f in os.scandir(os.getcwd()) if f.is_dir() ]
#print("Available directories:")
#print(subfolders)
directory = ""#input("Input directory to pull data from: ")

filename = "load_balance/wl_lb_bins.dat".format(directory)
wl_lb_bins = nc.Dataset(filename)
wl_lb_bins = np.array(wl_lb_bins["grid data"][:], dtype=np.float64).T

filename = "load_balance/wl_lb_avg_time.dat".format(directory)
wl_lb_avg_time = nc.Dataset(filename)
wl_lb_avg_time = np.array(wl_lb_avg_time["grid data"][:], dtype=np.float64).T

filename = "load_balance/wl_window_time.dat".format(directory)
wl_window_time = nc.Dataset(filename)
wl_window_time = np.array(wl_window_time["grid data"][:], dtype=np.float64)

filename = "load_balance/wl_lb_mc_steps.dat".format(directory)
wl_lb_mc_steps = nc.Dataset(filename)
wl_lb_mc_steps = np.array(wl_lb_mc_steps["grid data"][:], dtype=np.int64)

print(wl_lb_mc_steps, sum(wl_lb_mc_steps))

try:
  iter = np.where(wl_lb_bins[:,0] == -1)[0][0]
except:
  iter = np.shape(wl_lb_bins)[0]

filename = "load_balance/wl_lb_max_time.dat".format(directory)
wl_lb_max_time = nc.Dataset(filename)
wl_lb_max_time = np.array(wl_lb_max_time["grid data"][:], dtype=np.float64).T

print(np.sum(np.max(wl_lb_max_time, axis=1))/60)

window = np.shape(wl_lb_bins)[1]

x_axis = np.arange(1,iter+1,1).astype(np.int32)

columns = 8

for i in range(window):
  plt.plot(x_axis, wl_lb_bins[0:iter,i], label=i+1)
#plt.title("Bins per Window")
plt.ylabel("# of Bins per Domain")
plt.xlabel("Wang-Landau Iteration")
plt.xticks(x_axis)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(flip(handles, columns), flip(labels, columns), title="Domain", loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=columns)
plt.tight_layout()
plt.savefig('load_balance/load_balance_1.pdf', bbox_inches='tight')
plt.close()

time_adjusted = np.zeros(np.shape(wl_lb_avg_time))
for i in range(iter):
  time_adjusted[i] = wl_lb_avg_time[i]/(np.sum(wl_lb_avg_time[i])/window)
for i in range(window):
  plt.plot(x_axis, time_adjusted[0:iter,i], label=i+1)
plt.title("Window Average Time Ratio")
plt.ylabel("Time Ratio")
plt.xlabel("W-L Iteration")
plt.axhline(0.75, linestyle='--', color='red', alpha=0.5)
plt.axhline(1, linestyle='--', color='red')
plt.axhline(1.25, linestyle='--', color='red', alpha=0.5)
plt.xticks(x_axis)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend(flip(handles, columns), flip(labels, columns), loc='upper center', bbox_to_anchor=(0.5, -0.25), ncol=columns)
plt.tight_layout()
plt.savefig('load_balance/load_balance_2.pdf', bbox_inches='tight')
plt.close()

time_std = np.zeros([iter])
for i in range(iter):
  time_std[i] = np.std(wl_lb_avg_time[i]/(np.sum(wl_lb_avg_time[i])/window))
plt.title("Window Average Time Ratio Standard Deviation")
plt.ylabel("Standard Deviation")
plt.xlabel("W-L Iteration")
plt.plot(x_axis, time_std)
plt.xticks(x_axis)
plt.tight_layout()
plt.savefig('load_balance/load_balance_3.pdf', bbox_inches='tight')
plt.close()
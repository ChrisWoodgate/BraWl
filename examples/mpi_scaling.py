import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m
import netCDF4 as nc
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import itertools

font_size = 12
np.set_printoptions(suppress=True)
#plt.rcParams.update({"text.usetex": True,
#                     "font.size": font_size})
plt.rcParams.update({"font.size": font_size})
plt.rc('font', family='serif')
#plt.rc('text', usetex=True)

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

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

subfolders = [ f.name for f in os.scandir(os.getcwd()) if f.is_dir() ]
print("Available directories:")
print(subfolders)
elements = input("Input elements to pull data from: ")

windows = np.array([1, 2, 4, 6, 8])
directories = np.empty_like(windows, dtype='<U32')
sums = np.zeros(np.shape(windows))

for i in range(len(windows)):
  directories[i] = "{}_{:02d}".format(elements, windows[i])

for i, directory in enumerate(directories):
    filename = "{}/wl_lb_max_time.dat".format(directory)
    wl_lb_max_time = nc.Dataset(filename)
    wl_lb_max_time = np.array(wl_lb_max_time["grid data"][:], dtype=np.float64).T
    
    row_sums = np.max(wl_lb_max_time, axis=1)
    total_sum = np.sum(row_sums)
            
    sums[i] = total_sum
sums = sums/60/60
plt.plot(windows, sums)
plt.xscale('log', base=2)
plt.xticks(windows, labels=windows)
plt.xlabel("Windows")
plt.ylabel("Time Taken (hrs)")
plt.show()
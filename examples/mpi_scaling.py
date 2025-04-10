import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m
import netCDF4 as nc
import os
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import itertools
from scipy.optimize import curve_fit

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

windows = np.array([1, 2, 3, 4, 5, 6, 7, 8])
walkers = np.array([1, 2, 4, 8])
repeats = 4
directories = np.empty_like(windows, dtype='<U32')
sums = np.zeros([len(windows), repeats])

for i in range(len(windows)):
  directories[i] = "{}_{:02d}".format(elements, windows[i])

for i in range(len(windows)):
    total_sum = 0
    for j in range(1,repeats+1):
        filename = "{}_{:02d}_{}/wl_lb_max_time.dat".format(elements, windows[i], j)
        wl_lb_max_time = nc.Dataset(filename)
        wl_lb_max_time = np.array(wl_lb_max_time["grid data"][:], dtype=np.float64).T

        row_sums = np.max(wl_lb_max_time, axis=1)
        print(windows[i], j, row_sums)
        sums[i, j-1] = np.sum(row_sums)

sums /= 60
sums_mean = np.mean(sums, axis=1)
sums_error = np.std(sums, axis=1)

plt.errorbar(windows, sums_mean, yerr=sums_error, capsize=3, ecolor = "red")
plt.xticks(windows, labels=windows)
plt.xlabel("Windows")
plt.ylabel("Time Taken (mins)")
plt.savefig('../figures/{}_time.pdf'.format(''.join(elements)), bbox_inches='tight')
plt.show()

windows *= 1

sums = np.mean(sums[0])/sums
sums_mean = np.mean(sums, axis=1)
sums_error = np.std(sums, axis=1)

def linear(x, m):
    x0 = x[0]
    return m*x + (1-m*x0)

params, covariance = curve_fit(linear, windows, sums_mean)

plt.plot(windows, linear(windows, params))
plt.errorbar(windows, sums_mean, yerr=sums_error, capsize=3, ecolor = "#D62728", ls='none', fmt='o', color = "#D62728")
plt.xticks(windows, labels=windows)
plt.xlabel("Cores")
plt.ylabel("Speed Up")
plt.text(0.02, 0.98, f"Slope: {params[0]:.2f}", transform=plt.gca().transAxes, verticalalignment='top', horizontalalignment='left')
plt.savefig('../figures/{}_speedup.pdf'.format(''.join(elements)), bbox_inches='tight')
plt.show()
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math as m
import netCDF4 as nc
import os
import re
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import itertools
from scipy.optimize import curve_fit

font_size = 16
np.set_printoptions(suppress=True)
#plt.rcParams.update({"text.usetex": True,
#                     "font.size": font_size})
plt.rcParams.update({"font.size": font_size})
plt.rcParams["figure.figsize"] = (16,8)
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

methods = np.array([0, 1, 2, 3, 4, 5])
methods_label = np.array(["NU-LB-R", "NU-LB", "NU-R", "NU", "R", "Domains"])
windows = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
walkers = np.array([1, 2, 3, 4, 5, 6])
overlaps = np.array([0, 10, 25, 50, 75])
repeats = 5
bins=512

walkers = np.array([1])
#methods = np.array([0, 2, 4])
#methods_label = np.array(["Non-Uniform + Balance + Replica", "Non-Uniform + Replica", "Replica"])
#methods = np.array([0,1])
#walkers = np.array([1])
#overlaps = np.array([50])

sums = np.zeros([len(methods), len(walkers), len(windows), len(overlaps), repeats])
mc_steps = np.zeros([len(methods), len(walkers), len(windows), len(overlaps), repeats])
one_sums = np.zeros([len(methods), len(walkers), repeats])
for method_id, method in enumerate(methods):
  for walker in range(len(walkers)):
    for window in range(len(windows)):
      for overlap in range(len(overlaps)):
        for k in range(1,repeats+1):
          filename = "{}_{:02d}_{:02d}_{:02d}_{}/load_balance/wl_lb_max_time.dat".format(method, walkers[walker], windows[window], overlaps[overlap], k)
          wl_lb_max_time = nc.Dataset(filename)
          wl_lb_max_time = np.array(wl_lb_max_time["grid data"][:], dtype=np.float64).T

          row_sums = np.max(wl_lb_max_time, axis=1)
          sums[method_id, walker, window, overlap, k-1] = np.sum(row_sums)

          filename = "{}_{:02d}_{:02d}_{:02d}_{}/load_balance/wl_lb_mc_steps.dat".format(method, walkers[walker], windows[window], overlaps[overlap], k)
          wl_lb_mc_steps = nc.Dataset(filename)
          wl_lb_mc_steps = np.array(wl_lb_mc_steps["grid data"][:], dtype=np.int64)
          mc_steps[method_id, walker, window, overlap, k-1] = np.sum(wl_lb_mc_steps)

# REMOVE THIS LATER
#sums[0, 5, 8, 0, 0] = 250

one_sums = sums[:, :, 0, :, :]

sums_mean = np.mean(sums, axis=4)
sums_min = np.min(sums, axis=4)
sums_err = np.std(sums, axis=4)/np.sqrt(repeats)

mc_steps_mean = np.mean(mc_steps, axis=4)
mc_steps_err = np.std(mc_steps, axis=4)/np.sqrt(repeats)

one_sums_err = np.zeros([len(walkers),len(overlaps)])
one_sums_mean = np.zeros([len(walkers),len(overlaps)])

for walker in range(len(walkers)):
  for overlap in range(len(overlaps)):
    one_sums_err[walker, overlap] = np.std(one_sums[:, walker, overlap, :])/np.sqrt(len(methods)*repeats)
    one_sums_mean[walker, overlap] = np.mean(one_sums[:, walker, overlap, :])

    sums_mean[:, walker, 0, overlap] = one_sums_mean[walker, overlap]
    sums_err[:, walker, 0, overlap] = one_sums_err[walker, overlap]

    mc_steps_mean[:, walker, 0, overlap] = np.mean(mc_steps[:, walker, 0, overlap, :])
    mc_steps_err[:, walker, 0, overlap] = np.std(mc_steps[:, walker, 0, overlap, :])/np.sqrt(len(methods)*repeats)

mc_steps_err = mc_steps_err/mc_steps_mean

for method in range(len(methods)):
  for walker in range(len(walkers)):
    for overlap in range(len(overlaps)):
      mc_steps_mean[method, walker, :, overlap] = mc_steps_mean[method, walker, :, overlap]/mc_steps_mean[method, walker, 0, overlap]

mc_steps_err = mc_steps_err*mc_steps_mean

efficiency = np.empty_like(sums_mean)
efficiency_err = np.empty_like(sums_err)

# Window Efficiency
for overlap in range(len(overlaps)):
  max_y = 0
  for method_id, method in enumerate(methods):
    for walker in range(len(walkers)):
      efficiency[method_id, walker, :, overlap] = one_sums_mean[walker, overlap]/(sums_mean[method_id, walker, :, overlap]*windows)
      efficiency_err[method_id, walker, :, overlap] = efficiency[method_id, walker, :, overlap]*np.sqrt((one_sums_err[walker, overlap]/one_sums_mean[walker, overlap])**2+(sums_err[method_id, walker, :, overlap]/sums_mean[method_id, walker, :, overlap])**2)
      if max_y < np.max(efficiency[method_id, walker, :, overlap]):
        max_y = np.max(efficiency[method_id, walker, :, overlap])
  for method_id, method in enumerate(methods):
    plt.axhline(y=1, linestyle='-')
    for walker in range(len(walkers)):
      plt.errorbar(windows, efficiency[method_id, walker, :, overlap], yerr=efficiency_err[method_id, walker, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(walker+1))
    plt.xticks(windows, labels=windows)
    plt.gca().set_box_aspect(1)
    plt.xlabel("Domains")
    plt.ylabel("Domain Efficiency")
    plt.title("Method {}".format(methods_label[method_id]))
    plt.legend(loc="upper right", title="Walkers")
    plt.tight_layout()

    plt.ylim(0, max_y*1.25)
    plt.savefig('{}_{:02d}_window_efficiency.pdf'.format(method, overlaps[overlap]), bbox_inches='tight')
    plt.close()

# Speed Up
for overlap in range(len(overlaps)):
  for method in range(len(methods)):
    for walker in range(len(walkers)):
      efficiency[method, walker, :, overlap] = one_sums_mean[0, overlap]/(sums_mean[method, walker, :, overlap])
      efficiency_err[method, walker, :, overlap] = efficiency[method, walker, :, overlap]*np.sqrt((one_sums_err[0, overlap]/one_sums_mean[0, overlap])**2+(sums_err[method, walker, :, overlap]/sums_mean[method, walker, :, overlap])**2)

    plt.errorbar(windows*1, efficiency[method, 0, :, overlap], yerr=efficiency_err[method, 0, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(methods_label[method]))

  plt.plot(windows, windows, color="#D62728")
  plt.plot(windows, np.arange(1, np.max(windows)*0.75+np.max(windows)*0.75/len(windows), np.max(windows)*0.75/len(windows)), linestyle="--", color="#FF9896")
  plt.xticks(windows, labels=windows)
  plt.xlabel("Cores")
  plt.ylabel("Speed-up")
  plt.title("Speed-up for 1 walker per domain")
  plt.legend(loc="upper left", title="Methods")
  plt.tight_layout()
  plt.savefig('{:02d}_method_speedup.pdf'.format(overlaps[overlap]), bbox_inches='tight')
  plt.close()

# MC Steps
for overlap in range(len(overlaps)):
  overlap = 1
  max_y = 0
  for method in range(len(methods)):
  #for method in [0, 2, 4]:
    if max_y < np.max(mc_steps_mean[method, 0, :, overlap]):
      max_y = np.max(mc_steps_mean[method, 0, :, overlap])

for overlap in range(len(overlaps)):
  for method in range(len(methods)):
    plt.axhline(y=1, linestyle='-')
    plt.errorbar(windows, mc_steps_mean[method, 0, :, overlap], yerr=mc_steps_err[method, 0, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(methods_label[method]))
    plt.xticks(windows, labels=windows)
    plt.xlabel("Cores")
    plt.ylabel("Relative MC Steps")
    plt.title("Relative MC Steps")
    plt.tight_layout()
    plt.gca().set_box_aspect(1)
    plt.ylim(0, max_y*1.25)
    plt.savefig('{}_{:02d}_mc_steps.pdf'.format(methods[method], overlaps[overlap]), bbox_inches='tight')
    plt.close()

# Speed Up
for method in range(len(methods)):
  for overlap in range(len(overlaps)):
    for walker in range(len(walkers)):
      efficiency[method, walker, :, overlap] = one_sums_mean[0, overlap]/(sums_mean[method, walker, :, overlap])
      efficiency_err[method, walker, :, overlap] = efficiency[method, walker, :, overlap]*np.sqrt((one_sums_err[0, overlap]/one_sums_mean[0, overlap])**2+(sums_err[method, walker, :, overlap]/sums_mean[method, walker, :, overlap])**2)

    plt.errorbar(windows*1, efficiency[method, 0, :, overlap], yerr=efficiency_err[method, 0, :, overlap], capsize=3, ls='none', fmt='o', label="{}%".format(overlaps[overlap]))

  plt.plot(windows, windows, color="#D62728")
  plt.plot(windows, np.arange(1, np.max(windows)*0.75+np.max(windows)*0.75/len(windows), np.max(windows)*0.75/len(windows)), linestyle="--", color="#FF9896")
  plt.xticks(windows, labels=windows)
  plt.xlabel("Cores")
  plt.ylabel("Speed-up")
  plt.title("Speed-up for 1 walker per domain")
  plt.legend(loc="upper left", title="Overlaps")
  plt.tight_layout()
  plt.savefig('{}_overlap_speedup.pdf'.format(methods[method]), bbox_inches='tight')
  plt.close()

cores = np.linspace(1, int(np.max(windows)*np.max(walkers)), int(np.max(windows)*np.max(walkers)))
for overlap in range(len(overlaps)):
  max_y = 0
  for method in range(len(methods)):
    for walker in range(len(walkers)):
      plt.errorbar(windows*(walker+1), efficiency[method, walker, :, overlap], yerr=efficiency_err[method, walker, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(walkers[walker]))
      if max_y < np.max(efficiency[method, walker, :, overlap]):
        max_y = np.max(efficiency[method, walker, :, overlap])
    
    plt.plot(cores, cores, color="#D62728")
    plt.plot(cores, np.arange(1, np.max(cores)*0.75+np.max(cores)*0.75/len(cores), np.max(cores)*0.75/len(cores)), linestyle="--", color="#FF9896")
    #plt.xticks(cores, labels=cores)
    plt.ylim(0, max_y*1.25)
    plt.xlabel("Cores")
    plt.ylabel("Speed-up")
    plt.title("Speed-up for Method {}".format(methods_label[method]))
    plt.legend(loc="lower right", title="Walkers")
    plt.tight_layout()
    plt.savefig('{}_speedup.pdf'.format(methods[method]), bbox_inches='tight')
    plt.close()

one_sums = np.zeros([len(methods), len(windows), len(overlaps), repeats])
one_sums = sums[:, 0, :, :, :]
one_sums_err = np.zeros([len(windows),len(overlaps)])
one_sums_mean = np.zeros([len(windows),len(overlaps)])

for window in range(len(windows)):
  one_sums_err[window] = np.std(one_sums[:, window, :])/np.sqrt(len(methods)*repeats)
  one_sums_mean[window] = np.mean(one_sums[:, window, :])

  sums_mean[:, 0, window] = one_sums_mean[window]
  sums_err[:, 0, window] = one_sums_err[window]

for window in range(len(windows)):
  for overlap in range(len(overlaps)):
    one_sums_err[window, overlap] = np.std(one_sums[:, window, overlap, :])/np.sqrt(len(methods)*repeats)
    one_sums_mean[window, overlap] = np.mean(one_sums[:, window, overlap, :])

    sums_mean[:, 0, window, overlap] = one_sums_mean[window, overlap]
    sums_err[:, 0, window, overlap] = one_sums_err[window, overlap]

efficiency = np.empty_like(sums_mean)
efficiency_err = np.empty_like(sums_err)


# Walker Efficiency
for overlap in range(len(overlaps)):
  max_y = 0
  for method in range(len(methods)):
    plt.axhline(y=1, linestyle='-')
    for window in range(len(windows)):
      efficiency[method, :, window, overlap] = one_sums_mean[window, overlap]/(sums_mean[method, :, window, overlap]*walkers)
      efficiency_err[method, :, window, overlap] = efficiency[method, :, window, overlap]*np.sqrt((one_sums_err[window, overlap]/one_sums_mean[window, overlap])**2+(sums_err[method, :, window, overlap]/sums_mean[method, :, window, overlap])**2)

    for window in np.array([0, 1, 3, 7, 15]):
      plt.errorbar(walkers, efficiency[method, :, window, overlap], yerr=efficiency_err[method, :, window, overlap], capsize=3, ls='none', fmt='o', label="{}".format(window+1))
      if max_y < np.max(efficiency[method, :, window, overlap]):
        max_y = np.max(efficiency[method, :, window, overlap])

    plt.xticks(walkers, labels=walkers)
    plt.xlabel("Walkers")
    plt.ylabel("Walker Efficiency")
    plt.title("Method {}".format(methods[method]))
    plt.legend(loc="upper right", title="Windows")
    plt.tight_layout()
    plt.gca().set_box_aspect(1)
    plt.ylim(0, max_y*1.25)
    plt.savefig('{}_walker_efficiency.pdf'.format(methods[method]), bbox_inches='tight')
    plt.close()

exit()

wl_dos_one = np.zeros([repeats, bins])
window = 16
print("Window = ", window)
for method in methods:
  for k in range(1,repeats+1):
    filename = "{}_{:02d}_{:02d}_{}/data/wl_dos.dat".format(method, 6, window, k)
    wl_dos = nc.Dataset(filename)
    wl_dos_one[k-1] = np.array(wl_dos["grid data"][:], dtype=np.float64).T
  for k in range(0,repeats):
    wl_dos_one[k] = wl_dos_one[k] + np.mean(wl_dos_one[0] - wl_dos_one[k])
  print(method,np.mean(np.var(wl_dos_one, axis=0)))


#print(wl_dos_one[0]-wl_dos_one[1])
for k in range(0,repeats):
  wl_dos_one[k] = wl_dos_one[k] + np.mean(wl_dos_one[0] - wl_dos_one[k])
  plt.plot(np.linspace(1,512,512), wl_dos_one[k])
plt.savefig('xxxxx.pdf', bbox_inches='tight')



#for method in range(6):
#  for i in range(len(walkers)):
#    for j in range(len(windows)):
#        for k in range(1,repeats+1):
#          filename = "{}_{:02d}_{:02d}_{}/data/wl_dos.dat".format(method, walkers[i], windows[j], k)
#          wl_lb_max_time = nc.Dataset(filename)
#          wl_lb_max_time = np.array(wl_lb_max_time["grid data"][:], dtype=np.float64).T
#
#          row_sums = np.max(wl_lb_max_time, axis=1)
#          sums[i, j, k-1] = np.sum(row_sums)

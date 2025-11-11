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
import matplotlib.ticker as ticker
import matplotlib.ticker as mticker

def plot_domain_efficiency(selected_methods, replica=True):
  if replica:
    replica = "r"
  else:
    replica = "nr"
  
  for overlap in range(len(overlaps)):
    fig, axes = plt.subplots(1, 3, figsize=figsize_subplots)
    max_y = 0
    for method_id in selected_methods:
      for walker in range(len(walkers)):
        if max_y < np.max(efficiency[method_id, walker, :, overlap]):
          max_y = np.max(efficiency[method_id, walker, :, overlap])

    for ax_id, method_id in enumerate(selected_methods):
      ax = axes[ax_id]
      ax.axhline(y=1, linestyle='-')

      for walker in range(len(walkers)):
        ax.errorbar(windows, efficiency[method_id, walker, :, overlap], yerr=efficiency_err[method_id, walker, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(walker+1))
      ax.set_xticks(windows[::2])
      ax.set_xlabel(r"Sub-domains ($h$)")
      ax.tick_params(axis="y", direction='inout')
      if ax_id == 0:
        ax.set_ylabel("Sub-domain Efficiency")
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
      else:
        ax.set_yticklabels([])
      #plt.title("Method {}".format(methods_label[method_id]))
      #plt.legend(loc="upper right", title="Walkers")
      ax.set_ylim(0, max_y*1.25)
      ax.text(0.5, -0.225, titles[ax_id], transform=ax.transAxes, ha='center', va='top', fontweight="bold")
    
    
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('{:02d}_{}_domain_efficiency.pdf'.format(overlaps[overlap], replica), bbox_inches='tight')
    plt.close()

def plot_mc_step(selected_methods, replica=True):
  if replica:
    replica = "r"
  else:
    replica = "nr"
  
  for overlap in range(len(overlaps)):
    fig, axes = plt.subplots(1, 3, figsize=figsize_subplots)
    max_y = 0
    for method_id in selected_methods:
      for walker in range(len(walkers)):
        if max_y < np.max(mc_steps_mean[method_id, walker, :, overlap]):
          max_y = np.max(mc_steps_mean[method_id, walker, :, overlap])

    for ax_id, method_id in enumerate(selected_methods):
      ax = axes[ax_id]
      ax.axhline(y=1, linestyle='-')

      for walker in range(len(walkers)):
        ax.errorbar(windows, mc_steps_mean[method_id, walker, :, overlap], yerr=mc_steps_err[method_id, walker, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(methods_label[method_id]))
      ax.set_xticks(windows[::2])
      ax.set_xlabel("Wang-Landau Instances")
      ax.tick_params(axis="y", direction='inout')
      if ax_id == 0:
        ax.set_ylabel("Relative MC Steps")
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
      else:
        ax.set_yticklabels([])
      ax.set_ylim(0, max_y*1.25)
      ax.text(0.5, -0.225, titles[ax_id], transform=ax.transAxes, ha='center', va='top', fontweight="bold")
    
    
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('{:02d}_{}_mc_steps.pdf'.format(overlaps[overlap], replica), bbox_inches='tight')
    plt.close()

def plot_walker_efficiency(selected_methods, selected_overlaps, replica=True):
  if replica:
    replica = "r"
  else:
    replica = "nr"
  
  for overlap in selected_overlaps:
    fig, axes = plt.subplots(1, 3, figsize=figsize_subplots)
    max_y = 0
    for method_id in selected_methods:
      for window in range(len(windows)):
        if max_y < np.max(efficiency[method_id, :, window, overlap]):
          max_y = np.max(efficiency[method_id, :, window, overlap])

    for ax_id, method_id in enumerate(selected_methods):
      ax = axes[ax_id]
      ax.axhline(y=1, linestyle='-')

      for window in np.array([0, 1, 3, 7, 15]):
        ax.errorbar(walkers, efficiency[method, :, window, overlap], yerr=efficiency_err[method, :, window, overlap], capsize=3, ls='none', fmt='o', label="{}".format(window+1))
      ax.set_xticks(walkers)
      ax.set_xlabel(r"Walkers ($m$)")
      ax.tick_params(axis="y", direction='inout')
      if ax_id == 0:
        ax.set_ylabel("Walker Efficiency")
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
      else:
        ax.set_yticklabels([])
      ax.set_ylim(0, max_y*1.25)
      ax.text(0.5, -0.25, titles[ax_id], transform=ax.transAxes, ha='center', va='top', fontweight="bold")

      #if ax_id == 1:
        #ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.35), ncol=5, title="Sub-domains")
      if ax_id == 2:
        ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, title="Sub-domains")
    
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('{:02d}_{}_walker_efficiency.pdf'.format(overlaps[overlap], replica), bbox_inches='tight')
    plt.close()

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

font_size = 24
figsize = (16,8)
figsize_subplots = (18,6)
np.set_printoptions(suppress=True)
#plt.rcParams.update({"text.usetex": True,
#                     "font.size": font_size})
plt.rcParams.update({"font.size": font_size})
plt.rcParams["figure.figsize"] = figsize
#plt.rc('text', usetex=True)

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
methods_label = np.array(["1", "2", "3", "4", "5", "6"])
windows = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
walkers = np.array([1, 2, 3, 4, 5, 6])
overlaps = np.array([0, 10, 25, 50, 75])
repeats = 5
bins=512

titles = ["(a)", "(b)", "(c)"]
titles = ["Method 1", "Method 3", "Method 5"]

#methods = np.array([0, 2, 4])
#methods_label = np.array(["Non-Uniform + Balance + Replica", "Non-Uniform + Replica", "Replica"])
#methods = np.array([0,1])
walkers = np.array([1])
#overlaps = np.array([25])

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
  for method_id, method in enumerate(methods):
    for walker in range(len(walkers)):
      efficiency[method_id, walker, :, overlap] = one_sums_mean[walker, overlap]/(sums_mean[method_id, walker, :, overlap]*windows)
      efficiency_err[method_id, walker, :, overlap] = efficiency[method_id, walker, :, overlap]*np.sqrt((one_sums_err[walker, overlap]/one_sums_mean[walker, overlap])**2+(sums_err[method_id, walker, :, overlap]/sums_mean[method_id, walker, :, overlap])**2)

plot_domain_efficiency([0, 2, 4], True)

plot_domain_efficiency([1, 3, 5], False)

# Speed Up
for overlap in range(len(overlaps)):
  for method in range(len(methods)):
    for walker in range(len(walkers)):
      efficiency[method, walker, :, overlap] = one_sums_mean[0, overlap]/(sums_mean[method, walker, :, overlap])
      efficiency_err[method, walker, :, overlap] = efficiency[method, walker, :, overlap]*np.sqrt((one_sums_err[0, overlap]/one_sums_mean[0, overlap])**2+(sums_err[method, walker, :, overlap]/sums_mean[method, walker, :, overlap])**2)

    plt.errorbar(windows*1, efficiency[method, 0, :, overlap], yerr=efficiency_err[method, 0, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(methods_label[method]))

  plt.plot(windows, windows, color="#D62728")
  ymin, ymax = plt.ylim()
  plt.plot(windows, windows**2, linestyle="--", color="#FF9896")
  plt.ylim(ymin, ymax)
  plt.xticks(windows, labels=windows)
  plt.xlabel("Wang-Landau Instances")
  plt.ylabel("Speed-up")
  #plt.title("Speed-up for 1 walker per domain")
  #plt.legend(loc="upper left", title="Methods")
  #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(methods), title="Methods")
  plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, title="Methods")
  
  plt.savefig('{:02d}_method_speedup.pdf'.format(overlaps[overlap]), bbox_inches='tight')
  plt.close()

# MC Steps
plot_mc_step([0, 2, 4], True)

plot_mc_step([1, 3, 5], False)

# Speed Up walker comparison
for method in range(len(methods)):
  for overlap in range(len(overlaps)):
    for walker in range(len(walkers)):
      efficiency[method, walker, :, overlap] = one_sums_mean[0, overlap]/(sums_mean[method, walker, :, overlap])
      efficiency_err[method, walker, :, overlap] = efficiency[method, walker, :, overlap]*np.sqrt((one_sums_err[0, overlap]/one_sums_mean[0, overlap])**2+(sums_err[method, walker, :, overlap]/sums_mean[method, walker, :, overlap])**2)

    plt.errorbar(windows*1, efficiency[method, 0, :, overlap], yerr=efficiency_err[method, 0, :, overlap], capsize=3, ls='none', fmt='o', label="{}%".format(overlaps[overlap]))

  plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
  plt.plot(windows, windows, color="#D62728")
  ymin, ymax = plt.ylim()
  plt.plot(windows, windows**2, linestyle="--", color="#FF9896")
  plt.ylim(ymin, ymax)
  plt.xticks(windows, labels=windows)
  plt.xlabel("Wang-Landau Instances")
  plt.ylabel("Speed-up")
  #plt.title("Speed-up for 1 walker per domain")
  #plt.legend(loc="upper left", title="Overlap %")
  #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(overlaps), title="Overlap %")
  plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, title="Overlap %")
  
  plt.savefig('{}_overlap_speedup.pdf'.format(methods[method]), bbox_inches='tight')
  plt.close()

# Load Walker Data
methods = np.array([0, 2, 4])
methods_label = np.array(["1", "3", "5"])
windows = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16])
walkers = np.array([1, 2, 3, 4, 5, 6])
overlaps = np.array([25])

#methods = np.array([0, 1, 2, 3, 4, 5])
#methods_label = np.array(["1", "2", "3", "4", "5", "6"])

#methods = np.array([0, 2, 4])
#walkers = np.array([1])
#overlaps = np.array([0, 10, 25, 50, 75])

sums = np.zeros([len(methods), len(walkers), len(windows), len(overlaps), repeats])
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

sums_mean = np.mean(sums, axis=4)
sums_min = np.min(sums, axis=4)
sums_err = np.std(sums, axis=4)/np.sqrt(repeats)

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

for method in range(len(methods)):
  for overlap in range(len(overlaps)):
    for walker in range(len(walkers)):
      efficiency[method, walker, :, overlap] = one_sums_mean[0, overlap]/(sums_mean[method, walker, :, overlap])
      efficiency_err[method, walker, :, overlap] = efficiency[method, walker, :, overlap]*np.sqrt((one_sums_err[0, overlap]/one_sums_mean[0, overlap])**2+(sums_err[method, walker, :, overlap]/sums_mean[method, walker, :, overlap])**2)

cores = np.linspace(1, int(np.max(windows)*np.max(walkers)), int(np.max(windows)*np.max(walkers)))
for overlap in range(len(overlaps)):
  for method in range(len(methods)):
    max_y = 0
    for walker in range(len(walkers)):
      plt.errorbar(windows*(walker+1), efficiency[method, walker, :, overlap], yerr=efficiency_err[method, walker, :, overlap], capsize=3, ls='none', fmt='o', label="{}".format(walkers[walker]))
      if max_y < np.max(efficiency[method, walker, :, overlap]):
        max_y = np.max(efficiency[method, walker, :, overlap])
    
    plt.plot(cores, cores, color="#D62728")
    #plt.plot(cores, cores**2, linestyle="--", color="#FF9896")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    #plt.xticks(cores, labels=cores)
    plt.ylim(0, max_y*1.25)
    plt.xlabel("Wang-Landau Instances")
    plt.ylabel("Speed-up")
    #plt.title("Speed-up for Method {}".format(methods_label[method]))
    #plt.legend(loc="lower right", title="Walkers")
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(walkers), title="Walkers")
    plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, title="Walkers")
    
    plt.savefig('{:02d}_{}_walker_speedup.pdf'.format(overlaps[overlap], methods[method]), bbox_inches='tight')
    plt.close()

cores = np.linspace(1, int(np.max(windows)*np.max(walkers)), int(np.max(windows)*np.max(walkers)))
for overlap in range(len(overlaps)):
  for method in range(len(methods)):
    max_y = 0
    for walker in range(len(walkers)):
      plt.errorbar(windows*(walker+1), efficiency[method, walker, :, overlap]/(windows*(walker+1)), yerr=efficiency_err[method, walker, :, overlap]/(windows*(walker+1)), capsize=3, ls='none', fmt='o', label="{}".format(walkers[walker]))
      if max_y < np.max(efficiency[method, walker, :, overlap]/(windows*(walker+1))):
        max_y = np.max(efficiency[method, walker, :, overlap]/(windows*(walker+1)))
    
    plt.gca().axhline(y=1, linestyle='-')
    #plt.plot(cores, cores, color="#D62728")
    #plt.plot(cores, cores**2, linestyle="--", color="#FF9896")
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
    #plt.xticks(cores, labels=cores)
    plt.ylim(0, max_y*1.25)
    plt.xlabel("Wang-Landau Instances")
    plt.ylabel("Efficiency")
    #plt.title("Speed-up for Method {}".format(methods_label[method]))
    #plt.legend(loc="lower right", title="Walkers")
    #plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=len(walkers), title="Walkers")
    plt.legend(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, title="Walkers")
    
    plt.savefig('{:02d}_{}_walker_efficiency.pdf'.format(overlaps[overlap], methods[method]), bbox_inches='tight')
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

for overlap in range(len(overlaps)):
  for method in range(len(methods)):
    for window in range(len(windows)):
      efficiency[method, :, window, overlap] = one_sums_mean[window, overlap]/(sums_mean[method, :, window, overlap]*walkers)
      efficiency_err[method, :, window, overlap] = efficiency[method, :, window, overlap]*np.sqrt((one_sums_err[window, overlap]/one_sums_mean[window, overlap])**2+(sums_err[method, :, window, overlap]/sums_mean[method, :, window, overlap])**2)

# Walker Efficiency

plot_walker_efficiency([0, 1, 2], [0], replica=True)
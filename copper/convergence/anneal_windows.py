import argparse
import os
import shutil
import numpy as np

base_path = os.path.dirname(os.path.abspath(__file__))
ry_to_mev = 13605.693122
energy_max = 10          # meV/atom, as in the production runs
shrink = 0.05            # otherwise energy_min is the annealed energy made 5% smaller, truncated towards 0
# Production energy_min (meV/atom); used when it is closer to zero than the annealed energy
production = {"Cu0.00": -181, "Cu0.50": -110}

# Slow anneal from 8010 K down to 210 K (about 2 minutes for 2048 atoms on 4 ranks)
metropolis_inp = """mode = simulated_annealing
n_mc_steps = 1000000
burn_in_start = F
burn_in = T
n_burn_in_steps = 100000
n_sample_steps = 256
calculate_energies = T
n_sample_steps_trajectory = 256
write_trajectory_energy = F
write_trajectory_asro = F
calculate_asro = F
n_sample_steps_asro = 2560
calculate_alro = F
n_sample_steps_alro = 256
write_trajectory_xyz = F
write_final_config_xyz = F
write_final_config_nc = F
read_start_config_nc = F
T = 8010.0
T_steps = 40
delta_T = -200.0
nbr_swap = F
"""

def run_dirs():
    with open(os.path.join(base_path, "run_dirs.txt")) as f:
        return [line.strip() for line in f if line.strip()]

def setup():
    # Make anneal/ in every run directory: same system, but mode=301
    for d in run_dirs():
        path = os.path.join(base_path, d)
        anneal = os.path.join(path, "anneal")
        os.makedirs(anneal, exist_ok=True)
        with open(os.path.join(path, "brawl.inp")) as f:
            brawl = f.read().replace("mode=302", "mode=301")
        with open(os.path.join(anneal, "brawl.inp"), "w") as f:
            f.write(brawl)
        for name in os.listdir(path):
            if name.endswith("_V_ijs_Ry.txt"):
                shutil.copy(os.path.join(path, name), anneal)
        with open(os.path.join(anneal, "metropolis.inp"), "w") as f:
            f.write(metropolis_inp)
    print("anneal directories written")

def apply(write=True):
    # Set energy_min in wl_input.inp from the lowest temperature annealed energy
    # (write=False only prints, for visual inspection)
    for d in run_dirs():
        path = os.path.join(base_path, d)
        energy_file = os.path.join(path, "anneal", "energies", "av_energy_diagnostics.dat")
        if not os.path.exists(energy_file):
            print(f"{d}: no anneal output yet, skipping")
            continue
        energies = np.atleast_2d(np.loadtxt(energy_file, skiprows=1))
        e_low = energies[np.argmin(energies[:, 0]), 1]*ry_to_mev
        prod = production[d.split("/")[1]]
        e_min = prod if prod > e_low else int(e_low*(1.0 - shrink))
        if not write:
            print(f"{d}: annealed E = {e_low:.2f} meV/atom -> energy_min = {e_min}")
            continue
        with open(os.path.join(path, "wl_input.inp")) as f:
            lines = f.read().split("\n")
        with open(os.path.join(path, "wl_input.inp"), "w") as f:
            for line in lines:
                if line.startswith("energy_min"):
                    line = f"energy_min={e_min}"
                elif line.startswith("energy_max"):
                    line = f"energy_max={energy_max}"
                f.write(line + "\n")
        print(f"{d}: annealed E = {e_low:.2f} meV/atom -> energy_min = {e_min}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("step", choices=["setup", "show", "apply"])
    args = parser.parse_args()
    if args.step == "setup":
        setup()
    elif args.step == "show":
        apply(write=False)
    else:
        apply()

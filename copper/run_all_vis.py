import argparse, os, subprocess, sys
from glob import glob
from multiprocessing import Pool

dirs = sorted(d for d in glob("Cu*") if os.path.isdir(d))
script = os.path.abspath("copper_vis.py")

def run(d):
    subprocess.run([sys.executable, script], cwd=d)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cores", type=int, default=min(os.cpu_count(), len(dirs)))
    args = parser.parse_args()
    with Pool(args.cores) as pool:
        pool.map(run, dirs)

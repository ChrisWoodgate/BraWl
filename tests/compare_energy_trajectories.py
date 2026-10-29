#!/usr/bin/env python3
"""Compare two BraWl energy trajectory files step by step.

Usage: compare_energy_trajectories.py file_a file_b [tolerance_Ry]

The files must contain the same steps. Exit status 0 if every energy
agrees within the tolerance (default 1.0e-10 Ry), 1 otherwise.

Authors: C. D. Woodgate, Claude Fable 5.1 (2026)
"""
import sys

import numpy as np


def main():
    a = np.loadtxt(sys.argv[1], ndmin=2)
    b = np.loadtxt(sys.argv[2], ndmin=2)
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0e-10
    if a.shape != b.shape or not np.array_equal(a[:, 0], b[:, 0]):
        print("step columns differ:", a.shape, b.shape)
        sys.exit(1)
    diff = np.abs(a[:, 1] - b[:, 1]).max()
    print(f"{len(a)} steps compared; max |dE| = {diff:.3e} Ry "
          f"(tolerance {tol:.1e})")
    ok = diff <= tol
    print("OK" if ok else "MISMATCH")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()

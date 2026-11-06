import sys
import numpy as np

def species(copper, atoms):
    a = int(np.round(atoms * 0.25 * (1 - copper)))
    b = int(np.round(atoms * 0.75 * (1 - copper)))
    c = int(np.round(atoms * copper))
    print(a, b, c, a+b+c)

if __name__ == "__main__":
    copper = float(sys.argv[1])
    atoms = int(sys.argv[2])
    species(copper, atoms)


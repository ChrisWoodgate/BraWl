import numpy as np
import sys
# Need to
# Read in file
# Delete relevant bits
# Output file

def read_xyz(filename):
    """
    Function to read in an xyz file line by line to get atom names and positions

    C. D. Woodgate, Bristol, 2024
    """

    # Python list (will be mixed-type) where we will store the data
    xyz = []

    # Open the file
    with open(filename) as file:

        # Loop over the lines
        for i, line in enumerate(file):
            # Skip the first two lines as these are metadata in the .xyz format
            if i<2:
                continue

            # First portion of the line is the element name, so grab it and remove whitespace
            element = line[0:3].strip()

            # Positions are next
            positions = line[3:]

            # Might have to fine tune the amount of whitespace used to split
            # depending on your exact xyz format
            positions=positions.split('      ')

            # Grab the x, y, and z components and make them floats
            x = float(positions[0].strip())
            y = float(positions[1].strip())
            z = float(positions[2].strip())

            # Parcel it up into a Python list
            thing = [element, x, y, z]

            # Append it to the list of everything
            xyz.append(thing)
            
    return xyz

def write_xyz(array, filename):
    """
    Function to write an xyz file line by line based on atom names and positions

    C. D. Woodgate, Bristol, 2024
    """

    # How many atoms to write?
    natoms = len(array)

    # Write this to line one of the .xyz
    file = open(filename, 'w')
    file.write(f'{natoms}\n')
    file.write('\n')
    file.close()

    # Loop over the atoms and write line-by-line
    with open(filename, 'a') as file:
        for line in array:
            file.write(line[0] + f'{line[1]:12.8f} {line[2]:12.8f} {line[3]:12.8f}\n')
            
        file.close()
    return

def slice(filename, indices=[1,1,1], size=30):
    """
    Function to slice a named xyz file. 
    'Indices' tells you the normal vector to the slice plane.
    'size' tells you how big the slice is

    C. D. Woodgate, Bristol, 2024
    """

    # Read the xyz file
    positions = read_xyz(filename)

    # New array for the sliced positions
    fixed = []

    # Loop over and find if they are the right 'distance' 
    # away from the origin to include or not
    for atom in positions:
        x = atom[1]
        y = atom[2]
        z = atom[3]

        dist = x*indices[0] + y*indices[1] + z*indices[2]

        if ((dist) < size):
            fixed.append(atom)

    # Write this new array to a new file
    write_xyz(fixed, filename[0:-4]+ '_sliced.xyz') 
    
    return

param_size = int(sys.argv[1])
slice('AlTiVNb/grids/proc_000config_at_T_0050.0.xyz', indices=[1,0,1], size=param_size)
slice('AlTiVNb/grids/proc_000config_at_T_1200.0.xyz', indices=[1,0,1], size=param_size)
slice('AlTiVNb/grids/proc_000config_at_T_3000.0.xyz', indices=[1,0,1], size=param_size)

slice('AlTiCrMo/grids/proc_000config_at_T_0050.0.xyz', indices=[1,0,1], size=param_size)
slice('AlTiCrMo/grids/proc_000config_at_T_1200.0.xyz', indices=[1,0,1], size=param_size)
slice('AlTiCrMo/grids/proc_000config_at_T_3000.0.xyz', indices=[1,0,1], size=param_size)
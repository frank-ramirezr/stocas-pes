import math as m
import os
import glob
import numpy as np

np.random.seed(1234)

def radius():
    """
    Calculates a dynamic confinement sphere radius based on the maximum 
    intramolecular distance found across input XYZ files and the total 
    number of molecules (k), scaling gracefully for large systems.
    """
    current_dir = os.getcwd()
    scratch_dir = os.path.join(current_dir, "scratch")
    xyz_files = glob.glob(os.path.join(scratch_dir, "mol*xyz"))

    n_molecs = len(xyz_files)
    if n_molecs == 0:
        print("Warning: No .xyz files found in scratch directory. Returning default radius.")
        return 5.0

    global_max_diameter = 0.0

    for file_path in xyz_files:
        coords = []
        with open(file_path, 'r') as f:
            lines = f.readlines()
            
            # Safely skip typical XYZ headers (Line 1: atom count, Line 2: comment)
            # and parse lines that actually contain coordinate data
            for line in lines:
                tokens = [t.strip() for t in line.split() if t.strip()]
                if len(tokens) == 4:  # Element, X, Y, Z
                    try:
                        coords.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
                    except ValueError:
                        continue # Skip lines that don't map to numbers safely

        num_atoms = len(coords)
        if num_atoms < 2:
            continue

        # Find the maximum internal distance (diameter) of this specific molecule
        mol_max_dist = 0.0
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                dist = m.sqrt(
                    (coords[i][0] - coords[j][0])**2 +
                    (coords[i][1] - coords[j][1])**2 +
                    (coords[i][2] - coords[j][2])**2
                )
                if dist > mol_max_dist:
                    mol_max_dist = dist

        # Track the largest single molecule diameter in the cluster pool
        if mol_max_dist > global_max_diameter:
            global_max_diameter = mol_max_dist

    # Fallback if molecules are single atoms (like noble gases)
    if global_max_diameter == 0.0:
        global_max_diameter = 2.0  # Approximate van der Waals diameter

    k = float(n_molecs)
    
    base_radius = global_max_diameter * 1.2
    
    scaling_modifier = 1.15 - (0.3 / k)
    
    cluster_radius = base_radius * (k ** (1.0 / 3.0)) * scaling_modifier
    
    return round(cluster_radius, 3)

# Example output simulation
# If global_max_diameter for water is ~1.5 Angstroms:
# For k = 5, cluster_radius will evaluate closer to ~2.6 - 3.5 Angstroms, depending on base geometry.

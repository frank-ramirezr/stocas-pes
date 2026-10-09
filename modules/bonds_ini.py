import math as m
import numpy as np
import os as os
import glob as glob
import re as re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]

ptable = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
    'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
    'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
    'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53
}

covalent_radii = np.array([
    0.00, 0.32, 0.46, 1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.62, # 0-9
    1.55, 1.39, 1.30, 1.21, 1.11, 1.07, 1.05, 1.02, 2.03, 1.76, # 10-19
    1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, # 20-29
    1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 2.20, 1.95, 1.90, 1.75, # 30-39
    1.64, 1.56, 1.47, 1.46, 1.42, 1.39, 1.44, 1.44, 1.42, 1.39, # 40-49
    1.39, 1.39, 1.38, 1.39                                       # 50-53
])

multiple_bond_rules = {
    (6, 6): [(1.24, 2), (1.37, 1)], # C≡C gets +2 (total 3), C=C gets +1 (total 2), else single (1)
    (7, 7): [(1.12, 2), (1.28, 1)], # N≡N, N=N
    (8, 8): [(1.23, 1)],            # O=O
    (6, 7): [(1.18, 2), (1.32, 1)], # C≡N, C=N
    (6, 8): [(1.16, 2), (1.26, 1)], # C≡O, C=O
    (7, 8): [(1.18, 2), (1.26, 1)], # N≡O, N=O
    (5, 7): [(1.24, 2), (1.42, 1)]  # B≡N, B=N
}

def compute_bonds():

    # calling predefined parameters

    cdir = os.getcwd()
    scratch_dir = f'{cdir}/scratch'
    dirs = glob.glob(f'{scratch_dir}/mol*.xyz')
    dirs.sort(key=natural_keys)

    total_bonds = []

    for molecule in dirs:

        with open(molecule, 'r') as f:
            lines = f.readlines()

        if lines.__len__() == 1:
            total_bonds.append(0)

            continue

        xyz_coords = []
        atoms = []

        for line in lines:
            atoms.append(line.strip().split()[0])
            xyz_coords.append(line.strip().split()[1:4])

        znu = np.asarray([ptable[atom] for atom in atoms], dtype='int')
        coords = np.asarray(xyz_coords, dtype='float')
        x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

        n_bonds_s = 0

        for i in range(0, len(znu)):
            for j in range(i + 1, len(znu)):

                distance = m.sqrt((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2)

                z_min, z_max = sorted((znu[i], znu[j]))

                r_expected = covalent_radii[z_min] + covalent_radii[z_max]

                bond_cutoff = 1.1 * r_expected

                if 0.80 <= distance <= bond_cutoff:
                    bonds_to_add = 1

                    if (z_min, z_max) in multiple_bond_rules:
                        for max_d, extra_bonds in multiple_bond_rules[(z_min, z_max)]:
                            if distance <= max_d:
                                bonds_to_add += extra_bonds
                                break

                    n_bonds_s += bonds_to_add

        total_bonds.append(n_bonds_s)

    return total_bonds

# execute function compute bonds

if __name__ == "__main__":
    total_bonds = compute_bonds()


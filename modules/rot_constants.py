import re
import os as os
from numpy import linalg as la
import numpy as np

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]

def inertia():
    ptable = {
        'H': 1.00797,
        'He': 4.0026,
        'Li': 6.941,
        'Be': 9.01218,
        'B': 10.81,
        'C': 12.011,
        'N': 14.0067,
        'O': 15.9994,
        'F': 18.998403,
        'Ne': 20.179,
        'Na': 22.98977,
        'Mg': 24.305,
        'Al': 26.98154,
        'Si': 28.0855,
        'P': 30.97376,
        'S': 32.06,
        'Cl': 35.453,
        'Ar': 39.948,
        'K': 39.0983,
        'Ca': 40.08,
        'Sc': 44.9559,
        'Ti': 47.9,
        'V': 50.9415,
        'Cr': 51.996,
        'Mn': 54.938,
        'Fe': 55.847,
        'Co': 58.9332,
        'Ni': 58.7,
        'Cu': 63.546,
        'Zn': 65.38,
        'Ga': 69.72,
        'Ge': 72.59,
        'As': 74.9216,
        'Se': 78.96,
        'Br': 79.904,
        'Kr': 83.8,
        'Rb': 85.4678,
        'Sr': 87.62,
        'Y': 88.9059,
        'Zr': 91.22,
        'Nb': 92.9064,
        'Mo': 95.94,
        'Tc': 98.0,
        'Ru': 101.07,
        'Rh': 102.9055,
        'Pd': 106.4,
        'Ag': 107.868,
        'Cd': 112.41,
        'In': 114.82,
        'Sn': 118.69,
        'Sb': 121.75,
        'Te': 127.6,
        'I': 126.9045,
        'Xe': 131.3,
        'Cs': 132.9054,
        'Ba': 137.33,
        'La': 138.9055,
        'Ce': 140.12,
        'Pr': 140.9077,
        'Nd': 144.24,
        'Pm': 145.0,
        'Sm': 150.4,
        'Eu': 151.96,
        'Gd': 157.25,
        'Tb': 158.9254,
        'Dy': 162.5,
        'Ho': 164.9304,
        'Er': 167.26,
        'Tm': 168.9342,
        'Yb': 173.04,
        'Lu': 174.967,
        'Hf': 178.49,
        'Ta': 180.9479,
        'W': 183.85,
        'Re': 186.207,
        'Os': 190.2,
        'Ir': 192.22,
        'Pt': 195.09,
        'Au': 196.9665,
        'Hg': 200.59,
        'Tl': 204.37,
        'Pb': 207.2,
        'Bi': 208.9804,
        'Po': 209.0,
        'At': 210.0,
        'Rn': 222.0,
        'Fr': 223.0,
        'Ra': 226.0254,
        'Ac': 227.0278,
        'Th': 232.0381,
        'Pa': 231.0359,
        'U': 238.029,
        'Np': 237.0482,
        'Pu': 242.0,
        'Am': 243.0,
        'Cm': 247.0,
        'Bk': 247.0,
        'Cf': 251.0,
        'Es': 252.0,
        'Fm': 257.0,
        'Md': 258.0,
        'No': 250.0,
        'Lr': 260.0,
        'Rf': 261.0,
        'Db': 262.0,
        'Sg': 263.0,
        'Bh': 262.0,
        'Hs': 255.0,
        'Mt': 256.0,
        'Ds': 269.0,
        'Rg': 272.0,
        'Cn': 277.0,
        'Uub': 277.0,
    }

    cdir = os.getcwd()
    scratch_dir = f'{cdir}/scratch'
    opt_file = f'{cdir}/xtbopt.xyz'

    masses = []
    coords = []

    with open(opt_file, 'r') as f:
        num_atoms = int(f.readline().strip())
        f.readline()  # skip comment

        for _ in range(num_atoms):
            parts = f.readline().split()
            symbol = parts[0]
            masses.append(ptable.get(symbol.capitalize(), 12.011))  # default fallback
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    masses = np.array(masses)
    coords = np.array(coords)

    # 1. Shift origin to Center of Mass (COM)
    total_mass = np.sum(masses)
    com = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
    coords_com = coords - com  # Centered coordinates

    # 2. Build Inertia Tensor
    x = coords_com[:, 0]
    y = coords_com[:, 1]
    z = coords_com[:, 2]

    ixx = np.sum(masses * (y ** 2 + z ** 2))
    iyy = np.sum(masses * (x ** 2 + z ** 2))
    izz = np.sum(masses * (x ** 2 + y ** 2))

    ixy = -np.sum(masses * x * y)
    ixz = -np.sum(masses * x * z)
    iyz = -np.sum(masses * y * z)

    tensor = np.array([[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]])

    # 3. Calculate Principal Moments of Inertia
    ia, ib, ic = la.eigvalsh(tensor)  # eigvalsh is optimized for symmetric matrices

    # Sort eigenvalues ascending (I_A <= I_B <= I_C)
    ia, ib, ic = sorted([ia, ib, ic])

    return ia, ib, ic
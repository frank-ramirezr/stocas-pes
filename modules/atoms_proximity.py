import math
import os as os
from modules import minmax_rad

def proximity():
    atom = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
            'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
            'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
            'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
            'U', 'Np', 'Pu', 'Am']

    atom_radii = [0.53, 0.31, 1.67, 1.12, 0.87, 0.67, 0.56, 0.48, 0.42, 0.38, 1.9, 1.45, 1.18, 1.11, 0.98, 0.88, 0.79,
                  0.71, 2.43, 1.94, 1.84, 1.76, 1.71, 1.66, 1.61, 1.56, 1.52, 1.49, 1.45, 1.42, 1.36, 1.25, 1.14, 1.03,
                  0.94, 0.88, 2.65, 2.19, 2.12, 2.06, 1.98, 1.9, 1.83, 1.78, 1.73, 1.69, 1.65, 1.61, 1.56, 1.45, 1.33,
                  1.23, 1.15, 1.08, 2.98, 2.53, 1.95, 1.85, 2.47, 2.06, 2.05, 2.38, 2.31, 2.33, 2.25, 2.28, 2.26, 2.26,
                  2.22, 2.22, 2.17, 2.08, 2.0, 1.93, 1.88, 1.85, 1.8, 1.77, 1.74, 1.71, 1.56, 1.54, 1.43, 1.35, 1.27,
                  1.2, 2.7, 2.15, 1.95, 1.8, 1.8, 1.75, 1.75, 1.75, 1.75]

    pt_len = len(atom_radii)
    cdir = os.getcwd()
    scratch_dir = f'{cdir}/scratch'
    xyz_file = f'{scratch_dir}/xyz'

    with open(xyz_file, 'r') as data_file:
        nlines = 0
        cluster = []
        for line in data_file:
            nlines = nlines + 1
            line = [i.strip() for i in line.split(' ') if i]
            cluster.append(line)

    for i in range(0, nlines):
        for j in range(0, pt_len):
            if cluster[i][0].lower() == atom[j].lower():
                cluster[i][0] = cluster[i][0].replace(cluster[i][0], str(atom_radii[j]))

    mol_atoms = f'{scratch_dir}/comp.ini'

    with open(mol_atoms, 'r') as cluster_comp:
        file = []
        for line in cluster_comp:
            file.append(line.split(';'))

    comp_lines = sum(1 for item in file)

    natoms = []
    nmoles = []

    for i in range(0, comp_lines):
        natoms.append(file[i][0])
        nmoles.append(file[i][1])

    count1 = 0
    min_distance = 1000

    for m in range(0, comp_lines):
        nmoles[m] = int(nmoles[m])
        natoms[m] = int(natoms[m])
        if m == 0:
            count1 = 0
        if m != 0:
            count1 = count1 + natoms[m - 1] * nmoles[m - 1]
        for k in range(1, nmoles[m] + 1):
            for i in range((k - 1) * natoms[m] + count1, k * natoms[m] + count1):
                for j in range(k * natoms[m] + count1, nlines):
                    rad1 = float(cluster[i][0])
                    rad2 = float(cluster[j][0])
                    x1 = float(cluster[i][1])
                    x2 = float(cluster[j][1])
                    y1 = float(cluster[i][2])
                    y2 = float(cluster[j][2])
                    z1 = float(cluster[i][3])
                    z2 = float(cluster[j][3])
                    distance = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
                    atom_distance = rad1 + rad2
                    if distance < atom_distance:
                        if distance < min_distance:
                            min_distance = distance

    return min_distance




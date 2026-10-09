import math
import os as os
import random as rn
import multiprocessing as mp


def ellipsoid():
    # ##################################################################################################################
    atom = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
            'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb',
            'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs',
            'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta',
            'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',
            'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
            'Ds', 'Rg', 'Uub']

    atom_mass = [1.00797, 4.0026, 6.941, 9.01218, 10.81, 12.011, 14.0067, 15.9994, 18.998403, 20.179, 22.98977, 24.305,
                 26.98154, 28.0855, 30.97376, 32.06, 35.453, 39.948, 39.0983, 40.08, 44.9559, 47.9, 50.9415, 51.996,
                 54.938, 55.847, 58.9332, 58.7, 63.546, 65.38, 69.72, 72.59, 74.9216, 78.96, 79.904, 83.8, 85.4678,
                 87.62,
                 88.9059, 91.22, 92.9064, 95.94, 98, 101.07, 102.9055, 106.4, 107.868, 112.41, 114.82, 118.69, 121.75,
                 127.6, 126.9045, 131.3, 132.9054, 137.33, 138.9055, 140.12, 140.9077, 144.24, 145, 150.4, 151.96,
                 157.25,
                 158.9254, 162.5, 164.9304, 167.26, 168.9342, 173.04, 174.967, 178.49, 180.9479, 183.85, 186.207, 190.2,
                 192.22, 195.09, 196.9665, 200.59, 204.37, 207.2, 208.9804, 209, 210, 222, 223, 226.0254, 227.0278,
                 232.0381, 231.0359, 238.029, 237.0482, 242, 243, 247, 247, 251, 252, 257, 258, 250, 260, 261, 262, 263,
                 262, 255, 256, 269, 272, 277]

    # ######################################################################################################################

    cdir = os.getcwd()
    scratch_dir = f'{cdir}/scratch'

    xyz_check = f'{scratch_dir}/xyz'
    atoms_total = sum(1 for line in open(xyz_check))

    file = open(xyz_check, 'r')
    lines = file.readlines()

    xyz_file = []
    for i in range(0, int(atoms_total)):
        xyz_file.append(lines[i].split())

    p_table_len = len(atom)

    for i in range(0, p_table_len):
        for j in range(0, int(atoms_total)):
            if xyz_file[j][0].lower() == atom[i].lower():
                xyz_file[j][0] = xyz_file[j][0].lower()
                xyz_file[j][0] = xyz_file[j][0].replace(atom[i].lower(), str(atom_mass[i]), 1)

    mol_atoms = f'{scratch_dir}/comp.ini'

    try:
        os.remove(f'{scratch_dir}/ellipsoid_file.dat')
    except FileNotFoundError:
        pass

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

    for m in range(0, comp_lines):
        nmoles[m] = int(nmoles[m])
        natoms[m] = int(natoms[m])
        if m == 0:
            count1 = 0
        if m != 0:
            count1 = count1 + natoms[m - 1] * nmoles[m - 1]
        for k in range(1, nmoles[m] + 1):

            minx = 1000
            maxx = -1000
            miny = 1000
            maxy = -1000
            minz = 1000
            maxz = -1000

            for i in range((k - 1) * natoms[m] + count1, k * natoms[m] + count1):

                x = float(xyz_file[i][1])
                y = float(xyz_file[i][2])
                z = float(xyz_file[i][3])

                if x <= minx:
                    minx = x
                if x >= maxx:
                    maxx = x
                if y <= miny:
                    miny = y
                if y >= maxy:
                    maxy = y
                if z >= maxz:
                    maxz = z
                if z <= minz:
                    minz = z

            a = maxx - minx
            b = maxy - miny
            c = maxz - minz
            x0 = (maxx + minx) / 2.0
            y0 = (maxy + miny) / 2.0
            z0 = (maxz + minz) / 2.0

            x0 = '{0:.8f}'.format(x0)
            y0 = '{0:.8f}'.format(y0)
            z0 = '{0:.8f}'.format(z0)

            a = '{0:.8f}'.format(a)
            b = '{0:.8f}'.format(b)
            c = '{0:.8f}'.format(c)

            ellipsoid_file = f'{a}; {b}; {c}; {x0}; {y0}; {z0}\n'

            with open(f'{scratch_dir}/ellipsoid_file.dat', 'a') as ellipsoid_data:
                ellipsoid_data.write(ellipsoid_file)

    n_lines = 0
    ellipsoid_file = []

    with open(f'{scratch_dir}/ellipsoid_file.dat', 'r') as ellipsoid_data:
        for line in ellipsoid_data:
            n_lines = n_lines + 1
            ellipsoid_file.append(line.split(';'))

    iterations = 500

    min_distance = 1000

    for i in range(0, n_lines - 1):
        for j in range(i + 1, n_lines):

            a1 = ellipsoid_file[i][0]
            b1 = ellipsoid_file[i][1]
            c1 = ellipsoid_file[i][2]
            x01 = ellipsoid_file[i][3]
            y01 = ellipsoid_file[i][4]
            z01 = ellipsoid_file[i][5]

            a2 = ellipsoid_file[j][0]
            b2 = ellipsoid_file[j][1]
            c2 = ellipsoid_file[j][2]
            x02 = ellipsoid_file[j][3]
            y02 = ellipsoid_file[j][4]
            z02 = ellipsoid_file[j][5]

            a1 = float(a1)
            b1 = float(b1)
            c1 = float(c1)
            x01 = float(x01)
            y01 = float(y01)
            z01 = float(z01)

            a2 = float(a2)
            b2 = float(b2)
            c2 = float(c2)
            x02 = float(x02)
            y02 = float(y02)
            z02 = float(z02)

            for k in range(1, iterations + 1):
                theta1 = 2.0 * math.pi * rn.uniform(0, 1)
                chi1 = math.pi * rn.uniform(0, 1)
                x1 = x01 - a1 * math.cos(theta1) * math.sin(chi1)
                y1 = y01 - b1 * math.sin(theta1) * math.sin(chi1)
                z1 = z01 - c1 * math.cos(chi1)

                theta2 = 2.0 * math.pi * rn.uniform(0, 1)
                chi2 = math.pi * rn.uniform(0, 1)
                x2 = x02 - a2 * math.cos(theta2) * math.sin(chi2)
                y2 = y02 - b2 * math.sin(theta2) * math.sin(chi2)
                z2 = z02 - c2 * math.cos(chi2)

                distance = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
                if distance < min_distance:
                    min_distance = distance
    return min_distance
#
# data_points = 100
# max_dist = 0
# min_dist = 1000
#
# average_dist = 0
#
# for k in range(0, data_points):
#     distance = ellipsoid()
#     average_dist = average_dist + distance
#     if distance <= min_dist:
#         min_dist = distance
#     if distance >= max_dist:
#         max_dist = distance
#
# average_dist = average_dist / data_points
#
# print(f'average = {average_dist}; min_distance = {min_dist}; max_distance = {max_dist}')
#
#
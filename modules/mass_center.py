import re
import math
import os as os
import glob
import random as ran


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):

    return [atoi(c) for c in re.split(r'(\d+)', text)]


def mass_center():
    max_atm_z = 109

    atom = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K',
            'Ca',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr',
            'Y',
            'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La',
            'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re',
            'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np',
            'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg',
            'Uub']

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

    precision = 8
    current_dir = os.getcwd()
    scratch_dir = f'{current_dir}/scratch'

    dirs = glob.glob(f'{scratch_dir}/mol*xyz')
    dirs.sort(key=natural_keys)

    n_molecs = len(dirs)
    max_atomic_number = 112
    k = 0

    try:
        mass_center_loc = f'{scratch_dir}/mass_center.dat'
        os.remove(mass_center_loc)

    except:
        pass

    while k < n_molecs:

        with open(dirs[k]) as molecule:

            new_coords = []

            for line in molecule:
                line = [i.strip() for i in line.split(' ') if i]
                new_coords.append(line)

        num_lines = 0

        with open(file=dirs[k]) as file:
            for _ in file:
                num_lines = num_lines + 1

        mol_mass = 0.0
        mass_x = 0.0
        mass_y = 0.0
        mass_z = 0.0

        for i in range(0, num_lines):

            j = 0

            while j <= max_atomic_number - 1:

                if new_coords[i][0].lower() == atom[j].lower():

                    mass_x = mass_x + float(new_coords[i][1]) * float(atom_mass[j])
                    mass_y = mass_y + float(new_coords[i][2]) * float(atom_mass[j])
                    mass_z = mass_z + float(new_coords[i][3]) * float(atom_mass[j])
                    mol_mass = mol_mass + float(atom_mass[j])

                    # print(new_coords[i][0])

                j = j + 1

        center_x = mass_x / mol_mass
        center_y = mass_y / mol_mass
        center_z = mass_z / mol_mass

        center_x = "{0:.12f}".format(center_x, precision)
        center_y = "{0:.12f}".format(center_y, precision)
        center_z = "{0:.12f}".format(center_z, precision)

        mass_center = open(f'{scratch_dir}/mass_center.dat', 'a')
        mass_center_write = f'{center_x} {center_y} {center_z} \n'
        mass_center.write(mass_center_write)


        k = k + 1


# mass_center()


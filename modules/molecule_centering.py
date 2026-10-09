import re as re
import os as os
import glob as glob
from modules import mass_center

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]

def centering():

    cdir = os.getcwd()
    scratch_dir = f'{cdir}/scratch'

    mass_center.mass_center()

    dirs = glob.glob(f'{scratch_dir}/mol*xyz')
    dirs.sort(key=natural_keys)

    n_molecs = len(dirs)
    k = 0

    with open(f'{scratch_dir}/mass_center.dat', 'r') as center_mass:
        centers = []
        for line_center in center_mass:
            line_center = [i.strip() for i in line_center.split(' ') if i]
            centers.append(line_center)

    while k < n_molecs:
        with open(dirs[k]) as molecule:
            new_coords = []
            for line in molecule:
                line = [i.strip() for i in line.split(' ') if i]
                new_coords.append(line)

        num_lines = sum(1 for line in open(dirs[k]))

        mol_file = f'{scratch_dir}/mol{k + 1}.xyz'
        with open(mol_file, 'w') as file:

            for i in range(0, num_lines):
                atom = new_coords[i][0]
                x_trasl = float(new_coords[i][1]) - float(centers[k][0])
                y_trasl = float(new_coords[i][2]) - float(centers[k][1])
                z_trasl = float(new_coords[i][3]) - float(centers[k][2])

                x_trasl = '{0:.8f}'.format(x_trasl)
                y_trasl = '{0:.8f}'.format(y_trasl)
                z_trasl = '{0:.8f}'.format(z_trasl)

                fixed_cartesians = new_coords[i][4]

                file.write(f'{atom} {x_trasl} {y_trasl} {z_trasl} {fixed_cartesians}\n')
        k = k + 1

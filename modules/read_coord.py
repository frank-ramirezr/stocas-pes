import os as os
import shutil
from modules import molecule_centering

def mol_gen():

    current_dir = os.getcwd()

    try:
        shutil.rmtree(f'{current_dir}/scratch')
    except FileNotFoundError:
        pass

    os.makedirs('scratch')
    scratch_dir = f'{current_dir}/scratch'
    cluster = f'{scratch_dir}/comp.ini'

    with open("comp.xyz") as f:
        new_coords = []
        for line in f:
            line = [i.strip() for i in line.split(' ') if i]
            new_coords.append(line)

    num_lines = sum(1 for line in open('comp.xyz'))

    mol_count = 0
    precision = 8

    for i in range(0, num_lines):

        delim = new_coords[i][0]

        if delim == 'mol':

            natoms = int(new_coords[i][1])
            nmolecs = int(new_coords[i][2])
            nmolecs = int(nmolecs)
            fixed = str(new_coords[i][3])

            fixed_cartesian = False
            if fixed.lower() == 'true':
                fixed_cartesian = True

            with open(cluster, 'a') as atoms_molecules:
                atoms_molecules.write(f'{natoms}; {nmolecs}; {fixed_cartesian} \n')

            for j_molec in range(0, nmolecs):

                mol_count = mol_count + 1
                mol_index = f'mol{mol_count}.xyz'

                try:
                    mol_index = f'{scratch_dir}/{mol_index}'
                    os.remove(mol_index)

                except FileNotFoundError:
                    pass

                for k_atom in range(i + 1, i + natoms + 1):

                    atom = str(new_coords[k_atom][0])

                    x = float(new_coords[k_atom][1])
                    y = float(new_coords[k_atom][2])
                    z = float(new_coords[k_atom][3])

                    x = "{:.{}f}".format(x, precision)
                    y = "{:.{}f}".format(y, precision)
                    z = "{:.{}f}".format(z, precision)

                    coord_set = f'{atom} {x} {y} {z} {fixed_cartesian}\n'
                    orca_xyz = open(mol_index, 'a')

                    orca_xyz.write(coord_set)

        if delim == '*':
            pass

mol_gen()
molecule_centering.centering()





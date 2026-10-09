import re
import math as m
import os as os
import glob
import sys
import random as rn
import numpy as np
from modules.radius import radius
from modules.trasl_matrix import trasl_matrix
from modules.rot_matrix import rot_matrix

ptable = {
    'H': 1.00797, 'He': 4.0026, 'Li': 6.941, 'Be': 9.01218, 'B': 10.81, 'C': 12.011,
    'N': 14.0067, 'O': 15.9994, 'F': 18.998403, 'Ne': 20.179, 'Na': 22.98977, 'Mg': 24.305,
    'Al': 26.98154, 'Si': 28.0855, 'P': 30.97376, 'S': 32.06, 'Cl': 35.453, 'Ar': 39.948,
    'K': 39.0983, 'Ca': 40.08, 'Sc': 44.9559, 'Ti': 47.9, 'V': 50.9415, 'Cr': 51.996,
    'Mn': 54.938, 'Fe': 55.847, 'Co': 58.9332, 'Ni': 58.7, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.72, 'Ge': 72.59, 'As': 74.9216, 'Se': 78.96, 'Br': 79.904, 'Kr': 83.8,
    'Rb': 85.4678, 'Sr': 87.62, 'Y': 88.9059, 'Zr': 91.22, 'Nb': 92.9064, 'Mo': 95.94,
    'Tc': 98, 'Ru': 101.07, 'Rh': 102.9055, 'Pd': 106.4, 'Ag': 107.868, 'Cd': 112.41,
    'In': 114.82, 'Sn': 118.69, 'Sb': 121.75, 'Te': 127.6, 'I': 126.9045, 'Xe': 131.3,
    'Cs': 132.9054, 'Ba': 137.33, 'La': 138.9055, 'Ce': 140.12, 'Pr': 140.9077, 'Nd': 144.24,
    'Pm': 145, 'Sm': 150.4, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.9254, 'Dy': 162.5,
    'Ho': 164.9304, 'Er': 167.26, 'Tm': 168.9342, 'Yb': 173.04, 'Lu': 174.967, 'Hf': 178.49,
    'Ta': 180.9479, 'W': 183.85, 'Re': 186.207, 'Os': 190.2, 'Ir': 192.22, 'Pt': 195.09,
    'Au': 196.9665, 'Hg': 200.59, 'Tl': 204.37, 'Pb': 207.2, 'Bi': 208.9804, 'Po': 209,
    'At': 210, 'Rn': 222, 'Fr': 223, 'Ra': 226.0254, 'Ac': 227.0278, 'Th': 232.0381,
    'Pa': 231.0359, 'U': 238.029, 'Np': 237.0482, 'Pu': 242, 'Am': 243, 'Cm': 247,
    'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 250,
    'Lr': 260, 'Rf': 261, 'Db': 262, 'Sg': 263, 'Bh': 262, 'Hs': 255,
    'Mt': 256, 'Ds': 269, 'Rg': 272, 'Uub': 277
}

atomic_radii = {
    'H': 0.53, 'He': 0.31, 'Li': 1.67, 'Be': 1.12, 'B': 0.87, 'C': 0.67,
    'N': 0.56, 'O': 0.48, 'F': 0.42, 'Ne': 0.38, 'Na': 1.9, 'Mg': 1.45,
    'Al': 1.18, 'Si': 1.11, 'P': 0.98, 'S': 0.88, 'Cl': 0.79, 'Ar': 0.71,
    'K': 2.43, 'Ca': 1.94, 'Sc': 1.84, 'Ti': 1.76, 'V': 1.71, 'Cr': 1.66,
    'Mn': 1.61, 'Fe': 1.56, 'Co': 1.52, 'Ni': 1.49, 'Cu': 1.45, 'Zn': 1.42,
    'Ga': 1.36, 'Ge': 1.25, 'As': 1.14, 'Se': 1.03, 'Br': 0.94, 'Kr': 0.88,
    'Rb': 2.65, 'Sr': 2.19, 'Y': 2.12, 'Zr': 2.06, 'Nb': 1.98, 'Mo': 1.9,
    'Tc': 1.83, 'Ru': 1.78, 'Rh': 1.73, 'Pd': 1.69, 'Ag': 1.65, 'Cd': 1.61,
    'In': 1.56, 'Sn': 1.45, 'Sb': 1.33, 'Te': 1.23, 'I': 1.15, 'Xe': 1.08,
    'Cs': 2.98, 'Ba': 2.53, 'La': 1.95, 'Ce': 1.85, 'Pr': 2.47, 'Nd': 2.06,
    'Pm': 2.05, 'Sm': 2.38, 'Eu': 2.31, 'Gd': 2.33, 'Tb': 2.25, 'Dy': 2.28,
    'Ho': 2.26, 'Er': 2.26, 'Tm': 2.22, 'Yb': 2.22, 'Lu': 2.17, 'Hf': 2.08,
    'Ta': 2.0, 'W': 1.93, 'Re': 1.88, 'Os': 1.85, 'Ir': 1.8, 'Pt': 1.77,
    'Au': 1.74, 'Hg': 1.71, 'Tl': 1.56, 'Pb': 1.54, 'Bi': 1.43, 'Po': 1.35,
    'At': 1.27, 'Rn': 1.2, 'Fr': 2.7, 'Ra': 2.15, 'Ac': 1.95, 'Th': 1.8,
    'Pa': 1.8, 'U': 1.75, 'Np': 1.75, 'Pu': 1.75, 'Am': 1.75
}


def mass_center(coords, ptable=ptable):
    coords = np.asarray(coords)
    atoms = coords[:, 0]
    natoms, ncols = coords.shape
    coords = np.asarray(coords[:, 1:4], dtype='float')

    get_mass = np.vectorize(ptable.get, otypes=[float])
    masses = get_mass(atoms)
    com = np.average(coords, axis=0, weights=masses)
    return com


def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]


def coord_gen(global_iter,
              calc_rad,
              input_rad,
              max_rot,
              max_trasl,
              actual_temp,
              t0):

    current_dir = os.getcwd()
    scratch_dir = f'{current_dir}/scratch'
    dirs = glob.glob(f'{scratch_dir}/mol*.xyz')
    dirs.sort(key=natural_keys)

    try:
        os.remove(f'{scratch_dir}/xyz')
    except FileNotFoundError:
        pass

    n_molecs = len(dirs)
    total_atoms = 0

    if calc_rad:
        cluster_radius = radius()
    else:
        cluster_radius = input_rad

    for k in range(0, n_molecs):
        with open(dirs[k]) as molecule:
            init_xyz = []
            for line in molecule:
                line = [i.strip() for i in line.split(' ') if i]
                if line:
                    init_xyz.append(line)
            total_atoms = total_atoms + len(init_xyz)

    precision = 12
    max_iter = 10000  # Increased iterations since constraints are tighter now

    # Track coordinates of all accepted atoms across previous molecules
    # Format: (atom_symbol, x, y, z)
    accepted_atoms = []

    # Outer-loop: Molecules
    for k in range(0, n_molecs):
        success = False

        current_scale = actual_temp / t0

        max_trasl_step = max(0.05, max_trasl * current_scale)
        max_rot_step = max(0.05, max_rot * current_scale)

        # Inner-loop: Geometry variations
        for iteration in range(0, max_iter):

            if global_iter == 1:

                found_initial_displacement = False

                for _ in range(0, 1000):

                    phi = rn.uniform(a=0, b=2.0 * m.pi)
                    costheta = rn.uniform(a=-1.0, b=1.0)
                    theta = m.acos(costheta)

                    dx = cluster_radius * m.sin(theta) * m.cos(phi)
                    dy = cluster_radius * m.sin(theta) * m.sin(phi)
                    dz = cluster_radius * m.cos(theta)

                    displacement = m.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

                    if displacement < cluster_radius:
                        found_initial_displacement = True
                        break

                if not found_initial_displacement:
                    dx, dy, dz = 0.0, 0.0, 0.0

            elif global_iter > 1:
                dx, dy, dz = trasl_matrix(max_trasl=max_trasl_step)

            if global_iter == 1:
                r11, r12, r13, r21, r22, r23, r31, r32, r33 = rot_matrix(max_rot=max_rot)

            else:
                r11, r12, r13, r21, r22, r23, r31, r32, r33 = rot_matrix(max_rot=max_rot_step)

            with open(dirs[k]) as molecule:
                init_xyz = []
                for line in molecule:
                    line = [i.strip() for i in line.split(' ') if i]
                    if line:
                        init_xyz.append(line)

            num_lines = len(init_xyz)
            cmx, cmy, cmz = mass_center(init_xyz)

            new_x_coords = []
            new_y_coords = []
            new_z_coords = []
            atom_list = []
            fixed_cartesian_list = []

            rad_max_dist = -1.0

            # Transform each atom in molecule k
            for i in range(0, num_lines):

                x0, y0, z0 = np.asarray(init_xyz[i][1:4], dtype='float')
                atom = init_xyz[i][0]
                fixed_cartesian = str(init_xyz[i][4])

                x = x0 - cmx
                y = y0 - cmy
                z = z0 - cmz

                if fixed_cartesian.lower() == 'false':
                    x0, y0, z0 = x, y, z
                    x = x0 * r11 + y0 * r12 + z0 * r13
                    y = x0 * r21 + y0 * r22 + z0 * r23
                    z = x0 * r31 + y0 * r32 + z0 * r33

                    x = x + cmx + dx
                    y = y + cmy + dy
                    z = z + cmz + dz

                atom_list.append(atom)
                new_x_coords.append(x)
                new_y_coords.append(y)
                new_z_coords.append(z)

                fixed_cartesian_list.append(fixed_cartesian)

                distance = m.sqrt(x ** 2 + y ** 2 + z ** 2)
                if distance > rad_max_dist:
                    rad_max_dist = distance

            # 1. Boundary check
            out_of_sphere = rad_max_dist > cluster_radius

            if out_of_sphere:
                continue  # Out of bounds; discard step immediately

            # 2. Steric Proximity Check (Atom Space Overlap check)
            overlap_detected = False

            for i in range(num_lines):

                at_k = atom_list[i]
                rad_k = atomic_radii[at_k]
                xk, yk, zk = new_x_coords[i], new_y_coords[i], new_z_coords[i]

                # Check distances against all successfully accumulated atoms from molecules 0 to k-1
                for prev_atom, xp, yp, zp in accepted_atoms:
                    rad_p = atomic_radii[prev_atom]

                    # Euclidean metric calculation
                    spatial_dist = m.sqrt((xk - xp) ** 2 + (yk - yp) ** 2 + (zk - zp) ** 2)

                    # Rejection condition: shorter than combined covalent/VDW hard spheres
                    overlap_threshold = 0.7 * (rad_k + rad_p)

                    if spatial_dist < overlap_threshold:
                        overlap_detected = True
                        break
                if overlap_detected:
                    break

            if overlap_detected:
                continue  # Steric crash; discard step and try a new random coordinate set

            # Append molecule k's validated coordinates to our lookup pool
            for i in range(num_lines):
                accepted_atoms.append((atom_list[i], new_x_coords[i], new_y_coords[i], new_z_coords[i]))

            # Consolidated block writing for optimization
            with open(dirs[k], 'w') as molecule_file:
                for i in range(num_lines):
                    x_str = "{:.{}f}".format(new_x_coords[i], precision)
                    y_str = "{:.{}f}".format(new_y_coords[i], precision)
                    z_str = "{:.{}f}".format(new_z_coords[i], precision)
                    molecule_file.write(f'{atom_list[i]} {x_str} {y_str} {z_str} {fixed_cartesian_list[i]}\n')

            # Append coordinates into structural tracking logs
            mode = 'w' if k == 0 else 'a'
            with open(f'{scratch_dir}/xyz', mode) as xyz_file:
                for i in range(num_lines):
                    x_str = "{:.{}f}".format(new_x_coords[i], precision)
                    y_str = "{:.{}f}".format(new_y_coords[i], precision)
                    z_str = "{:.{}f}".format(new_z_coords[i], precision)
                    xyz_file.write(f'{atom_list[i]} {x_str} {y_str} {z_str}\n')

            success = True
            break  # Break inner loop, move to molecule k+1

    # Generate complete frame files
    with open(f'{current_dir}/anneal.xyz', 'w') as anneal_file:
        anneal_file.write(f'{total_atoms}\n\n')
        for atom, x, y, z in accepted_atoms:
            x_str = "{:.{}f}".format(x, precision)
            y_str = "{:.{}f}".format(y, precision)
            z_str = "{:.{}f}".format(z, precision)
            anneal_file.write(f'{atom} {x_str} {y_str} {z_str}\n')

    with open(f'{current_dir}/anneal_trj.xyz', 'a') as trj_file:
        trj_file.write(f'{total_atoms}\n\n')
        for atom, x, y, z in accepted_atoms:
            x_str = "{:.{}f}".format(x, precision)
            y_str = "{:.{}f}".format(y, precision)
            z_str = "{:.{}f}".format(z, precision)
            trj_file.write(f'{atom} {x_str} {y_str} {z_str}\n')
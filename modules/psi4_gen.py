import math as m
import os as os
import glob
from modules import radius
from modules import trasl_matrix
from modules import rot_matrix

def psi4_input_gen(charge,
                   mult,
                   psi_method,
                   psi_ref, basis,
                   psi_basis_jkfit,
                   psi_basis_ri,
                   calc_rad,
                   input_rad):

    current_dir = os.getcwd()
    scratch_dir = f'{current_dir}/scratch'
    
    try:
        os.remove(f'{current_dir}/anneal.inp')
    except FileNotFoundError:
        pass

    with open('anneal.inp', 'a') as psi4_input:
        psi4_input.write('molecule { \n')
        psi4_input.write(f'{charge} {mult} \n')

    dirs = glob.glob(f'{scratch_dir}/mol*xyz')

    try:
        os.remove(f'{scratch_dir}/xyz')
    except FileNotFoundError:
        pass

    k = 0
    n_molecs = len(dirs)
    total_atoms = 0
    
    if calc_rad:

        cluster_radius = radius.radius()
    else:

        cluster_radius = input_rad

    while (k < n_molecs):
    
        with open(dirs[k]) as molecule:
        
            new_coords = []
            
            for line in molecule:
            
                line = [i.strip() for i in line.split(' ') if i]
                new_coords.append(line)
                
            total_atoms = total_atoms + len(new_coords)
            
        k = k + 1

    precision = 12
    max_iter = 1000000
    iteration = 0

    while (iteration < max_iter):

        try:
            os.remove(f'{scratch_dir}/xyz')
            os.remove(f'{current_dir}/anneal.xyz')
        except FileNotFoundError:
            pass

        with open(f'{current_dir}/anneal.xyz', 'a') as anneal_xyz:
        
            anneal_xyz.write(f'{str(total_atoms)} \n')
            anneal_xyz.write(' \n')

        k = 0
        while k < n_molecs:

            with open(dirs[k]) as molecule:

                new_coords = []

                for line in molecule:
                    line = [i.strip() for i in line.split(' ') if i]
                    new_coords.append(line)

            num_lines = sum(1 for line in open(dirs[k]))

            tx, ty, tz = trasl_matrix.trasl_matrix(calc_rad, input_rad)

            tx = float(tx)
            ty = float(ty)
            tz = float(tz)

            r11, r12, r13, r21, r22, r23, r31, r32, r33 = rot_matrix.rot_matrix()
            for i in range(0, num_lines):

                atom = new_coords[i][0]
                
                x0 = float(new_coords[i][1])
                y0 = float(new_coords[i][2])
                z0 = float(new_coords[i][3])
                
                fixed_cartesian = str(new_coords[i][4])

                if fixed_cartesian.lower() == 'true':

                    x = x0
                    y = y0
                    z = z0

                elif fixed_cartesian.lower() == 'false':

                    x = x0 * r11 + y0 * r12 + z0 * r13 + tx
                    y = x0 * r21 + y0 * r22 + z0 * r23 + ty
                    z = x0 * r31 + y0 * r32 + z0 * r33 + tz

                x = "{:.{}f}".format(x, precision)
                y = "{:.{}f}".format(y, precision)
                z = "{:.{}f}".format(z, precision)

                xyz_coords = f'{atom}   {x}     {y}     {z}'

                with open(f'{scratch_dir}/xyz', 'a') as xyz_file:
                    xyz_file.write(f'{xyz_coords} \n')

            k = k + 1

        i = 0
        j = 0

        with open(f'{scratch_dir}/xyz', 'r') as xyz:
        
            xyz_check = []
            
            for line in xyz:
            
                line = [i.strip() for i in line.split(' ') if i]
                xyz_check.append(line)

        atoms_closed = False

        for i in range(0, total_atoms):
        
            for j in range(i + 1, total_atoms):
            
                atom1 = xyz_check[i][0]
                atom2 = xyz_check[j][0]
                
                x1 = xyz_check[i][1]
                x2 = xyz_check[j][1]
                y1 = xyz_check[i][2]
                y2 = xyz_check[j][2]
                z1 = xyz_check[i][3]
                z2 = xyz_check[j][3]
                
                x1 = float(x1)
                x2 = float(x2)
                y1 = float(y1)
                y2 = float(y2)
                z1 = float(z1)
                z2 = float(z2)

                distance_check = m.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)

                if atom1 == atom2 and distance_check < 1.1:
                    atoms_closed = True
                    break

                if atom1 != atom2 and distance_check < 0.9:
                    atoms_closed = True
                    break

            i = i + 1

        i = 0
        j = 0

        out_of_sphere = False

        while (i < total_atoms):

            x1 = xyz_check[i][1]
            y1 = xyz_check[i][2]
            z1 = xyz_check[i][3]
            
            x1 = float(x1)
            y1 = float(y1)
            z1 = float(z1)
            
            distance_check = m.sqrt(x1 ** 2 + y1 ** 2 + z1 ** 2)
            
            if distance_check > cluster_radius:
            
                out_of_sphere = True
                break
                
            i = i + 1

        if not out_of_sphere and not atoms_closed:
            break

        iteration = iteration + 1

    with open(f'{current_dir}/anneal_trj.xyz', 'a') as anneal_trj_xyz:
    
        anneal_trj_xyz.write(f'{total_atoms} \n')
        anneal_trj_xyz.write(f'iteration number = {iteration + 1}\n')

    with open(f'{scratch_dir}/xyz', 'r') as xyz_final:
    
        xyz_check = []
        
        for line in xyz_final:
        
            line = [i.strip() for i in line.split(' ') if i]
            xyz_check.append(line)

    i = 0

    while (i < total_atoms):
    
        atom = xyz_check[i][0]
        x1 = xyz_check[i][1]
        y1 = xyz_check[i][2]
        z1 = xyz_check[i][3]
        
        x1 = float(x1)
        y1 = float(y1)
        z1 = float(z1)
        
        x1 = "{0:.8f}".format(x1, precision)
        y1 = "{0:.8f}".format(y1, precision)
        z1 = "{0:.8f}".format(z1, precision)
        
        xyz_coords = f'{atom}   {x1}    {y1}    {z1}'
        
        with open(f'{current_dir}/anneal.inp', 'a') as psi4_input:
        
            psi4_input.write(f'{str(xyz_coords)} \n')
            
        with open(f'{current_dir}/anneal_trj.xyz', 'a') as anneal_trj_xyz:
        
            anneal_trj_xyz.write(f'{str(xyz_coords)} \n')
            
        i = i + 1

    with open(f'{current_dir}/anneal.inp', 'a') as psi4_input:
    
        psi4_input.write('} \n')
        psi4_input.write(f' \n')
        psi4_input.write('set { \n')
        psi4_input.write(f'reference {psi_ref} \n')
        psi4_input.write(f'basis {basis} \n')
        psi4_input.write(f'scf_type df \n')
        
        if str(psi_method).lower() == 'mp2':
        
            psi4_input.write(f'mp2_type df \n')
            
        psi4_input.write(f'df_basis_scf {psi_basis_jkfit} \n')
        
        if str(psi_method).lower() == 'mp2':
        
            psi4_input.write(f'df_basis_mp2 {psi_basis_ri} \n')
            psi4_input.write('freeze_core true \n')
            
        psi4_input.write('} \n')
        psi4_input.write(' \n')
        
        if str(psi_method).lower() == 'mp2':
        
            psi4_input.write('energy("mp2")')
            
        else:
        
            psi4_input.write('energy("scf")')
            
        psi4_input.write(' \n')
        psi4_input.write(' \n')
        psi4_input.write(' \n')




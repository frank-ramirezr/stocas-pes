#! /usr/bin/env python3

import glob as glob
import math as m
import os as os
import shutil as sh
import subprocess
import sys
import time
import random
from modules import read_coord
from modules import coord_gen
from modules import acceptance
from modules import bonds_ini
from modules import bonds_opt
from modules import orca_gen
from modules import psi4_gen
from modules import rot_constants
from modules import str_comparison
from modules import quenching_route
from modules import minmax_rad

########################################################################################################################
# Reading input file
########################################################################################################################

input_file = open('input.dat', 'r')

lines = input_file.readlines()
line_number = []

for line in lines:
    line_number.append(line.split(';')[0])

method = line_number[0]
e_eval = line_number[1]
e_eval = int(e_eval)
temp_i = line_number[2]
temp_i = float(temp_i)
temp_f = line_number[3]
temp_f = float(temp_f)
temp_n = line_number[4]
temp_n = float(temp_n)
q_route = line_number[5]
q_route = int(q_route)

calc_rad, input_rad = line_number[6].split(",")

if calc_rad.lower() == 'true':
    calc_rad = True
else:
    calc_rad = False

input_rad = float(input_rad)

orca_model, composite = line_number[7].split(",")

if composite.lower() == 'composite':
    composite = True

elif composite.lower() == 'non-composite':
    composite = False

semiempiric_model = line_number[8]
hamiltonian = line_number[9]
basis = line_number[10]
charge, mult = line_number[11].split(',')
charge = int(charge)
mult = int(mult)
scf_maxiter = line_number[12]
scf_maxiter = int(scf_maxiter)
runnings = line_number[13]
runnings = int(runnings)

parallel, ex_cores = line_number[14].split(',')

if parallel.lower() == 'true':
    parallel = True
else:
    parallel = False

ex_cores = int(ex_cores)

orca_path = line_number[15]
psi_method = line_number[16]
psi_ref = line_number[17]
psi_basis_jkfit = line_number[18]
psi_basis_ri = line_number[19]
max_rot, max_trasl = line_number[20].split(',')
max_rot = float(max_rot)
max_trasl = float(max_trasl)
rot_dev, min_delta_e = line_number[21].split(',')
rot_dev = float(rot_dev)
min_delta_e = float(min_delta_e)
xtb_theory = line_number[22]
xtb_path = line_number[23]

input_file.close()

########################################################################################################################
# Delete previous running data files
########################################################################################################################

cdir = os.getcwd()

try:
    old_files = glob.glob(f'{cdir}/anneal*')
    for line in old_files:
        os.remove(line)
except FileNotFoundError:
    pass

try:
    os.remove(f'{cdir}/structures.dat')
    os.remove(f'{cdir}/summary.dat')
except FileNotFoundError:
    pass

########################################################################################################################
# DONE READING DATA FILE: INPUT
########################################################################################################################

read_coord.mol_gen()

if method == 'orca':
    
    orca_gen.orca_file_gen(method,
                           orca_model,
                           composite,
                           charge,
                           mult,
                           scf_maxiter,
                           semiempiric_model,
                           ex_cores,
                           parallel,
                           basis,
                           hamiltonian)


# Start energy evaluation, cycles

########################################################################################################################
# ORCA energy evaluation
########################################################################################################################

xyz_files_old = glob.glob(f'{cdir}/s*.xyz')
for xyz_files in xyz_files_old:
    os.remove(xyz_files)

running_old = glob.glob(f'{cdir}/run*')
for running in running_old:
    sh.rmtree(running)

scratch_dir = f'{cdir}/scratch'

# ######################################################################################################################
# ORCA RUNNING SECTION
# ######################################################################################################################

if method == 'orca':

    orca_run = f'{orca_path} anneal.inp > anneal.out 2> .error'
    orca_normal = 'FINAL SINGLE'
    n_run = 1

    # ######################################################################################################################
    # simulation extend from total iteration
    with open('simm_adv.out', 'w') as simm_adv_file:
        simm_adv_file.write('E step, Temperature (K), Simulation advance (%), nrunls')

        while n_run <= runnings:

            os.mkdir(f'{cdir}/run{n_run}')

            e0 = 0.0
            accepted_structure = 0
            accepted_e = 0
            accepted_delta = 0
            dir_index = 0
            dir_index = dir_index + 1
            temp_count = 0

            annpy_data_structures = []
            annpy_data_summary = []

            time_i = time.time()

            # ######################################################################################################################
            # simulation extend from total iteration

            current_iteration = 0
            total_iterations = e_eval * (temp_n + 1)

            temp_count = 0

            while temp_count <= temp_n:
                if q_route == 1:
                
                    actual_temp = quenching_route.linear(temp_i, temp_f, temp_n, temp_count)
                    pass
                if q_route == 2:

                    actual_temp = quenching_route.exponential(temp_i, temp_f, temp_n, temp_count)
                    pass
                if q_route == 3:

                    actual_temp = quenching_route.decay(temp_i, temp_f, temp_n, temp_count)
                    pass

                energy_count = 1
                
                actual_temp_print = '{0:.2f}'.format(actual_temp)

                while energy_count <= e_eval:

                    coord_gen.coord_gen(calc_rad, input_rad)

                    subprocess.call(orca_run, shell=True)

                    # ######################################################################################################################
                    # simulation extend from total iteration

                    current_iteration = current_iteration + 1
                    simm_adv = 100.0 * current_iteration / total_iterations
                    temp_ij = float(actual_temp)

                    with open('simm_adv.out', 'a') as simm_adv_file:
                        simm_adv_file.write(f'E_step: {energy_count}, '
                                            f'Temp: {temp_ij:.2f}, '
                                            f'Completed: {simm_adv:.2f}, '
                                            f'Run number: {n_run} \n')

                    # ######################################################################################################################

                    with open('anneal.out', 'r') as orca_terminated:
                        if 'TERMINATED NORMALLY' in orca_terminated.read():
                            orca_terminated = True
                        else:
                            orca_terminated = False

                    if not orca_terminated:
                        print('ORCA single point crashed')

                        energy_count = energy_count + 1

                        continue

                    if orca_terminated:
                        orca_file = open('anneal.out', 'r')
                        for line in orca_file:
                            if orca_normal in line:
                                line = line.split(' ')
                                e_iter = line[-1]
                                e_iter = float(e_iter)

                    if float(e_iter) < float(e0):

                        accepted_structure = accepted_structure + 1
                        structure = f's{accepted_structure}.xyz'

                        accepted_e = accepted_e + 1
                        last_temp_e = actual_temp
                        last_conf_e = structure

                        e_iter = '{0:.6f}'.format(e_iter)

                        n_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

                        with open(structure, 'a') as mol_index:
                        
                            mol_index.write(f'{n_atoms} \n')
                            energy_accepted = f'Energy = {e_iter} \n'
                            mol_index.write(energy_accepted)

                        with open(f'{scratch_dir}/xyz', 'r') as accepted:
                        
                            new_system = []
                            for line in accepted:
                            
                                line = [i.strip() for i in line.split(' ') if i]
                                new_system.append(line)

                        for i in range(0, n_atoms):
                        
                            x_coord = float(new_system[i][1])
                            x_coord = '{0:.8f}'.format(x_coord)
                            y_coord = float(new_system[i][2])
                            y_coord = '{0:.8f}'.format(y_coord)
                            z_coord = float(new_system[i][3])
                            z_coord = '{0:.8f}'.format(z_coord)
                            
                            xyz_file = f'{new_system[i][0]} {x_coord} {y_coord} {z_coord} \n'
                            
                            with open(structure, 'a') as file:
                            
                                file.write(xyz_file)
                                
                            e0 = e_iter

                        output = f'Structure accepted = {structure}, temperature = {actual_temp_print} K, energy = {e_iter} H \n'
                        file_structures = open(f'{cdir}/structures.dat', 'a')
                        file_structures.write(output)
                        file.close()

                    if float(e_iter) > float(e0):
                    
                        random_prob = random.uniform(0.0, 1.0)
                        prob_accept = acceptance.boltz(actual_temp, e0, e_iter)
                        prob_accept = float(prob_accept)

                        if prob_accept > random_prob:
                        
                            last_temp_b = actual_temp
                            accepted_delta = accepted_delta + 1
                            
                            prob_accept = '{0:.6f}'.format(prob_accept)
                            random_prob = '{0:.6f}'.format(random_prob)
                            
                            accepted_structure = accepted_structure + 1
                            structure = f's{accepted_structure}.xyz'
                            last_conf_delta = structure

                            e_iter = '{0:.6f}'.format(e_iter)

                            n_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

                            with open(structure, 'a') as mol_index:
                            
                                mol_index.write(f'{n_atoms} \n')
                                energy_accepted = f'Energy = {e_iter} \n'
                                mol_index.write(energy_accepted)

                            with open(f'{scratch_dir}/xyz', 'r') as accepted:
                            
                                new_system = []
                                for line in accepted:
                                
                                    line = [i.strip() for i in line.split(' ') if i]
                                    new_system.append(line)

                            for i in range(0, n_atoms):
                            
                                x_coord = float(new_system[i][1])
                                x_coord = '{0:.8f}'.format(x_coord)
                                y_coord = float(new_system[i][2])
                                y_coord = '{0:.8f}'.format(y_coord)
                                z_coord = float(new_system[i][3])
                                z_coord = '{0:.8f}'.format(z_coord)
                                
                                xyz_file = f'{new_system[i][0]} {x_coord} {y_coord} {z_coord} \n'
                                
                                with open(structure, 'a') as file:
                                
                                    file.write(xyz_file)

                            output = f'Boltzmann accepted = {structure}, temperature = {actual_temp_print} K, energy = ' \
                                     f'{e_iter} H, Random_Prob = {random_prob}, P_Boltzmann = {prob_accept} \n'
                                     
                            file_structures = open(f'{cdir}/structures.dat', 'a')
                            file_structures.write(output)
                            
                            file.close()

                    energy_count = energy_count + 1
                temp_count = temp_count + 1

            file_summmary = open(f'{cdir}/summary.dat', 'w')
            output = f'Total structures accepted = {accepted_e + accepted_delta} \n'
            file_summmary.write(output)
            output = f'Total structures accepted by energy = {accepted_e} \n'
            file_summmary.write(output)
            output = f'Total structures accepted by Boltzmann = {accepted_delta} \n'
            file_summmary.write(output)
            output = f'Lowest energy reached = {e0} H, at {last_temp_e} K, configuration = {last_conf_e} \n'
            file_summmary.write(output)
            try:
                output = f'Latest structure accepted by Boltzmann = {last_conf_delta}, energy = {e_boltzman} H, ' \
                         f'at {last_temp_b} K \n'
                file_summmary.write(output)
            except NameError:
                pass
                
            time_f = time.time()
            time_total = time_f - time_i
            tot_min, seconds = divmod(time_total, 60)
            hours, minutes = divmod(tot_min, 60)
            h = int(hours)
            m_val = int(minutes)
            sec = '{:.3f}'.format(seconds)

            total_time = f'total CPU time = {h} hours {m_val} minutes {sec} seconds \n'
            
            file_summmary.write(total_time)
                
            sh.move('anneal_trj.xyz', f'{cdir}/run{n_run}')
            sh.move('structures.dat', f'{cdir}/run{n_run}')
            sh.move('summary.dat', f'{cdir}/run{n_run}')

            final_structures = glob.glob(f'{cdir}/s*.xyz')
            
            for finals in final_structures:
                sh.move(finals, f'{cdir}/run{n_run}')
                
            n_run = n_run + 1

# ######################################################################################################################
# Done: ORCA annealing
# ######################################################################################################################

# ######################################################################################################################
# PSI4: RUNNING SECTION
# ######################################################################################################################

if str(method).lower() == 'psi4':

    read_coord.mol_gen()

    psi4_run = f'psi4 -i anneal.inp -o anneal.out -n{ex_cores}'
    psi4_normal_scf = 'Total Energy ='
    psi4_normal_mp2 = 'Total Energy              ='

    n_run = 1
    while n_run <= runnings:

        os.mkdir(f'{cdir}/run{n_run}')

        e0 = 0.0
        accepted_structure = 0
        accepted_e = 0
        accepted_delta = 0
        dir_index = 0
        dir_index = dir_index + 1
        temp_count = 0

        annpy_data_structures = []
        annpy_data_summary = []
        time_i = time.time()

# ######################################################################################################################
# simulation extend from total iteration

        current_iteration = 0
        total_iterations = e_eval * (temp_n + 1)

        temp_count = 0

        while temp_count <= temp_n:
            if q_route == 1:
            
                actual_temp = quenching_route.linear(temp_i, temp_f, temp_n, temp_count)
                pass
            if q_route == 2:

                actual_temp = quenching_route.exponential(temp_i, temp_f, temp_n, temp_count)
                pass
            if q_route == 3:

                actual_temp = quenching_route.decay(temp_i, temp_f, temp_n, temp_count)
                pass

            energy_count = 1

            while energy_count <= e_eval:
            
                psi4_gen.psi4_input_gen(charge,
                                        mult,
                                        psi_method,
                                        psi_ref,
                                        basis,
                                        psi_basis_jkfit,
                                        psi_basis_ri,
                                        calc_rad,
                                        input_rad)

                os.system(psi4_run)

# ######################################################################################################################
# simulation extend from total iteration

                current_iteration = current_iteration + 1
                simm_adv = 100.0 * current_iteration / total_iterations
                temp_ij = float(actual_temp)

                with open('simm_adv.out', 'a') as simm_adv_file:
                    simm_adv_file.write(f'E_step: {energy_count}, '
                                        f'Temp: {temp_ij:.2f}, '
                                        f'Completed: {simm_adv:.2f}, '
                                        f'Run number: {n_run} \n')

# ######################################################################################################################

                with open('anneal.out', 'r') as psi4_terminated:
                    if 'Psi4 exiting successfully' in psi4_terminated.read():
                        psi4_terminated = True
                    else:
                        psi4_terminated = False

                if not psi4_terminated:
                    print('PSI4 single point crashed')

                    energy_count = energy_count + 1

                    continue

                if psi4_terminated:
                
                    psi4_file = open('anneal.out', 'r')
                    
                    for line in psi4_file:
                    
                        if str(psi_method).lower() == 'mp2':
                        
                            if psi4_normal_mp2 in line:
                            
                                line = line.split(' ')
                                e_iter = line[-2]
                                e_iter = float(e_iter)
                                
                        elif str(psi_method).lower() == 'scf':
                        
                            if psi4_normal_scf in line:
                            
                                line = line.split(' ')
                                e_iter = line[-1]
                                e_iter = float(e_iter)

                actual_temp_print = '{0:.2f}'.format(actual_temp)

                if float(e_iter) < float(e0):

                    accepted_structure = accepted_structure + 1
                    structure = f's{accepted_structure}.xyz'

                    accepted_e = accepted_e + 1
                    last_temp_e = actual_temp
                    last_conf_e = structure

                    e_iter = '{0:.6f}'.format(e_iter)

                    n_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

                    with open(structure, 'a') as mol_index:
                    
                        mol_index.write(f'{n_atoms} \n')
                        energy_accepted = f'Energy = {e_iter} \n'
                        mol_index.write(energy_accepted)

                    with open(f'{scratch_dir}/xyz', 'r') as accepted:
                    
                        new_system = []
                        for line in accepted:
                        
                            line = [i.strip() for i in line.split(' ') if i]
                            new_system.append(line)

                    for i in range(0, n_atoms):
                    
                        x_coord = float(new_system[i][1])
                        x_coord = '{0:.8f}'.format(x_coord)
                        y_coord = float(new_system[i][2])
                        y_coord = '{0:.8f}'.format(y_coord)
                        z_coord = float(new_system[i][3])
                        z_coord = '{0:.8f}'.format(z_coord)
                        
                        xyz_file = f'{new_system[i][0]} {x_coord} {y_coord} {z_coord} \n'
                        
                        with open(structure, 'a') as file:
                        
                            file.write(xyz_file)
                            
                        e0 = e_iter

                    output = f'Structure accepted = {structure}, temperature = {actual_temp_print} K, energy = {e_iter} H \n'
                    file_structures = open(f'{cdir}/structures.dat', 'a')
                    file_structures.write(output)
                    file.close()

                if float(e_iter) > float(e0):
                
                    random_prob = random.uniform(0.0, 1.0)
                    prob_accept = acceptance.boltz(actual_temp, e0, e_iter)
                    prob_accept = float(prob_accept)

                    if prob_accept > random_prob:
                    
                        last_temp_b = actual_temp
                        accepted_delta = accepted_delta + 1
                        
                        prob_accept = '{0:.6f}'.format(prob_accept)
                        random_prob = '{0:.6f}'.format(random_prob)
                        
                        accepted_structure = accepted_structure + 1
                        structure = f's{accepted_structure}.xyz'
                        last_conf_delta = structure

                        e_iter = '{0:.6f}'.format(e_iter)

                        n_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

                        with open(structure, 'a') as mol_index:
                        
                            mol_index.write(f'{n_atoms} \n')
                            energy_accepted = f'Energy = {e_iter} \n'
                            mol_index.write(energy_accepted)

                        with open(f'{scratch_dir}/xyz', 'r') as accepted:
                        
                            new_system = []
                            for line in accepted:
                            
                                line = [i.strip() for i in line.split(' ') if i]
                                new_system.append(line)

                        for i in range(0, n_atoms):
                        
                            x_coord = float(new_system[i][1])
                            x_coord = '{0:.8f}'.format(x_coord)
                            y_coord = float(new_system[i][2])
                            y_coord = '{0:.8f}'.format(y_coord)
                            z_coord = float(new_system[i][3])
                            z_coord = '{0:.8f}'.format(z_coord)
                            
                            xyz_file = f'{new_system[i][0]} {x_coord} {y_coord} {z_coord} \n'
                            
                            with open(structure, 'a') as file:
                            
                                file.write(xyz_file)

                        output = f'Boltzmann accepted = {structure}, temperature = {actual_temp_print} K, energy = ' \
                                 f'{e_iter} H, Random_Prob = {random_prob}, P_Boltzmann = {prob_accept} \n'
                                 
                        file_structures = open(f'{cdir}/structures.dat', 'a')
                        file_structures.write(output)
                        
                        file.close()

                energy_count = energy_count + 1
            temp_count = temp_count + 1

        file_summmary = open(f'{cdir}/summary.dat', 'w')
        output = f'Total structures accepted = {accepted_e + accepted_delta} \n'
        file_summmary.write(output)
        output = f'Total structures accepted by energy = {accepted_e} \n'
        file_summmary.write(output)
        output = f'Total structures accepted by Boltzmann = {accepted_delta} \n'
        file_summmary.write(output)
        output = f'Lowest energy reached = {e0} H, at {last_temp_e} K, configuration = {last_conf_e} \n'
        file_summmary.write(output)
        try:
            output = f'Latest structure accepted by Boltzmann = {last_conf_delta}, energy = {e_boltzman} H, ' \
                     f'at {last_temp_b} K \n'
            file_summmary.write(output)
        except NameError:
            pass
            
        time_f = time.time()
        time_total = time_f - time_i
        tot_min, seconds = divmod(time_total, 60)
        hours, minutes = divmod(tot_min, 60)
        h = int(hours)
        m_val = int(minutes)
        sec = '{:.3f}'.format(seconds)

        total_time = f'total CPU time = {h} hours {m_val} minutes {sec} seconds \n'
        
        file_summmary.write(total_time)
            
        sh.move('anneal_trj.xyz', f'{cdir}/run{n_run}')
        sh.move('structures.dat', f'{cdir}/run{n_run}')
        sh.move('summary.dat', f'{cdir}/run{n_run}')

        final_structures = glob.glob(f'{cdir}/s*.xyz')
        
        for finals in final_structures:
            sh.move(finals, f'{cdir}/run{n_run}')
            
        n_run = n_run + 1

# ######################################################################################################################
#  PSI4 running finished
# ######################################################################################################################

# ######################################################################################################################
#  XTB Modules start -- Single point calculation version
# ######################################################################################################################

if str(method).lower() == 'xtb-sp':

    import subprocess


    def call_xtb(program_path, *args, output_file='anneal.out', working_dir=None):
        # Enforce a specific directory if provided, otherwise use current
        cwd = working_dir if working_dir else os.getcwd()

        try:
            with open(output_file, 'w') as outfile:
                # Combine stdout and stderr into the same file so you don't miss panic crashes
                result = subprocess.run(
                    [program_path, *args],
                    stdout=outfile,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=True,
                    cwd=cwd
                )
        except subprocess.CalledProcessError as e:
            pass

    anneal_file = 'anneal.xyz'

    if parallel:
        arguments = [
            anneal_file,
            '--' + str(xtb_theory),
            '--charge',
            str(charge),
            '--uhf',
            str(int(mult - 1)),
            '--iterations',
            str(scf_maxiter),
            '-P',
            str(ex_cores),
        ]

    elif not parallel:
        arguments = [
            anneal_file,
            '--' + str(xtb_theory),
            '--charge',
            str(charge),
            '--uhf',
            str(int(mult - 1)),
            '--iterations',
            str(scf_maxiter)
        ]

    xtb_energy_line = 'TOTAL ENERGY'

# ######################################################################################################################
# simulation extend from total iteration
    with open('simm_adv.out', 'w') as simm_adv_file:
        simm_adv_file.write('E step, Temperature (K), Simulation advance (%), nrunls')

    n_run = 1
    while n_run <= runnings:
        os.mkdir(f'{cdir}/run{n_run}')

        e0 = 0.0
        accepted_structure = 0
        accepted_e = 0
        accepted_delta = 0
        dir_index = 0
        dir_index = dir_index + 1
        temp_count = 0

        time_i = time.time()

# ######################################################################################################################
# simulation extend from total iteration

        current_iteration = 0
        total_iterations = e_eval * (temp_n + 1)

        e_current = 0.0
        e_global_min = 0.0

        while temp_count <= temp_n:
            if q_route == 1:
            
                actual_temp = quenching_route.linear(temp_i, temp_f, temp_n, temp_count)
                pass
            if q_route == 2:

                actual_temp = quenching_route.exponential(temp_i, temp_f, temp_n, temp_count)
                pass
            if q_route == 3:

                actual_temp = quenching_route.decay(temp_i, temp_f, temp_n, temp_count)
                pass

            energy_count = 1
            
            actual_temp_print = '{0:.2f}'.format(actual_temp)

            while energy_count <= e_eval:

                current_iteration = current_iteration + 1
                coord_gen.coord_gen(current_iteration, calc_rad, input_rad, max_rot, max_trasl, actual_temp, temp_i)

                call_xtb(xtb_path, *arguments, output_file='anneal.out')

# ######################################################################################################################
# simulation extend from total iteration

                simm_adv = 100.0 * current_iteration / total_iterations

                with open('simm_adv.out', 'a') as simm_adv_file:
                
                    simm_adv_file.write(f'E_step: {energy_count}, '
                                        f'Temp: {actual_temp:.2f}, '
                                        f'Completed: {simm_adv:.2f}, '
                                        f'Run number: {n_run} \n')

                with open('anneal.out', 'r') as xtb_final_opt:
                
                    if 'normal termination of xtb' in xtb_final_opt.read():
                        xtb_final_opt = True
                        
                    else:
                        xtb_final_opt = False

                if not xtb_final_opt:
                    print('XTB single point crashed')

                    energy_count = energy_count + 1

                    continue

                if xtb_final_opt:
                
                    xtb_file = open('anneal.out', 'r')
                    
                    for line in xtb_file:
                        if xtb_energy_line in line:
                            line = line.split(' ')
                            e_iter = line[-5]
                            e_iter = float(e_iter)
                
                if float(e_iter) < float(e0):

                    accepted_structure = accepted_structure + 1
                    structure = f's{accepted_structure}.xyz'

                    accepted_e = accepted_e + 1
                    last_temp_e = actual_temp
                    last_conf_e = structure

                    e_iter = '{0:.6f}'.format(e_iter)

                    n_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

                    with open(structure, 'a') as mol_index:
                    
                        mol_index.write(f'{n_atoms} \n')
                        energy_accepted = f'Energy = {e_iter} \n'
                        mol_index.write(energy_accepted)

                    with open(f'{scratch_dir}/xyz', 'r') as accepted:
                    
                        new_system = []
                        for line in accepted:
                        
                            line = [i.strip() for i in line.split(' ') if i]
                            new_system.append(line)

                    for i in range(0, n_atoms):
                    
                        x_coord = float(new_system[i][1])
                        x_coord = '{0:.8f}'.format(x_coord)
                        y_coord = float(new_system[i][2])
                        y_coord = '{0:.8f}'.format(y_coord)
                        z_coord = float(new_system[i][3])
                        z_coord = '{0:.8f}'.format(z_coord)
                        
                        xyz_file = f'{new_system[i][0]} {x_coord} {y_coord} {z_coord} \n'
                        
                        with open(structure, 'a') as file:
                        
                            file.write(xyz_file)
                            
                    e0 = e_iter

                    output = f'Structure accepted = {structure}, temperature = {actual_temp_print} K, energy = {e_iter} H \n'
                    file_structures = open(f'{cdir}/structures.dat', 'a')
                    file_structures.write(output)
                    file.close()

                if float(e_iter) > float(e0):
                
                    random_prob = random.uniform(0.0, 1.0)

                    prob_accept = acceptance.tsallis_accept(actual_temp, e0, e_iter, q=1.2)
                    prob_accept = float(prob_accept)

                    if prob_accept > random_prob:
                    
                        last_temp_b = actual_temp
                        accepted_delta = accepted_delta + 1
                        
                        prob_accept = '{0:.6f}'.format(prob_accept)
                        random_prob = '{0:.6f}'.format(random_prob)
                        
                        accepted_structure = accepted_structure + 1
                        structure = f's{accepted_structure}.xyz'
                        last_conf_delta = structure

                        e_iter = '{0:.6f}'.format(e_iter)

                        n_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

                        with open(structure, 'a') as mol_index:
                        
                            mol_index.write(f'{n_atoms} \n')
                            energy_accepted = f'Energy = {e_iter} \n'
                            mol_index.write(energy_accepted)

                        with open(f'{scratch_dir}/xyz', 'r') as accepted:
                        
                            new_system = []
                            for line in accepted:
                            
                                line = [i.strip() for i in line.split(' ') if i]
                                new_system.append(line)

                        for i in range(0, n_atoms):
                        
                            x_coord = float(new_system[i][1])
                            x_coord = '{0:.8f}'.format(x_coord)
                            y_coord = float(new_system[i][2])
                            y_coord = '{0:.8f}'.format(y_coord)
                            z_coord = float(new_system[i][3])
                            z_coord = '{0:.8f}'.format(z_coord)
                            
                            xyz_file = f'{new_system[i][0]} {x_coord} {y_coord} {z_coord} \n'
                            
                            with open(structure, 'a') as file:
                            
                                file.write(xyz_file)

                        e0 = e_iter

                        output = f'Boltzmann accepted = {structure}, temperature = {actual_temp_print} K, energy = ' \
                                 f'{e_iter} H, Random_Prob = {random_prob}, P_Boltzmann = {prob_accept} \n'
                                 
                        file_structures = open(f'{cdir}/structures.dat', 'a')
                        file_structures.write(output)
                        
                        file.close()

                energy_count = energy_count + 1
            temp_count = temp_count + 1

        file_summmary = open(f'{cdir}/summary.dat', 'a')
        output = f'Total structures accepted = {accepted_e + accepted_delta} \n'
        file_summmary.write(output)
        output = f'Total structures accepted by energy = {accepted_e} \n'
        file_summmary.write(output)
        output = f'Total structures accepted by Boltzmann = {accepted_delta} \n'
        file_summmary.write(output)
        output = f'Lowest energy reached = {e0} H, at {last_temp_e} K, configuration = {last_conf_e} \n'
        file_summmary.write(output)
        try:
            output = f'Latest structure accepted by Boltzmann = {last_conf_delta}, energy = {e_boltzman} H, ' \
                     f'at {last_temp_b} K \n'
            file_summmary.write(output)
        except NameError:
            pass
            
        time_f = time.time()
        time_total = time_f - time_i
        tot_min, seconds = divmod(time_total, 60)
        hours, minutes = divmod(tot_min, 60)
        h = int(hours)
        m_val = int(minutes)
        sec = '{:.3f}'.format(seconds)

        total_time = f'total CPU time = {h} hours {m_val} minutes {sec} seconds \n'
        
        file_summmary.write(total_time)

        file_summmary.close()

        try:
            os.makedirs(f'{cdir}/run{n_run}/anneal_log/')

        except FileExistsError:
            pass
            
        sh.move('anneal_trj.xyz', f'{cdir}/run{n_run}/anneal_log/')
        sh.move('structures.dat', f'{cdir}/run{n_run}/anneal_log/')
        sh.move('summary.dat', f'{cdir}/run{n_run}/anneal_log/')

        final_structures = glob.glob(f'{cdir}/s*.xyz')
        
        for finals in final_structures:
            sh.move(finals, f'{cdir}/run{n_run}')
            
        n_run = n_run + 1

# #####################################################################################################################
#   XTB: OPTIMIZATION FOLLOWED BY CARTESIAN COORDINATES COMPARISON
# #####################################################################################################################

if str(method).lower() == 'xtb-opt':

    def call_xtb(program_path, *args, output_file='anneal.out'):

        try:

            with open(output_file, 'w') as outfile:
                result = subprocess.run([program_path, *args],
                                        stdout=outfile,
                                        stderr=subprocess.PIPE,
                                        text=True,
                                        check=True)

        except subprocess.CalledProcessError as e:
            pass

    anneal_file = 'anneal.xyz'

    arguments = [
        anneal_file,
        '--opt',  # Separate flag
        '--cycles',
        '999',  # Separate flag and value
        '--' + str(xtb_theory),
        '--charge',
        str(charge),
        '--uhf',
        str(int(mult - 1)),
        '--iterations',
        str(scf_maxiter),
        '-P',
        str(ex_cores),
    ]

    xtb_normal = 'GEOMETRY OPTIMIZATION CONVERGED'
    xtb_energy = 'TOTAL ENERGY'
    rot_all = 'rotational constants'
    rot_comp = f'{scratch_dir}/rot_comp.dat'

    with open('simm_adv.out', 'w') as simm_adv_file:
        simm_adv_file.write('E step, Temperature (K), Simulation advance (%), nrunls')

    n_run = 1

    energy_list = []
    str_list = []

    while n_run <= runnings:
        os.mkdir(f'{cdir}/run{n_run}')

        e0 = 0.0
        accepted_structure = 0
        accepted_e = 0
        dir_index = 0
        dir_index = dir_index + 1

        energy_count = 1

        time_i = time.time()

# ######################################################################################################################
# simulation extend from total iteration
# ######################################################################################################################

# save structure function:
        def save_structure(structure_name, energy_val):
            """Reads xtbopt.xyz and saves formatted XYZ file."""
            with open('xtbopt.xyz', 'r') as f:
                num_atoms = int(f.readline().strip())
                comment_line = f.readline()

                atoms = []
                coords = []
                for _ in range(num_atoms):
                    parts = f.readline().split()
                    atoms.append(parts[0])
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

            with open(structure_name, 'w') as opt_data:
                opt_data.write(f'{num_atoms}\n')
                opt_data.write(f'Energy = {energy_val} Hartree\n')
                for i in range(num_atoms):
                    atom = atoms[i]
                    x = '{0:.8f}'.format(coords[i][0])
                    y = '{0:.8f}'.format(coords[i][1])
                    z = '{0:.8f}'.format(coords[i][2])
                    opt_data.write(f'{atom} {x} {y} {z}\n')

            return num_atoms

        while energy_count <= e_eval:

            coord_gen.coord_gen(energy_count, calc_rad, input_rad, max_rot, max_trasl, actual_temp=1, t0=1)

            initial_bonds = bonds_ini.compute_bonds()

            call_xtb(xtb_path, *arguments, output_file='anneal.out')

            with open('anneal.out', 'r') as xtb_file:

                if xtb_normal in xtb_file.read():

                    xtb_final_opt = True
                else:

                    xtb_final_opt = False

            if not xtb_final_opt:

                print(f'Iteration = {energy_count}; XTB OPT crashed')
                energy_count = energy_count + 1
                continue

# continue with the next iteration if reactive PES is found

            opt_bonds = bonds_opt.compute_bonds()

            reaction = True

            if sorted(initial_bonds) != sorted(opt_bonds):
                reaction = True

            else:
                reaction = False

            if reaction:
                print(f'Iteration: {energy_count}: Structure rejected reactive PES (init bonds = {initial_bonds}, opt_bonds = {opt_bonds})')

                energy_count = energy_count + 1
                continue

# extract structure's energy from optimization

            e_iter = None
            if not reaction:
                with open('anneal.out', 'r') as f:

                    for line in f:
                        if xtb_energy in line:
                            line = line.split()
                            e_iter = float(line[-3])

# #####################################################################################################################
# Structures comparison non-reactive interaction -- Bonds in initial configuration as the optimized structure
# #####################################################################################################################

            if xtb_final_opt and energy_count == 1:

                accepted_structure = accepted_structure + 1
                structure = f's{accepted_structure}.xyz'

# read optimization xyfile and copy optimized structure into formatted file

                save_structure(structure, e_iter)

# variable update

                energy_accepted = e_iter
                accepted_e = accepted_e + 1
                e_accepted = '{0:.8f}'.format(e_iter)

# append structures and energy in list

                str_list.append(structure)
                energy_list.append(e_iter)

# output accepted structure to structure.dat data file

                output = f'Structure accepted = {structure}, energy = {e_accepted} H \n'

                with open(file=f'{cdir}/structures.dat', mode='a') as f:
                    f.write(output)

# update rotational constants file
                with open(rot_comp, 'a') as rot_data:

                    rot_a, rot_b, rot_c = rot_constants.inertia()
                    rot_abc = f'{rot_a} {rot_b} {rot_c} \n'
                    rot_data.write(rot_abc)

# current minimum is updated and continue to next iteration
                e0 = e_iter
                energy_count = energy_count + 1
                continue

# #####################################################################################################################

            delta_e = 627.51 * abs(float(e_iter) - float(e0))

            if energy_count > 1 and xtb_final_opt and delta_e > min_delta_e:

                noneq_str, ra_dev_min, rb_dev_min, rc_dev_min = str_comparison.rot(rot_dev)

                if noneq_str:

                    accepted_structure = accepted_structure + 1
                    structure = f's{accepted_structure}.xyz'

                    dev_ia = '{0:.3f}'.format(ra_dev_min)
                    dev_ib = '{0:.3f}'.format(rb_dev_min)
                    dev_ic = '{0:.3f}'.format(rc_dev_min)

                    delta = '{0:.3f}'.format(delta_e)
                    e_accepted = '{0:.8f}'.format(e_iter)

                    output = f'Structure accepted = {structure}, energy = {e_accepted} H , delta_E = {delta} ' \
                             f'kcal/mol, % deviation (I_A, I_B, I_C) = ({dev_ia}, {dev_ib}, {dev_ic}) \n'

                    with open(file=f'{cdir}/structures.dat', mode='a') as f:
                        f.write(output)

# read optimization xyfile and copy optimized structure into formatted file

                    save_structure(structure, e_iter)

# variable update

                    energy_accepted = e_iter
                    accepted_e = accepted_e + 1
                    e_accepted = '{0:.8f}'.format(e_iter)

                    with open(rot_comp, 'a') as rot_data:

                        rot_a, rot_b, rot_c = rot_constants.inertia()
                        rot_abc = f'{rot_a} {rot_b} {rot_c} \n'
                        rot_data.write(rot_abc)

# current minimum is updated and continue to next iteration

                    accepted_e = accepted_e + 1
                    e0 = e_iter

                    str_list.append(structure)
                    energy_list.append(e_iter)

                elif not noneq_str:

                    print(
                        f'Iteration {energy_count}: Structure rejected (rotational'
                        ' duplicate).'
                    )

            elif delta_e <= min_delta_e:

                print(
                    f'Iteration {energy_count}: Structure rejected (delta_E = {delta_e:.3f}'
                    f' <= {min_delta_e} kcal/mol).'
                )

            energy_count = energy_count + 1

# #####################################################################################################################
# Calculate delta_e if energy count is greater than one and the XTB optimization is completed successfully
# #####################################################################################################################

        e_min = +1.0e10
        e_max = -1.0e10

        for energy, structure in zip(energy_list, str_list):

            if energy <= e_min:
                e_min = energy
                str_min = structure

            if energy >= e_max:
                e_max = energy
                str_max = structure

        de_max = -627.51 * (float(e_min) - float(e_max))
        de_max = '{0:.3f}'.format(de_max)

        glob_min = f'Global minimum structure: {str_min}; E = {e_min} H \n'
        glob_max = f'Most energetic structure: {str_max}; E = {e_max} H \n'
        delta_e = f'Max relative binding energy: {de_max} kcal/mol \n'

        tot_str_accepted = f'Total structures accepted = {accepted_e} \n'

        with open(file=f'{cdir}/summary.dat', mode='a') as f:

            f.write(tot_str_accepted)
            f.write(glob_min)
            f.write(glob_max)
            f.write(delta_e)

        time_f = time.time()
        time_total = time_f - time_i
        tot_min, seconds = divmod(time_total, 60)
        hours, minutes = divmod(tot_min, 60)
        h = int(hours)
        m_val = int(minutes)
        sec = '{:.3f}'.format(seconds)

        total_time = f'total CPU time = {h} hours {m_val} minutes {sec} seconds \n'

        with open(file=f'{cdir}/summary.dat', mode='a') as f:
            f.write(total_time)

        final_structures = glob.glob(f'{cdir}/s*.xyz')
        
        sh.move('anneal_trj.xyz', f'{cdir}/run{n_run}')
        sh.move('structures.dat', f'{cdir}/run{n_run}')
        sh.move('summary.dat', f'{cdir}/run{n_run}')
        
        for finals in final_structures:
            sh.move(finals, f'{cdir}/run{n_run}')

        n_run = n_run + 1

# #####################################################################################################################
# TEST EVALUATION
# #####################################################################################################################

if str(method).lower() == 'test':

    cdir = os.getcwd()
    os.mkdir(f'{cdir}/run_test')

    try:
        sh.rmtree(f'{cdir}/anneal*')
        old_distances_dat = f'{cdir}/distances.dat'
        os.remove(old_distances_dat)
    except FileNotFoundError:
        pass

    scratch_dir = f'{cdir}/scratch'

    accepted_structure = 0
    energy_count = 1

    while energy_count <= e_eval:

        coord_gen.coord_gen(calc_rad, input_rad, max_rot, max_trasl)

        max_dist, avg_dist = minmax_rad.minmax()

        max_dist = '{0:.6f}'.format(max_dist)
        avg_dist = '{0:.6f}'.format(avg_dist)

        with open('distances.dat', 'a') as distance_file:
            distance_file.write(f'{max_dist} {avg_dist} \n')

        n_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

        accepted_structure = accepted_structure + 1
        structure = f's{accepted_structure}.xyz'

        with open(structure, 'a') as mol_index:
            mol_index.write(f'{n_atoms} \n')
            mol_index.write('\n')

        with open(f'{scratch_dir}/xyz', 'r') as accepted:
            new_system = []
            for line in accepted:
                line = [i.strip() for i in line.split(' ') if i]
                new_system.append(line)

        for i in range(0, n_atoms):
            x_coord = float(new_system[i][1])
            x_coord = '{0:.8f}'.format(x_coord)
            y_coord = float(new_system[i][2])
            y_coord = '{0:.8f}'.format(y_coord)
            z_coord = float(new_system[i][3])
            z_coord = '{0:.8f}'.format(z_coord)
            xyz_file = f'{new_system[i][0]} {x_coord} {y_coord} {z_coord} \n'
            with open(structure, 'a') as file:
                file.write(xyz_file)

        energy_count = energy_count + 1

    final_structures = glob.glob(f'{cdir}/s*.xyz')
    sh.move('distances.dat', f'{cdir}/run_test')
    sh.move('anneal_trj.xyz', f'{cdir}/run_test')
    for finals in final_structures:
        sh.move(finals, f'{cdir}/run_test')

# ######################################################################################################################
# Simulated annealing completed
# ######################################################################################################################

try:
    old_files = glob.glob(f'{cdir}/anneal*')
    for line in old_files:
        os.remove(line)
    scratch_files = glob.glob(f'{scratch_dir}/*')
    for line in scratch_files:
        os.remove(line)
    os.rmdir(scratch_dir)
    xtb_files = glob.glob(f'{cdir}/xtb*')
    for line in xtb_files:
        os.remove(line)
    os.remove(f'{cdir}/wbo')
    os.remove(f'{cdir}/charges')
except FileNotFoundError:
    pass





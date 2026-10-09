# Simply ORCA file generator
import sys

def orca_file_gen(method,
                  orca_model,
                  composite,
                  charge, mult,
                  scf_maxiter,
                  semiempiric_model,
                  ex_cores,
                  parallel,
                  basis,
                  hamiltonian):

    if method == 'orca':

        orca_file = open('anneal.inp', 'w')

# generate orca input file for non-composite methods

        if not composite:

            if orca_model.lower() == 'sem':
                orca_file.write(f'! notrah {semiempiric_model.lower()} kdiis sp \n')
                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif orca_model.lower() == 'scf':
                orca_file.write(f'! notrah rhf hf kdiis sp \n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%basis \n')
                orca_file.write(f'basis {basis} \n')
                orca_file.write(f'auxj "autoaux" \n')
                orca_file.write(f'end \n')
                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif orca_model.lower() == 'dft':
                orca_file.write(f'! notrah rks {hamiltonian.lower()} kdiis sp \n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%basis \n')
                orca_file.write(f'basis {basis} \n')
                orca_file.write(f'auxj "autoaux" \n')
                orca_file.write(f'end \n')
                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif orca_model.lower() == 'dft-d':
                orca_file.write(f'! notrah rks {hamiltonian.lower()} kdiis d3 sp \n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%basis \n')
                orca_file.write(f'basis {basis} \n')
                orca_file.write(f'auxj "autoaux" \n')
                orca_file.write(f'end \n')
                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

# exit code if a Hamiltonian and non-composite method is selected
        composite_list = ['hf', 'b97', 'r2scan', 'pbeh', 'b3lyp', 'wb97x']
        non_composite_list = ['sem', 'dft', 'dft-b', 'scf']

        model_corectness = True

        if composite and hamiltonian.lower() not in composite_list:
            model_corectness = False
            print("*************************************************************************")
            print("*** ORCA NOT PROPERLY DEFINED. SELECT ONE OF THE ALLOWED COMBINATIONS ***")
            print("***                ANNEALING EXECUTION WILL BE ABORTED                ***")
            print("*************************************************************************")
            sys.exit(1)

        elif not composite and orca_model.lower() not in non_composite_list:
            print("*************************************************************************")
            print("*** ORCA NOT PROPERLY DEFINED. SELECT ONE OF THE ALLOWED COMBINATIONS ***")
            print("***                ANNEALING EXECUTION WILL BE ABORTED                ***")
            print("*************************************************************************")
            sys.exit(1)

        elif composite and hamiltonian.lower() in composite_list:
            model_corectness = True

# generate input file for orca selected composite methods

        if composite and model_corectness:

            if hamiltonian.lower() == 'hf':
                orca_file.write(f'! notrah hf-3c kdiis sp\n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif hamiltonian.lower() == 'b97':
                orca_file.write(f'! notrah b97-3c kdiis sp\n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif hamiltonian.lower() == 'r2scan':
                orca_file.write(f'! notrah r2scan-3c kdiis sp\n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif hamiltonian.lower() == 'pbeh':
                orca_file.write(f'! notrah pbeh-3c kdiis sp\n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif hamiltonian.lower() == 'b3lyp':
                orca_file.write(f'! notrah b3lyp-3c kdiis sp\n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()

            elif hamiltonian.lower() == 'wb97x':
                orca_file.write(f'! notrah wb97x-3c kdiis sp\n')
                if parallel:
                    orca_file.write(f'%pal nprocs {ex_cores} end \n')

                orca_file.write(f'%scf \n')
                orca_file.write(f'maxiter {scf_maxiter} \n')
                orca_file.write(f'end \n')
                orca_file.write(f'*xyzfile {charge} {mult} anneal.xyz \n')
                orca_file.close()


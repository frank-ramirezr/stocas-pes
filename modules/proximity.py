import math
import os as os


def proximity():
    # ##################################################################################################################
    global cube
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
        os.remove(f'{scratch_dir}/parall_file.dat')
        os.remove(f'{scratch_dir}/cube.dat')
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

            xmin = 1000
            xmax = -1000
            ymin = 1000
            ymax = -1000
            zmin = 1000
            zmax = -1000

            for i in range((k - 1) * natoms[m] + count1, k * natoms[m] + count1):

                x = float(xyz_file[i][1])
                y = float(xyz_file[i][2])
                z = float(xyz_file[i][3])

                if x < xmin:
                    xmin = x
                if x > xmax:
                    xmax = x
                if y < ymin:
                    ymin = y
                if y > ymax:
                    ymax = y
                if z > zmax:
                    zmax = z
                if z < zmin:
                    zmin = z

            xmin = xmin - 0.25
            xmax = xmax + 0.25
            ymin = ymin - 0.25
            ymax = ymax + 0.25
            zmin = zmin - 0.25
            zmax = zmax + 0.25

            xmin = '{0:.8f}'.format(xmin)
            xmax = '{0:.8f}'.format(xmax)
            ymin = '{0:.8f}'.format(ymin)
            ymax = '{0:.8f}'.format(ymax)
            zmin = '{0:.8f}'.format(zmin)
            zmax = '{0:.8f}'.format(zmax)

            parall_pip = f'{xmin}; {xmax}; {ymin}; {ymax}; {zmin}; {zmax}\n'

            with open(f'{scratch_dir}/parall_file.dat', 'a') as parall_data:
                parall_data.write(parall_pip)

    n_lines = 0
    parall_file = []

    with open(f'{scratch_dir}/parall_file.dat', 'r') as parall_data:
        for line in parall_data:
            n_lines = n_lines + 1
            parall_file.append(line.split(';'))

    min_distance = 1000

    cube = f'{scratch_dir}/cube.dat'

    space_violation = False
    for i in range(0, n_lines - 1):
        for j in range(0, n_lines):
            xspace = False
            yspace = False
            zspace = False
            xmin1 = float(parall_file[i][0])
            xmax1 = float(parall_file[i][1])
            ymin1 = float(parall_file[i][2])
            ymax1 = float(parall_file[i][3])
            zmin1 = float(parall_file[i][4])
            zmax1 = float(parall_file[i][5])
            xmin2 = float(parall_file[j][0])
            xmax2 = float(parall_file[j][1])
            ymin2 = float(parall_file[j][2])
            ymax2 = float(parall_file[j][3])
            zmin2 = float(parall_file[j][4])
            zmax2 = float(parall_file[j][5])
            if (xmin1 <= xmin2 <= xmax1) or (xmin1 <= xmax2 < xmax1):
                xspace = True
            if (ymin1 <= ymin2 <= ymax1) or (ymin1 <= ymax2 < ymax1):
                yspace = True
            if (zmin1 <= zmin2 <= zmax1) or (zmin1 <= zmax2 < zmax1):
                zspace = True
            if xspace and yspace and zspace and (i != j):
                space_violation = True
    for i in range(0, n_lines):
        xmin = float(parall_file[i][0])
        xmax = float(parall_file[i][1])
        ymin = float(parall_file[i][2])
        ymax = float(parall_file[i][3])
        zmin = float(parall_file[i][4])
        zmax = float(parall_file[i][5])

        data1 = f'{xmin}; {ymin}; {zmin} \n'
        data2 = f'{xmin}; {ymax}; {zmin} \n'
        data3 = f'{xmax}; {ymax}; {zmin} \n'
        data4 = f'{xmax}; {ymin}; {zmin} \n'
        data5 = f'{xmin}; {ymin}; {zmax} \n'
        data6 = f'{xmin}; {ymax}; {zmax} \n'
        data7 = f'{xmax}; {ymax}; {zmax} \n'
        data8 = f'{xmax}; {ymin}; {zmax} \n'
        with open(f'{cube}', 'a') as line:
            line.write(data1)
            line.write(data2)
            line.write(data3)
            line.write(data4)
            line.write(data5)
            line.write(data6)
            line.write(data7)
            line.write(data8)

    with open(cube, 'r') as data:
        nlines = 0
        cube_data = []
        for line in data:
            nlines = nlines + 1
            cube_data.append(line.split(';'))

    for line_number in range(0, n_lines - 1):
        i0 = line_number * 8
        i1 = (line_number + 1) * 8
        i2 = (line_number + 2) * 8
        for i in range(i0, i1):
            for j in range(i1, i2):
                x1 = float(cube_data[i][0])
                x2 = float(cube_data[j][0])
                y1 = float(cube_data[i][1])
                y2 = float(cube_data[j][1])
                z1 = float(cube_data[i][2])
                z2 = float(cube_data[j][2])
                distance = math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
                if distance < min_distance:
                    min_distance = distance

    # try:
    #     os.remove(f'{scratch_dir}/parall_file.dat')
    #     os.remove(f'{scratch_dir}/cube.dat')
    # except FileNotFoundError:
    #     pass
    #
    return min_distance, space_violation

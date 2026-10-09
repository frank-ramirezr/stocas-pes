import os as os
from modules import radius
import math as m


def out_of_sphere(calc_rad, input_rad):

    cdir = os.getcwd()
    scratch_dir = f'{cdir}/scratch'

    if calc_rad:
        cluster_radius = radius.radius()
    else:
        cluster_radius = input_rad

    with open(f'{scratch_dir}/xyz', 'r') as xyz:
        xyz_check = []
        for line in xyz:
            line = [i.strip() for i in line.split(' ') if i]
            xyz_check.append(line)

    total_atoms = sum(1 for line in open(f'{scratch_dir}/xyz'))

    out_of_sphere = False

    for i in range(0, total_atoms):
        x = xyz_check[i][1]
        y = xyz_check[i][2]
        z = xyz_check[i][3]
        x = float(x)
        y = float(y)
        z = float(z)
        distance_check = m.sqrt(x ** 2 + y ** 2 + z ** 2)
        if distance_check > cluster_radius:
            out_of_sphere = True
            return out_of_sphere





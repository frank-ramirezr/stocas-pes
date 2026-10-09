import math as m
import os as os


def minmax():
    anneal = []
    cdir = os.getcwd()

    # Read and parse the XYZ file safely
    with open(f'{cdir}/anneal.xyz', 'r') as opt_file:
        for line in opt_file:
            line_elements = line.split()
            anneal.append(line_elements)

    # Cleanly cast the atom count to integer
    n_atoms = int(anneal[0][0].strip())
    n_lines = n_atoms + 2

    max_distance = -1e10
    min_distance = +1e10
    sum_distance = 0.0
    sum_iter = 0

    # Double loop to calculate all unique pairwise interatomic distances
    for i in range(2, n_lines):
        # Extract and cast coords for atom 1 once per outer iteration (improves speed)
        x1 = float(anneal[i][1])
        y1 = float(anneal[i][2])
        z1 = float(anneal[i][3])

        for j in range(i + 1, n_lines):
            x2 = float(anneal[j][1])
            y2 = float(anneal[j][2])
            z2 = float(anneal[j][3])

            # Distance computation via Euclidean metric
            distance = m.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
            sum_iter += 1

            if distance > max_distance:
                max_distance = distance

            sum_distance += distance

    # Guard against dividing by zero if the system contains fewer than 2 atoms
    average_distance = sum_distance / float(sum_iter) if sum_iter > 0 else 0.0

    return max_distance, average_distance
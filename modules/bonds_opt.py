from collections import defaultdict, deque
import math as m
import numpy as np
import os as os
import glob as glob
import re as re
import sys

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [atoi(c) for c in re.split(r'(\d+)', text)]

ptable = {
    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
    'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
    'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
    'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
    'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
    'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53
}

covalent_radii = np.array([
    0.00, 0.32, 0.46, 1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.62, # 0-9
    1.55, 1.39, 1.30, 1.21, 1.11, 1.07, 1.05, 1.02, 2.03, 1.76, # 10-19
    1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, # 20-29
    1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 2.20, 1.95, 1.90, 1.75, # 30-39
    1.64, 1.56, 1.47, 1.46, 1.42, 1.39, 1.44, 1.44, 1.42, 1.39, # 40-49
    1.39, 1.39, 1.38, 1.39                                       # 50-53
])

multiple_bond_rules = {
    (6, 6): [(1.24, 2), (1.37, 1)], # C≡C gets +2 (total 3), C=C gets +1 (total 2), else single (1)
    (7, 7): [(1.12, 2), (1.28, 1)], # N≡N, N=N
    (8, 8): [(1.23, 1)],            # O=O
    (6, 7): [(1.18, 2), (1.32, 1)], # C≡N, C=N
    (6, 8): [(1.16, 2), (1.26, 1)], # C≡O, C=O
    (7, 8): [(1.18, 2), (1.26, 1)], # N≡O, N=O
    (5, 7): [(1.24, 2), (1.42, 1)]  # B≡N, B=N
}

# 2. Read XYZ File

def compute_bonds():
    atoms = []
    coords = []

    with open('xtbopt.xyz', 'r') as f:
        num_atoms = int(f.readline().strip())
        comment_line = f.readline()

        for _ in range(num_atoms):
            parts = f.readline().split()
            atoms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    coords = np.array(coords)
    znus = np.array([ptable[atom] for atom in atoms], dtype='int')

    # 3. Build Adjacency Graph with Bond Order
    adj_list = defaultdict(list)

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            dist = m.dist(coords[i], coords[j])

            z_min, z_max = sorted((znus[i], znus[j]))
            r_expected = covalent_radii[z_min] + covalent_radii[z_max]
            bond_cutoff = 1.1 * r_expected

            # Check if a single covalent bond exists
            if 0.80 <= dist <= bond_cutoff:
                bond_order = 1  # Base single bond

                # Check for multiple bond rules
                pair_key = (z_min, z_max)
                if pair_key in multiple_bond_rules:
                    for threshold, extra_order in multiple_bond_rules[pair_key]:
                        if dist <= threshold:
                            bond_order += extra_order
                            break  # Matched the shortest/strongest threshold rule

                adj_list[i].append((j, bond_order))
                adj_list[j].append((i, bond_order))

    # 4. Connected Components Search (BFS to Assign Fragment/Molecule IDs)
    visited = [False] * num_atoms
    mol_labels = [0] * num_atoms
    current_mol_id = 1

    for i in range(num_atoms):
        if not visited[i]:
            queue = deque([i])
            visited[i] = True

            while queue:
                curr = queue.popleft()
                mol_labels[curr] = current_mol_id

                for neighbor, order in adj_list[curr]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)

            current_mol_id += 1

    # 5. Count Total Bonds (Including Bond Order) per Molecular Fragment
    num_molecules = max(mol_labels)
    bonds_per_molecule = []

    for mol_id in range(1, num_molecules + 1):
        atom_indices = [
            idx for idx, label in enumerate(mol_labels) if label == mol_id
        ]

        # Sum the bond order for every bond connected to atoms in this fragment
        total_bonds = (
                sum(
                    sum(order for neighbor, order in adj_list[atom_idx])
                    for atom_idx in atom_indices
                )
                // 2
        )

        bonds_per_molecule.append(total_bonds)


    return bonds_per_molecule

if __name__ == '__main__':

    total_bonds = compute_bonds()

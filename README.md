# STOCAS-PES

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python Version](https://img.shields.io/badge/python-3.8%2B-blue)](https://www.python.org/)
[![QC Engines](https://img.shields.io/badge/QC_Engines-xTB%20%7C%20ORCA%20%7C%20Psi4-teal)](#quantum-chemistry-interfaces)

**STOCAS-PES** (**STO**chastic **C**luster **A**rchitecture & **S**ampling for **P**otential **E**nergy **S**urfaces) is a ground-up, open-source Python framework designed for the unbiased stochastic exploration of potential energy surfaces (PES) of atomic, molecular, and mixed host-guest/solvated clusters.

By combining rigid-body pseudo-random walk perturbations with spherical confinement boundaries, steric van der Waals collision filtering, and a generalized **Tsallis non-extensive statistical mechanics** acceptance protocol, STOCAS-PES dynamically hops over high-energy barriers to discover low-lying conformational motifs and global minima.

---

## Key Features

* **Universal Cluster Sampling:** Works seamlessly on pure atomic clusters, pure molecular water/solvent networks, and mixed ion–molecule aggregates (e.g., $[(\text{H}_2\text{O})_n \cdot\cdot\cdot X]^-$).
* **Spherical Boundary Confinement:** Constrains translations ($\Delta r_{\text{max}}$) and rotations ($\Delta \theta_{\text{max}}$) within spherical limits, preventing molecular evaporation and eliminating box-corner trapping artifacts.
* **Steric Collision Filtering:** Pre-screens candidate trial geometries using interatomic van der Waals radii sums ($r_{ij} < r_{i,\text{vdW}} + r_{j,\text{vdW}}$) to avoid unphysical overlapping before invoking expensive quantum calculations.
* **Tsallis Barrier Hopping:** Features dual-acceptance logic—direct energy minimization plus a Tsallis non-extensive distribution mechanism to escape deep local minima on complex multi-basin landscapes.
* **Custom Quenching Schedules:** Offers three distinct cooling routines targeting different interaction strengths:
  * Linear Decay ($T_{\text{lin}}$)
  * Exponential Decay ($T_{\text{exp}}$)
  * Power/Aggressive Decay ($T_{\text{dec}}$)
* **Quantum Engine Interoperability:** Interfaces directly with open-source and commercial quantum chemistry packages:
  * **xTB** (GFN0-xTB, GFN1-xTB, GFN2-xTB, GFN-FF)
  * **Psi4** (DFT, MP2, HF)
  * **ORCA** (DFT, Composite methods like $r^2\text{SCAN-3c}$, Coupled Cluster)
* **Automated Multi-Seed Workspaces:** Manages multi-run ensemble sampling ($P$ independent runs) with distinct pseudo-random seed generation and automated directory structuring (`run_k/`).

---

## Algorithm Overview

Master Loop (Run k = 1 ... P)

├── Initialize PRNG Seed & Output Directory (run_k/)

├── Center Molecules (Subtract Center of Mass)

└── Temperature Loop (T_step = 1 ... T_steps)

├── Compute Current Temperature via Quenching Function

└── Stochastic Walk Loop (Cycle = 1 ... N_eval)

├── Sequential Translation & Rotation Differential

├── Overlap Filter ($r_{ij} < r_{i,vdW} + r_{j,vdW}$) $\rightarrow$ [If Overlap] $\rightarrow$ Re-perturb

├── Execute Quantum Engine (xTB / ORCA / Psi4)

└── Accept / Reject Logic:

├── Direct Energy Accept (E_new < E_current)

└── Tsallis Thermal Accept (P_Tsallis > U(0,1))

## Installation

### Prerequisites

Ensure you have Python 3.8+ installed along with standard scientific libraries:

```bash
pip install numpy scipy matplotlib
```

xTB: conda install -c conda-forge xtb-python or via binary executable.

Psi4: conda install -c conda-forge psi4

ORCA: Ensure the orca binary path is exported in your environment (export PATH="/path/to/orca:$PATH").

git clone [https://github.com/your-username/stocas-pes.git](https://github.com/frank-ramirezr/stocas-pes.git)
cd stocas-pes

## How it works

### Input File Specifications

In the current release (**v1.0**), STOCAS-PES requires two specific input files located in your execution directory. **The file names are strictly fixed in v1.0 and must not be altered**:

1. `input.dat` — Controls execution parameters, quantum chemistry engine selection, and search constraints.
2. `comp.xyz` — Defines the molecular composition, Cartesian coordinates, and mobility flags for each cluster component.

### Program Configuration File: `input.dat`

The `input.dat` file contains positional parameter lines (comments separated by `!` are supported). Below is a complete annotated template:

```text
xtb-sp;          ! 1. Engine & Mode: orca, xtb-sp (single point), xtb-opt (optimization), or psi4
1000;            ! 2. Energy evaluations per temperature step (Integer)
600.0;           ! 3. Initial temperature T0 in Kelvin (Real)
100.0;           ! 4. Final target temperature in Kelvin (Real)
100;             ! 5. Total temperature steps nT (Integer)
1;               ! 6. Quenching schedule: 1 = Linear, 2 = Exponential, 3 = Power Decay
false, 4.5;      ! 7. Estimate sphere radius (Boolean), override radius value (Real in Å)
scf, non-composite;  ! 8. ORCA Theory: dft, dft-d, sem, scf (HF), or defined composite method
pm3;             ! 9. Semiempirical Method (ORCA): AM1, MNDO, PM3
hf;              ! 10. Functional / Hamiltonian: hf, b97, r2scan, pbeh, b3lyp, wb97x
"6-31+g*";       ! 11. Orbital Basis Set (and general auxiliary basis)
-1, 1;           ! 12. Total Cluster Charge (Integer), Spin Multiplicity 2S+1 (Integer)
300;             ! 13. SCF Maximum Iterations (Integer)
3;               ! 14. Total independent runs / sampling repetitions (Integer)
false, 4;        ! 15. Parallel execution flag (Boolean), execution cores (Integer)
/home/frank/bin/orca/orca;  ! 16. Full path to ORCA binary executable
scf;             ! 17. Psi4 Method: scf or mp2
rhf;             ! 18. Psi4 Reference: rhf or uhf
def2-universal-jkfit;  ! 19. Psi4 JK auxiliary basis set for SCF
def2-svp-ri;     ! 20. Psi4 RI auxiliary basis set for MP2
1.0, 1.0;        ! 21. Max rotation angle (radians), Max translation step (Å)
2.5, 0.5;        ! 22. xtb-opt constraints: Rotational constant dev (%), Min energy dev (kcal/mol)
gxtb;            ! 23. xTB Theory level: gxtb, gfn1, gfn2, gfn-ff
xtb;             ! 24. Full path to xTB binary executable
```

### Cluster composition

This is an example of `mol.xyz` file of the water hexamer with bromide anion.

```text
mol 3 6 false
O      0.285672    0.905654   -0.000003
H      0.995949    0.685137    0.591079
H      0.572091    1.584266   -0.599961
*
mol 1 1 false
Br     0.572091    1.584266   -0.599961
*
```

**Header Syntax (mol <num_atoms> <num_mol> <frozen_flag>)**

**mol**: Keyword tag identifying a molecular block.

**<num_atoms>**: Number of atoms in a single molecule unit (e.g., 3 for $\text{H}_2\text{O}$, 1 for $\text{Br}^-$).

**<num_mol>**: Number of identical instances of this molecule in the cluster (e.g., 6 for six water molecules).

**<frozen_flag>**: Boolean flag (true/false). If set to true, this molecular species remains structurally frozen/anchored in space during the random-walk stochastic exploration.

*: Required block terminator line. Must be placed immediately after the Cartesian coordinates of each molecular species. If omitted, subsequent molecules in the file will not be parsed.

## Execution

*To launch STOCAS-PES, ensure **input.dat** and **comp.xyz** reside in your current working directory*:
```bash
python run_stocas.py
```

### Citation

If you use STOCAS-PES in your research, please cite our paper:
*Once the article is available online the full name and citation will be posted here*. 

### License 

Distributed under the GNU General Public License v3.0 (GPL-3.0). See LICENSE for full legal text.


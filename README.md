# Vibrational Frequency and Thermodynamics Calculations Using Quantum ESPRESSO

Python scripts for calculating vibrational frequencies using the **finite-displacement method** with **Quantum ESPRESSO (`pw.x`)**.

The scripts automate the generation of displaced atomic structures, preparation of Quantum ESPRESSO input files, extraction of atomic forces from `pw.x` output, and construction of the numerical Hessian, vibrational frequencies and thermodynamics corrections.

> **Note:** This workflow is different from the standard Quantum ESPRESSO `ph.x` phonon calculation. `ph.x` uses Density Functional Perturbation Theory (DFPT), whereas this repository uses explicit atomic displacements and numerical differentiation of forces.

 **Disclaimer**

This manual focuses on the **technical use of the Python scripts and their integration with Quantum ESPRESSO**. It does not provide a comprehensive discussion of the scientific principles underlying vibrational frequency calculations, finite-displacement methods, density functional theory (DFT), or the interpretation of calculated results.

The scripts are provided as computational tools and **the scientific validity, accuracy, and interpretation of the results are not guaranteed**. Users are responsible for ensuring that the computational settings, numerical parameters, and resulting data are appropriate for their specific research problem.

Users are expected to have a **sufficient background in computational chemistry, quantum chemistry, density functional theory, and Quantum ESPRESSO** to understand the calculations being performed and to critically evaluate the resulting data.

In particular, users should independently verify:

* the suitability of the chosen DFT method and pseudopotentials;
* convergence with respect to computational parameters;
* the choice of displacement magnitude;
* the quality and convergence of the reference structure;
* the numerical accuracy of the calculated forces and Hessian;
* the presence and origin of imaginary frequencies;
* the treatment of translational and rotational modes;
* and the physical interpretation of the calculated vibrational frequencies and thermodynamics correction.


---

## 1. Overview

The workflow implemented in this repository is:

```text
                    Optimized structure
                           │
                           ▼
                        displace.py
                           │
                           ▼
             Displaced QE input files
                           │
                           ▼
                     pw.x calculations
                           │
                           ▼
                    Atomic forces
                           │
                           ▼
                     qe_to_vib.py
                           │
                           ▼
                 Numerical Hessian
                           │
                           ▼
                  Mass-weighted Hessian
                           │
                           ▼
                 Vibrational frequencies and thermodynamics corrections
```

The method is based on numerical differentiation of the atomic forces.

For a displacement along coordinate (j),


$H_{ij} = -\frac{\partial F_i}{\partial x_j}$

where $(H_{ij})$ is an element of the Hessian matrix and ($F_i$) is the force acting on coordinate ($i$).

Using central finite differences,


$H_{ij} \approx -\frac{F_i(x_j+\Delta x)-F_i(x_j-\Delta x)}{2\Delta x}$


The resulting Hessian is mass-weighted and diagonalized to obtain the vibrational frequencies.

---

# 2. Requirements

## 2.1 Python with NumPy

Python 3 is required.

Check your Python version:

```bash
python --version
```

or:

```bash
python3 --version
```

Recommended:

```text
Python >= 3.9
```

---

## 2.2 Quantum ESPRESSO


The main executable required by this workflow is:

```text
pw.x
```

For reference, Quantum ESPRESSO documents the `pw.x` input variables in its official input documentation.

---

# 3. Repository Structure



```text
qe-finite-displacement/
│
├── README.md
│
├── scripts/
│   ├── displacement.py
│   └── qe_to_vib.py
│
└── examples/
    ├── molecule-H2/
    └── adsorbed-HNO3/

```

The main scripts are:

```text
displace.py
qe_to_vib.py
```

---

# 4. Calculation Workflow

The calculation consists of four main steps:

### Step 1

Obtain a fully optimized structure.

### Step 2

Generate positive and negative atomic displacements.

### Step 3

Run `pw.x` for every displaced structure.

### Step 4

Extract forces and calculate the numerical Hessian and vibrational frequencies.

---

# 5. Step 1 — Optimize the Reference Structure

Before generating displacements, the reference structure should be optimized.

For example:

```text
calculation='relax'
```

or, when optimizing both atomic positions and lattice parameters:

```text
calculation='vc-relax'
```

A typical relaxation input may look like:

```text
&CONTROL
    calculation = 'relax',
    prefix      = 'system',
    outdir      = './tmp/',
    pseudo_dir  = './pseudo/',
    tprnfor     = .true.,
/

&SYSTEM
    ibrav       = 0,
    nat         = 3,
    ntyp        = 2,
    ecutwfc     = 60.0,
    ecutrho     = 480.0,
/

&ELECTRONS
    conv_thr    = 1.0d-10,
/

&IONS
/

ATOMIC_SPECIES
H  1.008   H.pbe-rrkjus.UPF
O  15.999  O.pbe-rrkjus.UPF

ATOMIC_POSITIONS angstrom
O   ...
H   ...
H   ...

K_POINTS gamma
```

For a relaxation, the force convergence criterion is controlled by `forc_conv_thr`. Quantum ESPRESSO defines this as the threshold for the maximum force component during ionic optimization.

Run:

```bash
pw.x < relax.in > relax.out
```


After convergence, obtain the optimized geometry.

---

# 6. Step 2 — Generate Atomic Displacements

The script `displace.py` generates the displaced structures required for the numerical Hessian.

```bash
python3 displace.py
```

The initial QE input containing optimized structure without displacement must be present in the same directory with the `displace.py` script. The geometry must be given in angstrom unit.

```text
ATOMIC_POSITIONS angstrom
Zn   ...
O    ...
O    ...

```
 
Modify the `displace.py` script accordingly under `User settings`


```text
# =========================
# User settings
# =========================
qe_input = "disp_000.in" # name of the initial QE input containing optimized structure without displacement
disp = 0.02  # Displacement magnitude in Angstrom
atoms_to_move = [54, 55] # index of moved atoms, the ordering must be consistent with QE input

```

For each Cartesian coordinate, two structures are generated:

```text
+Δx
-Δx
```


For a system containing $N$ (moved) atoms, there are $3N$ Cartesian coordinates.

Therefore, a conventional central-difference calculation requires:

$2(3N)=6N$

displaced calculations.

For example, for a 10-(moved) atom(s):

```text
3 × 10 = 30 coordinates

30 × 2 = 60 displaced calculations
```

For example :

```text
Δx = 0.02 Å
```

For atom (i), the generated structures may contain:

```text
x + 0.02 Å
x - 0.02 Å
```

while all other coordinates remain unchanged.

There will be 6N generated QE input files with atomic displacement with the name formatting

```text
disp_001_atom(atoms_to_move[0])_+1x.in
disp_002_atom(atoms_to_move[0])_-1x.in
...
disp_6N_atom(atoms_to_move[-1])_-1z.in
```

---

# 7. Quantum ESPRESSO Input Files

The important point is that **all electronic-structure parameters should remain identical between displaced structures**.

These include:

* exchange-correlation functional;
* pseudopotentials;
* plane-wave cutoff;
* charge-density cutoff;
* k-point mesh;
* smearing parameters;
* spin configuration;
* DFT+U parameters;
* van der Waals correction;
* convergence criteria.

Only the atomic coordinates should change.

The displaced structures should normally be evaluated using:

```text
calculation='scf'
```

and **not** re-optimized.

The purpose of the displacement is to evaluate the force at a specific displaced geometry.

Therefore:

Do not use:

```text
calculation='relax'
```

for the displaced structures unless there is a specific reason to do so.

Otherwise, the imposed displacement would be removed by the geometry optimization.

---

# 8. Example Displacement Input

A displaced calculation may look like:

```text
&CONTROL
    calculation = 'scf',
    prefix      = 'disp_001',
    outdir      = './tmp/',
    pseudo_dir  = './pseudo/',
    tprnfor     = .true.,
/

&SYSTEM
    ibrav       = 0,
    nat         = 3,
    ntyp        = 2,
    ecutwfc     = 60.0,
    ecutrho     = 480.0,
/

&ELECTRONS
    conv_thr    = 1.0d-10,
/

ATOMIC_SPECIES
H  1.008   H.pbe-rrkjus.UPF
O  15.999  O.pbe-rrkjus.UPF

ATOMIC_POSITIONS angstrom
O   ...
H   ...
H   ...

K_POINTS gamma
```

The `tprnfor` option requests force calculation in `pw.x`

---

# 9. Running the Displacement Calculations

Run all displaced geometries

```bash
pw.x < disp_001_atom(atoms_to_move[0])_+1x.in > disp_001_atom(atoms_to_move[0])_+1x.out
...
pw.x < disp_6N_atom(atoms_to_move[-1])_-1z.in > disp_6N_atom(atoms_to_move[-1])_-1z.out
```


---

# 10. Checking Quantum ESPRESSO Output

Before calculating frequencies, verify that every `pw.x` calculation completed successfully.

A simple check is:

```bash
grep "JOB DONE" *.out
```

You should obtain one successful completion message for every displacement.

Also check that forces are present:

```bash
grep "Forces acting on atoms" */pw.out
```

---

# 11. Calculating Frequencies and Thermodynamics Correction

After all QE calculations have finished:

Modify the `qe_to_vib.py` script accordingly under `User input`

Example:

```text
disp_ang = 0.02 # displaced magnitude. Must be consistent with the displaced calculations. Unit in angstrom
nat = 56 # Total number of atoms 
ads_atoms = [54, 55, 56] Displaced atom line. The order must be consistent with the order in QE input
moved_atoms = [N, H, H] # Type of displaced atoms. The order must be consistent with the order in QE input
T = 298.15  # Temperature in K
P = 1.0 # Pressure in atm
freq_cutoff_cm = 50.0   # Cutoff frequency. All positive frequencies below cutoff will be replaced by the cutoff value for thermodynamic corrections. Unit in cm^-1

```
Run the script

```bash
python qe_to_vib.py
```

A typical output may be:

```text
Vibrational frequencies (cm^-1):
Mode  1:      38.21 cm^-1
Mode  2:      52.81 cm^-1
Mode  3:     247.11 cm^-1
Mode  4:     289.71 cm^-1
Mode  5:     474.01 cm^-1
Mode  6:     579.86 cm^-1
Mode  7:    1470.47 cm^-1
Mode  8:    3431.09 cm^-1
Mode  9:    3535.60 cm^-1

ZPE               = 0.628021 eV
E_vib(298.15 K)      = 0.081867 eV
S_vib(298.15 K)      = 6.061751e-04 eV/K
T*S_vib(298.15 K)    = 0.180731 eV

```

---


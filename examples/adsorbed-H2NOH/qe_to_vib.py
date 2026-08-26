import numpy as np
import glob
import re

# =========================
# Physical constants
# =========================
RY_TO_EV = 13.605693
BOHR_TO_M = 0.529177e-10
EV_TO_J = 1.602176634e-19
AMU_TO_KG = 1.66053906660e-27
ATM_TO_PA = 101325
HBAR = 1.054571817e-34
KB = 1.380649e-23
C_CM_S = 2.99792458e10

# =========================
# Atomic masses (amu)
# =========================
H  = 1.008
He = 4.002602

Li = 6.94
Be = 9.0121831
B  = 10.81
C  = 12.011
N  = 14.007
O  = 15.999
F  = 18.998403163
Ne = 20.1797

Na = 22.98976928
Mg = 24.305
Al = 26.9815385
Si = 28.085
P  = 30.973761998
S  = 32.06
Cl = 35.45
Ar = 39.948

K  = 39.0983
Ca = 40.078

# =========================
# User input
# =========================
disp_ang = 0.02 # must be consistent with the displaced atoms in angstrom
nat = 58
ads_atoms = [54, 55, 56, 57, 58]   
moved_atoms = [N, O, H, H, H] #the order must be consistent with the order in QE input
T = 298.15  # temperature (K)
P = 1.0 # Pressure (atm)
freq_cutoff_cm = 50.0   # cm^-1
# =========================
# System definition
# =========================
disp_m = disp_ang * 1e-10
masses = np.array(moved_atoms) * AMU_TO_KG
outfiles = sorted(glob.glob("disp_*.out"))
ads_atoms_py = [a - 1 for a in ads_atoms]
moved_atom = len(ads_atoms)
total_mass = np.sum(masses)
molecule = nat == moved_atom
P = P * ATM_TO_PA

# =========================
# QE force parser
# =========================
def read_forces(outfile):

    forces = []

    with open(outfile) as f:
        lines = f.readlines()

    read = False

    for line in lines:

        if "Forces acting on atoms" in line:
            read = True
            continue

        if read:

            if "atom" in line and "force =" in line:

                m = re.search(
                    r"force =\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)",
                    line
                )

                if m:
                    fx, fy, fz = map(float, m.groups())
                    forces.append([fx, fy, fz])

            # stop after nat atoms
            if len(forces) == nat:
                break

    return np.array(forces)

# =========================
# Build force-constant matrix
# =========================
ndof = 3 * moved_atom
Phi = np.zeros((ndof, ndof))

for out in outfiles:
    m = re.search(r"atom(\d+)_([+-])\d*([xyz])", out)
    atom = int(m.group(1)) - 1
    sign = 1 if m.group(2) == "+" else -1
    direction = "xyz".index(m.group(3))

    if atom not in ads_atoms_py:
        continue

    col = 3 * ads_atoms_py.index(atom) + direction
    forces = read_forces(out)

    # Ry/au → N
    forces_si = forces * RY_TO_EV * EV_TO_J / BOHR_TO_M

    for i, ai in enumerate(ads_atoms_py):
        for a in range(3):
            row = 3*i + a
            Phi[row, col] += -sign * forces_si[ai, a] / disp_m

# ==================================
# central finite difference correction
# ==================================
Phi *= 0.5

# =========================
# Mass-weighted Hessian
# =========================
D = np.zeros_like(Phi)
for i in range(moved_atom):
    for j in range(moved_atom):
        for a in range(3):
            for b in range(3):
                D[3*i+a, 3*j+b] = Phi[3*i+a, 3*j+b] / np.sqrt(masses[i]*masses[j])

# =========================
# Diagonalize
# =========================
eigvals, eigvecs = np.linalg.eigh(D)

# ==================================
# Frequencies including imaginary
# ==================================
freq_cm = []

for val in eigvals:

    if val >= 0:
        omega = np.sqrt(val)
        freq = omega / (2*np.pi*C_CM_S)
        freq_cm.append(freq)

    else:
        omega = np.sqrt(abs(val))
        freq = omega / (2*np.pi*C_CM_S)
        freq_cm.append(-freq)   # negative = imaginary

freq_cm = np.array(freq_cm)

# =========================
# Thermodynamic quantities
# =========================
all_omega = []

for val in eigvals:
    
    if val >= 0:
        omega = np.sqrt(val)
        all_omega.append(omega)
        
    else:
        omega = np.sqrt(abs(val))
        all_omega.append(-omega)
        
all_omega = sorted(all_omega)
all_omega = np.array(all_omega)
    
real_omega = []

for val in eigvals:

    if val > 0:

        omega = np.sqrt(val)

        freqs_cm = omega / (2*np.pi*C_CM_S)

        if freqs_cm < freq_cutoff_cm:
            freqs_cm = freq_cutoff_cm
            omega = 2*np.pi*C_CM_S*freqs_cm

        real_omega.append(omega)

real_omega = np.array(real_omega)

# =========================
# Translational energy
# =========================

Etrans_J = (3 / 2) * KB * T

Etrans_eV = Etrans_J / EV_TO_J

# =========================
# Translational entropy
# =========================

trans_part_func = (
    ((total_mass * KB * T) /
     (2 * np.pi * (HBAR ** 2))) ** (3/2)
) * (KB * T / P)

Strans_JK = (np.log(trans_part_func) + (5 / 2)) * KB

Strans_eVK = Strans_JK / EV_TO_J

# =========================
# Rotational energy
# =========================

if molecule:
    if moved_atom == 1:
        Erot_J = 0.0000
        
    elif moved_atom == 2:
        Erot_J = KB * T
        
    else:
        Erot_J = (3 / 2) * KB * T
        
else:
    Erot_J = 0.000
Erot_eV = Erot_J / EV_TO_J


# =========================
# ZPE
# =========================

if molecule:
    if moved_atom == 1:
        ZPE_J = 0.000
        
    elif moved_atom == 2:
        ZPE_J = 0.5 * HBAR * all_omega[5]
        
    else:
        ZPE_J = 0.5 * np.sum(HBAR * all_omega[6:])


else:
    ZPE_J = 0.5 * np.sum(HBAR * real_omega)

ZPE_eV = ZPE_J / EV_TO_J

# =========================
# Vibrational energy correction
# =========================
Evib_J = np.sum(
#    0.5 * HBAR * real_omega +
    (HBAR * real_omega) /
    (np.exp(HBAR * real_omega / (KB * T)) - 1)
)

Evib_eV = Evib_J / EV_TO_J


# =========================
# Translational Gibbs free energy
# =========================
#Gtrans_J = Etrans_J + KB * T - T * Strans_JK

#Gtrans_eV = Gtrans_J / EV_TO_J

# =========================
# Vibrational entropy
# =========================
if molecule:
    if moved_atom == 1:
        Svib_JK = 0.0000
        
    elif moved_atom == 2:
        Svib_JK = KB * (
            (HBAR * all_omega[5] / (KB * T)) /
            (np.exp(HBAR * all_omega[5] / (KB * T)) - 1)
            -
            np.log(1 - np.exp(-HBAR * all_omega[5] / (KB * T)))
        ) 
            
        
    else:
        Svib_JK = np.sum(
            KB * (
                (HBAR * all_omega[6:] / (KB * T)) /
                (np.exp(HBAR * all_omega[6:] / (KB * T)) - 1)
                -
                np.log(1 - np.exp(-HBAR * all_omega[6:] / (KB * T)))
            )
        )

else:
    Svib_JK = np.sum(
        KB * (
            (HBAR * real_omega / (KB * T)) /
            (np.exp(HBAR * real_omega / (KB * T)) - 1)
            -
            np.log(1 - np.exp(-HBAR * real_omega / (KB * T)))
        )
    )


Svib_eVK = Svib_JK / EV_TO_J

# =========================
# Gibbs free energy correction
# =========================
#Gvib_J = ZPE_J + Evib_J - T * Svib_JK
#Gvib_eV = Gvib_J / EV_TO_J

# =========================
# Output
# =========================
print(f"There is/are {moved_atom} atom(s) moved out of {nat} total atom(s)")
if molecule:
    print("\nCaution: the number of moved atom(s) is/are similar to the total atom in the unit cell")
    print("your system may be an isolated atom/molecule, i.e., gas")
  #  print(f"\nI have detected {moved_atom} moved atom(s)")    
    print("\nThere are also contributions from translational and rotational motion")
    print("\nTranslational contribution")
    print(f"E_trans({T} K)      = {Etrans_eV:.6e} eV")
    print(f"S_trans({T} K)      = {Strans_eVK:.6e} eV/K")
   
    
    
  
    if moved_atom == 1:
        print("Your system may be a monoatomic gas")
        print("\nNo rotational contribution")
        
    elif moved_atom == 2:
        
        print("\nRotational contribution")
        print(f"E_rot({T} K)      = {Erot_eV:.6e} eV")
        print("\nCalculations for entropy rotational contribution has not been implemented. Please use molecular program like Gaussian, ORCA, etc.")
        print("\nYour system may be a linear diatomic molecule")
        print("\nIf that is correct, only the 6th vibrational mode will be used for vibrational thermodynamics correction")
        
    else:
        
        print("your system may be a polyatomic molecule")
        print("I assume it is nonlinear")
        print("\nIf that is correct, only the 7th and after vibrational modes will be used for vibrational thermodynamics correction")
        print("\nCalculations for entropy rotational contribution has not been implemented. Please use molecular program like Gaussian, ORCA, etc.")
        
    print("\nOtherwise, use the printed vibrational frequency below for correction")

else:
    print("\nYour system is atom(s) bound to other system, i.e., adsorbed")
    
        


print("\nVibrational frequencies (cm^-1):")

for i, f in enumerate(freq_cm):

    if f >= 0:
        print(f"Mode {i+1:2d}: {f:10.2f} cm^-1")

    else:
        print(f"Mode {i+1:2d}: {abs(f):10.2f} i cm^-1")

print(f"\nZPE               = {ZPE_eV:.6f} eV")
print(f"E_vib({T} K)      = {Evib_eV:.6f} eV")
print(f"S_vib({T} K)      = {Svib_eVK:.6e} eV/K")
print(f"T*S_vib({T} K)    = {T*Svib_eVK:.6f} eV")
#print(f"G_vib({T} K)      = {Gvib_eV:.6f} eV")


#print("Atomic displacement")
#for mode in range(len(eigvals)):
#
#    print(f"\nMode {mode+1}")
#
#    vec = eigvecs[:, mode]
#
#    for i, atom in enumerate(ads_atoms):
#
#        dx = vec[3*i]
#        dy = vec[3*i+1]
#        dz = vec[3*i+2]
#
#        print(
#            f"Atom {atom:3d}: "
#            f"{dx:10.4f} "
#            f"{dy:10.4f} "
#            f"{dz:10.4f}"
#        )

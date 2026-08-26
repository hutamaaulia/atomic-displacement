import copy

# =========================
# User settings
# =========================
qe_input = "NO3-Ag-WS2.scf.in"
disp = 0.02  # Angstrom
atoms_to_move = [54, 55, 56, 57]  

# =========================
# Read QE input
# =========================
with open(qe_input) as f:
    lines = f.readlines()

# Locate ATOMIC_POSITIONS
pos_start = None
for i, line in enumerate(lines):
    if line.strip().startswith("ATOMIC_POSITIONS"):
        pos_start = i + 1
        break

if pos_start is None:
    raise RuntimeError("ATOMIC_POSITIONS not found")

# Read atomic positions
atoms = []
positions = []

i = pos_start
while i < len(lines) and lines[i].strip():
    fields = lines[i].split()
    atoms.append(fields[0])
    positions.append([float(x) for x in fields[1:4]])
    i += 1

# =========================
# Generate displacements
# =========================
directions = {
    "x": [1, 0, 0],
    "y": [0, 1, 0],
    "z": [0, 0, 1]
}

count = 1
for atom_number in atoms_to_move:

    atom = atom_number - 1   # convert to Python index
    for dname, vec in directions.items():
        for sign in [+1, -1]:
            new_pos = copy.deepcopy(positions)
            for k in range(3):
                new_pos[atom][k] += sign * disp * vec[k]

            # Write new QE input
            out = lines[:]
            for j, pos in enumerate(new_pos):
                out[pos_start + j] = (
                    f"{atoms[j]:2s} "
                    f"{pos[0]:15.10f} {pos[1]:15.10f} {pos[2]:15.10f}\n"
                )

            fname = f"disp_{count:03d}_atom{atom_number}_{sign:+d}{dname}.in"
            with open(fname, "w") as f:
                f.writelines(out)

            count += 1

print(f"Generated {count-1} displaced structures.")

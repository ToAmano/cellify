#!/usr/bin/env bash
set -e

echo "=== 3C-SiC (Quantum ESPRESSO) Examples ==="

# 1. Automatically convert primitive 3C-SiC cell (2 atoms) to conventional cubic cell (8 atoms)
# This preserves all QE control parameters and comments, but updates CELL_PARAMETERS and ATOMIC_POSITIONS to cubic.
echo "1. Automatically converting primitive 3C-SiC cell to conventional cubic cell..."
cellify -i 3csic.in -o 3csic_conventional.in --conventional

# 2. Directly generate a 2x2x2 conventional supercell (64 atoms: 32 Si, 32 C)
echo "2. Generating 2x2x2 conventional supercell directly from primitive cell..."
cellify -i 3csic.in -o 3csic_supercell_222.in --conventional --dim 2 2 2

# 3. Create a Silicon-Carbon double vacancy in the 2x2x2 conventional supercell (removes Si:0 and C:33)
# This will result in 62 atoms (31 Si, 31 C) and update nat=62 automatically.
echo "3. Creating a Silicon-Carbon double vacancy (Si at index 0 and its nearest C at index 33) in the supercell..."
cellify -i 3csic_supercell_222.in -o 3csic_vacancy.in --vacancy-index "Si:0" --vacancy-index "C:33"

# 4. Extract relaxed structure from a PWSCF relaxation output log and generate a new SCF input
# This reads the final relaxed structure from 3csic_relaxed.out (62 atoms), merges it with the parameters of 3csic_vacancy.in,
# and writes a new 3csic_scf.in file with calculation = 'scf' and updated nat = 62.
echo "4. Extracting relaxed structure from PWSCF output and creating scf input..."
cellify -i 3csic_relaxed.out --template 3csic_vacancy.in -o 3csic_scf.in --calc scf

echo "3C-SiC examples completed. Output files generated:"
echo "  - 3csic_conventional.in    (8 atoms, cubic cell)"
echo "  - 3csic_supercell_222.in   (64 atoms, conventional 2x2x2 supercell)"
echo "  - 3csic_vacancy.in         (62 atoms, 1 Si-C double vacancy, nat=62)"
echo "  - 3csic_scf.in             (62 atoms, extracted from relax.out using vacancy template, calculation='scf')"

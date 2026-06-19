#!/usr/bin/env bash
set -e

echo "=== Quantum ESPRESSO Examples ==="

# 1. Generate 2x2x2 supercell (updates nat and CELL_PARAMETERS/ATOMIC_POSITIONS)
echo "1. Generating 2x2x2 supercell while preserving input parameters..."
cellify -i qe.in -o qe_222.in --dim 2 2 2

# 2. n-type Doping (replaces Si at index 0 with P, which increments ntyp and appends ATOMIC_SPECIES)
echo "2. Doping Si supercell with P..."
cellify -i qe_222.in -o qe_doped.in --substitute "Si:P:0"

echo "Quantum ESPRESSO examples completed. Output files generated:"
echo "  - qe_222.in   (64 atoms, nat=64, preserves comments & control block)"
echo "  - qe_doped.in (63 Si, 1 P atoms, nat=64, ntyp=2, adds P atomic species)"

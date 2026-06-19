#!/usr/bin/env bash
set -e

echo "=== VASP Examples ==="

# 1. Simple 2x2x2 supercell
echo "1. Generating 2x2x2 supercell..."
cellify -i POSCAR -o POSCAR_222 --dim 2 2 2

# 2. n-type Doping (Replace 1 Si with P at index 0 in the 2x2x2 supercell)
echo "2. Doping 2x2x2 supercell with P at index 0..."
cellify -i POSCAR_222 -o POSCAR_doped --substitute "Si:P:0"

# 3. Minimum Periodic Image Distance Scaling
echo "3. Auto-scaling supercell to ensure minimum image distance >= 12.0 Å..."
cellify -i POSCAR -o POSCAR_scaled --min-dist 12.0

echo "VASP examples completed. Output files generated:"
echo "  - POSCAR_222   (64 Si atoms)"
echo "  - POSCAR_doped (63 Si, 1 P atoms)"
echo "  - POSCAR_scaled (scaled to keep defects >= 12.0 Å apart)"

#!/usr/bin/env bash
set -e

echo "=== Primitive to Conventional Cell Conversion ==="

# 1. Convert Silicon primitive cell (2 atoms) to conventional cubic cell (8 atoms)
echo "1. Automatically converting primitive cell to conventional cubic cell..."
cellify -i POSCAR_primitive -o POSCAR_conventional --conventional

# 2. Directly generate a 2x2x2 conventional supercell (64 atoms) from the primitive cell
echo "2. Generating 2x2x2 conventional supercell directly from primitive cell..."
cellify -i POSCAR_primitive -o POSCAR_conventional_222 --conventional --dim 2 2 2

echo "Conversion completed. Output files generated:"
echo "  - POSCAR_conventional     (8 Si atoms, cubic lattice)"
echo "  - POSCAR_conventional_222 (64 Si atoms, cubic lattice)"

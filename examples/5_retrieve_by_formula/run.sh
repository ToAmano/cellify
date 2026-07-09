#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

echo "=== Example 5: Retrieve Crystal Structures by Chemical Formula ==="

echo "1. Querying and downloading Silicon (Si) stable phase:"
cellify -i Si --select 1 -o Si_stable.cif

echo -e "\n2. Querying and downloading Titanium Dioxide (TiO2) stable phase:"
cellify -i TiO2 --select 1 -o TiO2_stable.cif

echo -e "\n3. Querying and downloading Trihydrogen Sulfide (H3S):"
cellify -i H3S --select 1 -o H3S_stable.cif

echo -e "\n4. Querying, downloading, and generating a 2x2x2 Supercell of Silicon (Si):"
cellify -i Si --select 1 -s 2 2 2 -o Si_supercell.cif

echo -e "\n=== Example Completed Successfully ==="

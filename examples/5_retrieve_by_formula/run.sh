#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

echo "=== Example 5: Retrieve Crystal Structures by Chemical Formula ==="

echo "1. Querying Silicon (Si):"
cellify -i Si

echo -e "\n2. Querying Titanium Dioxide (TiO2):"
cellify -i TiO2

echo -e "\n3. Querying Trihydrogen Sulfide (H3S):"
cellify -i H3S

echo -e "\n=== Example Completed Successfully ==="

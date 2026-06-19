# cellify Examples

This directory contains practical examples to demonstrate how to use `cellify` for common density functional theory (DFT) simulation prep workflows.

## Table of Contents
1. [VASP Supercell & Doping (POSCAR)](./1_vasp_supercell)
2. [Quantum ESPRESSO Parameter Preservation (qe.in)](./2_espresso_doping)
3. [Primitive to Conventional Cell Conversion (POSCAR_primitive)](./3_primitive_to_conventional)
4. [3C-SiC Primitive-to-Conventional Supercell (3csic.in)](./4_espresso_supercell)

---

## How to Run the Examples

First, ensure `cellify` is installed in your python environment (e.g., using `pip install -e .` from the root directory).

Then, you can navigate to any of the subdirectories and run the provided shell scripts (`run.sh`):

```bash
# Example 1: VASP
cd 1_vasp_supercell
bash run.sh

# Example 2: Quantum ESPRESSO
cd ../2_espresso_doping
bash run.sh

# Example 3: Primitive to Conventional Conversion
cd ../3_primitive_to_conventional
bash run.sh

# Example 4: 3C-SiC QE Conversion & Supercell
cd ../4_espresso_supercell
bash run.sh
```

Inspect the output files to see how the coordinates and calculation parameters (like `nat` and `ntyp` for Quantum ESPRESSO) are generated and preserved.

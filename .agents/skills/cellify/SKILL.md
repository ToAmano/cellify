---
name: cellify
description: Use the cellify CLI to generate supercells, conventional cells, vacancies, defects, and surface slabs from crystal structure files (POSCAR, CIF, Quantum ESPRESSO in).
---

# Cellify Agentic Skill Guide

This skill guides you in using `cellify` to manipulate crystal structure files (VASP POSCAR, CIF, and Quantum ESPRESSO input formats) for DFT calculations.

## 1. When to Use This Skill
Use this skill when you need to perform the following structure-building tasks:
- **Conventionalization**: Convert primitive/arbitrary unit cells to standard conventional cells (`--conventional`).
- **Supercell Generation**: Build larger periodic structures (`--supercell` or `-s`) from a unit cell.
- **Vacancy / Defect Builder**: Remove (`--vacancy` or `-v`) or substitute (`--doping` or `-d`) atoms to model defects or dopants.
- **Slab Cutting**: Create 2D surface slabs (`--slab`) with custom Miller indices, thickness, and vacuum space.
- **Index Mapping Inspection**: Show mapping of absolute atomic indices to element types and coordinates (`--show-indices`).

---

## 2. Command Line Interface Reference

The CLI syntax is:
```bash
cellify -i <INPUT_FILE> [OPTIONS] -o <OUTPUT_FILE>
```

### Main Options
- `-i, --input PATH`: Path to the input structure file (formats: VASP `POSCAR`/`CONTCAR`, `*.cif`, or Quantum ESPRESSO `*.in`).
- `-o, --output PATH`: Path to write the output structure file (format auto-detected by extension/filename).
- `-w, --view`: Opens the 3D WebGL viewer in your default browser.
- `--show-indices`: Dumps a neat table of absolute 0-based atomic indices, elements, fractional, and Cartesian coordinates to stdout, then exits.

### Transformation Options
- `--conventional`: Converts the input structure to its standard conventional cell representation *before* applying supercell scaling or defects.
- `-s, --supercell DIM`: Generates a supercell. `DIM` can be:
  - Three integers (e.g. `2 2 2`) for simple diagonal scaling.
  - A matrix string of 9 comma-separated integers (e.g. `2,0,0,0,2,0,0,0,2` representing the scaling matrix rows).
  - *If omitted*, `cellify` automatically scales the supercell based on a target minimum periodic distance (default is 10.0 Å).
- `-v, --vacancy RULES`: Removes atoms. Rules can be:
  - Element and count: `<Element>:<Count>` (e.g., `Si:2` - deletes the first 2 Si atoms).
  - Element and absolute indices: `<Element>:<index1>,<index2>,...` (e.g., `Si:0,4` - deletes Si atoms at absolute indices 0 and 4).
- `-d, --doping RULES`: Substitutes atoms. Rules can be:
  - Element, target element, and count: `<Src>:<Dst>:<Count>` (e.g., `Si:Ge:2` - replaces first 2 Si with Ge).
  - Element, target element, and absolute indices: `<Src>:<Dst>:<index1>,...` (e.g., `Si:Ge:0,4` - replaces Si at indices 0 and 4 with Ge).
- `--slab MILLER_THICKNESS_VACUUM`: Cuts a surface slab. Spec: `h,k,l,thickness,vacuum` (e.g., `1,1,1,3,15` cuts a (1,1,1) slab with a thickness of 3 layers/atomic layers and a 15 Å vacuum region).

---

## 3. Critical Workflow Conventions (LLM Best Practices)

### ① Index Validation (`--show-indices` First)
Before applying vacancies (`-v`) or doping (`-d`) by index:
1. Run the scaling command first (or dry-run) and use the `--show-indices` flag.
2. Inspect stdout to find the exact 0-based indices of the atoms you want to target.
3. Apply the defect using the retrieved indices.
*Never guess index values.* Due to SPGLib and Pymatgen cell conventionalization, atoms are sorted by element and coordinates, which might change their absolute index ordering.

### ② Conventionalization Ordering
Always apply `--conventional` if you are cutting a slab (`--slab`) or inserting pinpoint defects. Conventionalizing first ensures that:
1. The Miller indices (`h,k,l`) align with standard crystallographic orientations.
2. The coordinate representation is intuitive for indexing defects.

### ③ Quantum ESPRESSO (`text_replace` Engine)
When working with Quantum ESPRESSO `*.in` files:
- `cellify` uses a plain-text substitution engine. It preserves comments, namelists, and formatting exactly.
- It automatically updates `nat` (number of atoms) and `ntyp` (number of element types) inside the `&SYSTEM` namelist.
- You do *not* need to manually edit these parameters when using `cellify` output.

---

## 4. Practical Examples

### Example 1: Create a 2x2x2 Supercell of Silicon
```bash
cellify -i POSCAR -s 2 2 2 -o POSCAR_222
```

### Example 2: Inspect Indices of a Slab Model
```bash
cellify -i qe.in --conventional --slab 1,1,1,3,15 --show-indices
```

### Example 3: Create a Doped Divacancy Model in 3C-SiC
```bash
# 1. Convert to conventional cell, scale, and print indices to locate C atoms
cellify -i POSCAR --conventional -s 2 2 2 --show-indices

# 2. Re-run to delete Si at index 0 and C at index 32
cellify -i POSCAR --conventional -s 2 2 2 -v Si:0 -v C:32 -o POSCAR_divacancy
```

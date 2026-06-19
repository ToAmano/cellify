# cellify

A user-friendly command-line interface (CLI) tool to quickly, intuitively, and advancedly generate supercells and calculation-ready inputs from unit cells in density functional theory (DFT) calculation workflows (VASP, Quantum ESPRESSO, OpenMX, CP2K, etc.).

---

## 1. Target Users and Pain Points

### Target Users
*   Researchers in materials science, physics, and chemistry simulating crystals, interfaces, surfaces, defects, and amorphous structures using DFT.

### Current Pain Points (Limitations of Existing Tools)
1.  **"ASE and Pymatgen are powerful, but writing Python scripts is tedious"**
    *   Writing scripts with `read`, `make_supercell`, and `write` just to create a quick supercell is annoying.
2.  **"cif2cell and other tools are prone to broken installations"**
    *   Older python dependencies or compilation issues often cause setup problems.
3.  **"Specifying non-diagonal transformation matrices (orthogonalization, etc.) is unintuitive"**
    *   Quickly redefining lattices or cutting specific orientations from a terminal is difficult.
4.  **"Calculating sizes to avoid periodic boundary interferences is tedious"**
    *   Manually finding the smallest cell configuration to keep defect-to-defect distances above a threshold (e.g., $15\ \text{Å}$) is time-consuming.
5.  **"Creating surface slab models and inserting vacuum layers in separate tools is prone to errors"**

---

## 2. Requirements & Features

### ① Format-Free Multi-Format Conversion
*   Automatically determines file formats from file extensions or headers.
*   **Supported Formats**:
    *   VASP (`POSCAR`, `CONTCAR`)
    *   Quantum ESPRESSO (`.in`, `.txt`, `.qe`)
    *   Crystallographic Information File (`.cif`)
    *   XCrysDen Structure Format (`.xsf`, `.axsf`)
    *   XYZ format (`.xyz`)
    *   FHI-aims (`geometry.in`)

### ② Flexible Cell Expansion (Supercell Generation)
*   **Diagonal Scaling**: Simplest integer multiplication along lattice axes (e.g., `2 2 2`).
*   **Matrix-Based Redefinition**: Redefine lattices using an arbitrary $3 \times 3$ transformation matrix. Ideal for orthogonalizing hexagonal cells or extracting specific crystal orientations.
*   **Minimum Distance (Cutoff) Automatic Scaling**:
    *   Automatically calculates and generates the smallest diagonal supercell (or specific axis dimensions) that guarantees the distance between periodic images of any atom is $\ge d\ \text{Å}$. Extremely useful for defect and phonon calculations.

### ③ Easy Defect & Doping Modeling
*   **Substitutions**: Replace specific atoms at a given index (e.g., replacing Si at index 0 with P) or randomly replace a specified percentage of atoms (e.g., replacing $5\%$ of Si atoms with Al).
*   **Vacancies**: Remove atoms at specific indices or randomly delete a specified count of a given element.

### ④ Surface Slab Generation
*   Cut a surface slab from bulk structures by specifying Miller indices $(h, k, l)$, slab thickness (in $\text{Å}$ or layers), and vacuum thickness (in $\text{Å}$).

### ⑤ Logging and Metadata Analysis
*   Outputs structure logs to stderr during execution:
    *   Initial volume, atom count, and reduced formula.
    *   Final supercell volume, lattice constants, lattice angles, and atom count.
    *   Applied transformation matrix.
    *   Minimum atomic distance under periodic boundary conditions.

### ⑥ Calculation-Ready Input Generation
*   For formats like Quantum ESPRESSO where calculation parameters and coordinates coexist in a single file, the original parameters (`&CONTROL`, `&SYSTEM`, etc.) and comments are completely preserved.
*   The following parameters are automatically updated to match the generated supercell structure:
    *   **Total number of atoms (`nat`)**: Automatically updated to the supercell atom count.
    *   **Number of atomic types (`ntyp`)**: Dynamically incremented if new elements are added via doping.
    *   **Atomic species definitions (`ATOMIC_SPECIES`)**: Automatically appends definitions (mass, pseudopotentials) for newly introduced elements.

---

## 3. Installation

You can install `cellify` from the local repository directory:

```bash
# Clone the repository
git clone https://github.com/ToAmano/cellify.git
cd cellify

# Install in editable mode for development
pip install -e .

# Or install with test dependencies
pip install -e ".[test]"
```

After installation, the `cellify` command will be registered and executable from anywhere in your shell environment.

---

## 4. CLI Design

### Command-Line Arguments

```bash
cellify -i <input_file> -o <output_file> [options]
```

#### Arguments List
*   `-i`, `--input` : Input structure file path (Required).
*   `-o`, `--output` : Output structure file path (Default: `<input_base>_supercell.<ext>`).
*   `-d`, `--dim` : Diagonal scaling factors. 3 integers separated by spaces (e.g., `--dim 2 2 2`).
*   `-m`, `--matrix` : $3 \times 3$ transformation matrix. Specify row values separated by spaces, rows separated by slashes/commas/semicolons (e.g., `--matrix "1 -1 0 / 1 1 0 / 0 0 2"`).
*   `--min-dist` : Automatically generate a supercell with minimum periodic image distance $\ge$ specified distance (in $\text{Å}$).
*   `--substitute` : Substitution rule. Format: `element:target_element:index_or_percentage` (e.g., `--substitute "Si:P:0"` or `--substitute "Si:Al:5%"`).
*   `--vacancy` : Vacancy rule. Format: `element:index_or_count` (e.g., `--vacancy "Si:0"`, `--vacancy "O:2"`_).
*   `--slab` : Miller indices $h\ k\ l$ for surface slab model creation (e.g., `--slab 1 1 1`).
*   `--thick` : Slab thickness in $\text{Å}$ or layers (e.g., `--thick 15.0`).
*   `--vacuum` : Vacuum layer thickness in $\text{Å}$ (e.g., `--vacuum 15.0`).

---

## 5. Use Cases

### 1. Create a simple $2 \times 2 \times 3$ supercell (VASP POSCAR)
```bash
cellify -i POSCAR -o POSCAR_223 --dim 2 2 3
```

### 2. Orthogonalize a hexagonal cell (Quantum ESPRESSO input)
```bash
# Preserves &CONTROL and &SYSTEM settings, and updates nat, CELL_PARAMETERS, and ATOMIC_POSITIONS
cellify -i qe.in -o qe_ortho.in --matrix "1 -1 0 / 1 1 0 / 0 0 1"
```

### 3. Generate the smallest supercell keeping defect distance $\ge 15\ \text{Å}$
```bash
cellify -i POSCAR -o POSCAR_defect_bulk --min-dist 15.0
```

### 4. Create a silicon supercell and replace 1 atom with Phosphorus (n-type doped model)
```bash
cellify -i Si_unit.cif -o Si_doped.POSCAR --dim 3 3 3 --substitute "Si:P:0"
```

### 5. Generate a $\text{SrTiO}_3$ (100) surface slab model with $15\ \text{Å}$ vacuum
```bash
cellify -i STO_bulk.cif -o STO_100_slab.POSCAR --slab 1 0 0 --thick 12.0 --vacuum 15.0
```

---

## 6. Directory Structure

This project uses the standard Python `src-layout`:

```text
cellify/
├── README.md
├── NAMES.md
├── pyproject.toml
└── src/
    └── cellify/
        ├── __init__.py
        ├── cli.py            # Command-line argument parsing and execution flow
        ├── core.py           # Pure geometric modeling (supercell, defect, slab creation)
        └── adapters/         # Software-specific file I/O and parameter-preservation adapters
            ├── __init__.py
            ├── base.py       # Abstract base class for I/O adapters
            ├── espresso.py   # Quantum ESPRESSO adapter
            └── standard.py   # VASP/CIF generic format adapter
```

---

## 7. Technical Stack & Development Approach

1.  **Language**: **Python 3** (High affinity with scientific and DFT software ecosystems).
2.  **Core Libraries**: **pymatgen** and **ASE (Atomic Simulation Environment)**.
    *   **pymatgen**: Used for symmetry determination, structure analysis, defect modulations, and advanced slab generations.
    *   **ASE**: Used for format-free structure loading/writing and robust file parsed operations.
    *   Conversion between both frameworks is done seamlessly via `pymatgen.io.ase.AseAtomsAdaptor`.
3.  **Packaging**:
    *   Managed via `pyproject.toml` using `hatchling` as the build backend.
    *   Installable in editable mode using `pip install -e ".[test]"`.
    *   Registers `cellify` command as an entry point upon installation.

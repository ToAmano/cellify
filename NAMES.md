# Name Candidates for the DFT Supercell Generator Tool

A list of command-line name candidates representing the tool's identity (ease of use, parameter preservation, automatic expansion based on periodic distances, etc.), which are easy to type and memorable.

---

## Proposed Name Candidates

| Name | Pronunciation | Command Example | Description & Origin |
| :--- | :--- | :--- | :--- |
| **`scel`** | s-cell | `scel -i POSCAR -o SPOSCAR -d 2 2 2` | Short for **S**uper**cel**l. Extremely short and optimal for CLI typing. Least prone to command conflicts. |
| **`cellify`** | cell-i-fy | `cellify -i qe.in --min-dist 15.0` | **【Selected】** A verb meaning "to convert into a cell" or "to supercell". Sounds modern, friendly, and easy to remember. |
| **`dftsc`** | d-f-t-s-c | `dftsc -i POSCAR -o SPOSCAR` | Short for **D**FT **S**uper**c**ell. Immediately conveys its DFT-specific usage and contains no fluff. |
| **`readycell`** | ready-cell | `readycell -i qe.in -o qe_super.in` | Highlights the tool's core feature of creating **Calculation-ready** inputs (preserving calculations and updating headers automatically). |
| **`scgen`** | s-c-gen | `scgen -i input.cif -d 2 2 2` | Short for **S**uper**c**ell **Gen**erator. An orthodox name that clearly defines its purpose. |
| **`autosuper`** | auto-super | `autosuper -i POSCAR --min-dist 12.0` | Focuses on the smart feature of **automatically (Auto)** generating the optimal supercell based on a specified distance constraint. |

---

## Selection Criteria

1.  **CLI Ergonomics**:
    *   `scel` and `dftsc` are only 4–5 characters, making them highly efficient to type.
    *   `cellify` and `scgen` are easy to type without spelling errors.
2.  **PyPI Package Availability**:
    *   Ensures that when publishing the tool to PyPI, it does not conflict with existing packages.
    *   `cellify` is clear, clean, and stands out as a unique and descriptive project name.

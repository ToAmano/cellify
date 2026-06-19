# GEMINI.md - Development History and Code Conventions

This document records the development history, technical decisions, and code quality conventions of the `cellify` project, built collaboratively with the AI assistant "Gemini".

---

## 1. Project Background and Decisions

### Purpose of Development
To provide a command-line (CLI) tool that allows researchers to "easily" generate supercells and modify structures (doping, vacancies, surface slabs) for DFT calculations (VASP, Quantum ESPRESSO, etc.) while "completely preserving calculation parameters" in mixed format files.

### Decision Process
1.  **Survey of Existing Tools**: Analyzed the pros and cons of `cif2cell`, `phonopy`, `ASE`, and `pymatgen` (issues like brittle dependencies, manual header copy-pasting, etc.).
2.  **Requirements Definition (README.md)**: Defined specifications including "Calculation-ready input files (especially QE `nat` auto-updates)" and "Automatic scaling based on minimum periodic distances".
3.  **Naming Decision**: Chose **`cellify`** over other candidates like `scel` to support future general structure-modeling extensions.
4.  **GitHub Connection**: Created and connected a public repository: `ToAmano/cellify`.
5.  **Conventional Cell Conversion (Issue #5)**: Added `--conventional` flag to automatically convert primitive cells (commonly downloaded from materials databases) into conventional cells before performing other operations. This ensures intuitive defect index modeling and surface slab cutting.

---

## 2. Development & Coding Conventions

### ① Complete Type Hinting 【Strictly Required】
All Python modules under `src/cellify/` must have **complete type annotations** to ensure maintainability, static verification, and readability.

*   **Function Signatures**: Specify parameter and return types using the `typing` module or built-in generics.
    ```python
    from typing import Tuple, Dict, Any
    from pymatgen.core import Structure

    def load_structure_file(filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        ...
    ```
*   **Local Variables**: Annotate local variables when they hold complex types.
*   **Static Analysis**: Ensure all changes are validated and raise no errors when checked with static type checkers like `mypy`.

### ② Parameter Preservation (Plain-Text Substitution)
To prevent parse errors or accidental default overwrites (such as `K_POINTS None` introduction) in mixed-format calculations like Quantum ESPRESSO, adapters must keep calculation parameters (namelists, comments, formatting) in plain text and only modify parameters like `nat` and `ntyp` using regex matching. This is referred to as the `text_replace` engine.

### ③ Modular Architecture (Separation of I/O Adapters and Core)
To keep the codebase maintainable and open for future DFT software additions:
*   **`adapters` Package**: Any format-specific file reading/writing (and parameter-preserving text substitutions) must be encapsulated under adapters inheriting the abstract `BaseAdapter` class (e.g., `EspressoAdapter`, `StandardAdapter`).
*   **`core` Module**: Structural manipulation logic (supercells, defects, slabs) must be written as modular, pure-like functions that receive a `Structure` and return a modified `Structure` without being aware of input/output files.
*   **Future Modulations**: The core logic remains in a single `core.py` for simplicity now, but is kept decoupled so it can be easily split into a `core/` package when the number of features increases.

### ④ Branch Management and PR Lifecycle
To maintain clean git history and keep features isolated:
*   **One Branch Per Task/Feature**: Always create a new, dedicated branch for each feature or task. Do not implement unrelated requests or new tasks on an existing, open feature branch.
*   **Cleanup After Merge**: Once a branch's Pull Request is merged into `main`, delete the branch locally and remotely to prevent reuse.

---

## 3. Testing & Linting Conventions and Execution
To guarantee code quality and stability, unit/integration tests are managed using `pytest` inside the [tests/](file:///Users/amano/works/research/supercell/tests/) directory, and static analysis/formatting are managed by `pre-commit`.

### ① Execution Before Commit 【Strictly Required】
You **must** run and pass both the test suite and pre-commit checks locally before executing any `git commit`:

1.  **Run Pytest**:
    ```bash
    # Install dependencies including test tools
    pip install -e ".[test]"

    # Run tests
    pytest
    ```
2.  **Run Pre-commit**:
    ```bash
    # Install pre-commit hook (first time only to enable commit-time validation)
    pre-commit install

    # Run checks manually on all files before committing
    pre-commit run --all-files
    ```

### ② Conventions
*   Whenever implementing new features or fixing bugs, corresponding tests **must** be added to the `tests/` directory.
*   Always maintain a fully passing (`passed`) state for all tests and linters in local verification and CI. Do not commit or push failing code.

---

## 4. Future Roadmap and Milestones
*   Further stabilization and validation of vacancy generation (`--vacancy`) logic.
*   Adding test data and verification for surface slab models (`--slab`).
*   Integrating a lattice mismatch auto-relaxation module for heterostructures.
*   Setting up linting and formatting configuration files (`mypy`, `flake8`, `black`, `ruff`).

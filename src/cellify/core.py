"""
Core modeling logic for cellify.
Handles structure loading, supercell generation, substitutions,
vacancies, slab generation, and file saving using pymatgen and ASE.
"""

import contextlib
import io
import math
import re
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from cellify.adapters import BaseAdapter, get_adapter


def load_structure_file(filepath: str) -> Tuple[Structure, Dict[str, Any]]:
    """
    Loads a file and returns the structure object along with metadata.
    """
    adapter: BaseAdapter = get_adapter(filepath)
    return adapter.read(filepath)


def save_structure_file(
    filepath: str, structure: Structure, meta_data: Dict[str, Any]
) -> None:
    """
    Saves the structure to a file.
    """
    adapter: BaseAdapter = get_adapter(filepath)
    adapter.write(filepath, structure, meta_data)


def convert_to_conventional(structure: Structure) -> Structure:
    """
    Finds and returns the standard conventional cell of the structure.
    """
    sga = SpacegroupAnalyzer(structure)
    return sga.get_conventional_standard_structure()


def parse_matrix_string(matrix_str: str) -> np.ndarray:
    """
    Parses a matrix string like "1 -1 0 / 1 1 0 / 0 0 1" into a 3x3 numpy array.
    """
    # Split rows by slash, comma, or semicolon
    rows_raw: List[str] = re.split(r"[/,;]", matrix_str)
    if len(rows_raw) != 3:
        raise ValueError(
            "Matrix string must define exactly 3 rows (separated by /, , or ;)"
        )

    matrix: List[List[float]] = []
    for r in rows_raw:
        vals: List[float] = [float(x) for x in r.strip().split()]
        if len(vals) != 3:
            raise ValueError("Each row in the matrix must have exactly 3 elements")
        matrix.append(vals)

    return np.array(matrix)


def calculate_min_dist_scaling(
    structure: Structure, min_dist: float
) -> Tuple[int, int, int]:
    """
    Calculates the minimum diagonal scaling factors (nx, ny, nz) so that
    the perpendicular distance (plane-to-plane distance) along all lattice vectors
    is at least min_dist under periodic boundary conditions.
    """
    lattice = structure.lattice
    matrix = lattice.matrix
    a_vec, b_vec, c_vec = matrix[0], matrix[1], matrix[2]

    vol: float = lattice.volume

    # Perpendicular distance along each lattice vector (plane-to-plane distance d_i)
    # d_a = V / |b x c|
    # d_b = V / |c x a|
    # d_c = V / |a x b|
    d_a: float = vol / np.linalg.norm(np.cross(b_vec, c_vec))
    d_b: float = vol / np.linalg.norm(np.cross(c_vec, a_vec))
    d_c: float = vol / np.linalg.norm(np.cross(a_vec, b_vec))

    # Calculate required scaling factors
    nx: int = int(math.ceil(min_dist / d_a))
    ny: int = int(math.ceil(min_dist / d_b))
    nz: int = int(math.ceil(min_dist / d_c))

    return max(1, nx), max(1, ny), max(1, nz)


def apply_substitutions(structure: Structure, substitute_rules: List[str]) -> None:
    """
    Applies substitution rules to the structure.
    Rule formats:
        "Si:P:0" (replaces Si at absolute index 0 with P)
        "Si:Al:5%" (randomly replaces 5% of Si atoms with Al)
    """
    for rule in substitute_rules:
        _apply_single_substitution(structure, rule)


def _apply_single_substitution(structure: Structure, rule: str) -> None:
    """
    Applies a single substitution rule to the structure.
    """
    parts: List[str] = rule.split(":")
    if len(parts) != 3:
        raise ValueError(
            f"Invalid substitution rule: {rule}. Must be 'element:target_element:index_or_percentage'"
        )

    src_el, dest_el, target = parts[0], parts[1], parts[2]

    matching_indices: List[int] = [
        i for i, site in enumerate(structure) if site.specie.symbol == src_el
    ]
    if not matching_indices:
        print(f"Warning: No matching elements found for substitution source '{src_el}'")
        return

    if target.endswith("%"):
        _substitute_percentage(structure, src_el, dest_el, target, matching_indices)
    else:
        _substitute_index(structure, src_el, dest_el, target)


def _substitute_percentage(
    structure: Structure,
    src_el: str,
    dest_el: str,
    target: str,
    matching_indices: List[int],
) -> None:
    """
    Helper to apply substitution by percentage.
    """
    percentage: float = float(target[:-1]) / 100.0
    num_to_replace: int = int(round(len(matching_indices) * percentage))
    if num_to_replace == 0 and percentage > 0:
        num_to_replace = 1

    replace_indices = np.random.choice(matching_indices, num_to_replace, replace=False)
    for replace_idx in replace_indices:
        structure.replace(replace_idx, dest_el)
    print(f"Replaced {num_to_replace} of {src_el} with {dest_el} ({target})")


def _substitute_index(
    structure: Structure, src_el: str, dest_el: str, target: str
) -> None:
    """
    Helper to apply substitution by absolute index.
    """
    try:
        idx: int = int(target)
        if idx < 0 or idx >= len(structure):
            raise IndexError(
                f"Index {idx} out of range (structure size: {len(structure)})"
            )

        actual_symbol: str = structure[idx].specie.symbol
        if actual_symbol != src_el:
            print(
                f"Warning: Site index {idx} is '{actual_symbol}', not source element '{src_el}'. Replacing anyway."
            )

        structure.replace(idx, dest_el)
        print(f"Replaced site {idx} ({actual_symbol}) with {dest_el}")
    except ValueError as exc:
        raise ValueError(
            f"Invalid substitution target index or percentage: {target}"
        ) from exc


def apply_vacancies_by_index(structure: Structure, vacancy_rules: List[str]) -> None:
    """
    Applies vacancy rules to the structure by removing atoms at specific absolute indices.
    Rule format:
        "Si:0" (deletes Si atom at index 0)
    """
    indices_to_remove: List[int] = []

    for rule in vacancy_rules:
        parts: List[str] = rule.split(":")
        if len(parts) != 2:
            raise ValueError(
                f"Invalid vacancy index rule: {rule}. Must be 'element:index'"
            )

        src_el: str = parts[0]
        try:
            val: int = int(parts[1])
        except ValueError as exc:
            raise ValueError(
                f"Invalid vacancy index: {parts[1]} in rule {rule}"
            ) from exc

        if val < 0 or val >= len(structure):
            raise IndexError(
                f"Index {val} out of range (structure size: {len(structure)})"
            )

        actual_symbol: str = structure[val].specie.symbol
        if actual_symbol != src_el:
            print(
                f"Warning: Site index {val} is '{actual_symbol}', not vacancy element '{src_el}'. Removing anyway."
            )

        indices_to_remove.append(val)
        print(f"Removed site {val} ({actual_symbol}) to create vacancy")

    if indices_to_remove:
        # Sort indices in descending order to avoid shift errors when removing sites
        indices_to_remove = sorted(list(set(indices_to_remove)), reverse=True)
        structure.remove_sites(indices_to_remove)


def apply_vacancies_by_count(structure: Structure, vacancy_rules: List[str]) -> None:
    """
    Applies vacancy rules to the structure by randomly removing a specified count of atoms of a given element.
    Rule format:
        "O:2" (randomly deletes 2 oxygen atoms)
    """
    indices_to_remove: List[int] = []

    for rule in vacancy_rules:
        parts: List[str] = rule.split(":")
        if len(parts) != 2:
            raise ValueError(
                f"Invalid vacancy count rule: {rule}. Must be 'element:count'"
            )

        src_el: str = parts[0]
        try:
            val: int = int(parts[1])
        except ValueError as exc:
            raise ValueError(
                f"Invalid vacancy count: {parts[1]} in rule {rule}"
            ) from exc

        if val < 0:
            raise ValueError(f"Vacancy count cannot be negative: {val} in rule {rule}")

        matching_indices: List[int] = [
            i for i, site in enumerate(structure) if site.specie.symbol == src_el
        ]
        if not matching_indices:
            print(f"Warning: No matching elements found for vacancy source '{src_el}'")
            continue

        if val > len(matching_indices):
            raise ValueError(
                f"Requested vacancy count {val} exceeds available {src_el} atoms ({len(matching_indices)})"
            )

        remove_subset = np.random.choice(matching_indices, val, replace=False)
        indices_to_remove.extend(remove_subset)
        print(f"Created {val} vacancies of {src_el} (randomly selected)")

    if indices_to_remove:
        # Sort indices in descending order to avoid shift errors when removing sites
        indices_to_remove = sorted(list(set(indices_to_remove)), reverse=True)
        structure.remove_sites(indices_to_remove)


def generate_surface_slab(
    structure: Structure,
    miller_index: List[int],
    thick: Optional[float],
    vacuum: Optional[float],
) -> Structure:
    """
    Generates a surface slab model using pymatgen's SlabGenerator.
    """
    if all(h == 0 for h in miller_index):
        raise ValueError("Miller indices cannot all be zero.")

    slab_thick: float = thick if thick else 10.0
    vac_thick: float = vacuum if vacuum else 15.0

    gen = SlabGenerator(
        initial_structure=structure,
        miller_index=miller_index,
        min_slab_size=slab_thick,
        min_vacuum_size=vac_thick,
        center_slab=True,
    )

    slabs = gen.get_slabs()
    if not slabs:
        raise ValueError(f"Could not generate slab for Miller index {miller_index}")

    # Adopt the first generated slab model (often the most symmetric and stable one)
    slab = slabs[0]
    return slab


def get_structure_summary(structure: Structure, label: str = "") -> str:
    """
    Returns a formatted summary string of the structure.
    """
    lines: List[str] = []
    if label:
        lines.append(label)
    lines.append(f"  Formula: {structure.composition.reduced_formula}")
    lines.append(f"  Volume:  {structure.volume:.3f} A^3")
    lines.append(f"  Number of atoms: {len(structure)}")
    if label:
        lines.append("  Lattice constants:")
        lines.append(
            f"    a = {structure.lattice.a:.4f} A, b = {structure.lattice.b:.4f} A, c = {structure.lattice.c:.4f} A"
        )
        lines.append(
            f"    alpha = {structure.lattice.alpha:.2f} deg, beta = {structure.lattice.beta:.2f} deg, gamma = {structure.lattice.gamma:.2f} deg"
        )
    return "\n".join(lines)


def get_atomic_indices_table(structure: Structure) -> str:
    """
    Returns a formatted table string of all atomic indices and coordinates.
    """
    lines: List[str] = []
    lines.append("\nAbsolute Atomic Indices & Coordinates:")
    lines.append("-" * 78)
    header: str = (
        f"{'Index':<6} {'Element':<8} "
        f"{'Fractional Coordinates (a, b, c)':<36} "
        f"{'Cartesian (x, y, z)':<22}"
    )
    lines.append(header)
    lines.append("-" * 78)
    for idx, site in enumerate(structure):
        frac: str = (
            f"[{site.frac_coords[0]:.4f}, "
            f"{site.frac_coords[1]:.4f}, "
            f"{site.frac_coords[2]:.4f}]"
        )
        cart: str = (
            f"[{site.coords[0]:.3f}, "
            f"{site.coords[1]:.3f}, "
            f"{site.coords[2]:.3f}]"
        )
        lines.append(f"{idx:<6} {site.species_string:<8} {frac:<36} {cart:<22}")
    lines.append("-" * 78)
    lines.append(f"Total: {len(structure)} atoms ({structure.composition.formula})")
    return "\n".join(lines)


def run_cellify_pipeline(  # noqa: C901,CCR001 # pylint: disable=too-many-arguments,too-many-positional-arguments
    structure: Structure,
    conventional: bool = False,
    dim: Optional[Any] = None,
    min_dist: Optional[float] = None,
    substitute: Optional[List[str]] = None,
    vacancy_index: Optional[List[str]] = None,
    vacancy_count: Optional[List[str]] = None,
    slab: Optional[Any] = None,
    thick: Optional[float] = None,
    vacuum: Optional[float] = None,
) -> Tuple[Structure, str]:
    """
    Runs the entire modeling pipeline on the structure, capturing all console
    outputs and returning the final structure along with the captured log string.
    """
    log_stream: io.StringIO = io.StringIO()
    with contextlib.redirect_stdout(log_stream):
        # 1. Conventional cell conversion
        if conventional:
            print("Converting structure to standard conventional cell...")
            structure = convert_to_conventional(structure)

        # 2. Supercell generation
        if dim or min_dist is not None:
            if dim:
                if isinstance(dim, list):
                    print(f"Generating supercell with diagonal scaling: {dim}")
                    structure.make_supercell(dim)
                else:
                    clean_dim: str = str(dim).strip()
                    if "," in clean_dim or "/" in clean_dim or ";" in clean_dim:
                        matrix: np.ndarray = parse_matrix_string(clean_dim)
                        print(f"Generating supercell with matrix:\n{matrix}")
                        structure.make_supercell(matrix)
                    else:
                        factors: List[int] = [
                            int(x) for x in re.split(r"[\s,]+", clean_dim) if x
                        ]
                        if len(factors) != 3:
                            raise ValueError(
                                f"Diagonal scaling 'dim' must contain exactly 3 integers, got {factors}"
                            )
                        print(f"Generating supercell with diagonal scaling: {factors}")
                        structure.make_supercell(factors)
            elif min_dist is not None:
                nx: int
                ny: int
                nz: int
                nx, ny, nz = calculate_min_dist_scaling(structure, min_dist)
                print(
                    f"Calculated scaling for minimum distance >= {min_dist} A: [{nx}, {ny}, {nz}]"
                )
                structure.make_supercell([nx, ny, nz])

        # 3. Substitutions
        if substitute:
            apply_substitutions(structure, substitute)

        # 4. Vacancy index
        if vacancy_index:
            apply_vacancies_by_index(structure, vacancy_index)

        # 5. Vacancy count
        if vacancy_count:
            apply_vacancies_by_count(structure, vacancy_count)

        # 6. Slab generation
        if slab:
            if isinstance(slab, (list, tuple)):
                miller_indices: List[int] = [int(x) for x in slab]
            else:
                clean_slab: str = str(slab).strip()
                miller_indices = [int(x) for x in re.split(r"[\s,;]+", clean_slab) if x]
            if len(miller_indices) != 3:
                raise ValueError(
                    f"Miller indices must contain exactly 3 integers, got {miller_indices}"
                )
            print(f"Generating slab model for Miller indices: {miller_indices}")
            structure = generate_surface_slab(structure, miller_indices, thick, vacuum)

    return structure, log_stream.getvalue()

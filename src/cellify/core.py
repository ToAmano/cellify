"""
Core modeling logic for cellify.
Handles structure loading, supercell generation, substitutions,
vacancies, slab generation, and file saving using pymatgen and ASE.
"""

import contextlib
import io
import math
import os
import re
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from pymatgen.core import Lattice, Structure
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


def determine_output_path(input_path: str, output_path: Optional[str] = None) -> str:
    """
    Determines the output file path.
    """
    if output_path:
        return output_path

    base, ext = os.path.splitext(input_path)
    # Special case: VASP files like POSCAR or CONTCAR with no extension
    if not ext and os.path.basename(base) in ["POSCAR", "CONTCAR"]:
        return f"{base}_supercell"
    return f"{base}_supercell{ext}"


def process_template_and_validation(
    meta_data: Dict[str, Any],
    output_path: Optional[str],
    template_path: Optional[str] = None,
    calc: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Handles template loading, calculation overrides, and QE I/O format validations.
    """
    if template_path:
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template file '{template_path}' not found.")
        _, template_meta = load_structure_file(template_path)
        meta_data = template_meta

    if calc:
        meta_data["calculation"] = calc

    if output_path:
        is_input_qe_output = meta_data.get("mode") == "espresso_out"
        lower_out_path = output_path.lower()
        is_output_qe_input = (
            any(lower_out_path.endswith(ext) for ext in [".in", ".qe", ".pwi"])
            or "qe" in lower_out_path
            or "espresso" in lower_out_path
        )
        if is_input_qe_output and is_output_qe_input and not template_path:
            raise ValueError(
                "A template QE input file must be specified when reading from a QE output log file and writing to a QE input file."
            )

    return meta_data


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
    if not slabs:  # pragma: no cover
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


def apply_supercell(
    structure: Structure,
    dim: Optional[List[int]] = None,
    matrix: Optional[str] = None,
    min_dist: Optional[float] = None,
) -> Structure:
    """
    Applies supercell generation options to the structure.
    """
    if dim:
        print(f"Generating supercell with diagonal scaling: {dim}")
        structure.make_supercell(dim)
    elif matrix:
        try:
            mat: np.ndarray = parse_matrix_string(matrix)
            print(f"Generating supercell with matrix:\n{mat}")
            structure.make_supercell(mat)
        except Exception as e:
            raise e
    elif min_dist:
        nx, ny, nz = calculate_min_dist_scaling(structure, min_dist)
        print(
            f"Calculated scaling for minimum distance >= {min_dist} A: [{nx}, {ny}, {nz}]"
        )
        structure.make_supercell([nx, ny, nz])

    return structure


def scale_structure_volume(structure: Structure, factor: float) -> Structure:
    """
    Scales the structure's volume by a given factor while preserving length proportions and angles.
    """
    if factor <= 0:
        raise ValueError("Scaling factor must be positive.")
    struct_copy: Structure = structure.copy()
    struct_copy.scale_lattice(structure.volume * factor)
    return struct_copy


def scale_structure_lattice(structure: Structure, factor: float) -> Structure:
    """
    Scales the structure's lattice constants (lattice vectors) by a given factor.
    Fractional coordinates are kept unchanged.
    """
    if factor <= 0:
        raise ValueError("Scaling factor must be positive.")
    struct_copy: Structure = structure.copy()
    new_lattice = Lattice(struct_copy.lattice.matrix * factor)
    struct_copy.lattice = new_lattice
    return struct_copy


def scale_structure_axes(
    structure: Structure, fa: float, fb: float, fc: float
) -> Structure:
    """
    Scales the individual lattice vectors (axes) of the structure by factors fa, fb, and fc.
    Fractional coordinates are kept unchanged.
    """
    if fa <= 0 or fb <= 0 or fc <= 0:
        raise ValueError("All scaling factors must be positive.")
    struct_copy: Structure = structure.copy()
    matrix = np.array(struct_copy.lattice.matrix)
    matrix[0] *= fa
    matrix[1] *= fb
    matrix[2] *= fc
    struct_copy.lattice = Lattice(matrix)
    return struct_copy


def apply_defects_and_slab(  # noqa: C901,CCR001 # pylint: disable=too-many-arguments,too-many-positional-arguments
    structure: Structure,
    substitute: Optional[List[str]] = None,
    vacancy_index: Optional[List[str]] = None,
    vacancy_count: Optional[List[str]] = None,
    slab: Optional[List[int]] = None,
    thick: Optional[float] = None,
    vacuum: Optional[float] = None,
) -> Structure:
    """
    Applies substitutions, vacancies, and surface slab options to the structure.
    """
    if substitute:
        apply_substitutions(structure, substitute)

    if vacancy_index:
        apply_vacancies_by_index(structure, vacancy_index)

    if vacancy_count:
        apply_vacancies_by_count(structure, vacancy_count)

    if slab:
        print(f"Generating slab model for Miller indices: {slab}")
        structure = generate_surface_slab(structure, slab, thick, vacuum)

    return structure


def run_cellify_pipeline(  # noqa: C901,CCR001 # pylint: disable=too-many-arguments,too-many-positional-arguments
    structure: Structure,
    conventional: bool = False,
    dim: Optional[List[int]] = None,
    matrix: Optional[str] = None,
    min_dist: Optional[float] = None,
    scale_vol: Optional[float] = None,
    scale_lat: Optional[float] = None,
    scale_axes: Optional[List[float]] = None,
    substitute: Optional[List[str]] = None,
    vacancy_index: Optional[List[str]] = None,
    vacancy_count: Optional[List[str]] = None,
    slab: Optional[List[int]] = None,
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
        structure = apply_supercell(
            structure, dim=dim, matrix=matrix, min_dist=min_dist
        )

        # 2.5. Lattice scaling
        if scale_vol is not None:
            print(f"Scaling structure volume by factor: {scale_vol}")
            structure = scale_structure_volume(structure, scale_vol)
        elif scale_lat is not None:
            print(f"Scaling structure lattice constants by factor: {scale_lat}")
            structure = scale_structure_lattice(structure, scale_lat)
        elif scale_axes is not None:
            fa, fb, fc = scale_axes
            print(f"Scaling structure lattice axes by factors: a={fa}, b={fb}, c={fc}")
            structure = scale_structure_axes(structure, fa, fb, fc)

        # 3. Defects and Slab generation
        structure = apply_defects_and_slab(
            structure,
            substitute=substitute,
            vacancy_index=vacancy_index,
            vacancy_count=vacancy_count,
            slab=slab,
            thick=thick,
            vacuum=vacuum,
        )

    return structure, log_stream.getvalue()

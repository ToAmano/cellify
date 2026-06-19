"""
Core modeling logic for cellify.
Handles structure loading, supercell generation, substitutions, 
vacancies, slab generation, and file saving using pymatgen and ASE.
"""

import math
import re
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.ase import AseAtomsAdaptor

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
            print(
                f"Warning: No matching elements found for substitution source '{src_el}'"
            )
            continue

        if target.endswith("%"):
            # Random replacement based on percentage
            percentage: float = float(target[:-1]) / 100.0
            num_to_replace: int = int(round(len(matching_indices) * percentage))
            if num_to_replace == 0 and percentage > 0:
                num_to_replace = 1

            replace_indices = np.random.choice(
                matching_indices, num_to_replace, replace=False
            )
            for idx in replace_indices:
                structure.replace(idx, dest_el)
            print(f"Replaced {num_to_replace} of {src_el} with {dest_el} ({target})")
        else:
            # Absolute index replacement
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
            except ValueError:
                raise ValueError(
                    f"Invalid substitution target index or percentage: {target}"
                )


def apply_vacancies(structure: Structure, vacancy_rules: List[str]) -> None:
    """
    Applies vacancy rules to the structure (deletes specified atoms).
    Rule formats:
        "Si:0" (deletes Si atom at index 0)
        "O:2" (randomly deletes 2 oxygen atoms)
    """
    indices_to_remove: List[int] = []

    for rule in vacancy_rules:
        parts: List[str] = rule.split(":")
        if len(parts) != 3 and len(parts) != 2:
            raise ValueError(
                f"Invalid vacancy rule: {rule}. Must be 'element:index' or 'element:count'"
            )

        src_el: str = parts[0]
        target: str = parts[1]

        matching_indices: List[int] = [
            i for i, site in enumerate(structure) if site.specie.symbol == src_el
        ]
        if not matching_indices:
            print(f"Warning: No matching elements found for vacancy source '{src_el}'")
            continue

        try:
            val: int = int(target)

            # If the value is less than or equal to the count of matching elements, treat as count-based vacancy creation
            if val <= len(matching_indices) and val > 0 and len(structure) > 20:
                remove_subset = np.random.choice(matching_indices, val, replace=False)
                indices_to_remove.extend(remove_subset)
                print(f"Created {val} vacancies of {src_el} (randomly selected)")
            else:
                # Treat as index-based vacancy creation
                if val < 0 or val >= len(structure):
                    raise IndexError(f"Index {val} out of range")

                actual_symbol: str = structure[val].specie.symbol
                if actual_symbol != src_el:
                    print(
                        f"Warning: Site index {val} is '{actual_symbol}', not vacancy element '{src_el}'. Removing anyway."
                    )

                indices_to_remove.append(val)
                print(f"Removed site {val} ({actual_symbol}) to create vacancy")
        except ValueError:
            raise ValueError(f"Invalid vacancy target: {target}")

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
    return slab.generate_unique_slab_structs()[0]

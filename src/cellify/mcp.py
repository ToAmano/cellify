"""
MCP (Model Context Protocol) server implementation for cellify.
Exposes structure manipulation tools to external LLM agents.
"""

import os
import re
import sys
import uuid
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from pymatgen.core import Structure

from mcp.server.fastmcp import FastMCP

from cellify.core import (
    apply_substitutions,
    apply_vacancies_by_count,
    apply_vacancies_by_index,
    calculate_min_dist_scaling,
    convert_to_conventional,
    generate_surface_slab,
    load_structure_file,
    parse_matrix_string,
    save_structure_file,
)

# Initialize the FastMCP server
mcp: FastMCP = FastMCP("cellify")


def _get_structure_summary_str(structure: Structure, label: str = "") -> str:
    """
    Returns a formatted summary of the structure.
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


def _get_atomic_indices_str(structure: Structure) -> str:
    """
    Returns a formatted table of all atomic indices and coordinates.
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


def _get_structure_string(
    structure: Structure, meta_data: Dict[str, Any], input_filepath: str
) -> str:
    """
    Saves the structure to a temporary file to capture its formatted output string
    and returns it, then cleans up the temporary file.
    """
    base: str
    ext: str
    base, ext = os.path.splitext(input_filepath)
    if not ext:
        if meta_data.get("mode") == "standard" or os.path.basename(base) in ["POSCAR", "CONTCAR"]:
            ext = "_POSCAR"
        elif meta_data.get("mode") in ["espresso_text_replace", "espresso_out"]:
            ext = ".in"

    temp_filename: str = f".__tmp_cellify_mcp_{uuid.uuid4().hex}{ext}"
    try:
        save_structure_file(temp_filename, structure, meta_data)
        with open(temp_filename, "r", encoding="utf-8") as f:
            content: str = f.read()
        return content
    finally:
        if os.path.exists(temp_filename):
            os.remove(temp_filename)


def _apply_scaling(
    structure: Structure, dim: Optional[str], min_dist: Optional[float]
) -> None:
    """
    Applies dimension scaling or min-distance scaling to the structure.
    """
    if dim:
        clean_dim: str = dim.strip()
        if "," in clean_dim or "/" in clean_dim or ";" in clean_dim:
            matrix: np.ndarray = parse_matrix_string(clean_dim)
            structure.make_supercell(matrix)
        else:
            factors: List[int] = [int(x) for x in re.split(r"[\s,]+", clean_dim) if x]
            if len(factors) != 3:
                raise ValueError(
                    f"Diagonal scaling 'dim' must contain exactly 3 integers, got {factors}"
                )
            structure.make_supercell(factors)
    elif min_dist is not None:
        nx: int
        ny: int
        nz: int
        nx, ny, nz = calculate_min_dist_scaling(structure, min_dist)
        structure.make_supercell([nx, ny, nz])


def _normalize_rules(rules: Optional[List[str]]) -> List[str]:
    """
    Ensures the rules parameter is returned as a list of strings.
    """
    if not rules:
        return []
    if isinstance(rules, str):
        return [rules]
    return rules


@mcp.tool()
def cellify_info(filepath: str) -> str:
    """
    Get detailed information about a crystal structure, including lattice parameters,
    atomic species, and a table of 0-based absolute atomic indices and coordinates.
    This is highly recommended to call before generating defects to locate atom indices.

    Args:
        filepath: Path to the input structure file (e.g. POSCAR, cif, qe.in, qe.out).
    """
    if not os.path.exists(filepath):
        return f"Error: File '{filepath}' not found."

    try:
        structure: Structure
        meta_data: Dict[str, Any]
        structure, meta_data = load_structure_file(filepath)
        summary: str = _get_structure_summary_str(
            structure, f"Structure Info for {os.path.basename(filepath)}"
        )
        indices: str = _get_atomic_indices_str(structure)
        return f"{summary}\n{indices}"
    except Exception as e:
        return f"Error loading file details: {str(e)}"


@mcp.tool()
def cellify_conventional(filepath: str, output_path: Optional[str] = None) -> str:
    """
    Convert a crystal structure (e.g., primitive cell downloaded from materials databases)
    into its conventional standard cell representation.

    Args:
        filepath: Path to the input structure file.
        output_path: Optional path to save the resulting conventional cell. If omitted,
                     the formatted structure content is returned as text.
    """
    if not os.path.exists(filepath):
        return f"Error: File '{filepath}' not found."

    try:
        structure: Structure
        meta_data: Dict[str, Any]
        structure, meta_data = load_structure_file(filepath)
        conv_structure: Structure = convert_to_conventional(structure)

        if output_path:
            save_structure_file(output_path, conv_structure, meta_data)
            return (
                f"Successfully converted to conventional cell and saved to {output_path}.\n"
                f"{_get_structure_summary_str(conv_structure)}"
            )
        else:
            return _get_structure_string(conv_structure, meta_data, filepath)
    except Exception as e:
        return f"Error converting to conventional cell: {str(e)}"


@mcp.tool()
def cellify_supercell(
    filepath: str,
    dim: Optional[str] = None,
    min_dist: Optional[float] = None,
    conventional: bool = False,
    output_path: Optional[str] = None,
) -> str:
    """
    Generate a supercell from a structure file. Supports diagonal scaling, transformation matrix,
    or automatic scaling based on minimum periodic distance.

    Args:
        filepath: Path to the input structure file.
        dim: Scaling dimensions. Specify either as diagonal factors (e.g., '2 2 2') or
             a 3x3 transformation matrix (e.g., '2 0 0 / 0 2 0 / 0 0 2').
        min_dist: Target minimum periodic distance in Angstroms for automatic scaling.
        conventional: If True, converts the structure to its conventional cell before scaling.
        output_path: Optional path to save the resulting supercell. If omitted,
                     the formatted structure content is returned as text.
    """
    if not os.path.exists(filepath):
        return f"Error: File '{filepath}' not found."

    try:
        structure: Structure
        meta_data: Dict[str, Any]
        structure, meta_data = load_structure_file(filepath)

        if conventional:
            structure = convert_to_conventional(structure)

        _apply_scaling(structure, dim, min_dist)

        if output_path:
            save_structure_file(output_path, structure, meta_data)
            return (
                f"Successfully generated supercell and saved to {output_path}.\n"
                f"{_get_structure_summary_str(structure)}"
            )
        else:
            return _get_structure_string(structure, meta_data, filepath)
    except Exception as e:
        return f"Error generating supercell: {str(e)}"


@mcp.tool()
def cellify_defect(
    filepath: str,
    substitute: Optional[List[str]] = None,
    vacancy_index: Optional[List[str]] = None,
    vacancy_count: Optional[List[str]] = None,
    dim: Optional[str] = None,
    min_dist: Optional[float] = None,
    conventional: bool = False,
    output_path: Optional[str] = None,
) -> str:
    """
    Introduce defects (vacancies and/or doping/substitutions) into a crystal structure.
    Optionally scales the structure to a supercell before applying defects.

    Args:
        filepath: Path to the input structure file.
        substitute: List of substitution rules, e.g. ["Si:Ge:0", "Si:Al:12%"].
        vacancy_index: List of vacancy rules by index, e.g. ["Si:0,4"].
        vacancy_count: List of vacancy rules by count, e.g. ["Si:2"].
        dim: Scaling dimensions for supercell generation (diagonal factors or matrix).
        min_dist: Target minimum periodic distance in Angstroms for automatic scaling.
        conventional: If True, converts the structure to its conventional cell before scaling/defects.
        output_path: Optional path to save the resulting defect structure. If omitted,
                     the formatted structure content is returned as text.
    """
    if not os.path.exists(filepath):
        return f"Error: File '{filepath}' not found."

    try:
        structure: Structure
        meta_data: Dict[str, Any]
        structure, meta_data = load_structure_file(filepath)

        if conventional:
            structure = convert_to_conventional(structure)

        _apply_scaling(structure, dim, min_dist)

        # Apply substitutions/doping
        sub_list: List[str] = _normalize_rules(substitute)
        if sub_list:
            apply_substitutions(structure, sub_list)

        # Apply vacancies by index
        vac_idx_list: List[str] = _normalize_rules(vacancy_index)
        if vac_idx_list:
            apply_vacancies_by_index(structure, vac_idx_list)

        # Apply vacancies by count
        vac_cnt_list: List[str] = _normalize_rules(vacancy_count)
        if vac_cnt_list:
            apply_vacancies_by_count(structure, vac_cnt_list)

        if output_path:
            save_structure_file(output_path, structure, meta_data)
            return (
                f"Successfully applied defects and saved to {output_path}.\n"
                f"{_get_structure_summary_str(structure)}"
            )
        else:
            return _get_structure_string(structure, meta_data, filepath)
    except Exception as e:
        return f"Error applying defects: {str(e)}"


@mcp.tool()
def cellify_slab(
    filepath: str,
    miller: str,
    thick: float,
    vacuum: float,
    dim: Optional[str] = None,
    min_dist: Optional[float] = None,
    conventional: bool = True,
    output_path: Optional[str] = None,
) -> str:
    """
    Generate a surface slab model with a specified Miller indices and vacuum thickness.
    Highly recommended to conventionalize the input cell first to ensure intuitive slab cutting.

    Args:
        filepath: Path to the input structure file.
        miller: Miller indices of the surface plane. Format as e.g. '1 1 1' or '1,1,1'.
        thick: Slab thickness in Angstroms or layers (minimum layer/atomic layers).
        vacuum: Vacuum layer thickness in Angstroms.
        dim: Scaling dimensions for final supercell generation (diagonal factors or matrix).
        min_dist: Target minimum periodic distance in Angstroms for automatic scaling.
        conventional: If True (default), converts the structure to its conventional cell first.
        output_path: Optional path to save the resulting slab model. If omitted,
                     the formatted structure content is returned as text.
    """
    if not os.path.exists(filepath):
        return f"Error: File '{filepath}' not found."

    try:
        structure: Structure
        meta_data: Dict[str, Any]
        structure, meta_data = load_structure_file(filepath)

        if conventional:
            structure = convert_to_conventional(structure)

        # Parse Miller indices
        miller_indices: List[int] = [
            int(x) for x in re.split(r"[\s,;]+", miller.strip()) if x
        ]
        if len(miller_indices) != 3:
            raise ValueError(
                f"Miller indices 'miller' must contain exactly 3 integers, got {miller_indices}"
            )

        # Generate surface slab
        slab_structure: Structure = generate_surface_slab(
            structure, miller_indices, thick, vacuum
        )

        _apply_scaling(slab_structure, dim, min_dist)

        if output_path:
            save_structure_file(output_path, slab_structure, meta_data)
            return (
                f"Successfully generated slab model and saved to {output_path}.\n"
                f"{_get_structure_summary_str(slab_structure)}"
            )
        else:
            return _get_structure_string(slab_structure, meta_data, filepath)
    except Exception as e:
        return f"Error generating surface slab: {str(e)}"


def main() -> None:
    """
    Main entry point for running the cellify MCP server.
    """
    mcp.run()


if __name__ == "__main__":
    main()

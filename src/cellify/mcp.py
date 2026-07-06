"""
MCP (Model Context Protocol) server implementation for cellify.
Exposes a single structured crystal structure modeling tool to external LLM agents.
"""

import os
import re
import uuid
from typing import Any, Dict, List, Optional

import numpy as np
from mcp.server.fastmcp import FastMCP  # pylint: disable=import-error
from pymatgen.core import Structure

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
        if meta_data.get("mode") == "standard" or os.path.basename(base) in [
            "POSCAR",
            "CONTCAR",
        ]:
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


@mcp.tool()  # type: ignore[misc]
def cellify(  # noqa: C901,CCR001 # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-return-statements,broad-exception-caught
    input_path: str,
    output_path: Optional[str] = None,
    template: Optional[str] = None,
    calc: Optional[str] = None,
    dim: Optional[str] = None,
    min_dist: Optional[float] = None,
    conventional: bool = False,
    substitute: Optional[List[str]] = None,
    vacancy_index: Optional[List[str]] = None,
    vacancy_count: Optional[List[str]] = None,
    slab: Optional[str] = None,
    thick: Optional[float] = None,
    vacuum: Optional[float] = None,
    show_indices: bool = False,
) -> str:
    """
    All-in-one crystal structure modeling tool for cellify.
    Generate supercells, conventional cells, slabs, vacancies, and substitutions.

    Args:
        input_path: Path to the input structure file (e.g. POSCAR, cif, qe.in, qe.out).
        output_path: Optional path to save the resulting structure. If omitted,
                     the formatted structure content is returned as text.
        template: Optional template QE input file to preserve computational parameters.
        calc: Optional override for the QE calculation parameter (e.g. scf, nscf).
        dim: Scaling dimensions (diagonal factors '2 2 2' or 3x3 matrix).
        min_dist: Target minimum periodic distance in Angstroms for automatic scaling.
        conventional: Convert input structure to conventional standard cell first.
        substitute: List of substitution rules (e.g. ['Si:Ge:0', 'Si:Al:12%']).
        vacancy_index: List of vacancy rules by index (e.g. ['Si:0,4']).
        vacancy_count: List of vacancy rules by count (e.g. ['Si:2']).
        slab: Miller indices for surface slab generation (e.g. '1 1 1' or '1,1,1').
        thick: Slab thickness in Angstroms or layers (required if slab is specified).
        vacuum: Vacuum layer thickness in Angstroms (required if slab is specified).
        show_indices: Print absolute atomic indices and coordinate mapping.
    """
    if not os.path.exists(input_path):
        return f"Error: Input file '{input_path}' not found."

    log: List[str] = []
    log.append(f"Loading structure from: {input_path}")

    try:
        structure: Structure
        meta_data: Dict[str, Any]
        structure, meta_data = load_structure_file(input_path)
    except Exception as e:
        return f"Error loading file: {str(e)}"

    log.append(_get_structure_summary_str(structure))

    # 1. Template & calculation override handling
    if template:
        if not os.path.exists(template):
            return f"Error: Template file '{template}' not found."
        try:
            template_meta: Dict[str, Any]
            _, template_meta = load_structure_file(template)
            meta_data = template_meta
        except Exception as e:
            return f"Error loading template file: {str(e)}"

    if calc:
        meta_data["calculation"] = calc

    # 2. Conventional cell conversion
    if conventional:
        log.append("Converting structure to standard conventional cell...")
        structure = convert_to_conventional(structure)

    # 3. Supercell generation
    if dim or min_dist is not None:
        try:
            _apply_scaling(structure, dim, min_dist)
            if dim:
                log.append(f"Applied scaling: {dim}")
            else:
                log.append(f"Applied min-dist scaling (>= {min_dist} A)")
        except Exception as e:
            return f"Error applying scaling: {str(e)}"

    # 4. Substitutions/defects
    sub_list: List[str] = _normalize_rules(substitute)
    if sub_list:
        try:
            apply_substitutions(structure, sub_list)
            log.append(f"Applied substitutions: {sub_list}")
        except Exception as e:
            return f"Error applying substitutions: {str(e)}"

    vac_idx_list: List[str] = _normalize_rules(vacancy_index)
    if vac_idx_list:
        try:
            apply_vacancies_by_index(structure, vac_idx_list)
            log.append(f"Applied vacancy index rules: {vac_idx_list}")
        except Exception as e:
            return f"Error applying vacancy index: {str(e)}"

    vac_cnt_list: List[str] = _normalize_rules(vacancy_count)
    if vac_cnt_list:
        try:
            apply_vacancies_by_count(structure, vac_cnt_list)
            log.append(f"Applied vacancy count rules: {vac_cnt_list}")
        except Exception as e:
            return f"Error applying vacancy count: {str(e)}"

    # 5. Slab generation
    if slab:
        if thick is None or vacuum is None:
            return (
                "Error: Both 'thick' and 'vacuum' must be specified when 'slab' is set."
            )
        try:
            miller_indices: List[int] = [
                int(x) for x in re.split(r"[\s,;]+", slab.strip()) if x
            ]
            if len(miller_indices) != 3:
                return f"Error: Miller indices must contain exactly 3 integers, got {miller_indices}"
            structure = generate_surface_slab(structure, miller_indices, thick, vacuum)
            log.append(f"Generated slab model for Miller indices: {miller_indices}")
        except Exception as e:
            return f"Error generating slab: {str(e)}"

    # Summary of final structure
    log.append(_get_structure_summary_str(structure, label="Final structure summary:"))

    if show_indices:
        log.append(_get_atomic_indices_str(structure))

    # Output handling
    if output_path:
        try:
            save_structure_file(output_path, structure, meta_data)
            log.append(f"Successfully saved final structure to: {output_path}")
            return "\n".join(log)
        except Exception as e:
            return f"Error saving file: {str(e)}"
    else:
        # Return formatted structure content directly
        try:
            struct_str: str = _get_structure_string(structure, meta_data, input_path)
            log_str: str = "\n".join(log)
            # Separate diagnostics from structure content
            return f"{log_str}\n\n=== STRUCTURE CONTENT ===\n{struct_str}"
        except Exception as e:
            return f"Error generating structure text output: {str(e)}"


def main() -> None:
    """
    Main entry point for running the cellify MCP server.
    """
    mcp.run()


if __name__ == "__main__":
    main()

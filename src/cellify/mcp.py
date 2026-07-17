"""
MCP (Model Context Protocol) server implementation for cellify.
Exposes a single structured crystal structure modeling tool to external LLM agents.
"""

# pylint: disable=duplicate-code
import os
from typing import Any, Dict, List, Optional

from mcp.server.fastmcp import FastMCP  # pylint: disable=import-error
from pymatgen.core import Structure

from cellify.core import (
    determine_output_path,
    get_atomic_indices_table,
    get_structure_summary,
    load_structure_file,
    process_template_and_validation,
    run_cellify_pipeline,
    save_structure_file,
)

# Initialize the FastMCP server
mcp: FastMCP = FastMCP("cellify")


@mcp.tool()  # type: ignore[misc]
def cellify(  # noqa: C901,CCR001 # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-return-statements,broad-exception-caught
    input_path: str,
    output_path: Optional[str] = None,
    template: Optional[str] = None,
    calc: Optional[str] = None,
    dim: Optional[List[int]] = None,
    matrix: Optional[str] = None,
    min_dist: Optional[float] = None,
    conventional: bool = False,
    scale_vol: Optional[float] = None,
    scale_lat: Optional[float] = None,
    scale_axes: Optional[List[float]] = None,
    substitute: Optional[List[str]] = None,
    vacancy_index: Optional[List[str]] = None,
    vacancy_count: Optional[List[str]] = None,
    slab: Optional[List[int]] = None,
    thick: Optional[float] = None,
    vacuum: Optional[float] = None,
    show_indices: bool = False,
    select: Optional[int] = None,
) -> str:
    """
    All-in-one crystal structure modeling tool for cellify.
    Generate supercells, conventional cells, slabs, vacancies, and substitutions.

    Args:
        input_path: Path to the input structure file (e.g. POSCAR, cif, qe.in, qe.out).
        output_path: Optional path to save the resulting structure. If omitted,
                     automatically determines the filename (e.g. <input_base>_supercell.<ext>).
        template: Optional template QE input file to preserve computational parameters.
        calc: Optional override for the QE calculation parameter (e.g. scf, nscf).
        dim: Diagonal scaling factors for the supercell (e.g., [2, 2, 2]).
        matrix: 3x3 transformation matrix (e.g. '1 0 0 / 0 1 0 / 0 0 2').
        min_dist: Target minimum periodic distance in Angstroms for automatic scaling.
        conventional: Convert input structure to conventional standard cell first.
        scale_vol: Scale the structure's volume by a given factor.
        scale_lat: Scale the structure's lattice constants by a given factor.
        scale_axes: Scale the individual lattice vectors (axes) by three factors (fa, fb, fc).
        substitute: List of substitution rules (e.g. ['Si:Ge:0', 'Si:Al:12%']).
        vacancy_index: List of vacancy rules by index (e.g. ['Si:0', 'Si:4']).
        vacancy_count: List of vacancy rules by count (e.g. ['Si:2']).
        slab: Miller indices for surface slab generation (e.g. [1, 1, 1]).
        thick: Slab thickness in Angstroms or layers (required if slab is specified).
        vacuum: Vacuum layer thickness in Angstroms (required if slab is specified).
        show_indices: Print absolute atomic indices and coordinate mapping.
        select: 1-based index to select a structure from query results non-interactively.
    """
    log: List[str] = []
    structure: Structure
    meta_data: Dict[str, Any]

    if not os.path.exists(input_path):
        if "/" not in input_path and "\\" not in input_path and "." not in input_path:
            from cellify.optimade import (
                SelectionError,
                retrieve_cif_by_formula,
                select_and_download_structure,
            )

            if select is None:
                try:
                    return retrieve_cif_by_formula(input_path)
                except Exception as e:
                    return f"Error querying formula: {str(e)}"

            try:
                db_name: str
                entry: Dict[str, Any]
                structure, db_name, entry, _ = select_and_download_structure(
                    input_path, select=select
                )
                meta_data = {}
                entry_id = entry.get("id", "unknown")
                log.append(f"Downloading structure from {db_name} (ID: {entry_id})...")
                input_path = f"{input_path}_{entry_id}.cif"
            except SelectionError as e:
                return f"Error: {str(e)}\n\n{e.summary}"
            except ValueError as e:
                return f"Error: {str(e)}"
            except Exception as e:
                return f"Error downloading structure: {str(e)}"
        else:
            return f"Error: Input file '{input_path}' not found."
    else:
        log.append(f"Loading structure from: {input_path}")
        try:
            structure, meta_data = load_structure_file(input_path)
        except Exception as e:
            return f"Error loading file: {str(e)}"

    log.append(get_structure_summary(structure))

    output_path = determine_output_path(input_path, output_path)
    try:
        meta_data = process_template_and_validation(
            meta_data, output_path, template, calc
        )
    except Exception as e:
        return f"Error: {str(e)}"

    # Execute modeling pipeline
    try:
        structure, pipeline_log = run_cellify_pipeline(
            structure,
            conventional=conventional,
            dim=dim,
            matrix=matrix,
            min_dist=min_dist,
            scale_vol=scale_vol,
            scale_lat=scale_lat,
            scale_axes=scale_axes,
            substitute=substitute,
            vacancy_index=vacancy_index,
            vacancy_count=vacancy_count,
            slab=slab,
            thick=thick,
            vacuum=vacuum,
        )
        if pipeline_log:
            log.append(pipeline_log.strip())
    except Exception as e:
        return f"Error processing structure: {str(e)}"

    # Summary of final structure
    log.append(get_structure_summary(structure, label="Final structure summary:"))

    if show_indices:
        log.append(get_atomic_indices_table(structure))

    log.append(f"\nSaving final structure to: {output_path}")
    try:
        save_structure_file(output_path, structure, meta_data)
        log.append("Success!")
        return "\n".join(log)
    except Exception as e:
        return f"Error saving file: {str(e)}"


def main() -> None:
    """
    Main entry point for running the cellify MCP server.
    """
    mcp.run()


if __name__ == "__main__":  # pragma: no cover
    main()

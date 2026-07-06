"""
MCP (Model Context Protocol) server implementation for cellify.
Exposes a single structured crystal structure modeling tool to external LLM agents.
"""

# pylint: disable=duplicate-code
import os
import uuid
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
    substitute: Optional[List[str]] = None,
    vacancy_index: Optional[List[str]] = None,
    vacancy_count: Optional[List[str]] = None,
    slab: Optional[List[int]] = None,
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
        dim: Diagonal scaling factors for the supercell (e.g., [2, 2, 2]).
        matrix: 3x3 transformation matrix (e.g. '1 0 0 / 0 1 0 / 0 0 2').
        min_dist: Target minimum periodic distance in Angstroms for automatic scaling.
        conventional: Convert input structure to conventional standard cell first.
        substitute: List of substitution rules (e.g. ['Si:Ge:0', 'Si:Al:12%']).
        vacancy_index: List of vacancy rules by index (e.g. ['Si:0', 'Si:4']).
        vacancy_count: List of vacancy rules by count (e.g. ['Si:2']).
        slab: Miller indices for surface slab generation (e.g. [1, 1, 1]).
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

    log.append(get_structure_summary(structure))

    # 1. Template & calculation validation/processing
    try:
        # Determine output path for validation checks
        validation_output_path = determine_output_path(input_path, output_path)
        meta_data = process_template_and_validation(
            meta_data, validation_output_path, template, calc
        )
    except Exception as e:
        return f"Error: {str(e)}"

    if slab and (thick is None or vacuum is None):
        return "Error: Both 'thick' and 'vacuum' must be specified when 'slab' is set."

    # Execute modeling pipeline
    try:
        structure, pipeline_log = run_cellify_pipeline(
            structure,
            conventional=conventional,
            dim=dim,
            matrix=matrix,
            min_dist=min_dist,
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

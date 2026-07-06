"""
Command-line interface (CLI) for cellify.
Handles arg parsing, workflow orchestration, and user output reporting.
"""

import argparse
import os
import sys
from typing import Any, Dict, List, Optional, cast

from cellify import __version__
from cellify.core import (
    get_atomic_indices_table,
    get_structure_summary,
    load_structure_file,
    run_cellify_pipeline,
    save_structure_file,
)


def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="cellify: A friendly DFT helper CLI for generating supercells and calculation-ready inputs."
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"cellify {__version__}"
    )

    # I/O options
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input structure file path (e.g. POSCAR, input.cif, qe.in, qe.out)",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output structure file path (default: <input_base>_supercell.<ext>)",
    )
    parser.add_argument(
        "--template",
        type=str,
        help="Template QE input file to preserve computational parameters and comments when generating output.",
    )
    parser.add_argument(
        "--calc",
        "--calculation",
        dest="calc",
        type=str,
        help="Override the calculation parameter in the QE input file (e.g. scf, nscf, bands).",
    )

    # Supercell options
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-d",
        "--dim",
        nargs=3,
        type=int,
        metavar=("nx", "ny", "nz"),
        help="Diagonal scaling factors for the supercell (e.g., -d 2 2 2)",
    )
    group.add_argument(
        "-m",
        "--matrix",
        help="3x3 transformation matrix. Specify as 'r11 r12 r13 / r21 r22 r23 / r31 r32 r33'",
    )
    group.add_argument(
        "--min-dist",
        type=float,
        metavar="DISTANCE",
        help="Automatically generate a supercell where the minimum distance between periodic images is >= DISTANCE (in Angstroms)",
    )
    parser.add_argument(
        "--conventional",
        action="store_true",
        help="Automatically convert the input structure to its standard conventional representation before applying other operations.",
    )

    # Doping / Defect options
    parser.add_argument(
        "--substitute",
        action="append",
        default=[],
        help="Substitution rule: 'element:target_element:index_or_percentage' (e.g., 'Si:P:0' or 'Si:Al:5%%')",
    )
    parser.add_argument(
        "--vacancy-index",
        action="append",
        default=[],
        help="Vacancy index rule: 'element:index' (e.g., 'Si:0' or 'C:33')",
    )
    # Keep --vacancy as an alias for backward compatibility
    parser.add_argument(
        "--vacancy",
        dest="vacancy_index",
        action="append",
        default=[],
        help="Deprecated alias for --vacancy-index",
    )
    parser.add_argument(
        "--vacancy-count",
        action="append",
        default=[],
        help="Vacancy count rule: 'element:count' (e.g., 'O:2')",
    )

    # Slab options
    parser.add_argument(
        "--slab",
        nargs=3,
        type=int,
        metavar=("h", "k", "l"),
        help="Miller indices for surface slab generation (e.g., --slab 1 0 0)",
    )
    parser.add_argument(
        "--thick", type=float, help="Slab thickness (in Angstroms or layers)"
    )
    parser.add_argument(
        "--vacuum", type=float, help="Vacuum layer thickness (in Angstroms)"
    )
    parser.add_argument(
        "-w",
        "--view",
        action="store_true",
        help="Quickly visualize the generated structure in 3D using ASE (requires GUI environment).",
    )
    parser.add_argument(
        "--show-indices",
        action="store_true",
        help="Print absolute atomic indices and coordinate mapping of the final structure.",
    )

    return parser.parse_args(args)


def _determine_output_path(args: argparse.Namespace) -> str:
    """
    Determines the output file path.
    """
    if args.output:
        return cast(str, args.output)

    base, ext = os.path.splitext(args.input)
    # Special case: VASP files like POSCAR or CONTCAR with no extension
    if not ext and os.path.basename(base) in ["POSCAR", "CONTCAR"]:
        return f"{base}_supercell"
    return f"{base}_supercell{ext}"


def _process_template_and_validation(
    args: argparse.Namespace, meta_data: Dict[str, Any], output_path: str
) -> Dict[str, Any]:
    """
    Handles template loading, calculation overrides, and QE I/O format validations.
    """
    if args.template:
        if not os.path.exists(args.template):
            print(f"Error: Template file '{args.template}' not found.", file=sys.stderr)
            sys.exit(1)
        print(f"Loading calculation parameters template from: {args.template}")
        try:
            _, template_meta = load_structure_file(args.template)
            # Retain the original file content/formatting from the template
            meta_data = template_meta
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error loading template file: {e}", file=sys.stderr)
            sys.exit(1)

    if args.calc:
        meta_data["calculation"] = args.calc

    # Validate output format if reading from QE output log
    is_input_qe_output = meta_data.get("mode") == "espresso_out"
    lower_out_path = output_path.lower()
    is_output_qe_input = (
        any(lower_out_path.endswith(ext) for ext in [".in", ".qe", ".pwi"])
        or "qe" in lower_out_path
        or "espresso" in lower_out_path
    )

    if is_input_qe_output and is_output_qe_input and not args.template:
        print(
            "Error: A template QE input file must be specified via --template when reading from a QE output log file.",
            file=sys.stderr,
        )
        sys.exit(1)

    return meta_data


def main() -> None:  # noqa: C901,CCR001
    """
    Main entry point for the cellify CLI utility.
    """
    args: argparse.Namespace = parse_args()

    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading structure from: {args.input}")
    try:
        structure, meta_data = load_structure_file(args.input)
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error loading file: {e}", file=sys.stderr)
        sys.exit(1)

    print(get_structure_summary(structure))

    output_path: str = _determine_output_path(args)
    meta_data = _process_template_and_validation(args, meta_data, output_path)

    # Execute modeling pipeline
    try:
        structure, pipeline_log = run_cellify_pipeline(
            structure,
            conventional=args.conventional,
            dim=args.dim or args.matrix,
            min_dist=args.min_dist,
            substitute=args.substitute,
            vacancy_index=args.vacancy_index,
            vacancy_count=args.vacancy_count,
            slab=args.slab,
            thick=args.thick,
            vacuum=args.vacuum,
        )
        if pipeline_log:
            print(pipeline_log.strip())
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error processing structure: {e}", file=sys.stderr)
        sys.exit(1)

    # Print final structure summary
    print(get_structure_summary(structure, label="Final structure summary:"))

    if args.show_indices:
        print(get_atomic_indices_table(structure))

    print(f"\nSaving final structure to: {output_path}")
    try:
        save_structure_file(output_path, structure, meta_data)
        print("Success!")
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error saving file: {e}", file=sys.stderr)
        sys.exit(1)

    # 3. Optional visualization
    if args.view:
        print("\nOpening structure WebGL viewer...")
        try:
            # pylint: disable=import-outside-toplevel
            from cellify.viewer import open_browser_viewer

            open_browser_viewer(structure)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error launching WebGL viewer: {e}", file=sys.stderr)
            try:
                print("Attempting to fall back to ASE native GUI viewer...")
                import _tkinter  # noqa: F401 # pylint: disable=unused-import
                from ase.visualize import view
                from pymatgen.io.ase import AseAtomsAdaptor

                atoms = AseAtomsAdaptor.get_atoms(structure)
                view(atoms)
            except Exception as ase_err:  # pylint: disable=broad-exception-caught
                print(f"ASE GUI viewer not available: {ase_err}", file=sys.stderr)
                try:
                    print("Falling back to matplotlib 2D projection viewer.")
                    import matplotlib.pyplot as plt
                    from ase.visualize.plot import plot_atoms
                    from pymatgen.io.ase import AseAtomsAdaptor

                    atoms = AseAtomsAdaptor.get_atoms(structure)
                    _, ax = plt.subplots(figsize=(6, 6))
                    plot_atoms(atoms, ax, rotation="10x,10y,0z")
                    ax.set_axis_off()
                    plt.tight_layout()
                    print("Close the matplotlib window to continue.")
                    plt.show()
                except (
                    Exception
                ) as fallback_err:  # pylint: disable=broad-exception-caught
                    print(
                        f"Error launching matplotlib viewer fallback: {fallback_err}",
                        file=sys.stderr,
                    )
                    sys.exit(1)


if __name__ == "__main__":
    main()

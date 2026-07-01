"""
Command-line interface (CLI) for cellify.
Handles arg parsing, workflow orchestration, and user output reporting.
"""

import argparse
import os
import sys
from typing import Any, Dict, List, Optional, cast

import numpy as np
from pymatgen.core import Structure

from cellify import __version__
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


def _print_structure_summary(structure: Structure, label: str = "") -> None:
    """
    Prints a formatted summary of the structure.
    """
    if label:
        print(f"\n{label}")
    print(f"  Formula: {structure.composition.reduced_formula}")
    print(f"  Volume:  {structure.volume:.3f} A^3")
    print(f"  Number of atoms: {len(structure)}")
    if label:
        print("  Lattice constants:")
        print(
            f"    a = {structure.lattice.a:.4f} A, b = {structure.lattice.b:.4f} A, c = {structure.lattice.c:.4f} A"
        )
        print(
            f"    alpha = {structure.lattice.alpha:.2f} deg, beta = {structure.lattice.beta:.2f} deg, gamma = {structure.lattice.gamma:.2f} deg"
        )


def _print_atomic_indices(structure: Structure) -> None:
    """
    Prints a formatted table of all atomic indices and coordinates.
    """
    print("\nAbsolute Atomic Indices & Coordinates:")
    print("-" * 78)
    header = (
        f"{'Index':<6} {'Element':<8} "
        f"{'Fractional Coordinates (a, b, c)':<36} "
        f"{'Cartesian (x, y, z)':<22}"
    )
    print(header)
    print("-" * 78)
    for idx, site in enumerate(structure):
        frac = (
            f"[{site.frac_coords[0]:.4f}, "
            f"{site.frac_coords[1]:.4f}, "
            f"{site.frac_coords[2]:.4f}]"
        )
        cart = (
            f"[{site.coords[0]:.3f}, "
            f"{site.coords[1]:.3f}, "
            f"{site.coords[2]:.3f}]"
        )
        print(f"{idx:<6} {site.species_string:<8} {frac:<36} {cart:<22}")
    print("-" * 78)
    print(f"Total: {len(structure)} atoms ({structure.composition.formula})")


def _apply_supercell(structure: Structure, args: argparse.Namespace) -> None:
    """
    Applies supercell generation options to the structure.
    """
    if args.dim:
        print(f"Generating supercell with diagonal scaling: {args.dim}")
        structure.make_supercell(args.dim)
    elif args.matrix:
        try:
            matrix: np.ndarray = parse_matrix_string(args.matrix)
            print(f"Generating supercell with matrix:\n{matrix}")
            structure.make_supercell(matrix)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error parsing matrix: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.min_dist:
        nx, ny, nz = calculate_min_dist_scaling(structure, args.min_dist)
        print(
            f"Calculated scaling for minimum distance >= {args.min_dist} A: [{nx}, {ny}, {nz}]"
        )
        structure.make_supercell([nx, ny, nz])


def _apply_defects_and_slab(
    structure: Structure, args: argparse.Namespace
) -> Structure:
    """
    Applies substitutions, vacancies, and surface slab options to the structure.
    """
    if args.substitute:
        try:
            apply_substitutions(structure, args.substitute)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error applying substitutions: {e}", file=sys.stderr)
            sys.exit(1)

    if args.vacancy_index:
        try:
            apply_vacancies_by_index(structure, args.vacancy_index)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error applying vacancy index: {e}", file=sys.stderr)
            sys.exit(1)

    if args.vacancy_count:
        try:
            apply_vacancies_by_count(structure, args.vacancy_count)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error applying vacancy count: {e}", file=sys.stderr)
            sys.exit(1)

    if args.slab:
        print(f"Generating slab model for Miller indices: {args.slab}")
        try:
            structure = generate_surface_slab(
                structure, args.slab, args.thick, args.vacuum
            )
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error generating slab: {e}", file=sys.stderr)
            sys.exit(1)

    return structure


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


def main() -> None:
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

    _print_structure_summary(structure)

    output_path: str = _determine_output_path(args)
    meta_data = _process_template_and_validation(args, meta_data, output_path)

    # 0. Conventional cell conversion
    if args.conventional:
        print("Converting structure to standard conventional cell...")
        structure = convert_to_conventional(structure)

    # 1. Supercell generation
    _apply_supercell(structure, args)

    # 2. Defect and slab modifications
    structure = _apply_defects_and_slab(structure, args)

    # Print final structure summary
    _print_structure_summary(structure, label="Final structure summary:")

    if args.show_indices:
        _print_atomic_indices(structure)

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

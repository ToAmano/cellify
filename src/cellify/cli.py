"""
Command-line interface (CLI) for cellify.
Handles arg parsing, workflow orchestration, and user output reporting.
"""

import argparse
import os
import sys
import time
from typing import Any, Dict, List, Optional

from pymatgen.core import Structure
from rich.align import Align
from rich.console import Console
from rich.panel import Panel
from rich.text import Text

from cellify import __version__
from cellify.core import (
    determine_output_path,
    get_atomic_indices_table,
    get_structure_summary,
    load_structure_file,
    process_template_and_validation,
    run_cellify_pipeline,
    save_structure_file,
)


def animate_print(text: str, delay: float = 0.03) -> None:
    """
    Prints text line-by-line with a small delay if stdout is a TTY.
    """
    if not text:
        return
    lines = text.splitlines()
    if sys.stdout.isatty():
        actual_delay = delay if len(lines) <= 30 else 0.0
        for line in lines:
            print(line, flush=True)
            time.sleep(actual_delay)
    else:
        print(text.rstrip("\r\n"), flush=True)


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
        "--select",
        type=int,
        help="1-based index to select a structure from query results (non-interactive).",
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


def main() -> None:  # noqa: C901,CCR001
    """
    Main entry point for the cellify CLI utility.
    """
    console: Console = Console()
    if sys.stdout.isatty():
        banner_text: Text = Text()
        banner_text.append("C E L L I F Y\n", style="bold cyan")
        banner_text.append(
            f"A friendly DFT helper for crystal structures | v{__version__}",
            style="dim italic",
        )
        console.print(
            Panel(Align.center(banner_text), border_style="cyan", expand=False)
        )

    args: argparse.Namespace = parse_args()

    structure: Structure
    meta_data: Dict[str, Any]

    if not os.path.exists(args.input):
        # Fallback detection: If the input file does not exist and has no path separator
        # or extension, we treat it as a chemical formula and query external databases
        # via OPTIMADE API.
        if "/" not in args.input and "\\" not in args.input and "." not in args.input:
            from cellify.optimade import select_and_download_structure

            formula = args.input
            status = None

            def cli_prompt(summary: str, limit: int) -> int:
                """Handles interactive user inputs from stdout/stdin and validates the choice index.

                Displays the search results summary and prompts the user to enter a
                valid 1-based index corresponding to a candidate structure.
                """
                if status is not None:
                    status.stop()
                animate_print(summary)
                try:
                    choice_str = input(f"Select a structure (1-{limit}): ").strip()
                    if not choice_str:
                        raise ValueError("Invalid selection.")
                    choice = int(choice_str)
                    if not 1 <= choice <= limit:
                        raise ValueError("Invalid selection.")
                    return choice
                except (KeyboardInterrupt, EOFError):
                    print()
                    raise ValueError("Invalid selection.") from None
                except ValueError:
                    raise ValueError("Invalid selection.") from None

            # Selection mode configuration: If the user specified a selection index via
            # --select, it directly determines the structure to download. Otherwise,
            # we trigger the interactive prompt selection only if we are in a TTY environment.
            interactive_prompt = cli_prompt if sys.stdin.isatty() else None

            try:
                db_name: str
                entry: Dict[str, Any]
                summary: str
                if sys.stdout.isatty():
                    status = console.status(
                        "[bold green]Querying databases...", spinner="dots"
                    )
                    status.start()
                try:
                    structure, db_name, entry, summary = select_and_download_structure(
                        formula,
                        select=args.select,
                        interactive_prompt=interactive_prompt,
                    )
                finally:
                    if status is not None:
                        status.stop()
            except ValueError as e:
                summary = getattr(e, "summary", "")
                if summary:
                    animate_print(summary)
                print(f"Error: {e}", file=sys.stderr)
                sys.exit(1)
            except Exception as e:  # pylint: disable=broad-exception-caught
                print(f"Error downloading structure: {e}", file=sys.stderr)
                sys.exit(1)

            if args.select is not None:
                animate_print(summary)

            print(f"Downloading structure from {db_name} (ID: {entry.get('id')})...")
            entry_id = entry.get("id", "unknown")
            args.input = f"{formula}_{entry_id}.cif"
            meta_data = {}
        else:
            print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Loading structure from: {args.input}")
        try:
            structure, meta_data = load_structure_file(args.input)
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(f"Error loading file: {e}", file=sys.stderr)
            sys.exit(1)

    animate_print(get_structure_summary(structure))

    output_path: str = determine_output_path(args.input, args.output)
    try:
        meta_data = process_template_and_validation(
            meta_data, output_path, args.template, args.calc
        )
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    # Execute modeling pipeline
    try:
        if sys.stdout.isatty():
            with console.status("[bold green]Processing structure...", spinner="dots"):
                structure, _ = run_cellify_pipeline(
                    structure,
                    conventional=args.conventional,
                    dim=args.dim,
                    matrix=args.matrix,
                    min_dist=args.min_dist,
                    substitute=args.substitute,
                    vacancy_index=args.vacancy_index,
                    vacancy_count=args.vacancy_count,
                    slab=args.slab,
                    thick=args.thick,
                    vacuum=args.vacuum,
                    capture_output=False,
                )
        else:
            structure, _ = run_cellify_pipeline(
                structure,
                conventional=args.conventional,
                dim=args.dim,
                matrix=args.matrix,
                min_dist=args.min_dist,
                substitute=args.substitute,
                vacancy_index=args.vacancy_index,
                vacancy_count=args.vacancy_count,
                slab=args.slab,
                thick=args.thick,
                vacuum=args.vacuum,
                capture_output=False,
            )
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error processing structure: {e}", file=sys.stderr)
        sys.exit(1)

    # Print final structure summary
    animate_print(get_structure_summary(structure, label="Final structure summary:"))

    if args.show_indices:
        animate_print(get_atomic_indices_table(structure))

    print(f"\nSaving final structure to: {output_path}")
    try:
        save_structure_file(output_path, structure, meta_data)
        animate_print("Success!")
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
                import _tkinter  # noqa: F401 # pylint: disable=unused-import,import-error
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


if __name__ == "__main__":  # pragma: no cover
    main()

import argparse
import sys
from cellify import __version__

def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description="cellify: A friendly DFT helper CLI for generating supercells and calculation-ready inputs."
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"cellify {__version__}"
    )
    
    # I/O
    parser.add_argument(
        "-i", "--input", required=True, help="Input structure file path (e.g. POSCAR, input.cif, qe.in)"
    )
    parser.add_argument(
        "-o", "--output", help="Output structure file path (default: input_supercell.<ext>)"
    )

    # Supercell options
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-d", "--dim", nargs=3, type=int, metavar=("nx", "ny", "nz"),
        help="Diagonal scaling factors for the supercell (e.g., -d 2 2 2)"
    )
    group.add_argument(
        "-m", "--matrix",
        help="3x3 transformation matrix. Specify as 'r11 r12 r13 / r21 r22 r23 / r31 r32 r33'"
    )
    group.add_argument(
        "--min-dist", type=float, metavar="DISTANCE",
        help="Automatically generate a supercell where the minimum distance between periodic images is >= DISTANCE (in Angstroms)"
    )

    # Doping / Defect options
    parser.add_argument(
        "--substitute", action="append", default=[],
        help="Substitution rule: 'element:target_element:index_or_percentage' (e.g., 'Si:P:0' or 'Si:Al:5%%')"
    )
    parser.add_argument(
        "--vacancy", action="append", default=[],
        help="Vacancy rule: 'element:index_or_count' (e.g., 'Si:0' or 'O:2')"
    )

    # Slab options
    parser.add_argument(
        "--slab", nargs=3, type=int, metavar=("h", "k", "l"),
        help="Miller indices for surface slab generation (e.g., --slab 1 0 0)"
    )
    parser.add_argument(
        "--thick", type=float, help="Slab thickness (in Angstroms or layers)"
    )
    parser.add_argument(
        "--vacuum", type=float, help="Vacuum layer thickness (in Angstroms)"
    )

    return parser.parse_args(args)

def main():
    args = parse_args()
    print("Welcome to cellify!")
    print(f"Input file: {args.input}")
    if args.output:
        print(f"Output file: {args.output}")
    else:
        print("Output file: (automatically generated)")

    if args.dim:
        print(f"Supercell dimension: {args.dim}")
    elif args.matrix:
        print(f"Transformation matrix: {args.matrix}")
    elif args.min_dist:
        print(f"Minimum distance constraint: {args.min_dist} A")

    if args.substitute:
        print(f"Substitutions: {args.substitute}")
    if args.vacancy:
        print(f"Vacancies: {args.vacancy}")

    if args.slab:
        print(f"Slab Miller indices: {args.slab}")
        print(f"Slab thickness: {args.thick}")
        print(f"Vacuum thickness: {args.vacuum}")

    # TODO: core.py / parser.py のロジックを呼び出す

if __name__ == "__main__":
    main()

import argparse
import sys
import os
from typing import List, Optional, Tuple, Dict, Any
import numpy as np
from pymatgen.core import Structure
from cellify import __version__
from cellify.core import (
    load_structure_file,
    save_structure_file,
    parse_matrix_string,
    calculate_min_dist_scaling,
    apply_substitutions,
    apply_vacancies,
    generate_surface_slab
)

def parse_args(args: Optional[List[str]] = None) -> argparse.Namespace:
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
        "-o", "--output", help="Output structure file path (default: <input_base>_supercell.<ext>)"
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

def main() -> None:
    args: argparse.Namespace = parse_args()
    
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
        
    print(f"Loading structure from: {args.input}")
    try:
        structure, meta_data = load_structure_file(args.input)
    except Exception as e:
        print(f"Error loading file: {e}", file=sys.stderr)
        sys.exit(1)
        
    print(f"  Formula: {structure.composition.reduced_formula}")
    print(f"  Volume:  {structure.volume:.3f} A^3")
    print(f"  Number of atoms: {len(structure)}")
    
    # 1. Supercell generation
    if args.dim:
        print(f"Generating supercell with diagonal scaling: {args.dim}")
        structure.make_supercell(args.dim)
    elif args.matrix:
        try:
            matrix: np.ndarray = parse_matrix_string(args.matrix)
            print(f"Generating supercell with matrix:\n{matrix}")
            structure.make_supercell(matrix)
        except Exception as e:
            print(f"Error parsing matrix: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.min_dist:
        nx, ny, nz = calculate_min_dist_scaling(structure, args.min_dist)
        print(f"Calculated scaling for minimum distance >= {args.min_dist} A: [{nx}, {ny}, {nz}]")
        structure.make_supercell([nx, ny, nz])
        
    # 2. Defect and doping modifications
    if args.substitute:
        try:
            apply_substitutions(structure, args.substitute)
        except Exception as e:
            print(f"Error applying substitutions: {e}", file=sys.stderr)
            sys.exit(1)
            
    if args.vacancy:
        try:
            apply_vacancies(structure, args.vacancy)
        except Exception as e:
            print(f"Error applying vacancies: {e}", file=sys.stderr)
            sys.exit(1)

    # 3. Slab generation
    if args.slab:
        print(f"Generating slab model for Miller indices: {args.slab}")
        try:
            structure = generate_surface_slab(structure, args.slab, args.thick, args.vacuum)
        except Exception as e:
            print(f"Error generating slab: {e}", file=sys.stderr)
            sys.exit(1)
            
    # Print final structure summary
    print("\nFinal structure summary:")
    print(f"  Formula: {structure.composition.reduced_formula}")
    print(f"  Volume:  {structure.volume:.3f} A^3")
    print(f"  Number of atoms: {len(structure)}")
    print("  Lattice constants:")
    print(f"    a = {structure.lattice.a:.4f} A, b = {structure.lattice.b:.4f} A, c = {structure.lattice.c:.4f} A")
    print(f"    alpha = {structure.lattice.alpha:.2f} deg, beta = {structure.lattice.beta:.2f} deg, gamma = {structure.lattice.gamma:.2f} deg")

    # Determine default output filename if not specified
    if not args.output:
        base, ext = os.path.splitext(args.input)
        if not ext and base in ["POSCAR", "CONTCAR"]:
            args.output = f"{base}_supercell"
        else:
            args.output = f"{base}_supercell{ext}"
            
    print(f"\nSaving final structure to: {args.output}")
    try:
        save_structure_file(args.output, structure, meta_data)
        print("Success!")
    except Exception as e:
        print(f"Error saving file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()

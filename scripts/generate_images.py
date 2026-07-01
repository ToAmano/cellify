import os
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read as ase_read
from ase.visualize.plot import plot_atoms
from ase.build import bulk
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor

from cellify.core import (
    load_structure_file,
    convert_to_conventional,
    calculate_min_dist_scaling,
    apply_substitutions,
    apply_vacancies_by_index,
    generate_surface_slab,
)


def save_render(atoms, filepath, figsize=(6, 6), rotation="15x,30y,0z", show_unit_cell=2):
    fig, ax = plt.subplots(figsize=figsize)
    # Use show_unit_cell=show_unit_cell to draw cell borders
    plot_atoms(atoms, ax, rotation=rotation, show_unit_cell=show_unit_cell)
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(filepath, bbox_inches="tight", transparent=True, dpi=200)
    plt.close(fig)
    print(f"Saved: {filepath}")


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    root_dir = os.path.abspath(os.path.join(script_dir, ".."))
    output_dir = os.path.join(root_dir, "docs/images")
    os.makedirs(output_dir, exist_ok=True)

    poscar_path = os.path.join(root_dir, "tests/POSCAR")
    qe_in_path = os.path.join(root_dir, "tests/qe.in")

    # ==========================================
    # Example 1: POSCAR -> 2x2x3 Supercell
    # ==========================================
    ex1_in, _ = load_structure_file(poscar_path)
    ex1_out = ex1_in.copy()
    ex1_out.make_supercell([2, 2, 3])

    save_render(AseAtomsAdaptor.get_atoms(ex1_in), os.path.join(output_dir, "ex1_input.png"), rotation="15x,45y,0z")
    save_render(AseAtomsAdaptor.get_atoms(ex1_out), os.path.join(output_dir, "ex1_output.png"), rotation="15x,45y,0z")

    # ==========================================
    # Example 2: Primitive Si -> Conventional 2x2x2 Supercell
    # ==========================================
    ex2_in, _ = load_structure_file(poscar_path)
    # 1. Convert to conventional standard structure
    ex2_conv = convert_to_conventional(ex2_in)
    # 2. Scale conventional cell to 2x2x2 supercell
    ex2_out = ex2_conv.copy()
    ex2_out.make_supercell([2, 2, 2])

    save_render(AseAtomsAdaptor.get_atoms(ex2_in), os.path.join(output_dir, "ex2_input.png"), rotation="15x,45y,0z")
    save_render(AseAtomsAdaptor.get_atoms(ex2_out), os.path.join(output_dir, "ex2_output.png"), rotation="15x,45y,0z")

    # ==========================================
    # Example 3: Hexagonal qe.in -> Orthogonal matrix cell
    # ==========================================
    ex3_in, _ = load_structure_file(qe_in_path)
    ex3_out = ex3_in.copy()
    # Apply matrix transformation: 1 -1 0 / 1 1 0 / 0 0 1
    matrix = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])
    ex3_out.make_supercell(matrix)

    save_render(AseAtomsAdaptor.get_atoms(ex3_in), os.path.join(output_dir, "ex3_input.png"), rotation="15x,45y,0z")
    save_render(AseAtomsAdaptor.get_atoms(ex3_out), os.path.join(output_dir, "ex3_output.png"), rotation="15x,45y,0z")

    # ==========================================
    # Example 4: POSCAR -> min-dist 15.0 supercell
    # ==========================================
    ex4_in, _ = load_structure_file(poscar_path)
    nx, ny, nz = calculate_min_dist_scaling(ex4_in, 15.0)
    ex4_out = ex4_in.copy()
    ex4_out.make_supercell([nx, ny, nz])

    save_render(AseAtomsAdaptor.get_atoms(ex4_in), os.path.join(output_dir, "ex4_input.png"), rotation="15x,45y,0z")
    save_render(AseAtomsAdaptor.get_atoms(ex4_out), os.path.join(output_dir, "ex4_output.png"), rotation="15x,45y,0z")

    # ==========================================
    # Example 5: Si_unit.cif -> 3x3x3 Doped Supercell (substitute Si:P:0)
    # ==========================================
    # Use conventional standard silicon cell as unit cell input
    si_prim_atoms = bulk("Si")
    si_prim_pmg = AseAtomsAdaptor.get_structure(si_prim_atoms)
    ex5_in = convert_to_conventional(si_prim_pmg)

    # Generate 3x3x3 supercell
    ex5_out = ex5_in.copy()
    ex5_out.make_supercell([3, 3, 3])
    # Dope: substitute Silicon atom at index 0 with Phosphorus
    apply_substitutions(ex5_out, ["Si:P:0"])

    save_render(AseAtomsAdaptor.get_atoms(ex5_in), os.path.join(output_dir, "ex5_input.png"), rotation="15x,45y,0z")
    save_render(AseAtomsAdaptor.get_atoms(ex5_out), os.path.join(output_dir, "ex5_output.png"), rotation="15x,45y,0z")

    # ==========================================
    # Example 6: Si_supercell -> vacancy-index "Si:0"
    # ==========================================
    # Setup Si supercell (conventional 2x2x2)
    ex6_in = convert_to_conventional(si_prim_pmg)
    ex6_in.make_supercell([2, 2, 2])

    # Deletes Silicon atom at absolute index 0
    ex6_out = ex6_in.copy()
    apply_vacancies_by_index(ex6_out, ["Si:0"])

    # Render input first and get position of atom index 0
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_atoms(
        AseAtomsAdaptor.get_atoms(ex6_in),
        ax,
        rotation="15x,45y,0z",
        show_unit_cell=2,
    )
    from matplotlib.patches import Circle

    circles = [p for p in ax.patches if isinstance(p, Circle)]
    vacancy_pos = circles[0].center
    vacancy_radius = circles[0].radius

    ax.set_axis_off()
    plt.tight_layout()
    ex6_in_path = os.path.join(output_dir, "ex6_input.png")
    plt.savefig(ex6_in_path, bbox_inches="tight", transparent=True, dpi=200)
    plt.close(fig)
    print(f"Saved: {ex6_in_path}")

    # Render output and overlay the ghost atom + annotation
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_atoms(
        AseAtomsAdaptor.get_atoms(ex6_out),
        ax,
        rotation="15x,45y,0z",
        show_unit_cell=2,
    )

    # Draw ghost atom (red dashed circle)
    ghost = Circle(
        vacancy_pos,
        vacancy_radius,
        fill=False,
        edgecolor="red",
        linestyle="--",
        linewidth=2.0,
    )
    ax.add_patch(ghost)

    # Draw arrow and label pointing to it
    ax.annotate(
        "Vacancy (Si:0)",
        xy=vacancy_pos,
        xytext=(vacancy_pos[0] - 2.5, vacancy_pos[1] + 2.5),
        arrowprops=dict(
            facecolor="red", shrink=0.1, width=1.5, headwidth=6, headlength=6
        ),
        color="red",
        fontweight="bold",
        fontsize=11,
    )

    ax.set_axis_off()
    plt.tight_layout()
    ex6_out_path = os.path.join(output_dir, "ex6_output.png")
    plt.savefig(ex6_out_path, bbox_inches="tight", transparent=True, dpi=200)
    plt.close(fig)
    print(f"Saved: {ex6_out_path}")

    # ==========================================
    # Example 7: STO_bulk -> STO_100_slab
    # ==========================================
    # Build SrTiO3 unit cell
    lattice = Lattice.cubic(3.905)
    species = ["Sr", "Ti", "O", "O", "O"]
    coords = [[0, 0, 0], [0.5, 0.5, 0.5], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]
    ex7_in = Structure(lattice, species, coords)

    # Generate slab model (100) miller index, 12A thick, 15A vacuum
    ex7_out = generate_surface_slab(ex7_in, [1, 0, 0], 12.0, 15.0)

    save_render(AseAtomsAdaptor.get_atoms(ex7_in), os.path.join(output_dir, "ex7_input.png"), rotation="15x,45y,0z")
    save_render(
        AseAtomsAdaptor.get_atoms(ex7_out),
        os.path.join(output_dir, "ex7_output.png"),
        figsize=(3, 9),
        rotation="10x,45y,0z",
    )




if __name__ == "__main__":
    main()

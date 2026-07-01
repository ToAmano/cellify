import os
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
from ase.build import bulk
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor


def save_render(atoms, filepath, figsize=(5, 5), rotation="15x,30y,0z"):
    fig, ax = plt.subplots(figsize=figsize)
    plot_atoms(atoms, ax, rotation=rotation)
    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(filepath, bbox_inches="tight", transparent=True, dpi=150)
    plt.close(fig)
    print(f"Saved: {filepath}")


def main():
    output_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "../docs/images")
    )
    os.makedirs(output_dir, exist_ok=True)

    # 1. Conventional Cell Conversion (Example 3)
    # Get primitive Silicon
    si_prim_atoms = bulk("Si")  # 2 atoms, rhombohedral
    si_prim_pmg = AseAtomsAdaptor.get_structure(si_prim_atoms)

    # Convert to conventional
    sga = SpacegroupAnalyzer(si_prim_pmg)
    si_conv_pmg = sga.get_conventional_standard_structure()
    si_conv_atoms = AseAtomsAdaptor.get_atoms(si_conv_pmg)

    save_render(
        si_prim_atoms,
        os.path.join(output_dir, "ex3_primitive.png"),
        rotation="15x,45y,0z",
    )
    save_render(
        si_conv_atoms,
        os.path.join(output_dir, "ex3_conventional.png"),
        rotation="15x,45y,0z",
    )

    # 2. Doping / Substitution (Example 5)
    # Create clean Si 2x2x2 supercell
    si_bulk_pmg = sga.get_conventional_standard_structure()
    si_bulk_pmg.make_supercell([2, 2, 2])
    save_render(
        AseAtomsAdaptor.get_atoms(si_bulk_pmg),
        os.path.join(output_dir, "ex5_bulk.png"),
    )

    # Substitute index 0 with P
    si_doped_pmg = si_bulk_pmg.copy()
    si_doped_pmg.replace(0, "P")
    save_render(
        AseAtomsAdaptor.get_atoms(si_doped_pmg),
        os.path.join(output_dir, "ex5_doped.png"),
    )

    # 3. Vacancy (Example 6)
    # Use the same bulk Si 2x2x2 supercell
    # Remove index 0 to create a vacancy
    si_vac_pmg = si_bulk_pmg.copy()
    si_vac_pmg.remove_sites([0])
    save_render(
        AseAtomsAdaptor.get_atoms(si_vac_pmg),
        os.path.join(output_dir, "ex6_vacancy.png"),
    )

    # 4. Surface Slab (Example 7)
    # Create SrTiO3 bulk unit cell
    lattice = Lattice.cubic(3.905)
    species = ["Sr", "Ti", "O", "O", "O"]
    coords = [
        [0, 0, 0],
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 0],
        [0.5, 0, 0.5],
        [0, 0.5, 0.5],
    ]
    sto_bulk = Structure(lattice, species, coords)
    save_render(
        AseAtomsAdaptor.get_atoms(sto_bulk),
        os.path.join(output_dir, "ex7_bulk.png"),
    )

    # Create Slab (100) miller index, 12A thick, 15A vacuum
    slab_gen = SlabGenerator(
        sto_bulk,
        miller_index=(1, 0, 0),
        min_slab_size=12.0,
        min_vacuum_size=15.0,
    )
    slab = slab_gen.get_slabs()[0]
    save_render(
        AseAtomsAdaptor.get_atoms(slab),
        os.path.join(output_dir, "ex7_slab.png"),
        figsize=(4, 7),
        rotation="10x,45y,0z",
    )


if __name__ == "__main__":
    main()

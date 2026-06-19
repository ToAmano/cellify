import os
import re

import numpy as np
import pytest
from pymatgen.core import Structure

from cellify.core import (
    apply_substitutions,
    apply_vacancies,
    calculate_min_dist_scaling,
    convert_to_conventional,
    generate_surface_slab,
    load_structure_file,
    parse_matrix_string,
    save_structure_file,
)


@pytest.fixture
def poscar_path():
    return os.path.join(os.path.dirname(__file__), "POSCAR")


@pytest.fixture
def qe_path():
    return os.path.join(os.path.dirname(__file__), "qe.in")


def test_load_and_save_poscar(poscar_path, tmp_path):
    # Test loading structure
    structure, meta_data = load_structure_file(poscar_path)
    assert isinstance(structure, Structure)
    assert meta_data["mode"] == "standard"
    assert len(structure) == 2
    assert structure.composition.reduced_formula == "Si"

    # Test saving structure
    out_path = tmp_path / "POSCAR_out"
    save_structure_file(str(out_path), structure, meta_data)
    assert out_path.exists()

    # Re-load and verify
    struct_new, _ = load_structure_file(str(out_path))
    assert len(struct_new) == 2


def test_load_and_save_qe(qe_path, tmp_path):
    # Test loading structure
    structure, meta_data = load_structure_file(qe_path)
    assert isinstance(structure, Structure)
    assert meta_data["mode"] == "espresso_text_replace"
    assert len(structure) == 2

    # Test saving and dynamic parameter updates (save as a 2x2x2 supercell)
    structure.make_supercell([2, 2, 2])
    out_path = tmp_path / "qe_out.in"
    save_structure_file(str(out_path), structure, meta_data)

    assert out_path.exists()
    with open(out_path, "r") as f:
        content = f.read()

    # Check if nat is updated to 16
    assert re.search(r"nat\s*=\s*16", content) is not None
    # Check if ntyp remains 1
    assert re.search(r"ntyp\s*=\s*1", content) is not None
    # Check if original calculation parameter is preserved
    assert "calculation = 'scf'" in content


def test_parse_matrix_string():
    matrix_str = "1 -1 0 / 1 1 0 / 0 0 2"
    matrix = parse_matrix_string(matrix_str)
    assert isinstance(matrix, np.ndarray)
    assert matrix.shape == (3, 3)
    assert np.allclose(matrix, [[1, -1, 0], [1, 1, 0], [0, 0, 2]])

    with pytest.raises(ValueError):
        parse_matrix_string("1 0 0 / 0 1 0")  # Insufficient number of rows


def test_calculate_min_dist_scaling(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    # Get required scaling for min-dist >= 10.0 A
    nx, ny, nz = calculate_min_dist_scaling(structure, 10.0)
    assert nx == 4
    assert ny == 4
    assert nz == 4


def test_apply_substitutions(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    structure.make_supercell([2, 2, 2])  # 16 atoms

    # Substitute the 0th atom with P
    apply_substitutions(structure, ["Si:P:0"])
    assert structure.composition.reduced_formula == "Si15P"
    assert structure[0].specie.symbol == "P"


def test_apply_vacancies(poscar_path):
    structure, _ = load_structure_file(poscar_path)

    # Create a vacancy by removing the 0th atom
    apply_vacancies(structure, ["Si:0"])
    assert len(structure) == 1


def test_convert_to_conventional(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    assert len(structure) == 2  # Primitive cell has 2 atoms

    conv_structure = convert_to_conventional(structure)
    assert len(conv_structure) == 8  # Conventional cell has 8 atoms
    assert conv_structure.composition.reduced_formula == "Si"

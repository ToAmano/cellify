import os
import re
from unittest.mock import patch

import numpy as np
import pytest
from pymatgen.core import Structure

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
    apply_vacancies_by_index(structure, ["Si:0"])
    assert len(structure) == 1


def test_convert_to_conventional(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    assert len(structure) == 2  # Primitive cell has 2 atoms

    conv_structure = convert_to_conventional(structure)
    assert len(conv_structure) == 8  # Conventional cell has 8 atoms
    assert conv_structure.composition.reduced_formula == "Si"


def test_parse_matrix_string_errors():
    with pytest.raises(ValueError, match="Matrix string must define exactly 3 rows"):
        parse_matrix_string("1 0 0 / 0 1 0")
    with pytest.raises(ValueError, match="Each row in the matrix must have exactly 3 elements"):
        parse_matrix_string("1 0 0 / 0 1 0 / 0 0")


def test_apply_substitutions_errors(poscar_path):
    structure, _ = load_structure_file(poscar_path)

    # Rule split error
    with pytest.raises(ValueError, match="Invalid substitution rule"):
        apply_substitutions(structure, ["Si:P"])

    # Matching elements not found (warning path)
    apply_substitutions(structure, ["H:P:0"])

    # Index out of range
    with pytest.raises(IndexError, match="out of range"):
        apply_substitutions(structure, ["Si:P:999"])

    # Invalid index target
    with pytest.raises(ValueError, match="Invalid substitution target index or percentage"):
        apply_substitutions(structure, ["Si:P:abc"])

    # Warning path: actual element does not match src_el
    structure2, _ = load_structure_file(poscar_path)
    apply_substitutions(structure2, ["Si:P:0"])
    apply_substitutions(structure2, ["Si:Al:0"])


def test_apply_substitutions_percentage(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    structure.make_supercell([2, 2, 2]) # 16 atoms

    # Percentage substitution
    apply_substitutions(structure, ["Si:P:50%"])
    assert structure.composition["Si"] == 8
    assert structure.composition["P"] == 8

    # Small percentage resulting in at least 1 atom replaced
    structure2, _ = load_structure_file(poscar_path)
    apply_substitutions(structure2, ["Si:P:0.01%"])
    assert structure2.composition.reduced_formula == "SiP"


def test_apply_vacancies_errors(poscar_path):
    structure, _ = load_structure_file(poscar_path)

    # 1. Index-based vacancy errors
    with pytest.raises(ValueError, match="Invalid vacancy index rule"):
        apply_vacancies_by_index(structure, ["Si"])

    with pytest.raises(ValueError, match="Invalid vacancy index"):
        apply_vacancies_by_index(structure, ["Si:abc"])

    with pytest.raises(IndexError, match="out of range"):
        apply_vacancies_by_index(structure, ["Si:999"])

    # Warning path: actual element does not match vacancy element
    structure2, _ = load_structure_file(poscar_path)
    apply_substitutions(structure2, ["Si:P:0"])
    apply_vacancies_by_index(structure2, ["Si:0"])

    # 2. Count-based vacancy errors
    with pytest.raises(ValueError, match="Invalid vacancy count rule"):
        apply_vacancies_by_count(structure, ["Si"])

    with pytest.raises(ValueError, match="Invalid vacancy count"):
        apply_vacancies_by_count(structure, ["Si:abc"])

    with pytest.raises(ValueError, match="cannot be negative"):
        apply_vacancies_by_count(structure, ["Si:-1"])

    with pytest.raises(ValueError, match="exceeds available"):
        apply_vacancies_by_count(structure, ["Si:999"])

    # Warning path: matching elements not found
    apply_vacancies_by_count(structure, ["H:2"])


def test_apply_vacancies_random(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    # Scale to 64 atoms to trigger the random vacancy branch
    structure.make_supercell([2, 4, 4])
    assert len(structure) == 64

    # Apply count-based vacancies (e.g. remove 4 Si atoms)
    apply_vacancies_by_count(structure, ["Si:4"])
    assert len(structure) == 60


def test_generate_surface_slab_errors(poscar_path):
    structure, _ = load_structure_file(poscar_path)

    # Invalid Miller index
    with pytest.raises(ValueError):
        generate_surface_slab(structure, [0, 0, 0], 5.0, 10.0)


# CLI main integration tests
def test_cli_main_simple(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "-d", "2", "2", "2"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()
    structure, _ = load_structure_file(str(out_file))
    assert len(structure) == 16


def test_cli_main_conventional(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--conventional"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()
    structure, _ = load_structure_file(str(out_file))
    assert len(structure) == 8


def test_cli_main_min_dist(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--min-dist", "10.0"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()
    structure, _ = load_structure_file(str(out_file))
    assert len(structure) == 128


def test_cli_main_substitute_and_vacancy(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = [
        "cellify", "-i", poscar_path, "-o", str(out_file),
        "--substitute", "Si:P:0", "--vacancy-index", "Si:1"
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()
    structure, _ = load_structure_file(str(out_file))
    assert len(structure) == 1
    assert structure[0].specie.symbol == "P"


def test_cli_main_vacancy_deprecated_alias(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = [
        "cellify", "-i", poscar_path, "-o", str(out_file),
        "--vacancy", "Si:0"
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()
    structure, _ = load_structure_file(str(out_file))
    assert len(structure) == 1


def test_cli_main_vacancy_count(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = [
        "cellify", "-i", poscar_path, "-o", str(out_file),
        "--vacancy-count", "Si:1"
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()
    structure, _ = load_structure_file(str(out_file))
    assert len(structure) == 1


def test_cli_main_slab(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = [
        "cellify", "-i", poscar_path, "-o", str(out_file),
        "--slab", "1", "0", "0", "--thick", "5.0", "--vacuum", "10.0"
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()


def test_cli_main_matrix(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--matrix", "1 0 0 / 0 1 0 / 0 0 2"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()
    structure, _ = load_structure_file(str(out_file))
    assert len(structure) == 4


def test_cli_main_qe(qe_path, tmp_path):
    out_file = tmp_path / "qe_out.in"
    test_args = ["cellify", "-i", qe_path, "-o", str(out_file), "-d", "2", "2", "2"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    assert out_file.exists()


def test_cli_main_missing_file():
    test_args = ["cellify", "-i", "nonexistent_file_path.POSCAR"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_cli_main_invalid_matrix(poscar_path):
    test_args = ["cellify", "-i", poscar_path, "--matrix", "1 0 / 0 1"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_cli_main_corrupt_file(tmp_path):
    corrupt_file = tmp_path / "corrupt.POSCAR"
    corrupt_file.write_text("corrupt contents")
    test_args = ["cellify", "-i", str(corrupt_file)]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


# EspressoAdapter errors
def test_espresso_adapter_errors(tmp_path):
    from cellify.adapters.espresso import EspressoAdapter
    adapter = EspressoAdapter()

    # 1. File not found
    with pytest.raises(FileNotFoundError):
        adapter.read(str(tmp_path / "nonexistent_file_path.in"))

    # 2. Corrupt/parse error
    corrupt_file = tmp_path / "corrupt_espresso.in"
    corrupt_file.write_text("invalid contents")
    with pytest.raises(ValueError, match="Failed to parse structure from Quantum ESPRESSO file"):
        adapter.read(str(corrupt_file))


# CLI error paths for substitutions, vacancies, and slabs
def test_cli_substitute_error(poscar_path):
    test_args = ["cellify", "-i", poscar_path, "--substitute", "Si:P"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_cli_vacancy_index_error(poscar_path):
    test_args = ["cellify", "-i", poscar_path, "--vacancy-index", "Si"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_cli_vacancy_count_error(poscar_path):
    test_args = ["cellify", "-i", poscar_path, "--vacancy-count", "Si"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_cli_slab_error(poscar_path):
    test_args = ["cellify", "-i", poscar_path, "--slab", "0", "0", "0"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


# CLI default output path determination
def test_cli_default_output_poscar(poscar_path, tmp_path):
    import shutil
    temp_poscar = tmp_path / "POSCAR"
    shutil.copy(poscar_path, temp_poscar)

    test_args = ["cellify", "-i", str(temp_poscar), "-d", "2", "2", "2"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()

    expected_out = tmp_path / "POSCAR_supercell"
    assert expected_out.exists()


def test_cli_default_output_with_ext(poscar_path, tmp_path):
    import shutil
    temp_poscar_ext = tmp_path / "POSCAR.vasp"
    shutil.copy(poscar_path, temp_poscar_ext)

    test_args = ["cellify", "-i", str(temp_poscar_ext), "-d", "2", "2", "2"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()

    expected_out = tmp_path / "POSCAR_supercell.vasp"
    assert expected_out.exists()


# CLI --view tests
def test_cli_main_view_browser(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--view"]
    with patch("sys.argv", test_args):
        with patch("webbrowser.open") as mock_open:
            from cellify.cli import main
            main()
            mock_open.assert_called_once()
            args, _ = mock_open.call_args
            assert args[0].startswith("file://")
    assert out_file.exists()


def test_cli_main_view_browser_fallback_to_ase(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--view"]
    from unittest.mock import MagicMock
    with patch("sys.argv", test_args):
        with patch("webbrowser.open", side_effect=RuntimeError("No browser available")):
            with patch.dict("sys.modules", {"_tkinter": MagicMock()}):
                with patch("ase.visualize.view") as mock_view:
                    from cellify.cli import main
                    main()
                    mock_view.assert_called_once()
    assert out_file.exists()


def test_cli_main_view_browser_fallback_to_matplotlib(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--view"]
    import ase.visualize.plot
    with patch("sys.argv", test_args):
        with patch("webbrowser.open", side_effect=RuntimeError("No browser available")):
            with patch("ase.visualize.view", side_effect=RuntimeError("No display available")):
                with patch("matplotlib.pyplot.show") as mock_show, \
                     patch("ase.visualize.plot.plot_atoms") as mock_plot_atoms:
                    from cellify.cli import main
                    main()
                    mock_show.assert_called_once()
                    mock_plot_atoms.assert_called_once()
    assert out_file.exists()


def test_cli_main_view_browser_fallback_error(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--view"]
    with patch("sys.argv", test_args):
        with patch("webbrowser.open", side_effect=RuntimeError("No browser available")):
            with patch("ase.visualize.view", side_effect=RuntimeError("No display available")):
                with patch("matplotlib.pyplot.show", side_effect=RuntimeError("No display")) as mock_show:
                    from cellify.cli import main
                    with pytest.raises(SystemExit) as excinfo:
                        main()
                    assert excinfo.value.code == 1
                    mock_show.assert_called_once()


def test_open_browser_viewer(poscar_path):
    from cellify.viewer import open_browser_viewer
    from cellify.cli import load_structure_file
    structure, _ = load_structure_file(poscar_path)
    with patch("webbrowser.open") as mock_open:
        open_browser_viewer(structure)
        mock_open.assert_called_once()
        args, _ = mock_open.call_args
        assert args[0].startswith("file://")


def test_cli_extract_relaxation(tmp_path):
    # 1. Create a mock PWSCF output log file
    qe_relax_out = tmp_path / "relax.out"
    content = """
Program PWSCF
celldm(1)    8.0
number of atoms/cell = 2
number of atomic types = 2
crystal axes:
  a(1) = (  1.000000   0.000000   0.000000 )
  a(2) = (  0.000000   1.000000   0.000000 )
  a(3) = (  0.000000   0.000000   1.000000 )
positions (alat units)
    1           Si  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
    2           C   tau(   2) = (   0.2500000   0.2500000   0.2500000  )

!    total energy              =     -120.00000000 Ry

CELL_PARAMETERS (angstrom)
   4.5000  0.0000  0.0000
   0.0000  4.5000  0.0000
   0.0000  0.0000  4.5000

ATOMIC_POSITIONS (crystal)
  Si  0.100000  0.100000  0.100000
  C   0.300000  0.300000  0.300000

!    total energy              =     -120.10000000 Ry
"""
    qe_relax_out.write_text(content)

    # 2. Create a template QE input file
    qe_relax_in = tmp_path / "relax.in"
    template_content = """&CONTROL
  calculation = 'vc-relax'
  restart_mode = 'from_scratch'
  pseudo_dir = './'
  outdir = './'
/
&SYSTEM
  ibrav = 0
  nat = 2
  ntyp = 2
/
&ELECTRONS
/

ATOMIC_SPECIES
  Si  28.085  Si.UPF
  C   12.011  C.UPF

CELL_PARAMETERS angstrom
  4.0 0.0 0.0
  0.0 4.0 0.0
  0.0 0.0 4.0

ATOMIC_POSITIONS crystal
  Si  0.0 0.0 0.0
  C   0.25 0.25 0.25
"""
    qe_relax_in.write_text(template_content)

    # 3. Test running cellify with --template and --calc
    out_file = tmp_path / "scf.in"
    test_args = [
        "cellify", "-i", str(qe_relax_out),
        "--template", str(qe_relax_in),
        "-o", str(out_file),
        "--calc", "scf"
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()

    assert out_file.exists()
    out_content = out_file.read_text()
    assert "calculation = 'scf'" in out_content or 'calculation = "scf"' in out_content
    # Check that relaxed coordinates and cell parameters are written
    assert "4.500000" in out_content
    assert "0.100000" in out_content


def test_cli_extract_relaxation_missing_template(tmp_path):
    qe_relax_out = tmp_path / "relax.out"
    content = """
Program PWSCF
celldm(1)    8.0
number of atoms/cell = 2
number of atomic types = 2
crystal axes:
  a(1) = (  1.000000   0.000000   0.000000 )
  a(2) = (  0.000000   1.000000   0.000000 )
  a(3) = (  0.000000   0.000000   1.000000 )
positions (alat units)
    1           Si  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
    2           C   tau(   2) = (   0.2500000   0.2500000   0.2500000  )

!    total energy              =     -120.00000000 Ry

CELL_PARAMETERS (angstrom)
   4.5000  0.0000  0.0000
   0.0000  4.5000  0.0000
   0.0000  0.0000  4.5000

ATOMIC_POSITIONS (crystal)
  Si  0.100000  0.100000  0.100000
  C   0.300000  0.300000  0.300000

!    total energy              =     -120.10000000 Ry
"""
    qe_relax_out.write_text(content)

    out_file = tmp_path / "scf.in"
    test_args = [
        "cellify", "-i", str(qe_relax_out),
        "-o", str(out_file),
        "--calc", "scf"
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_cli_extract_relaxation_template_not_found(tmp_path):
    qe_relax_out = tmp_path / "relax.out"
    content = """
Program PWSCF
celldm(1)    8.0
number of atoms/cell = 2
number of atomic types = 2
crystal axes:
  a(1) = (  1.000000   0.000000   0.000000 )
  a(2) = (  0.000000   1.000000   0.000000 )
  a(3) = (  0.000000   0.000000   1.000000 )
positions (alat units)
    1           Si  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
    2           C   tau(   2) = (   0.2500000   0.2500000   0.2500000  )

!    total energy              =     -120.00000000 Ry

CELL_PARAMETERS (angstrom)
   4.5000  0.0000  0.0000
   0.0000  4.5000  0.0000
   0.0000  0.0000  4.5000

ATOMIC_POSITIONS (crystal)
  Si  0.100000  0.100000  0.100000
  C   0.300000  0.300000  0.300000

!    total energy              =     -120.10000000 Ry
"""
    qe_relax_out.write_text(content)
    test_args = [
        "cellify", "-i", str(qe_relax_out),
        "--template", "non_existent_template.in",
        "-o", str(tmp_path / "scf.in")
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_cli_extract_relaxation_template_corrupt(tmp_path):
    qe_relax_out = tmp_path / "relax.out"
    content = """
Program PWSCF
celldm(1)    8.0
number of atoms/cell = 2
number of atomic types = 2
crystal axes:
  a(1) = (  1.000000   0.000000   0.000000 )
  a(2) = (  0.000000   1.000000   0.000000 )
  a(3) = (  0.000000   0.000000   1.000000 )
positions (alat units)
    1           Si  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
    2           C   tau(   2) = (   0.2500000   0.2500000   0.2500000  )

!    total energy              =     -120.00000000 Ry

CELL_PARAMETERS (angstrom)
   4.5000  0.0000  0.0000
   0.0000  4.5000  0.0000
   0.0000  0.0000  4.5000

ATOMIC_POSITIONS (crystal)
  Si  0.100000  0.100000  0.100000
  C   0.300000  0.300000  0.300000

!    total energy              =     -120.10000000 Ry
"""
    qe_relax_out.write_text(content)
    corrupt_template = tmp_path / "corrupt.in"
    corrupt_template.write_text("corrupt contents")
    test_args = [
        "cellify", "-i", str(qe_relax_out),
        "--template", str(corrupt_template),
        "-o", str(tmp_path / "scf.in")
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        with pytest.raises(SystemExit) as excinfo:
            main()
        assert excinfo.value.code == 1


def test_espresso_adapter_write_error_out_mode(poscar_path):
    from cellify.adapters.espresso import EspressoAdapter
    structure, _ = load_structure_file(poscar_path)
    adapter = EspressoAdapter()
    with pytest.raises(ValueError, match="Cannot write a QE input file using a QE output log"):
        adapter.write("dummy.in", structure, {"mode": "espresso_out"})


def test_espresso_adapter_read_fallback(tmp_path):
    # Triggers is_output=False (due to &control comment) but espresso-in fails, falling back to espresso-out
    fake_qe_in = tmp_path / "fake_qe.in"
    content = """
# This is a comment containing &control namelist to trigger is_output=False
Program PWSCF
celldm(1)    8.0
number of atoms/cell = 2
number of atomic types = 2
crystal axes:
  a(1) = (  1.000000   0.000000   0.000000 )
  a(2) = (  0.000000   1.000000   0.000000 )
  a(3) = (  0.000000   0.000000   1.000000 )
positions (alat units)
    1           Si  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
    2           C   tau(   2) = (   0.2500000   0.2500000   0.2500000  )

!    total energy              =     -120.00000000 Ry
"""
    fake_qe_in.write_text(content)
    from cellify.adapters.espresso import EspressoAdapter
    adapter = EspressoAdapter()
    struct, meta = adapter.read(str(fake_qe_in))
    assert len(struct) == 2
    assert meta["mode"] == "espresso_out"


def test_espresso_adapter_read_fallback_input(tmp_path):
    # Triggers is_output=True (no &control) but espresso-out fails, falling back to espresso-in
    fake_qe_in = tmp_path / "fake_qe.in"
    content = """&SYSTEM
  ibrav = 0
  nat = 1
  ntyp = 1
/
ATOMIC_SPECIES
  Si  28.085  Si.UPF
CELL_PARAMETERS angstrom
  4.0 0.0 0.0
  0.0 4.0 0.0
  0.0 0.0 4.0
ATOMIC_POSITIONS crystal
  Si 0.0 0.0 0.0
"""
    fake_qe_in.write_text(content)
    from cellify.adapters.espresso import EspressoAdapter
    adapter = EspressoAdapter()
    struct, meta = adapter.read(str(fake_qe_in))
    assert len(struct) == 1
    assert meta["mode"] == "espresso_text_replace"


def test_cli_extract_relaxation_real_file(tmp_path):
    # Test with the real divacancy_relax_gamma.out file
    import shutil
    from pathlib import Path
    real_out_src = Path(__file__).parent / "divacancy_relax_gamma.out"
    real_out = tmp_path / "divacancy_relax_gamma.out"
    shutil.copy(real_out_src, real_out)

    # Create a template
    template_in = tmp_path / "template.in"
    template_content = """&CONTROL
  calculation = 'vc-relax'
/
&SYSTEM
  ibrav = 0
  nat = 2
  ntyp = 2
/
&ELECTRONS
/
ATOMIC_SPECIES
  Si  28.085  Si.UPF
  C   12.011  C.UPF
CELL_PARAMETERS angstrom
  8.0 0.0 0.0
  0.0 8.0 0.0
  0.0 0.0 8.0
ATOMIC_POSITIONS crystal
  Si 0.0 0.0 0.0
  C 0.25 0.25 0.25
"""
    template_in.write_text(template_content)

    out_file = tmp_path / "scf.in"
    test_args = [
        "cellify", "-i", str(real_out),
        "--template", str(template_in),
        "-o", str(out_file),
        "--calc", "scf"
    ]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()

    assert out_file.exists()
    out_content = out_file.read_text()
    assert "calculation = 'scf'" in out_content or 'calculation = "scf"' in out_content
    assert "nat = 62" in out_content


def test_espresso_adapter_write_default_template(tmp_path):
    from cellify.adapters.espresso import EspressoAdapter
    from pymatgen.core import Structure
    structure = Structure([[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]], ["Si"], [[0.0, 0.0, 0.0]])
    adapter = EspressoAdapter()
    out_file = tmp_path / "default_scf.in"
    adapter.write(str(out_file), structure, {})
    assert out_file.exists()
    out_content = out_file.read_text()
    assert "calculation = 'scf'" in out_content or 'calculation = "scf"' in out_content
    assert "nat = 1" in out_content


def test_cli_default_output_no_ext_poscar(tmp_path):
    import shutil
    from pathlib import Path
    poscar_src = Path(__file__).parent / "POSCAR"
    temp_poscar = tmp_path / "POSCAR"
    shutil.copy(poscar_src, temp_poscar)
    test_args = ["cellify", "-i", str(temp_poscar), "-d", "2", "2", "2"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()
    expected_out = tmp_path / "POSCAR_supercell"
    assert expected_out.exists()


def test_cli_save_file_error(poscar_path, tmp_path):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "-d", "2", "2", "2"]
    with patch("sys.argv", test_args):
        with patch("cellify.cli.save_structure_file", side_effect=RuntimeError("Save failed")) as mock_save:
            from cellify.cli import main
            with pytest.raises(SystemExit) as excinfo:
                main()
            assert excinfo.value.code == 1
            mock_save.assert_called_once()


def test_get_adapter_exception(tmp_path):
    from cellify.adapters import get_adapter
    from cellify.adapters.standard import StandardAdapter
    # A directory path will raise IsADirectoryError when read, triggering the content check exception block
    adapter = get_adapter(str(tmp_path))
    assert isinstance(adapter, StandardAdapter)

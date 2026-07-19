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


def test_cli_main_formula_query(capsys, tmp_path):
    mock_mp_response = {
        "data": [
            {
                "id": "mp-165",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {
                        "gga_gga+u": {"energy_above_hull": 0.0}
                    },
                    "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                    "species_at_sites": ["Si", "Si"],
                },
            },
            {
                "id": "mp-9999",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"energy_above_hull": "N/A"},
                    "lattice_vectors": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0]],
                    "species_at_sites": ["Si"],
                },
            },
            {
                "id": "mp-1000",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"energy_above_hull": 0.05},
                },
            },
            {
                "id": "mp-9998",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": "invalid-stability-type",
                },
            },
            {
                "id": "mp-2000",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"energy_above_hull": 0.10},
                },
            },
            {
                "id": "mp-3000",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"energy_above_hull": 0.15},
                },
            },
        ]
    }
    mock_cod_response = {
        "data": [
            {
                "id": "1526655",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_cod_sg": "F d -3 m :1",
                    "_cod_vol": 160.0,
                    "_cod_a": 5.43,
                    "_cod_b": 5.43,
                    "_cod_c": 5.43,
                    "_cod_commonname": "Silicon",
                    "_cod_chemname": "Silicon",
                    "_cod_mineral": "Silicon mineral",
                },
            },
            {
                "id": "1526656",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_cod_sg": "F d -3 m",
                    "_cod_vol": "invalid-vol",
                    "_cod_a": "invalid-a",
                    "_cod_b": "invalid-b",
                    "_cod_c": "invalid-c",
                },
            },
            {
                "id": "1526657",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_cod_chemname": "Silicon Chem",
                },
            },
        ]
    }

    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code

            def json(self):
                return self.json_data

        if "materialsproject.org" in url:
            return MockResponse(mock_mp_response, 200)
        if "crystallography.net" in url:
            return MockResponse(mock_cod_response, 200)
        return MockResponse({}, 404)

    out_file = str(tmp_path / "Si_supercell.cif")
    with patch("requests.get", side_effect=mock_get):
        test_args = ["cellify", "-i", "Si", "--select", "1", "-o", out_file]
        with patch("sys.argv", test_args):
            from cellify.cli import main
            main()

    captured = capsys.readouterr()
    assert "Found 6 structures in Materials Project" in captured.out
    assert "Warning: Only the top 5 most relevant structures are shown. There are 1 more structures in Materials Project." in captured.out
    assert "mp-165" in captured.out
    assert "Space Group: R-3m" in captured.out
    assert "Volume: 157.46 A^3" in captured.out
    assert "Lattice: a=5.40, b=5.40, c=5.40 A" in captured.out
    assert "0.0000 eV/atom [Stable ★]" in captured.out
    assert "mp-9999" in captured.out
    assert "N/A eV/atom" in captured.out
    assert "mp-1000" in captured.out
    assert "0.0500 eV/atom" in captured.out
    assert "Querying Crystallography Open Database (COD) OPTIMADE" in captured.out
    assert "Found 3 structures in Crystallography Open Database (COD)" in captured.out
    assert "1526655" in captured.out
    assert "Space Group: F d -3 m :1" in captured.out
    assert "Volume: 160.00 A^3" in captured.out
    assert "Lattice: a=5.43, b=5.43, c=5.43 A" in captured.out
    assert "Name: Silicon, Silicon mineral" in captured.out
    assert "1526656" in captured.out
    assert "Volume: invalid-vol A^3" in captured.out
    assert "Lattice: a=invalid-a, b=invalid-b, c=invalid-c A" in captured.out


def test_cli_main_formula_query_sg_error(capsys, tmp_path):
    mock_mp_response = {
        "data": [
            {
                "id": "mp-165",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"energy_above_hull": 0.0},
                    "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                    "species_at_sites": ["Si", "Si"],
                },
            }
        ]
    }

    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code

            def json(self):
                return self.json_data

        return MockResponse(mock_mp_response, 200)

    out_file = str(tmp_path / "Si_supercell.cif")
    with patch("requests.get", side_effect=mock_get):
        with patch("cellify.optimade.SpacegroupAnalyzer", side_effect=RuntimeError("spglib error")):
            test_args = ["cellify", "-i", "Si", "--select", "1", "-o", out_file]
            with patch("sys.argv", test_args):
                from cellify.cli import main
                main()

    captured = capsys.readouterr()
    assert "Querying Materials Project OPTIMADE" in captured.out
    assert "Space Group" not in captured.out
    assert "Volume: 157.46 A^3" in captured.out


def test_cli_main_formula_query_status_non_200(capsys):
    def mock_get_error(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, status_code):
                self.status_code = status_code

            def json(self):
                return {}

        return MockResponse(500)

    with patch("requests.get", side_effect=mock_get_error):
        test_args = ["cellify", "-i", "Si"]
        with patch("sys.argv", test_args):
            from cellify.cli import main
            with pytest.raises(SystemExit) as excinfo:
                main()
            assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Found 0 structures in Materials Project" in captured.out
    assert "Found 0 structures in Crystallography Open Database (COD)" in captured.out
    assert "Error: No structures found for formula 'Si'." in captured.err


def test_cli_main_formula_query_error(capsys):
    def mock_get_error(url, *args, **kwargs):
        raise RuntimeError("Connection timed out")

    with patch("requests.get", side_effect=mock_get_error):
        test_args = ["cellify", "-i", "Si"]
        with patch("sys.argv", test_args):
            from cellify.cli import main
            with pytest.raises(SystemExit) as excinfo:
                main()
            assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Found 0 structures in Materials Project" in captured.out
    assert "Found 0 structures in Crystallography Open Database (COD)" in captured.out
    assert "Error: No structures found for formula 'Si'." in captured.err


def test_cli_main_formula_query_invalid_formula(capsys):
    def mock_get_error(url, *args, **kwargs):
        raise RuntimeError("early exit")

    with patch("requests.get", side_effect=mock_get_error):
        test_args = ["cellify", "-i", "invalid-formula-123!"]
        with patch("sys.argv", test_args):
            from cellify.cli import main
            with pytest.raises(SystemExit) as excinfo:
                main()
            assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Found 0 structures in Materials Project" in captured.out
    assert "Found 0 structures in Crystallography Open Database (COD)" in captured.out
    assert "Error: No structures found for formula 'invalid-formula-123!'." in captured.err


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


def test_cli_show_indices(poscar_path, tmp_path, capsys):
    out_file = tmp_path / "POSCAR_out"
    test_args = ["cellify", "-i", poscar_path, "-o", str(out_file), "--show-indices"]
    with patch("sys.argv", test_args):
        from cellify.cli import main
        main()

    captured = capsys.readouterr()
    assert "Absolute Atomic Indices & Coordinates:" in captured.out
    assert "Index" in captured.out
    assert "Fractional Coordinates (a, b, c)" in captured.out
    assert "Cartesian (x, y, z)" in captured.out
    assert "0      Si" in captured.out
    assert "Total: 2 atoms" in captured.out


def test_granular_optimade_query_and_download():
    # 1. Test imports
    from cellify.optimade import (
        query_formula_structures,
        format_formula_structures,
        download_structure_from_entry,
    )

    # Mock data for Materials Project (MP)
    mock_mp_data = [
        {
            "id": "mp-123",
            "attributes": {
                "chemical_formula_descriptive": "Si",
                "_mp_stability": {"energy_above_hull": 0.0},
                "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                "species_at_sites": ["Si", "Si"],
            },
        }
    ]

    # Mock data for Crystallography Open Database (COD)
    mock_cod_data = [
        {
            "id": "1000000",
            "attributes": {
                "chemical_formula_descriptive": "Si",
                "_cod_sg": "Fd-3m",
                "_cod_vol": 160.0,
                "_cod_a": 5.43,
                "_cod_b": 5.43,
                "_cod_c": 5.43,
                "_cod_commonname": "Silicon",
            },
        }
    ]

    # Mock requests.get
    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, text_or_json, status_code):
                self.text_or_json = text_or_json
                self.status_code = status_code

            def json(self):
                return self.text_or_json

            @property
            def text(self):
                return self.text_or_json

        if "materialsproject.org" in url:
            return MockResponse({"data": mock_mp_data}, 200)
        elif "crystallography.net/cod/optimade" in url:
            return MockResponse({"data": mock_cod_data}, 200)
        elif "crystallography.net/cod/1000000.cif" in url:
            # Return dummy CIF text
            cif_text = """data_global
_cell_length_a 5.43
_cell_length_b 5.43
_cell_length_c 5.43
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
loop_
_space_group_symop_operation_xyz
'x,y,z'
loop_
_atom_site_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si 0 0 0
"""
            return MockResponse(cif_text, 200)
        return MockResponse("Not Found", 404)

    with patch("requests.get", side_effect=mock_get):
        # 2. Test query
        mp_res, cod_res = query_formula_structures("Si")
        assert len(mp_res) == 1
        assert mp_res[0]["id"] == "mp-123"
        assert len(cod_res) == 1
        assert cod_res[0]["id"] == "1000000"

        # 3. Test formatting
        formatted_str, selection_map = format_formula_structures(mp_res, cod_res, "Si")
        assert "Found 1 structures in Materials Project" in formatted_str
        assert "Found 1 structures in Crystallography Open Database (COD)" in formatted_str
        assert "[1] ID: mp-123" in formatted_str
        assert "[2] ID: 1000000" in formatted_str

        assert 1 in selection_map
        assert selection_map[1] == ("Materials Project", mp_res[0])
        assert 2 in selection_map
        assert selection_map[2] == ("Crystallography Open Database (COD)", cod_res[0])

        # 4. Test download structure
        struct_mp = download_structure_from_entry("Materials Project", mp_res[0])
        assert isinstance(struct_mp, Structure)
        assert len(struct_mp) == 2
        assert set(struct_mp.symbol_set) == {"Si"}

        struct_cod = download_structure_from_entry("Crystallography Open Database (COD)", cod_res[0])
        assert isinstance(struct_cod, Structure)
        assert len(struct_cod) == 1
        assert set(struct_cod.symbol_set) == {"Si"}

        # 5. Test download structure with unsupported database name
        import pytest
        with pytest.raises(ValueError, match="Unsupported database name: UnknownDB"):
            download_structure_from_entry("UnknownDB", mp_res[0])


def test_select_and_download_structure():
    from cellify.optimade import select_and_download_structure, SelectionError
    from unittest.mock import patch, MagicMock
    from pymatgen.core import Structure

    # Dummy Structure
    dummy_structure = Structure([[1,0,0],[0,1,0],[0,0,1]], ["Si"], [[0,0,0]])
    dummy_entry = {"id": "mp-1"}

    mock_mp = [dummy_entry]
    mock_cod = []

    # 1. Test success with select
    with patch("cellify.optimade.query_formula_structures", return_value=(mock_mp, mock_cod)) as mock_query, \
         patch("cellify.optimade.format_formula_structures", return_value=("summary", {1: ("Materials Project", dummy_entry)})) as mock_format, \
         patch("cellify.optimade.download_structure_from_entry", return_value=dummy_structure) as mock_download:

        struct, db, entry, summary = select_and_download_structure("Si", select=1)
        assert struct == dummy_structure
        assert db == "Materials Project"
        assert entry == dummy_entry
        assert summary == "summary"
        mock_query.assert_called_once_with("Si")
        mock_format.assert_called_once_with(mock_mp, mock_cod, "Si")
        mock_download.assert_called_once_with("Materials Project", dummy_entry)

    # 2. Test success with interactive callback
    mock_prompt = MagicMock(return_value=1)
    with patch("cellify.optimade.query_formula_structures", return_value=(mock_mp, mock_cod)), \
         patch("cellify.optimade.format_formula_structures", return_value=("summary", {1: ("Materials Project", dummy_entry)})), \
         patch("cellify.optimade.download_structure_from_entry", return_value=dummy_structure):

        struct, db, entry, summary = select_and_download_structure("Si", interactive_prompt=mock_prompt)
        assert struct == dummy_structure
        assert db == "Materials Project"
        assert entry == dummy_entry
        assert summary == "summary"
        mock_prompt.assert_called_once_with("summary", 1)

    # 3. Test no structures found
    with patch("cellify.optimade.query_formula_structures", return_value=([], [])), \
         patch("cellify.optimade.format_formula_structures", return_value=("empty_summary", {})):
        with pytest.raises(SelectionError, match="No structures found for formula 'Si'.") as excinfo:
            select_and_download_structure("Si", select=1)
        assert excinfo.value.summary == "empty_summary"

    # 4. Test no select or interactive_prompt provided
    with patch("cellify.optimade.query_formula_structures", return_value=(mock_mp, mock_cod)), \
         patch("cellify.optimade.format_formula_structures", return_value=("summary", {1: ("Materials Project", dummy_entry)})):
        with pytest.raises(SelectionError, match="Chemical formula specified, but the session is non-interactive") as excinfo:
            select_and_download_structure("Si")
        assert excinfo.value.summary == "summary"

    # 5. Test selection index out of range
    with patch("cellify.optimade.query_formula_structures", return_value=(mock_mp, mock_cod)), \
         patch("cellify.optimade.format_formula_structures", return_value=("summary", {1: ("Materials Project", dummy_entry)})):
        with pytest.raises(SelectionError, match="Selection index 2 is out of range.") as excinfo:
            select_and_download_structure("Si", select=2)
        assert excinfo.value.summary == "summary"


def test_cli_main_formula_query_select_out_of_range(capsys):
    mock_mp_response = {
        "data": [
            {
                "id": "mp-165",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"gga_gga+u": {"energy_above_hull": 0.0}},
                    "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                    "species_at_sites": ["Si", "Si"],
                },
            }
        ]
    }

    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code
            def json(self):
                return self.json_data
        return MockResponse(mock_mp_response, 200)

    with patch("requests.get", side_effect=mock_get):
        test_args = ["cellify", "-i", "Si", "--select", "99"]
        with patch("sys.argv", test_args):
            from cellify.cli import main
            with pytest.raises(SystemExit) as excinfo:
                main()
            assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Error: Selection index 99 is out of range." in captured.err


def test_cli_main_formula_query_non_interactive_no_select(capsys):
    mock_mp_response = {
        "data": [
            {
                "id": "mp-165",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"gga_gga+u": {"energy_above_hull": 0.0}},
                    "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                    "species_at_sites": ["Si", "Si"],
                },
            }
        ]
    }

    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code
            def json(self):
                return self.json_data
        return MockResponse(mock_mp_response, 200)

    with patch("requests.get", side_effect=mock_get):
        test_args = ["cellify", "-i", "Si"]
        with patch("sys.argv", test_args):
            from cellify.cli import main
            with pytest.raises(SystemExit) as excinfo:
                main()
            assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Error: Chemical formula specified, but the session is non-interactive" in captured.err


def test_cli_main_formula_query_interactive_select(capsys, tmp_path):
    mock_mp_response = {
        "data": [
            {
                "id": "mp-165",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"gga_gga+u": {"energy_above_hull": 0.0}},
                    "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                    "species_at_sites": ["Si", "Si"],
                },
            }
        ]
    }

    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code
            def json(self):
                return self.json_data
        return MockResponse(mock_mp_response, 200)

    out_file = str(tmp_path / "Si_supercell.cif")
    with patch("requests.get", side_effect=mock_get):
        with patch("sys.stdin.isatty", return_value=True):
            with patch("builtins.input", return_value="1"):
                test_args = ["cellify", "-i", "Si", "-o", out_file]
                with patch("sys.argv", test_args):
                    from cellify.cli import main
                    main()

    captured = capsys.readouterr()
    assert "Downloading structure from Materials Project (ID: mp-165)..." in captured.out
    assert os.path.exists(out_file)


def test_cli_main_formula_query_interactive_select_invalid(capsys):
    mock_mp_response = {
        "data": [
            {
                "id": "mp-165",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"gga_gga+u": {"energy_above_hull": 0.0}},
                    "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                    "species_at_sites": ["Si", "Si"],
                },
            }
        ]
    }

    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code
            def json(self):
                return self.json_data
        return MockResponse(mock_mp_response, 200)

    # test invalid choice input
    with patch("requests.get", side_effect=mock_get):
        with patch("sys.stdin.isatty", return_value=True):
            with patch("builtins.input", return_value="invalid_choice"):
                test_args = ["cellify", "-i", "Si"]
                with patch("sys.argv", test_args):
                    from cellify.cli import main
                    with pytest.raises(SystemExit) as excinfo:
                        main()
                    assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Error: Invalid selection." in captured.err

    # test out of range choice input
    with patch("requests.get", side_effect=mock_get):
        with patch("sys.stdin.isatty", return_value=True):
            with patch("builtins.input", return_value="5"):
                test_args = ["cellify", "-i", "Si"]
                with patch("sys.argv", test_args):
                    from cellify.cli import main
                    with pytest.raises(SystemExit) as excinfo:
                        main()
                    assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Error: Invalid selection." in captured.err

    # test empty choice input
    with patch("requests.get", side_effect=mock_get):
        with patch("sys.stdin.isatty", return_value=True):
            with patch("builtins.input", return_value=""):
                test_args = ["cellify", "-i", "Si"]
                with patch("sys.argv", test_args):
                    from cellify.cli import main
                    with pytest.raises(SystemExit) as excinfo:
                        main()
                    assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Error: Invalid selection." in captured.err


def test_cli_main_formula_query_interactive_keyboard_interrupt(capsys):
    mock_mp_response = {
        "data": [
            {
                "id": "mp-165",
                "attributes": {
                    "chemical_formula_descriptive": "Si",
                    "_mp_stability": {"gga_gga+u": {"energy_above_hull": 0.0}},
                    "lattice_vectors": [[5.4, 0.0, 0.0], [0.0, 5.4, 0.0], [0.0, 0.0, 5.4]],
                    "cartesian_site_positions": [[0.0, 0.0, 0.0], [1.35, 1.35, 1.35]],
                    "species_at_sites": ["Si", "Si"],
                },
            }
        ]
    }
    def mock_get(url, *args, **kwargs):
        class MockResponse:
            def __init__(self, json_data, status_code):
                self.json_data = json_data
                self.status_code = status_code
            def json(self):
                return self.json_data
        return MockResponse(mock_mp_response, 200)

    with patch("requests.get", side_effect=mock_get):
        with patch("sys.stdin.isatty", return_value=True):
            with patch("builtins.input", side_effect=KeyboardInterrupt):
                test_args = ["cellify", "-i", "Si"]
                with patch("sys.argv", test_args):
                    from cellify.cli import main
                    with pytest.raises(SystemExit) as excinfo:
                        main()
                    assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Error: Invalid selection." in captured.err


def test_cli_main_formula_query_general_download_error(capsys):
    with patch("cellify.optimade.select_and_download_structure", side_effect=RuntimeError("Some network error")):
        test_args = ["cellify", "-i", "Si", "--select", "1"]
        with patch("sys.argv", test_args):
            from cellify.cli import main
            with pytest.raises(SystemExit) as excinfo:
                main()
            assert excinfo.value.code == 1

    captured = capsys.readouterr()
    assert "Error downloading structure: Some network error" in captured.err


def test_optimade_cod_data_more_than_5_warning():
    from cellify.optimade import format_formula_structures
    mp_data = []
    cod_data = [
        {"id": f"cod-{i}", "attributes": {"chemical_formula_descriptive": "Si"}}
        for i in range(1, 7)
    ]
    summary, selection_map = format_formula_structures(mp_data, cod_data, "Si")
    assert "Warning: Only the top 5 most relevant structures are shown." in summary
    assert "There are 1 more structures in Crystallography Open Database (COD)." in summary
    assert len(selection_map) == 5


def test_optimade_download_structure_failures():
    from cellify.optimade import download_structure_from_entry
    import pytest
    from unittest.mock import patch

    with pytest.raises(ValueError, match="Failed to parse structure from Materials Project attributes."):
        download_structure_from_entry("Materials Project", {"id": "mp-1", "attributes": {}})

    with pytest.raises(ValueError, match="COD entry has no ID."):
        download_structure_from_entry("Crystallography Open Database (COD)", {"attributes": {}})

    class MockResponse:
        def __init__(self, status_code):
            self.status_code = status_code
    with patch("requests.get", return_value=MockResponse(404)):
        with pytest.raises(RuntimeError, match="Failed to download CIF from COD: status 404"):
            download_structure_from_entry("Crystallography Open Database (COD)", {"id": "12345"})


def test_animate_print(capsys):
    from cellify.cli import animate_print
    import time

    # Test isatty=False case
    with patch("sys.stdout.isatty", return_value=False):
        animate_print("line1\nline2", delay=0.01)
        captured = capsys.readouterr()
        assert captured.out == "line1\nline2\n"

    # Test isatty=True case
    start_time = time.time()
    with patch("sys.stdout.isatty", return_value=True):
        animate_print("line1\nline2", delay=0.05)
        duration = time.time() - start_time
        captured = capsys.readouterr()
        # Since there are 2 lines, it should sleep at least 0.05s * 2 = 0.1s
        assert duration >= 0.08
        assert captured.out == "line1\nline2\n"

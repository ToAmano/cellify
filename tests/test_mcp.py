"""
Unit and integration tests for the cellify MCP server tool.
"""

import os
import shutil
from unittest.mock import patch

import pytest

pytest.importorskip("mcp")

from pymatgen.core import Structure

import cellify.mcp as mcp_mod
from cellify.core import load_structure_file
from cellify.mcp import cellify


@pytest.fixture
def poscar_path() -> str:
    """
    Returns the path to the test POSCAR file.
    """
    return os.path.join(os.path.dirname(__file__), "POSCAR")


@pytest.fixture
def qe_path() -> str:
    """
    Returns the path to the test qe.in file.
    """
    return os.path.join(os.path.dirname(__file__), "qe.in")


def test_mcp_show_indices(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool with show_indices=True.
    """
    # Copy POSCAR to tmp_path to avoid contaminating test folder
    tmp_poscar = os.path.join(tmp_path, "POSCAR")
    shutil.copy(poscar_path, tmp_poscar)

    res: str = cellify(tmp_poscar, show_indices=True)
    assert "Formula: Si" in res
    assert "Number of atoms: 2" in res
    assert "Absolute Atomic Indices & Coordinates" in res
    assert "0      Si" in res
    assert "Saving final structure to" in res
    assert os.path.exists(os.path.join(tmp_path, "POSCAR_supercell"))


def test_mcp_conventional(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool conventional conversion.
    """
    out_file = os.path.join(tmp_path, "POSCAR_conv")
    res_msg: str = cellify(poscar_path, conventional=True, output_path=out_file)
    assert "Saving final structure to" in res_msg
    assert os.path.exists(out_file)


def test_mcp_supercell(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool supercell scaling.
    """
    # Diagonal scaling output path
    out_file_diag = os.path.join(tmp_path, "POSCAR_super_diag")
    res_msg_diag: str = cellify(poscar_path, dim=[2, 2, 2], output_path=out_file_diag)
    assert "diagonal scaling" in res_msg_diag
    assert "Number of atoms: 16" in res_msg_diag
    assert os.path.exists(out_file_diag)

    # Matrix scaling output path
    out_file_mat = os.path.join(tmp_path, "POSCAR_super_mat")
    res_msg_mat: str = cellify(
        poscar_path, matrix="2 0 0 / 0 2 0 / 0 0 2", output_path=out_file_mat
    )
    assert "with matrix" in res_msg_mat
    assert "Number of atoms: 16" in res_msg_mat
    assert os.path.exists(out_file_mat)

    # min_dist scaling output path
    out_file_dist = os.path.join(tmp_path, "POSCAR_super_dist")
    res_msg_dist: str = cellify(poscar_path, min_dist=12.0, output_path=out_file_dist)
    assert "minimum distance" in res_msg_dist
    assert os.path.exists(out_file_dist)


def test_mcp_defect(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool defect application (doping, vacancies).
    """
    out_file = os.path.join(tmp_path, "POSCAR_defect")
    # Scaling to 2x2x2 (16 atoms), replace Si at index 0 with Ge, remove Si at index 1
    res_msg: str = cellify(
        poscar_path,
        substitute=["Si:Ge:0"],
        vacancy_index=["Si:1"],
        dim=[2, 2, 2],
        output_path=out_file,
    )
    assert "Replaced site 0" in res_msg
    assert "Removed site 1" in res_msg
    assert "Number of atoms: 15" in res_msg
    assert os.path.exists(out_file)

    # Re-load and verify composition contains Ge
    struct_new: Structure
    struct_new, _ = load_structure_file(out_file)
    assert struct_new.composition.reduced_formula == "Si14Ge"


def test_mcp_slab(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool slab generation.
    """
    out_file = os.path.join(tmp_path, "POSCAR_slab")
    res_msg: str = cellify(
        poscar_path, slab=[1, 1, 1], thick=4.0, vacuum=10.0, output_path=out_file
    )
    assert "Generating slab model for Miller indices" in res_msg
    assert os.path.exists(out_file)


def test_mcp_errors(poscar_path: str) -> None:
    """
    Tests error handling for invalid files or options.
    """
    # Non-existent file error
    res_info: str = cellify("nonexistent_file.POSCAR")
    assert "Error:" in res_info

    # Invalid scaling dim format
    res_scale_err: str = cellify(poscar_path, dim=[1, 2])
    assert "Error" in res_scale_err


def test_mcp_load_error(tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests error handling when the structure file exists but is empty/invalid.
    """
    empty_file = os.path.join(tmp_path, "empty_POSCAR")
    with open(empty_file, "w", encoding="utf-8") as f:
        f.write("")
    res_msg: str = cellify(empty_file)
    assert "Error" in res_msg


def test_mcp_validation_error(poscar_path: str) -> None:
    """
    Tests error handling when process_template_and_validation fails.
    """
    res_msg: str = cellify(poscar_path, template="nonexistent_template_file")
    assert "Error:" in res_msg


def test_mcp_save_error(poscar_path: str) -> None:
    """
    Tests error handling when saving the output file fails.
    """
    res_msg: str = cellify(poscar_path, output_path="/nonexistent_dir/POSCAR_out")
    assert "Error saving file" in res_msg or "Error" in res_msg


def test_mcp_main() -> None:
    """
    Tests running the MCP server main loop.
    """
    with patch("cellify.mcp.mcp.run") as mock_run:
        mcp_mod.main()
        mock_run.assert_called_once()


def test_mcp_formula_query() -> None:
    """
    Tests the cellify tool formula query fallback.
    """
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
            }
        ]
    }
    mock_cod_response = {
        "data": []
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

    with patch("requests.get", side_effect=mock_get):
        res: str = cellify("Si")
        assert "Querying Materials Project OPTIMADE" in res
        assert "Found 1 structures in Materials Project" in res
        assert "mp-165" in res
        assert "Space Group: R-3m" in res
        assert "Volume: 157.46 A^3" in res


def test_mcp_formula_query_error() -> None:
    """
    Tests error handling during the formula query inside the cellify tool.
    """
    def mock_get_error(url, *args, **kwargs):
        raise RuntimeError("Connection failed")

    with patch("requests.get", side_effect=mock_get_error):
        res: str = cellify("Si")
        assert "Error querying Materials Project: Connection failed" in res
        assert "Error querying Crystallography Open Database (COD): Connection failed" in res


def test_mcp_formula_query_raise_error() -> None:
    """
    Tests error handling when retrieve_cif_by_formula raises a direct exception.
    """
    with patch("cellify.optimade.retrieve_cif_by_formula", side_effect=RuntimeError("Generic error")):
        res: str = cellify("Si")
        assert "Error querying formula: Generic error" in res


def test_mcp_process_error(poscar_path: str) -> None:
    """
    Tests error handling when run_cellify_pipeline raises an exception.
    """
    with patch("cellify.mcp.run_cellify_pipeline", side_effect=RuntimeError("Process error")):
        res: str = cellify(poscar_path)
        assert "Error processing structure: Process error" in res


def test_mcp_formula_query_with_select(tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool formula query fallback with --select.
    """
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
            }
        ]
    }
    mock_cod_response = {
        "data": []
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

    # Output file path
    out_file = os.path.join(tmp_path, "POSCAR_select")

    with patch("requests.get", side_effect=mock_get):
        # 1. Valid select (1)
        res: str = cellify("Si", select=1, output_path=out_file)
        assert "Downloading structure from Materials Project (ID: mp-165)..." in res
        assert "Saving final structure to" in res
        assert os.path.exists(out_file)

        # 2. Out of range select (99)
        res_err: str = cellify("Si", select=99)
        assert "Error: Selection index 99 is out of range." in res_err

    # 3. Download exception
    with patch("requests.get", side_effect=mock_get), \
         patch("cellify.optimade.download_structure_from_entry", side_effect=RuntimeError("Download failed")):
        res_dl_err: str = cellify("Si", select=1)
        assert "Error downloading structure: Download failed" in res_dl_err

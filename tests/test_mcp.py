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
    res_info: str = cellify("nonexistent_file")
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

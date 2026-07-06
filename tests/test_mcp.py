"""
Unit and integration tests for the cellify MCP server tool.
"""

import os
from unittest.mock import patch

import pytest
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


def test_mcp_show_indices(poscar_path: str) -> None:
    """
    Tests the cellify tool with show_indices=True.
    """
    res: str = cellify(poscar_path, show_indices=True)
    assert "Formula: Si" in res
    assert "Number of atoms: 2" in res
    assert "Absolute Atomic Indices & Coordinates" in res
    assert "0      Si" in res


def test_mcp_conventional(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool conventional conversion.
    """
    # Test returning string content
    res_str: str = cellify(poscar_path, conventional=True)
    assert "=== STRUCTURE CONTENT ===" in res_str
    assert "Si" in res_str

    # Test writing to file
    out_file = os.path.join(tmp_path, "POSCAR_conv")
    res_msg: str = cellify(poscar_path, conventional=True, output_path=out_file)
    assert "Successfully saved final structure to" in res_msg
    assert os.path.exists(out_file)


def test_mcp_supercell(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify tool supercell scaling.
    """
    # Diagonal scaling text output
    res_str: str = cellify(poscar_path, dim="2 2 2")
    assert "=== STRUCTURE CONTENT ===" in res_str

    # Diagonal scaling output path
    out_file_diag = os.path.join(tmp_path, "POSCAR_super_diag")
    res_msg_diag: str = cellify(poscar_path, dim="2 2 2", output_path=out_file_diag)
    assert "diagonal scaling" in res_msg_diag
    assert "Number of atoms: 16" in res_msg_diag
    assert os.path.exists(out_file_diag)

    # Matrix scaling output path
    out_file_mat = os.path.join(tmp_path, "POSCAR_super_mat")
    res_msg_mat: str = cellify(
        poscar_path, dim="2 0 0 / 0 2 0 / 0 0 2", output_path=out_file_mat
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
        dim="2 2 2",
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
        poscar_path, slab="1 1 1", thick=4.0, vacuum=10.0, output_path=out_file
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

    # Slab without thickness/vacuum
    res_slab_err: str = cellify(poscar_path, slab="1 1 1")
    assert "Error:" in res_slab_err

    # Invalid scaling dim format
    res_scale_err: str = cellify(poscar_path, dim="1 2")
    assert "Error" in res_scale_err


def test_mcp_main() -> None:
    """
    Tests running the MCP server main loop.
    """
    with patch("cellify.mcp.mcp.run") as mock_run:
        mcp_mod.main()
        mock_run.assert_called_once()

"""
Unit and integration tests for the cellify MCP server tools.
"""

import os
from unittest.mock import patch

import pytest
from pymatgen.core import Structure

import cellify.mcp
from cellify.mcp import (
    cellify_conventional,
    cellify_defect,
    cellify_info,
    cellify_slab,
    cellify_supercell,
)


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


def test_mcp_info(poscar_path: str) -> None:
    """
    Tests the cellify_info tool.
    """
    res: str = cellify_info(poscar_path)
    assert "Structure Info for POSCAR" in res
    assert "Formula: Si" in res
    assert "Number of atoms: 2" in res
    assert "Absolute Atomic Indices & Coordinates" in res
    assert "0      Si" in res


def test_mcp_conventional(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify_conventional tool.
    """
    # Test returning string content
    res_str: str = cellify_conventional(poscar_path)
    assert "Si" in res_str

    # Test writing to file
    out_file = os.path.join(tmp_path, "POSCAR_conv")
    res_msg: str = cellify_conventional(poscar_path, output_path=out_file)
    assert "Successfully converted" in res_msg
    assert os.path.exists(out_file)


def test_mcp_supercell(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify_supercell tool with various scaling methods.
    """
    # Test returning string content
    res_str: str = cellify_supercell(poscar_path, dim="2 2 2")
    assert "Si" in res_str

    # Diagonal scaling output path
    out_file_diag = os.path.join(tmp_path, "POSCAR_super_diag")
    res_msg_diag: str = cellify_supercell(
        poscar_path, dim="2 2 2", output_path=out_file_diag
    )
    assert "Successfully generated supercell" in res_msg_diag
    assert "Number of atoms: 16" in res_msg_diag
    assert os.path.exists(out_file_diag)

    # Matrix scaling output path
    out_file_mat = os.path.join(tmp_path, "POSCAR_super_mat")
    res_msg_mat: str = cellify_supercell(
        poscar_path, dim="2 0 0 / 0 2 0 / 0 0 2", output_path=out_file_mat
    )
    assert "Successfully generated supercell" in res_msg_mat
    assert "Number of atoms: 16" in res_msg_mat
    assert os.path.exists(out_file_mat)

    # min_dist scaling output path
    out_file_dist = os.path.join(tmp_path, "POSCAR_super_dist")
    res_msg_dist: str = cellify_supercell(
        poscar_path, min_dist=12.0, output_path=out_file_dist
    )
    assert "Successfully generated supercell" in res_msg_dist
    assert os.path.exists(out_file_dist)


def test_mcp_defect(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify_defect tool.
    """
    out_file = os.path.join(tmp_path, "POSCAR_defect")
    # Scaling to 2x2x2 (16 atoms), replace Si at index 0 with Ge, remove Si at index 1
    res_msg: str = cellify_defect(
        poscar_path,
        substitute=["Si:Ge:0"],
        vacancy_index=["Si:1"],
        dim="2 2 2",
        output_path=out_file,
    )
    assert "Successfully applied defects" in res_msg
    assert "Number of atoms: 15" in res_msg
    assert os.path.exists(out_file)

    # Re-load and verify composition contains Ge
    struct_new: Structure
    struct_new, _ = cellify.mcp.load_structure_file(out_file)
    assert struct_new.composition.reduced_formula == "Si14Ge"


def test_mcp_slab(poscar_path: str, tmp_path: pytest.TempPathFactory) -> None:
    """
    Tests the cellify_slab tool.
    """
    out_file = os.path.join(tmp_path, "POSCAR_slab")
    res_msg: str = cellify_slab(
        poscar_path, miller="1 1 1", thick=4.0, vacuum=10.0, output_path=out_file
    )
    assert "Successfully generated slab model" in res_msg
    assert os.path.exists(out_file)


def test_mcp_errors() -> None:
    """
    Tests error handling for invalid files or options.
    """
    # Non-existent file error
    res_info: str = cellify_info("nonexistent_file")
    assert "Error:" in res_info

    res_conv: str = cellify_conventional("nonexistent_file")
    assert "Error:" in res_conv

    res_super: str = cellify_supercell("nonexistent_file")
    assert "Error:" in res_super

    res_defect: str = cellify_defect("nonexistent_file")
    assert "Error:" in res_defect

    res_slab: str = cellify_slab("nonexistent_file", miller="1 1 1", thick=4.0, vacuum=10.0)
    assert "Error:" in res_slab


def test_mcp_main() -> None:
    """
    Tests running the MCP server main loop.
    """
    with patch("cellify.mcp.mcp.run") as mock_run:
        cellify.mcp.main()
        mock_run.assert_called_once()

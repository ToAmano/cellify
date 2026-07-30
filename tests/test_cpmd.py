"""
Unit tests for CPMD adapter, templates, and JSON index mapping.
"""

import json
import os
from typing import Any, Dict

import pytest
from pymatgen.core import Lattice, Structure

from cellify.adapters.cpmd import CpmdAdapter
from cellify.core import process_template_and_validation, run_cellify_pipeline, save_structure_file


def test_cpmd_adapter_read_write(tmp_path: Any) -> None:
    # Create an unsorted structure with mixed atoms (O, Ti, O, Ti)
    lattice: Lattice = Lattice.cubic(4.0)
    species = ["O", "Ti", "O", "Ti"]
    coords = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.75, 0.75, 0.75],
    ]
    structure: Structure = Structure(lattice, species, coords, coords_are_cartesian=True)

    out_inp: str = str(tmp_path / "georelax.inp")
    adapter: CpmdAdapter = CpmdAdapter()
    adapter.write(out_inp, structure, {})

    # 1. Check output file exists
    assert os.path.exists(out_inp)

    # 2. Check JSON index map sidecar file exists and contains valid mapping
    json_map_path: str = f"{out_inp}.index_map.json"
    assert os.path.exists(json_map_path)

    with open(json_map_path, "r", encoding="utf-8") as f_json:
        data: Dict[str, Any] = json.load(f_json)

    assert data["num_atoms"] == 4
    assert data["sorted_symbols"] == ["O", "O", "Ti", "Ti"]
    assert "sort_indices" in data
    assert "inverse_indices" in data

    # 3. Check embedded ATOM_MAP_INVERSE comment in .inp file
    with open(out_inp, "r", encoding="utf-8") as f_inp:
        content: str = f_inp.read()

    assert "! ATOM_MAP_INVERSE:" in content
    assert "&ATOMS" in content
    assert "*pseudo/O_BLYP.psp" in content
    assert "*pseudo/Ti_BLYP.psp" in content

    # 4. Check reading written CPMD file back
    read_struct, meta = adapter.read(out_inp)
    assert len(read_struct) == 4
    assert meta["mode"] == "cpmd_text_replace"


def test_cpmd_template_aliases() -> None:
    meta_data: Dict[str, Any] = {"mode": "standard"}

    # Test template_path alias resolution
    resolved_meta = process_template_and_validation(
        meta_data, "output.inp", template_path="cpmd_georelax"
    )
    assert resolved_meta["mode"] == "cpmd_text_replace"

    # Test calc alias resolution
    resolved_calc_meta = process_template_and_validation(
        meta_data, "output.inp", calc="cpmd_bomd_relax"
    )
    assert resolved_calc_meta["mode"] == "cpmd_text_replace"


def test_cpmd_pipeline_integration(tmp_path: Any) -> None:
    lattice: Lattice = Lattice.cubic(5.0)
    structure: Structure = Structure(
        lattice, ["H", "O", "H"], [[0, 0, 0], [0.5, 0.5, 0.5], [0.2, 0.2, 0.2]], coords_are_cartesian=True
    )

    out_inp: str = str(tmp_path / "cpmd_run.inp")

    # Run pipeline with supercell scaling
    mod_struct, _ = run_cellify_pipeline(structure, dim=[2, 2, 2])
    meta_data: Dict[str, Any] = process_template_and_validation(
        {}, out_inp, template_path="cpmd_georelax"
    )
    save_structure_file(out_inp, mod_struct, meta_data)

    assert os.path.exists(out_inp)
    assert os.path.exists(f"{out_inp}.index_map.json")

    with open(out_inp, "r", encoding="utf-8") as f:
        text: str = f.read()

    assert "&ATOMS" in text
    assert "*pseudo/H_BLYP.psp" in text
    assert "*pseudo/O_BLYP.psp" in text
    assert "LMAX=S" in text
    assert "LMAX=P" in text


def test_cpmd_edge_cases(tmp_path: Any) -> None:
    adapter: CpmdAdapter = CpmdAdapter()

    # Test FileNotFoundError
    with pytest.raises(FileNotFoundError):
        adapter.read("non_existent_cpmd.inp")

    # Test parsing single-float cubic cell and Bohr coordinates
    bohr_cpmd_content = """&SYSTEM
   CELL
      18.897
&END

&ATOMS
*pseudo/O_BLYP.psp KLEINMAN-BYLANDER
   LMAX=P
   1
   0.0 0.0 0.0
invalid line to ignore
&END
"""
    cpmd_file = str(tmp_path / "bohr_sample.inp")
    with open(cpmd_file, "w", encoding="utf-8") as f:
        f.write(bohr_cpmd_content)

    struct, _ = adapter.read(cpmd_file)
    assert len(struct) == 1
    assert abs(struct.lattice.a - 10.0) < 1e-2

    # Test writing without &INFO block and without CELL block
    no_info_content = """&SYSTEM
   ANGSTROM
&END

&ATOMS
&END
"""
    no_info_file = str(tmp_path / "no_info.inp")
    with open(no_info_file, "w", encoding="utf-8") as f:
        f.write(no_info_content)

    adapter.write(no_info_file, struct, {"content": no_info_content})
    with open(no_info_file, "r", encoding="utf-8") as f:
        written_text = f.read()

    assert "! ATOM_MAP_INVERSE:" in written_text
    assert "CELL" in written_text

    # Test fallback minimal template generation
    min_tpl_file = str(tmp_path / "min_tpl.inp")
    adapter.write(min_tpl_file, struct, {"content": ""})
    assert os.path.exists(min_tpl_file)

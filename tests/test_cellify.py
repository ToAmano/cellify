import os
import re
import pytest
import numpy as np
from pymatgen.core import Structure
from cellify.core import (
    load_structure_file,
    save_structure_file,
    parse_matrix_string,
    calculate_min_dist_scaling,
    apply_substitutions,
    apply_vacancies,
    generate_surface_slab
)

@pytest.fixture
def poscar_path():
    return os.path.join(os.path.dirname(__file__), "POSCAR")

@pytest.fixture
def qe_path():
    return os.path.join(os.path.dirname(__file__), "qe.in")

def test_load_and_save_poscar(poscar_path, tmp_path):
    # ロードテスト
    structure, meta_data = load_structure_file(poscar_path)
    assert isinstance(structure, Structure)
    assert meta_data["mode"] == "standard"
    assert len(structure) == 2
    assert structure.composition.reduced_formula == "Si"
    
    # 保存テスト
    out_path = tmp_path / "POSCAR_out"
    save_structure_file(str(out_path), structure, meta_data)
    assert out_path.exists()
    
    # 再ロードして確認
    struct_new, _ = load_structure_file(str(out_path))
    assert len(struct_new) == 2

def test_load_and_save_qe(qe_path, tmp_path):
    # ロードテスト
    structure, meta_data = load_structure_file(qe_path)
    assert isinstance(structure, Structure)
    assert meta_data["mode"] == "espresso_text_replace"
    assert len(structure) == 2
    
    # 保存とパラメータ自動更新テスト (2x2x2スーパーセル化して保存)
    structure.make_supercell([2, 2, 2])
    out_path = tmp_path / "qe_out.in"
    save_structure_file(str(out_path), structure, meta_data)
    
    assert out_path.exists()
    with open(out_path, 'r') as f:
        content = f.read()
        
    # nat = 16 への更新チェック
    assert re.search(r'nat\s*=\s*16', content) is not None
    # ntyp = 1 の維持チェック
    assert re.search(r'ntyp\s*=\s*1', content) is not None
    # calculation = 'scf' の維持チェック
    assert "calculation = 'scf'" in content

def test_parse_matrix_string():
    matrix_str = "1 -1 0 / 1 1 0 / 0 0 2"
    matrix = parse_matrix_string(matrix_str)
    assert isinstance(matrix, np.ndarray)
    assert matrix.shape == (3, 3)
    assert np.allclose(matrix, [[1, -1, 0], [1, 1, 0], [0, 0, 2]])
    
    with pytest.raises(ValueError):
        parse_matrix_string("1 0 0 / 0 1 0") # 行が足りない

def test_calculate_min_dist_scaling(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    # min-dist 10.0 A に必要なスケーリング
    nx, ny, nz = calculate_min_dist_scaling(structure, 10.0)
    assert nx == 4
    assert ny == 4
    assert nz == 4

def test_apply_substitutions(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    structure.make_supercell([2, 2, 2]) # 16原子化
    
    # 0番目の原子を P に置換
    apply_substitutions(structure, ["Si:P:0"])
    assert structure.composition.reduced_formula == "Si15P"
    assert structure[0].specie.symbol == "P"

def test_apply_vacancies(poscar_path):
    structure, _ = load_structure_file(poscar_path)
    
    # 0番目の原子を削除
    apply_vacancies(structure, ["Si:0"])
    assert len(structure) == 1

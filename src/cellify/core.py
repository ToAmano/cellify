"""
cellify のコアモデリングロジック。
pymatgen と ASE を適宜使い分けて構造の読み込み、スーパーセル生成、
置換、空孔生成、スラブ作成、書き出しなどを行います。
"""
import math
import re
from typing import Tuple, Dict, Any, List, Optional
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator

from cellify.parser import parse_qe_input, write_qe_input

def load_structure_file(filepath: str) -> Tuple[Structure, Dict[str, Any]]:
    """
    ファイルをロードして structure オブジェクトとパースメタデータを返します。
    QE入力ファイルは特別扱いし、それ以外はpymatgenで自動判別読み込みします。
    """
    lower_path: str = filepath.lower()
    is_qe: bool = any(lower_path.endswith(ext) for ext in [".in", ".qe", ".pwi"]) or "qe" in lower_path or "espresso" in lower_path
    
    if is_qe:
        qe_data: Dict[str, Any] = parse_qe_input(filepath)
        return qe_data["structure"], qe_data
    else:
        struct: Structure = Structure.from_file(filepath)
        return struct, {"mode": "standard", "filepath": filepath}

def save_structure_file(filepath: str, structure: Structure, meta_data: Dict[str, Any]) -> None:
    """
    構造ファイルを保存します。QEの場合は元のパラメータを引き継ぎます。
    """
    if meta_data["mode"] != "standard":
        write_qe_input(filepath, structure, meta_data)
    else:
        structure.to(filename=filepath)

def parse_matrix_string(matrix_str: str) -> np.ndarray:
    """
    "1 -1 0 / 1 1 0 / 0 0 1" などの文字列を 3x3 の numpy 行列にパースします。
    """
    rows_raw: List[str] = re.split(r'[/,;]', matrix_str)
    if len(rows_raw) != 3:
        raise ValueError("Matrix string must define exactly 3 rows (separated by /, , or ;)")
    
    matrix: List[List[float]] = []
    for r in rows_raw:
        vals: List[float] = [float(x) for x in r.strip().split()]
        if len(vals) != 3:
            raise ValueError("Each row in the matrix must have exactly 3 elements")
        matrix.append(vals)
    
    return np.array(matrix)

def calculate_min_dist_scaling(structure: Structure, min_dist: float) -> Tuple[int, int, int]:
    """
    すべての格子ベクトルの周期境界における面間距離（法線距離）が
    min_dist 以上となるための最小の対角スケーリング因子 (nx, ny, nz) を計算します。
    """
    lattice = structure.lattice
    matrix = lattice.matrix
    a_vec, b_vec, c_vec = matrix[0], matrix[1], matrix[2]
    
    vol: float = lattice.volume
    
    d_a: float = vol / np.linalg.norm(np.cross(b_vec, c_vec))
    d_b: float = vol / np.linalg.norm(np.cross(c_vec, a_vec))
    d_c: float = vol / np.linalg.norm(np.cross(a_vec, b_vec))
    
    nx: int = int(math.ceil(min_dist / d_a))
    ny: int = int(math.ceil(min_dist / d_b))
    nz: int = int(math.ceil(min_dist / d_c))
    
    return max(1, nx), max(1, ny), max(1, nz)

def apply_substitutions(structure: Structure, substitute_rules: List[str]) -> None:
    """
    置換ルールを構造に適用します。
    """
    for rule in substitute_rules:
        parts: List[str] = rule.split(':')
        if len(parts) != 3:
            raise ValueError(f"Invalid substitution rule: {rule}. Must be 'element:target_element:index_or_percentage'")
        
        src_el, dest_el, target = parts[0], parts[1], parts[2]
        
        matching_indices: List[int] = [i for i, site in enumerate(structure) if site.specie.symbol == src_el]
        if not matching_indices:
            print(f"Warning: No matching elements found for substitution source '{src_el}'")
            continue
            
        if target.endswith('%'):
            percentage: float = float(target[:-1]) / 100.0
            num_to_replace: int = int(round(len(matching_indices) * percentage))
            if num_to_replace == 0 and percentage > 0:
                num_to_replace = 1
                
            replace_indices = np.random.choice(matching_indices, num_to_replace, replace=False)
            for idx in replace_indices:
                structure.replace(idx, dest_el)
            print(f"Replaced {num_to_replace} of {src_el} with {dest_el} ({target})")
        else:
            try:
                idx: int = int(target)
                if idx < 0 or idx >= len(structure):
                    raise IndexError(f"Index {idx} out of range (structure size: {len(structure)})")
                
                actual_symbol: str = structure[idx].specie.symbol
                if actual_symbol != src_el:
                    print(f"Warning: Site index {idx} is '{actual_symbol}', not source element '{src_el}'. Replacing anyway.")
                
                structure.replace(idx, dest_el)
                print(f"Replaced site {idx} ({actual_symbol}) with {dest_el}")
            except ValueError:
                raise ValueError(f"Invalid substitution target index or percentage: {target}")

def apply_vacancies(structure: Structure, vacancy_rules: List[str]) -> None:
    """
    空孔ルールを適用します。指定された原子を削除します。
    """
    indices_to_remove: List[int] = []
    
    for rule in vacancy_rules:
        parts: List[str] = rule.split(':')
        if len(parts) != 3 and len(parts) != 2:
            raise ValueError(f"Invalid vacancy rule: {rule}. Must be 'element:index' or 'element:count'")
        
        src_el: str = parts[0]
        target: str = parts[1]
        
        matching_indices: List[int] = [i for i, site in enumerate(structure) if site.specie.symbol == src_el]
        if not matching_indices:
            print(f"Warning: No matching elements found for vacancy source '{src_el}'")
            continue
            
        try:
            val: int = int(target)
            
            if val <= len(matching_indices) and val > 0 and len(structure) > 20:
                remove_subset = np.random.choice(matching_indices, val, replace=False)
                indices_to_remove.extend(remove_subset)
                print(f"Created {val} vacancies of {src_el} (randomly selected)")
            else:
                if val < 0 or val >= len(structure):
                    raise IndexError(f"Index {val} out of range")
                
                actual_symbol: str = structure[val].specie.symbol
                if actual_symbol != src_el:
                    print(f"Warning: Site index {val} is '{actual_symbol}', not vacancy element '{src_el}'. Removing anyway.")
                
                indices_to_remove.append(val)
                print(f"Removed site {val} ({actual_symbol}) to create vacancy")
        except ValueError:
            raise ValueError(f"Invalid vacancy target: {target}")
            
    if indices_to_remove:
        indices_to_remove = sorted(list(set(indices_to_remove)), reverse=True)
        structure.remove_sites(indices_to_remove)

def generate_surface_slab(
    structure: Structure, 
    miller_index: List[int], 
    thick: Optional[float], 
    vacuum: Optional[float]
) -> Structure:
    """
    pymatgen の SlabGenerator を使用して表面スラブモデルを構築します。
    """
    slab_thick: float = thick if thick else 10.0
    vac_thick: float = vacuum if vacuum else 15.0
    
    gen = SlabGenerator(
        initial_structure=structure,
        miller_index=miller_index,
        min_slab_size=slab_thick,
        min_vacuum_size=vac_thick,
        center_slab=True
    )
    
    slabs = gen.get_slabs()
    if not slabs:
        raise ValueError(f"Could not generate slab for Miller index {miller_index}")
    
    slab = slabs[0]
    return slab.generate_unique_slab_structs()[0]

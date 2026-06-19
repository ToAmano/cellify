"""
cellify のコアモデリングロジック。
pymatgen と ASE を適宜使い分けて構造の読み込み、スーパーセル生成、
置換、空孔生成、スラブ作成、書き出しなどを行います。
"""
import math
import re
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.core.surface import SlabGenerator

from cellify.parser import parse_qe_input, write_qe_input

def load_structure_file(filepath):
    """
    ファイルをロードして structure オブジェクトとパースメタデータを返します。
    QE入力ファイルは特別扱いし、それ以外はpymatgenで自動判別読み込みします。
    """
    # 拡張子やファイル名からQE入力ファイルを判定
    # 例: *.in, *.qe, もしくはファイル名が 'qe.in' や 'espresso.in'
    lower_path = filepath.lower()
    is_qe = any(lower_path.endswith(ext) for ext in [".in", ".qe", ".pwi"]) or "qe" in lower_path or "espresso" in lower_path
    
    if is_qe:
        qe_data = parse_qe_input(filepath)
        return qe_data["structure"], qe_data
    else:
        # VASP, CIF, XYZ などは pymatgen の自動判別ロードを使用
        struct = Structure.from_file(filepath)
        return struct, {"mode": "standard", "filepath": filepath}

def save_structure_file(filepath, structure, meta_data):
    """
    構造ファイルを保存します。QEの場合は元のパラメータを引き継ぎます。
    """
    if meta_data["mode"] in ["pymatgen", "fallback"]:
        write_qe_input(filepath, structure, meta_data)
    else:
        # standard モードの場合
        structure.to(filename=filepath)

def parse_matrix_string(matrix_str):
    """
    "1 -1 0 / 1 1 0 / 0 0 1" などの文字列を 3x3 の numpy 行列にパースします。
    """
    # スラッシュ、カンマ、またはセミコロンで各行を分割
    rows_raw = re.split(r'[/,;]', matrix_str)
    if len(rows_raw) != 3:
        raise ValueError("Matrix string must define exactly 3 rows (separated by /, , or ;)")
    
    matrix = []
    for r in rows_raw:
        vals = [float(x) for x in r.strip().split()]
        if len(vals) != 3:
            raise ValueError("Each row in the matrix must have exactly 3 elements")
        matrix.append(vals)
    
    return np.array(matrix)

def calculate_min_dist_scaling(structure, min_dist):
    """
    すべての格子ベクトルの周期境界における面間距離（法線距離）が
    min_dist 以上となるための最小の対角スケーリング因子 (nx, ny, nz) を計算します。
    """
    lattice = structure.lattice
    matrix = lattice.matrix
    a_vec, b_vec, c_vec = matrix[0], matrix[1], matrix[2]
    
    # 体積 V
    vol = lattice.volume
    
    # 各面の法線方向の間隔（面間距離） d_i
    # d_a = V / |b x c|
    # d_b = V / |c x a|
    # d_c = V / |a x b|
    d_a = vol / np.linalg.norm(np.cross(b_vec, c_vec))
    d_b = vol / np.linalg.norm(np.cross(c_vec, a_vec))
    d_c = vol / np.linalg.norm(np.cross(a_vec, b_vec))
    
    # 必要なスケーリング因子
    nx = int(math.ceil(min_dist / d_a))
    ny = int(math.ceil(min_dist / d_b))
    nz = int(math.ceil(min_dist / d_c))
    
    # 最低でも 1x1x1 にする
    return max(1, nx), max(1, ny), max(1, nz)

def apply_substitutions(structure, substitute_rules):
    """
    置換ルールを構造に適用します。
    ルール形式例: "Si:P:0" (インデックス0のSiをPに置換)
                 "Si:Al:5%" (Siの5%をAlにランダム置換)
    """
    for rule in substitute_rules:
        parts = rule.split(':')
        if len(parts) != 3:
            raise ValueError(f"Invalid substitution rule: {rule}. Must be 'element:target_element:index_or_percentage'")
        
        src_el, dest_el, target = parts[0], parts[1], parts[2]
        
        # 指定されたソース元素に該当するサイトのインデックスを抽出
        matching_indices = [i for i, site in enumerate(structure) if site.specie.symbol == src_el]
        if not matching_indices:
            print(f"Warning: No matching elements found for substitution source '{src_el}'")
            continue
            
        if target.endswith('%'):
            # パーセント比率に基づくランダム置換
            percentage = float(target[:-1]) / 100.0
            num_to_replace = int(round(len(matching_indices) * percentage))
            # 最低1つは置換するようにする（0%でない限り）
            if num_to_replace == 0 and percentage > 0:
                num_to_replace = 1
                
            # ランダムに選択
            replace_indices = np.random.choice(matching_indices, num_to_replace, replace=False)
            for idx in replace_indices:
                structure.replace(idx, dest_el)
            print(f"Replaced {num_to_replace} of {src_el} with {dest_el} ({target})")
        else:
            # インデックス指定
            try:
                # 該当する元素のリストにおける相対インデックスか、絶対インデックスか？
                # ここでは使いやすさのため、ファイル全体の絶対インデックスとして解釈
                idx = int(target)
                if idx < 0 or idx >= len(structure):
                    raise IndexError(f"Index {idx} out of range (structure size: {len(structure)})")
                
                actual_symbol = structure[idx].specie.symbol
                if actual_symbol != src_el:
                    print(f"Warning: Site index {idx} is '{actual_symbol}', not source element '{src_el}'. Replacing anyway.")
                
                structure.replace(idx, dest_el)
                print(f"Replaced site {idx} ({actual_symbol}) with {dest_el}")
            except ValueError:
                raise ValueError(f"Invalid substitution target index or percentage: {target}")

def apply_vacancies(structure, vacancy_rules):
    """
    空孔ルールを適用します。指定された原子を削除します。
    ルール形式例: "Si:0" (インデックス0のSiを削除)
                 "O:2" (O原子をランダムに2個削除)
    """
    indices_to_remove = []
    
    for rule in vacancy_rules:
        parts = rule.split(':')
        if len(parts) != 3 and len(parts) != 2:
            raise ValueError(f"Invalid vacancy rule: {rule}. Must be 'element:index' or 'element:count'")
        
        src_el = parts[0]
        target = parts[1]
        
        matching_indices = [i for i, site in enumerate(structure) if site.specie.symbol == src_el]
        if not matching_indices:
            print(f"Warning: No matching elements found for vacancy source '{src_el}'")
            continue
            
        try:
            # 数値のパース。もし指定された数値が総数以下であれば「個数」としてランダム削除、
            # 指定された元素の絶対インデックスに該当する場合は「その位置」を削除
            val = int(target)
            
            # 判別：該当元素の個数以下の場合は「個数指定」として扱う
            # ただし、インデックス指定と曖昧になるのを避けるため、ルール指定を明確にする
            if val <= len(matching_indices) and val > 0 and len(structure) > 20: # 簡易判定
                # 個数指定としてランダム削除
                remove_subset = np.random.choice(matching_indices, val, replace=False)
                indices_to_remove.extend(remove_subset)
                print(f"Created {val} vacancies of {src_el} (randomly selected)")
            else:
                # インデックス指定
                if val < 0 or val >= len(structure):
                    raise IndexError(f"Index {val} out of range")
                
                actual_symbol = structure[val].specie.symbol
                if actual_symbol != src_el:
                    print(f"Warning: Site index {val} is '{actual_symbol}', not vacancy element '{src_el}'. Removing anyway.")
                
                indices_to_remove.append(val)
                print(f"Removed site {val} ({actual_symbol}) to create vacancy")
        except ValueError:
            raise ValueError(f"Invalid vacancy target: {target}")
            
    if indices_to_remove:
        # 重複を排除して降順でソート（インデックスのズレを防ぐため）
        indices_to_remove = sorted(list(set(indices_to_remove)), reverse=True)
        structure.remove_sites(indices_to_remove)

def generate_surface_slab(structure, miller_index, thick, vacuum):
    """
    pymatgen の SlabGenerator を使用して表面スラブモデルを構築します。
    """
    # デフォルトの厚み値などを設定
    slab_thick = thick if thick else 10.0
    vac_thick = vacuum if vacuum else 15.0
    
    # SlabGenerator の初期化
    # center_slab=True でスラブをセルの中心に配置し、両側に真空を作る
    gen = SlabGenerator(
        initial_structure=structure,
        miller_index=miller_index,
        min_slab_size=slab_thick,
        min_vacuum_size=vac_thick,
        center_slab=True
    )
    
    # 可能な終端構造（terminations）のリストを取得し、最初のスラブを返す
    slabs = gen.get_slabs()
    if not slabs:
        raise ValueError(f"Could not generate slab for Miller index {miller_index}")
    
    # 最も対称性が高く安定な可能性がある最初のスラブモデルを採用
    slab = slabs[0]
    # Structure型に変換して返す
    return slab.generate_unique_slab_structs()[0]

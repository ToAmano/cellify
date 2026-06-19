"""
cellify のための入力・出力パーサモジュール。
主に Quantum ESPRESSO などの計算パラメータと構造データが混在するファイルフォーマットについて、
パラメータ、コメント、記述順を完全に維持しつつ、構造パラメータ（nat, ntyp等）を自動更新して出力する機能を提供します。
"""
import re
import os
from pymatgen.core import Structure

def parse_qe_input(filepath):
    """
    Quantum ESPRESSOの入力ファイルをテキストとして読み込み、構造データを抽出します。
    構造のパースには最も堅牢なASEのespresso-inパーサを利用します。
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")

    with open(filepath, 'r') as f:
        content = f.read()

    # 構造データをASE経由で安全にロードしてpymatgen Structureに変換
    try:
        from ase.io import read as ase_read
        from pymatgen.io.ase import AseAtomsAdaptor
        atoms = ase_read(filepath, format="espresso-in")
        structure = AseAtomsAdaptor.get_structure(atoms)
    except Exception as ase_err:
        raise ValueError(f"Failed to parse structure from Quantum ESPRESSO file: {ase_err}")

    return {
        "mode": "text_replace",
        "content": content,
        "structure": structure,
        "filepath": filepath
    }

def write_qe_input(filepath, structure, original_data):
    """
    元のパラメータ、コメント、その他のブロック（K_POINTS等）を完全に維持し、
    nat / ntyp および構造座標データのみを更新したQE入力ファイルを書き出します。
    """
    content = original_data["content"]
    
    # 1. 新しい原子数と原子種数を計算
    nat_new = len(structure)
    ntyp_new = len(structure.composition.elements)
    
    # 2. ネームリスト内の nat と ntyp を更新
    # 単語境界 (\b) を考慮して正確に置換
    content = re.sub(r'(\bnat\s*=\s*)\d+', r'\g<1>' + str(nat_new), content, flags=re.IGNORECASE)
    content = re.sub(r'(\bntyp\s*=\s*)\d+', r'\g<1>' + str(ntyp_new), content, flags=re.IGNORECASE)
    
    # 3. 元の構造関連ブロック（ATOMIC_SPECIES, CELL_PARAMETERS, ATOMIC_POSITIONS）をテキストから除去
    # 各ブロックはキーワードで始まり、次のブロック名、ネームリストの開始(&)、またはファイル末尾まで続く
    cleaned_content = content
    struct_keywords = ["ATOMIC_SPECIES", "CELL_PARAMETERS", "ATOMIC_POSITIONS"]
    for kw in struct_keywords:
        # キーワードから始まり、別のキーワード、ネームリスト、またはファイル末尾の手前までを削除
        pattern = r'(?i)^\s*' + kw + r'\b.*?(?=\n\s*(?:ATOMIC_SPECIES|CELL_PARAMETERS|ATOMIC_POSITIONS|K_POINTS|KPOINTS|&[A-Za-z]+)|\Z)'
        cleaned_content = re.sub(pattern, '', cleaned_content, flags=re.DOTALL | re.MULTILINE)
    
    # 前後の余分な改行を整理
    cleaned_content = cleaned_content.strip() + "\n\n"
    
    # 4. 元ファイルから擬ポテンシャル情報を抽出
    pseudos = {}
    species_match = re.search(r'(?i)ATOMIC_SPECIES\s*\n(.*?)(?=\n\s*(?:ATOMIC_|CELL_|K_POINTS|KPOINTS|&[A-Za-z]+)|\Z)', content, re.DOTALL)
    if species_match:
        for line in species_match.group(1).strip().split('\n'):
            parts = line.split()
            if len(parts) >= 3:
                pseudos[parts[0]] = (parts[1], parts[2])
    
    # 5. 各構造ブロックの再構築
    # ATOMIC_SPECIES
    species_str = "ATOMIC_SPECIES\n"
    for el in structure.composition.elements:
        el_symbol = el.symbol
        mass, pseudo = pseudos.get(el_symbol, (str(el.atomic_mass), f"{el_symbol}.UPF"))
        species_str += f"  {el_symbol}  {mass}  {pseudo}\n"
    
    # CELL_PARAMETERS (angstrom単位)
    cell_str = "\nCELL_PARAMETERS angstrom\n"
    for vec in structure.lattice.matrix:
        cell_str += f"  {vec[0]:.10f}  {vec[1]:.10f}  {vec[2]:.10f}\n"
    
    # ATOMIC_POSITIONS (crystal/fractional座標)
    pos_str = "\nATOMIC_POSITIONS crystal\n"
    for site in structure:
        pos_str += f"  {site.specie.symbol}  {site.a:.10f}  {site.b:.10f}  {site.c:.10f}\n"
    
    # 6. 保存
    with open(filepath, 'w') as f:
        f.write(cleaned_content)
        f.write(species_str)
        f.write(cell_str)
        f.write(pos_str)


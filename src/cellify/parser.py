"""
cellify のための入力・出力パーサモジュール。
主に Quantum ESPRESSO などの計算パラメータと構造データが混在するファイルフォーマットについて、
パラメータを維持しつつ構造パラメータ（nat, ntyp等）を自動更新して出力する機能を提供します。
"""
import re
import os
from pymatgen.io.qe.inputs import PWInput
from pymatgen.core import Structure

def parse_qe_input(filepath):
    """
    Quantum ESPRESSOの入力ファイルをパースし、構造データとパラメータメタデータを取得します。
    PymatgenのPWInputでの厳密なパースを試み、失敗した場合は正規表現によるフォールバック処理を行います。
    
    Returns:
        dict: パース結果オブジェクト
    """
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Input file not found: {filepath}")

    try:
        # Pymatgenのパーサで試行
        pw = PWInput.from_file(filepath)
        return {
            "mode": "pymatgen",
            "data": pw,
            "structure": pw.structure
        }
    except Exception as e:
        # Pymatgenが非標準パラメータや擬ポテンシャル定義不足等で失敗した場合のフォールバック
        with open(filepath, 'r') as f:
            content = f.read()
        
        # 構造をパースするために最低限の情報を一時的に抽出してpymatgenに渡すか、
        # またはASEを中介して構造のみを読み込む
        # ここではASEがQE入力ファイルの座標パースに対して頑健であるため、ASEで構造をロードする
        try:
            from ase.io import read as ase_read
            from pymatgen.io.ase import AseAtomsAdaptor
            atoms = ase_read(filepath, format="espresso-in")
            structure = AseAtomsAdaptor.get_structure(atoms)
        except Exception as ase_err:
            raise ValueError(f"Failed to parse structure from Quantum ESPRESSO file: {ase_err}")

        return {
            "mode": "fallback",
            "content": content,
            "structure": structure,
            "filepath": filepath
        }

def write_qe_input(filepath, structure, original_data):
    """
    スーパーセル化または修飾された structure を用いて、
    元のパラメータを引き継ぎつつ nat / ntyp / 構造部を自動更新したQE入力ファイルを書き出します。
    
    Args:
        filepath (str): 出力ファイルパス
        structure (Structure): スーパーセル化後の pymatgen Structure オブジェクト
        original_data (dict): parse_qe_input の返り値
    """
    if original_data["mode"] == "pymatgen":
        pw = original_data["data"]
        pw.structure = structure
        # PWInputは write_file 時に nat, ntyp などを自動更新してくれる
        pw.write_file(filepath)
    
    else:
        # フォールバック処理（プレーンテキストのネームリスト部を抽出し、構造部を差し替える）
        content = original_data["content"]
        
        # 1. 新しい原子数と原子種数を計算
        nat_new = len(structure)
        ntyp_new = len(structure.composition.elements)
        
        # 2. ネームリスト内の nat と ntyp を置換（大文字小文字対応、コメント行を除外する配慮）
        # 簡単な正規表現置換（インラインコメントや空白を考慮）
        content = re.sub(r'(nat\s*=\s*)\d+', r'\g<1>' + str(nat_new), content, flags=re.IGNORECASE)
        content = re.sub(r'(ntyp\s*=\s*)\d+', r'\g<1>' + str(ntyp_new), content, flags=re.IGNORECASE)
        
        # 3. &で始まるネームリストブロック（&CONTROL, &SYSTEM, &ELECTRONSなど）だけを抽出してヘッダー化
        # 各ブロックは '&name ... /' という形式
        namelists = []
        for match in re.finditer(r'&\w+.*?/', content, re.DOTALL | re.IGNORECASE):
            namelists.append(match.group(0))
        
        header = "\n\n".join(namelists) + "\n\n"
        
        # 4. 元ファイルから擬ポテンシャル情報 (ATOMIC_SPECIES) を抽出
        # 新しい元素がドーピング等で追加された場合に備え、既存の定義を辞書化する
        species_match = re.search(r'ATOMIC_SPECIES\s*\n(.*?)(?=\n\S|\Z)', content, re.DOTALL | re.IGNORECASE)
        pseudos = {}
        if species_match:
            for line in species_match.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 3:
                    pseudos[parts[0]] = (parts[1], parts[2])  # { 'Element': ('Mass', 'Pseudo.UPF') }
        
        # 5. ATOMIC_SPECIES ブロックの再構築
        species_str = "ATOMIC_SPECIES\n"
        for el in structure.composition.elements:
            el_symbol = el.symbol
            mass, pseudo = pseudos.get(el_symbol, (str(el.atomic_mass), f"{el_symbol}.UPF"))
            species_str += f"  {el_symbol}  {mass}  {pseudo}\n"
        
        # 6. CELL_PARAMETERS ブロックの再構築
        cell_str = "CELL_PARAMETERS angstrom\n"
        for vec in structure.lattice.matrix:
            cell_str += f"  {vec[0]:.10f}  {vec[1]:.10f}  {vec[2]:.10f}\n"
        
        # 7. ATOMIC_POSITIONS ブロックの再構築（結晶座標 fractional で統一）
        pos_str = "ATOMIC_POSITIONS crystal\n"
        for site in structure:
            pos_str += f"  {site.specie.symbol}  {site.a:.10f}  {site.b:.10f}  {site.c:.10f}\n"
        
        # 8. 結合してファイル保存
        with open(filepath, 'w') as f:
            f.write(header)
            f.write(species_str + "\n")
            f.write(cell_str + "\n")
            f.write(pos_str)

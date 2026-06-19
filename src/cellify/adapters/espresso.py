import re
import os
from typing import Dict, Any, Tuple
from pymatgen.core import Structure
from cellify.adapters.base import BaseAdapter

class EspressoAdapter(BaseAdapter):
    """
    Quantum ESPRESSO 入力ファイル用の I/O アダプター。
    元の制御・計算パラメータ（&CONTROL, &SYSTEM等）やコメント行を完全に保護したまま、
    nat, ntyp を自動更新し、構造座標ブロックのみを置換して書き出します。
    """
    
    def read(self, filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Input file not found: {filepath}")

        with open(filepath, 'r') as f:
            content: str = f.read()

        # 構造データをASE経由で安全にロードしてpymatgen Structureに変換
        try:
            from ase.io import read as ase_read
            from pymatgen.io.ase import AseAtomsAdaptor
            atoms = ase_read(filepath, format="espresso-in")
            structure: Structure = AseAtomsAdaptor.get_structure(atoms)
        except Exception as ase_err:
            raise ValueError(f"Failed to parse structure from Quantum ESPRESSO file: {ase_err}")

        meta_data: Dict[str, Any] = {
            "mode": "espresso_text_replace",
            "content": content,
            "filepath": filepath
        }
        return structure, meta_data

    def write(self, filepath: str, structure: Structure, meta_data: Dict[str, Any]) -> None:
        content: str = meta_data["content"]
        
        # 1. 新しい原子数と原子種数を計算
        nat_new: int = len(structure)
        ntyp_new: int = len(structure.composition.elements)
        
        # 2. ネームリスト内の nat と ntyp を更新
        content = re.sub(r'(\bnat\s*=\s*)\d+', r'\g<1>' + str(nat_new), content, flags=re.IGNORECASE)
        content = re.sub(r'(\bntyp\s*=\s*)\d+', r'\g<1>' + str(ntyp_new), content, flags=re.IGNORECASE)
        
        # 3. 元の構造関連ブロックをテキストから除去
        cleaned_content: str = content
        struct_keywords = ["ATOMIC_SPECIES", "CELL_PARAMETERS", "ATOMIC_POSITIONS"]
        for kw in struct_keywords:
            pattern = r'(?i)^\s*' + kw + r'\b.*?(?=\n\s*(?:ATOMIC_SPECIES|CELL_PARAMETERS|ATOMIC_POSITIONS|K_POINTS|KPOINTS|&[A-Za-z]+)|\Z)'
            cleaned_content = re.sub(pattern, '', cleaned_content, flags=re.DOTALL | re.MULTILINE)
        
        # 前後の余分な改行を整理
        cleaned_content = cleaned_content.strip() + "\n\n"
        
        # 4. 元ファイルから擬ポテンシャル情報を抽出
        pseudos: Dict[str, Any] = {}
        species_match = re.search(r'(?i)ATOMIC_SPECIES\s*\n(.*?)(?=\n\s*(?:ATOMIC_|CELL_|K_POINTS|KPOINTS|&[A-Za-z]+)|\Z)', content, re.DOTALL)
        if species_match:
            for line in species_match.group(1).strip().split('\n'):
                parts = line.split()
                if len(parts) >= 3:
                    pseudos[parts[0]] = (parts[1], parts[2])
        
        # 5. 各構造ブロックの再構築
        # ATOMIC_SPECIES
        species_str: str = "ATOMIC_SPECIES\n"
        for el in structure.composition.elements:
            el_symbol: str = el.symbol
            mass, pseudo = pseudos.get(el_symbol, (str(el.atomic_mass), f"{el_symbol}.UPF"))
            species_str += f"  {el_symbol}  {mass}  {pseudo}\n"
        
        # CELL_PARAMETERS
        cell_str: str = "\nCELL_PARAMETERS angstrom\n"
        for vec in structure.lattice.matrix:
            cell_str += f"  {vec[0]:.10f}  {vec[1]:.10f}  {vec[2]:.10f}\n"
        
        # ATOMIC_POSITIONS
        pos_str: str = "\nATOMIC_POSITIONS crystal\n"
        for site in structure:
            pos_str += f"  {site.specie.symbol}  {site.a:.10f}  {site.b:.10f}  {site.c:.10f}\n"
        
        # 6. 保存
        with open(filepath, 'w') as f:
            f.write(cleaned_content)
            f.write(species_str)
            f.write(cell_str)
            f.write(pos_str)

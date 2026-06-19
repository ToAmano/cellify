from typing import Dict, Any, Tuple
from pymatgen.core import Structure
from cellify.adapters.base import BaseAdapter

class StandardAdapter(BaseAdapter):
    """
    VASP (POSCAR), CIF, XYZ などの一般的な結晶構造ファイル用の I/O アダプター。
    特別なパラメータ維持処理は行わず、pymatgenの標準機能を用いて直接入出力を行います。
    """
    
    def read(self, filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        struct: Structure = Structure.from_file(filepath)
        meta_data: Dict[str, Any] = {
            "mode": "standard",
            "filepath": filepath
        }
        return struct, meta_data

    def write(self, filepath: str, structure: Structure, meta_data: Dict[str, Any]) -> None:
        structure.to(filename=filepath)

from cellify.adapters.base import BaseAdapter
from cellify.adapters.espresso import EspressoAdapter
from cellify.adapters.standard import StandardAdapter

def get_adapter(filepath: str) -> BaseAdapter:
    """
    ファイルパスや拡張子に基づいて、適切な I/O アダプターオブジェクトを返します。
    """
    lower_path: str = filepath.lower()
    # Quantum ESPRESSO 入力ファイルの判定
    is_qe: bool = any(lower_path.endswith(ext) for ext in [".in", ".qe", ".pwi"]) or "qe" in lower_path or "espresso" in lower_path
    
    if is_qe:
        return EspressoAdapter()
    else:
        return StandardAdapter()

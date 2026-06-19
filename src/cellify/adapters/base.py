from abc import ABC, abstractmethod
from typing import Dict, Any, Tuple
from pymatgen.core import Structure

class BaseAdapter(ABC):
    """
    cellify がサポートする結晶構造ファイルの入出力を抽象化する基底クラス。
    各DFT計算ソフトに特化したパラメータ引き継ぎロジックは、このクラスを継承して実装します。
    """
    
    @abstractmethod
    def read(self, filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        """
        結晶構造ファイルをロードし、構造オブジェクトとパースメタデータを返します。
        
        Args:
            filepath (str): 読み込み対象のファイルパス
            
        Returns:
            Tuple[Structure, Dict[str, Any]]: (Structureオブジェクト, メタデータ辞書)
        """
        pass

    @abstractmethod
    def write(self, filepath: str, structure: Structure, meta_data: Dict[str, Any]) -> None:
        """
        構造オブジェクトを指定されたパスに書き出します。
        メタデータ（元の計算パラメータやコメントなど）が存在する場合は引き継ぎます。
        
        Args:
            filepath (str): 書き出し先のファイルパス
            structure (Structure): スーパーセル化・修飾後の Structure オブジェクト
            meta_data (Dict[str, Any]): read 時に取得したメタデータ
        """
        pass

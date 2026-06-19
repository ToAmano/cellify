from abc import ABC, abstractmethod
from typing import Any, Dict, Tuple

from pymatgen.core import Structure


class BaseAdapter(ABC):
    """
    Abstract base class for structure file I/O supported by cellify.
    Parameter-preserving and software-specific output adapters should inherit this class.
    """

    @abstractmethod
    def read(self, filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        """
        Loads a structure file and returns the structure object along with metadata.

        Args:
            filepath (str): Path to the input file.

        Returns:
            Tuple[Structure, Dict[str, Any]]: A tuple of the pymatgen Structure object and a metadata dictionary.
        """
        pass

    @abstractmethod
    def write(
        self, filepath: str, structure: Structure, meta_data: Dict[str, Any]
    ) -> None:
        """
        Writes the structure to the specified path while preserving original metadata.

        Args:
            filepath (str): Path to the output file.
            structure (Structure): The modified/supercell Structure object.
            meta_data (Dict[str, Any]): Metadata retrieved during the read phase.
        """
        pass

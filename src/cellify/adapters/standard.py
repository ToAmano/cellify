from typing import Any, Dict, Tuple

from pymatgen.core import Structure

from cellify.adapters.base import BaseAdapter


class StandardAdapter(BaseAdapter):
    """
    Standard structure file adapter for formats like VASP (POSCAR), CIF, XYZ, etc.
    Does not perform parameter-preserving text replacements, and uses pymatgen's
    default I/O functionalities.
    """

    def read(self, filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        struct: Structure = Structure.from_file(filepath)
        meta_data: Dict[str, Any] = {"mode": "standard", "filepath": filepath}
        return struct, meta_data

    def write(
        self, filepath: str, structure: Structure, meta_data: Dict[str, Any]
    ) -> None:
        structure.to(filename=filepath)

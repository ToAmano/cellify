"""
I/O Adapters package for cellify.
"""

from cellify.adapters.base import BaseAdapter
from cellify.adapters.espresso import EspressoAdapter
from cellify.adapters.standard import StandardAdapter

__all__ = ["BaseAdapter", "EspressoAdapter", "StandardAdapter"]


def get_adapter(filepath: str) -> BaseAdapter:
    """
    Returns an appropriate I/O adapter object based on the filepath or extension.
    """
    lower_path: str = filepath.lower()
    # Check if the file is a Quantum ESPRESSO input file
    is_qe: bool = (
        any(lower_path.endswith(ext) for ext in [".in", ".qe", ".pwi"])
        or "qe" in lower_path
        or "espresso" in lower_path
    )

    if is_qe:
        return EspressoAdapter()
    return StandardAdapter()

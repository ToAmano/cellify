"""
I/O Adapters package for cellify.
"""

from cellify.adapters.base import BaseAdapter
from cellify.adapters.cpmd import CpmdAdapter
from cellify.adapters.espresso import EspressoAdapter
from cellify.adapters.standard import StandardAdapter

__all__ = ["BaseAdapter", "CpmdAdapter", "EspressoAdapter", "StandardAdapter"]


def get_adapter(filepath: str) -> BaseAdapter:
    """
    Returns an appropriate I/O adapter object based on the filepath, extension, or content.
    """
    import os

    lower_path: str = filepath.lower()
    is_cpmd: bool = lower_path.endswith(".cpmd") or "cpmd" in lower_path
    is_qe: bool = (
        any(lower_path.endswith(ext) for ext in [".in", ".qe", ".pwi", ".pwo"])
        or "qe" in lower_path
        or "espresso" in lower_path
        or "pwscf" in lower_path
    )

    # Content-based detection if file exists and name checks aren't conclusive
    if os.path.exists(filepath):
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                head: str = f.read(1000).lower()
                if "&cpmd" in head or "&atoms" in head:
                    is_cpmd = True
                elif "&control" in head or "program pwscf" in head or "pwscf" in head:
                    is_qe = True
        except Exception:  # pylint: disable=broad-exception-caught
            pass

    if is_cpmd:
        return CpmdAdapter()
    if is_qe:
        return EspressoAdapter()
    return StandardAdapter()

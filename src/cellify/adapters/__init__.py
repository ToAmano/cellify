"""
I/O Adapters package for cellify.
"""

from cellify.adapters.base import BaseAdapter
from cellify.adapters.espresso import EspressoAdapter
from cellify.adapters.standard import StandardAdapter

__all__ = ["BaseAdapter", "EspressoAdapter", "StandardAdapter"]


def get_adapter(filepath: str) -> BaseAdapter:
    """
    Returns an appropriate I/O adapter object based on the filepath, extension, or content.
    """
    import os

    lower_path: str = filepath.lower()
    is_qe: bool = (
        any(lower_path.endswith(ext) for ext in [".in", ".qe", ".pwi", ".pwo"])
        or "qe" in lower_path
        or "espresso" in lower_path
        or "pwscf" in lower_path
    )

    # Content-based detection if file exists and name checks aren't conclusive
    if not is_qe and os.path.exists(filepath):
        try:
            with open(filepath, "r", encoding="utf-8") as f:
                head: str = f.read(1000)
                if (
                    "&control" in head.lower()
                    or "program pwscf" in head.lower()
                    or "pwscf" in head.upper()
                ):
                    is_qe = True
        except Exception:  # pylint: disable=broad-exception-caught
            pass

    if is_qe:
        return EspressoAdapter()
    return StandardAdapter()

"""
OPTIMADE API client and query handlers for cellify.
Allows querying materials databases by chemical formula.
"""

from typing import Any, Dict, List, Optional, Tuple

import requests  # type: ignore[import-untyped]
from pymatgen.core import Composition, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def parse_optimade_entry_to_structure(entry: Dict[str, Any]) -> Optional[Structure]:
    """Converts an OPTIMADE entry dictionary to a pymatgen Structure object."""
    try:
        attrs = entry.get("attributes", {})
        lattice = attrs.get("lattice_vectors")
        species = attrs.get("species_at_sites")
        coords = attrs.get("cartesian_site_positions")
        if lattice and species and coords:
            return Structure(lattice, species, coords, coords_are_cartesian=True)
    except Exception:  # pylint: disable=broad-exception-caught
        pass
    return None


def _format_cod_fallback_info(attrs: Dict[str, Any], info_parts: List[str]) -> None:
    """Helper to extract and format fallback attributes from COD metadata."""
    sg_symbol = attrs.get("_cod_sg")
    if sg_symbol:
        info_parts.append(f"Space Group: {sg_symbol}")

    vol = attrs.get("_cod_vol") or attrs.get("_cod_volume")
    if vol:
        try:
            info_parts.append(f"Volume: {float(vol):.2f} A^3")
        except ValueError:
            info_parts.append(f"Volume: {vol} A^3")

    a = attrs.get("_cod_a")
    b = attrs.get("_cod_b")
    c = attrs.get("_cod_c")
    if a and b and c:
        try:
            info_parts.append(
                f"Lattice: a={float(a):.2f}, b={float(b):.2f}, c={float(c):.2f} A"
            )
        except ValueError:
            info_parts.append(f"Lattice: a={a}, b={b}, c={c} A")


def _extract_energy_above_hull(attrs: Dict[str, Any]) -> Optional[Any]:
    """Helper to extract energy_above_hull from flat or nested stability dicts."""
    stability = attrs.get("_mp_stability", {})
    if not isinstance(stability, dict):
        return None
    if "energy_above_hull" in stability:
        return stability.get("energy_above_hull")
    for val in stability.values():
        if isinstance(val, dict) and "energy_above_hull" in val:
            return val.get("energy_above_hull")
    return None


def _format_stability_info(attrs: Dict[str, Any], info_parts: List[str]) -> None:
    """Helper to extract and format stability information from Materials Project metadata."""
    e_above_hull = _extract_energy_above_hull(attrs)
    if e_above_hull is None:
        return

    if isinstance(e_above_hull, float):
        if abs(e_above_hull) < 1e-6:
            info_parts.append("Energy Above Hull: 0.0000 eV/atom [Stable ★]")
        else:
            info_parts.append(f"Energy Above Hull: {e_above_hull:.4f} eV/atom")
    else:
        info_parts.append(f"Energy Above Hull: {e_above_hull} eV/atom")


def _format_cod_names_info(attrs: Dict[str, Any], info_parts: List[str]) -> None:
    """Helper to extract, de-duplicate, and format names from COD metadata."""
    names = []
    seen = set()
    for name_key in ("_cod_commonname", "_cod_chemname", "_cod_mineral"):
        val = attrs.get(name_key)
        if val:
            clean_val = str(val).strip()
            if clean_val and clean_val.lower() not in seen:
                seen.add(clean_val.lower())
                names.append(clean_val)
    if names:
        info_parts.append(f"Name: {', '.join(names)}")


def _format_structure_info(idx: int, entry: Dict[str, Any], formula: str) -> str:
    """Helper to format a single structure OPTIMADE entry info."""
    attrs = entry.get("attributes", {})
    desc = attrs.get("chemical_formula_descriptive", formula)
    info_parts = [f"  [{idx+1}] ID: {entry.get('id')}", f"Formula: {desc}"]

    structure = parse_optimade_entry_to_structure(entry)
    if structure is not None:
        try:
            sga = SpacegroupAnalyzer(structure)
            sg_symbol = sga.get_space_group_symbol()
            info_parts.append(f"Space Group: {sg_symbol}")
        except Exception:  # pylint: disable=broad-exception-caught
            pass

        info_parts.append(f"Volume: {structure.volume:.2f} A^3")
        lat = structure.lattice
        info_parts.append(f"Lattice: a={lat.a:.2f}, b={lat.b:.2f}, c={lat.c:.2f} A")
    else:
        _format_cod_fallback_info(attrs, info_parts)

    _format_stability_info(attrs, info_parts)
    _format_cod_names_info(attrs, info_parts)

    return " | ".join(info_parts)


def _sort_materials_project_data(data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sorts Materials Project entries by thermodynamic stability."""

    def get_sort_key(entry: Dict[str, Any]) -> float:
        attrs = entry.get("attributes", {})
        val = _extract_energy_above_hull(attrs)
        if isinstance(val, (int, float)):
            return float(val)
        return float("inf")

    return sorted(data, key=get_sort_key)


def _sort_cod_data(data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sorts COD entries to prioritize named structures."""

    def get_cod_sort_key(entry: Dict[str, Any]) -> Tuple[int, int]:
        attrs = entry.get("attributes", {})
        if attrs.get("_cod_commonname"):
            priority = 0
        elif attrs.get("_cod_mineral") or attrs.get("_cod_chemname"):
            priority = 1
        else:
            priority = 2
        try:
            entry_id = int(entry.get("id") or 0)
        except ValueError:
            entry_id = 0
        return (priority, entry_id)

    return sorted(data, key=get_cod_sort_key)


def _query_single_database(
    db_name: str, base_url: str, formula: str, hill_formula: str
) -> str:
    """Helper to query a single OPTIMADE database and return sorted results as a string."""
    out: List[str] = [f"Querying {db_name} OPTIMADE for '{formula}'..."]
    try:
        if db_name == "Materials Project":
            url = f"{base_url}?filter=chemical_formula_reduced=%22{hill_formula}%22&page_limit=100"
        else:
            url = f"{base_url}?filter=chemical_formula_hill=%22{hill_formula}%22&page_limit=100"
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            out.append(f"{db_name} returned status code: {r.status_code}")
            return "\n".join(out)

        data = r.json().get("data", [])
        if db_name == "Materials Project":
            data = _sort_materials_project_data(data)
        elif db_name == "Crystallography Open Database (COD)":
            data = _sort_cod_data(data)

        out.append(f"Found {len(data)} structures in {db_name}:")
        for idx, entry in enumerate(data[:5]):
            out.append(_format_structure_info(idx, entry, formula))
    except Exception as e:  # pylint: disable=broad-exception-caught
        out.append(f"Error querying {db_name}: {e}")
    return "\n".join(out)


def retrieve_cif_by_formula(formula: str) -> str:
    """
    Queries OPTIMADE servers for the given chemical formula and returns the results.
    """
    try:
        hill_formula = Composition(formula).hill_formula.replace(" ", "")
    except Exception:  # pylint: disable=broad-exception-caught
        hill_formula = formula

    databases = {
        "Materials Project": "https://optimade.materialsproject.org/v1/structures",
        "Crystallography Open Database (COD)": "https://www.crystallography.net/cod/optimade/v1/structures",
    }

    results: List[str] = []
    for db_name, base_url in databases.items():
        results.append(_query_single_database(db_name, base_url, formula, hill_formula))
    return "\n\n".join(results)

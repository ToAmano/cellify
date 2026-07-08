"""
OPTIMADE API client and query handlers for cellify.
Allows querying materials databases by chemical formula.
"""

from typing import Any, Dict, List, Optional

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


def _query_single_database(
    db_name: str, base_url: str, formula: str, hill_formula: str
) -> None:
    """Helper to query a single OPTIMADE database and display sorted results."""
    print(f"Querying {db_name} OPTIMADE for '{formula}'...")
    try:
        if db_name == "Materials Project":
            url = f"{base_url}?filter=chemical_formula_reduced=%22{hill_formula}%22"
        else:
            url = f"{base_url}?filter=chemical_formula_hill=%22{hill_formula}%22&page_limit=5"
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            print(f"{db_name} returned status code: {r.status_code}")
            return

        data = r.json().get("data", [])
        if db_name == "Materials Project":
            # We sort Materials Project structures by thermodynamic stability (energy_above_hull)
            # in ascending order. This guarantees that the most stable ground-state phases
            # are positioned at the beginning of the list and not sliced out by the top 5 limit.
            def get_sort_key(entry: Dict[str, Any]) -> float:
                attrs = entry.get("attributes", {})
                val = _extract_energy_above_hull(attrs)
                if isinstance(val, (int, float)):
                    return float(val)
                return float("inf")

            data = sorted(data, key=get_sort_key)

        print(f"Found {len(data)} structures in {db_name}:")
        for idx, entry in enumerate(data[:5]):
            print(_format_structure_info(idx, entry, formula))
    except Exception as e:  # pylint: disable=broad-exception-caught
        print(f"Error querying {db_name}: {e}")


def retrieve_cif_by_formula(formula: str) -> None:
    """
    Queries OPTIMADE servers for the given chemical formula and prints the results.
    """
    try:
        hill_formula = Composition(formula).hill_formula.replace(" ", "")
    except Exception:  # pylint: disable=broad-exception-caught
        hill_formula = formula

    databases = {
        "Materials Project": "https://optimade.materialsproject.org/v1/structures",
        "Crystallography Open Database (COD)": "https://www.crystallography.net/cod/optimade/v1/structures",
    }

    for db_name, base_url in databases.items():
        _query_single_database(db_name, base_url, formula, hill_formula)

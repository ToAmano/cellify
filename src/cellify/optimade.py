"""
OPTIMADE API client and query handlers for cellify.
Allows querying materials databases by chemical formula.
"""

from typing import Any, Callable, Dict, List, Optional, Tuple

import requests  # type: ignore[import-untyped]
from pymatgen.core import Composition, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class SelectionError(ValueError):
    """Exception raised when structure selection fails during OPTIMADE query."""

    @property
    def summary(self) -> str:
        """Returns the search summary table."""
        return str(self.args[1])

    def __str__(self) -> str:
        return str(self.args[0])


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
    """Sorts Materials Project entries by thermodynamic stability.

    We sort Materials Project structures by their energy above hull (most stable first).
    Unstable or unmeasured structures with missing or invalid values are assigned a
    fallback value of float("inf") so that they are pushed to the bottom of the list.
    """

    def get_sort_key(entry: Dict[str, Any]) -> float:
        attrs = entry.get("attributes", {})
        val = _extract_energy_above_hull(attrs)
        if isinstance(val, (int, float)):
            return float(val)
        # Fallback to float("inf") for entries without a valid energy_above_hull
        # to push unstable/unmeasured structures to the bottom of the list.
        return float("inf")

    return sorted(data, key=get_sort_key)


def _sort_cod_data(data: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Sorts COD entries to prioritize named structures.

    We sort COD structures to prioritize named ones (commonname, mineral, or chemname)
    over unnamed ones so that users get descriptive crystal names (e.g., Graphite,
    Diamond) first. Within each priority group, entries are sorted by their database
    entry ID in ascending order.
    """

    def get_cod_sort_key(entry: Dict[str, Any]) -> Tuple[int, int]:
        attrs = entry.get("attributes", {})
        # We sort COD structures to prioritize named ones so that users get
        # descriptive crystal names (Graphite, Diamond, Lonsdaleite) first.
        # Priority: commonname (0) > mineral or chemname (1) > unnamed (2).
        # Sub-sorted by ascending database ID (older, more standard entries first).
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


def retrieve_cif_by_formula(formula: str) -> str:
    mp_data, cod_data = query_formula_structures(formula)
    summary, _ = format_formula_structures(mp_data, cod_data, formula)
    return summary


def query_formula_structures(
    formula: str,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """Queries OPTIMADE databases and returns the raw entries for Materials Project and COD.

    Currently, we only support querying the Materials Project and Crystallography
    Open Database (COD) via their respective OPTIMADE endpoints.
    """
    try:
        hill_formula = Composition(formula).hill_formula.replace(" ", "")
    except Exception:  # pylint: disable=broad-exception-caught
        hill_formula = formula

    mp_url = f"https://optimade.materialsproject.org/v1/structures?filter=chemical_formula_reduced=%22{hill_formula}%22&page_limit=100"
    cod_url = f"https://www.crystallography.net/cod/optimade/v1/structures?filter=chemical_formula_hill=%22{hill_formula}%22&page_limit=100"

    mp_data: List[Dict[str, Any]] = []
    try:
        r = requests.get(mp_url, timeout=10)
        if r.status_code == 200:
            mp_data = _sort_materials_project_data(r.json().get("data", []))
    except Exception:  # pylint: disable=broad-exception-caught
        pass

    cod_data: List[Dict[str, Any]] = []
    try:
        r = requests.get(cod_url, timeout=10)
        if r.status_code == 200:
            cod_data = _sort_cod_data(r.json().get("data", []))
    except Exception:  # pylint: disable=broad-exception-caught
        pass

    return mp_data, cod_data


def format_formula_structures(
    mp_data: List[Dict[str, Any]], cod_data: List[Dict[str, Any]], formula: str
) -> Tuple[str, Dict[int, Tuple[str, Dict[str, Any]]]]:
    """Formats the structures query results and constructs a selection mapping.

    Note that we limit the displayed list to the top 5 candidate entries from each
    database (Materials Project and COD) to prevent terminal clutter and keep the
    interactive choice selection manageable for the user.
    """
    out: List[str] = []
    selection_map: Dict[int, Tuple[str, Dict[str, Any]]] = {}
    current_idx = 1

    out.append(f"Querying Materials Project OPTIMADE for '{formula}'...")
    out.append(f"Found {len(mp_data)} structures in Materials Project:")
    for entry in mp_data[:5]:
        desc = _format_structure_info(current_idx - 1, entry, formula)
        out.append(desc)
        selection_map[current_idx] = ("Materials Project", entry)
        current_idx += 1
    if len(mp_data) > 5:
        out.append(
            f"  Warning: Only the top 5 most relevant structures are shown. "
            f"There are {len(mp_data) - 5} more structures in Materials Project."
        )

    out.append(
        f"\nQuerying Crystallography Open Database (COD) OPTIMADE for '{formula}'..."
    )
    out.append(
        f"Found {len(cod_data)} structures in Crystallography Open Database (COD):"
    )
    for entry in cod_data[:5]:
        desc = _format_structure_info(current_idx - 1, entry, formula)
        out.append(desc)
        selection_map[current_idx] = ("Crystallography Open Database (COD)", entry)
        current_idx += 1
    if len(cod_data) > 5:
        out.append(
            f"  Warning: Only the top 5 most relevant structures are shown. "
            f"There are {len(cod_data) - 5} more structures in Crystallography Open Database (COD)."
        )

    return "\n".join(out), selection_map


def download_structure_from_entry(db_name: str, entry: Dict[str, Any]) -> Structure:
    """Downloads and parses the crystal structure corresponding to the selection entry."""
    # pylint: disable=no-else-return
    if db_name == "Materials Project":
        struct = parse_optimade_entry_to_structure(entry)
        if struct is None:
            raise ValueError(
                "Failed to parse structure from Materials Project attributes."
            )
        return struct
    elif db_name == "Crystallography Open Database (COD)":
        entry_id = entry.get("id")
        if not entry_id:
            raise ValueError("COD entry has no ID.")
        url = f"https://www.crystallography.net/cod/{entry_id}.cif"
        r = requests.get(url, timeout=10)
        if r.status_code != 200:
            raise RuntimeError(
                f"Failed to download CIF from COD: status {r.status_code}"
            )
        return Structure.from_str(r.text, fmt="cif")
    else:
        raise ValueError("Unsupported database name: " + db_name)


def select_and_download_structure(
    formula: str,
    select: Optional[int] = None,
    interactive_prompt: Optional[Callable[[str, int], int]] = None,
) -> Tuple[Structure, str, Dict[str, Any], str]:
    """Unified function to query OPTIMADE databases for a formula, format/index the results,

    handle interactive/non-interactive selection, and download the chosen structure.

    Args:
        formula: The chemical formula of the target crystal structure (e.g., "C").
        select: An optional 1-indexed selection index to bypass the interactive prompt.
        interactive_prompt: An optional callable that accepts the formatted summary text
            and the maximum selection choice, returning the user's selected choice.

    Returns:
        A tuple of:
            - structure: The downloaded pymatgen Structure object.
            - db_name: The name of the database the structure was sourced from.
            - entry: The raw OPTIMADE response dictionary for the selected entry.
            - summary: The text table summarizing all candidate structures.
    """
    mp_data, cod_data = query_formula_structures(formula)
    # The summary is a formatted text table of candidate structures from Materials Project and COD
    summary, selection_map = format_formula_structures(mp_data, cod_data, formula)

    if not selection_map:
        raise SelectionError(f"No structures found for formula '{formula}'.", summary)

    if select is not None:
        choice = select
    elif interactive_prompt is not None:
        choice = interactive_prompt(summary, len(selection_map))
    else:
        raise SelectionError(
            "Chemical formula specified, but the session is non-interactive "
            "and --select was not provided.",
            summary,
        )

    if choice not in selection_map:
        raise SelectionError(f"Selection index {choice} is out of range.", summary)

    db_name, entry = selection_map[choice]
    structure = download_structure_from_entry(db_name, entry)
    return structure, db_name, entry, summary

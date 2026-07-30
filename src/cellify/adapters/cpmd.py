"""
CPMD (Car-Parrinello Molecular Dynamics) input adapter for cellify.
"""

import json
import os
import re
from typing import Any, Dict, List, Tuple

import numpy as np
from pymatgen.core import Lattice, Structure

from cellify.adapters.base import BaseAdapter


class CpmdAdapter(BaseAdapter):
    """
    CPMD input file adapter.
    Preserves calculation parameters (&INFO, &CPMD, &SYSTEM, &DFT, &VDW, etc.),
    automatically sorts atomic species to avoid CPMD 'TOO MANY ATOMIC SPECIES' errors,
    generates a JSON index map sidecar file (<filepath>.index_map.json),
    and embeds inverse index mapping into the input header.
    """

    def read(self, filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Input file not found: {filepath}")

        with open(filepath, "r", encoding="utf-8") as f:
            content: str = f.read()

        # Parse lattice and atoms from CPMD content
        lattice: Lattice = self._parse_cpmd_lattice(content)
        species, coords = self._parse_cpmd_atoms(content)

        if not species or len(species) == 0:
            # Fallback to dummy atom if reading a template file with empty &ATOMS section
            species = ["H"]
            coords = [[0.0, 0.0, 0.0]]

        structure: Structure = Structure(
            lattice, species, coords, coords_are_cartesian=True
        )

        meta_data: Dict[str, Any] = {
            "mode": "cpmd_text_replace",
            "content": content,
            "filepath": filepath,
        }
        return structure, meta_data

    def write(
        self, filepath: str, structure: Structure, meta_data: Dict[str, Any]
    ) -> None:
        content: str = meta_data.get("content", "")
        if not content:
            # Fallback to standard default CPMD georelax template
            template_path = os.path.join(
                os.path.dirname(__file__),
                "..",
                "templates",
                "cpmd",
                "georelax.inp.tpl",
            )
            if os.path.exists(template_path):
                with open(template_path, "r", encoding="utf-8") as f:
                    content = f.read()
            else:
                content = self._get_minimal_cpmd_template()

        # 1. Atom sorting by chemical species
        original_symbols: List[str] = [site.specie.symbol for site in structure]
        original_coords: List[np.ndarray] = [site.coords for site in structure]

        # Create enumerated tuples (orig_idx, symbol, coords)
        items = list(enumerate(zip(original_symbols, original_coords)))
        # Sort stably by chemical symbol
        sorted_items = sorted(items, key=lambda item: item[1][0])

        sort_indices: List[int] = [item[0] for item in sorted_items]
        sorted_symbols: List[str] = [item[1][0] for item in sorted_items]
        sorted_coords: List[np.ndarray] = [item[1][1] for item in sorted_items]

        # Compute inverse_indices mapping: sorted_index -> original_index
        inverse_indices: List[int] = list(range(len(structure)))
        for sorted_pos, (orig_pos, _) in enumerate(sorted_items):
            inverse_indices[sorted_pos] = orig_pos

        # 2. Write JSON sidecar index map file
        json_path: str = f"{filepath}.index_map.json"
        index_map_data: Dict[str, Any] = {
            "source_filepath": meta_data.get("filepath", ""),
            "output_filepath": filepath,
            "num_atoms": len(structure),
            "original_symbols": original_symbols,
            "sorted_symbols": sorted_symbols,
            "sort_indices": sort_indices,
            "inverse_indices": inverse_indices,
        }
        with open(json_path, "w", encoding="utf-8") as f_json:
            json.dump(index_map_data, f_json, indent=2)

        # 3. Update CELL parameters inside &SYSTEM
        content = self._update_cpmd_cell(content, structure.lattice)

        # 4. Embed ATOM_MAP_INVERSE metadata comment in &INFO or top header
        content = self._embed_atom_map_comment(content, inverse_indices)

        # 5. Build new &ATOMS section
        content = self._update_cpmd_atoms(
            content, sorted_symbols, sorted_coords
        )

        with open(filepath, "w", encoding="utf-8") as f:
            f.write(content)

    def _parse_cpmd_lattice(self, content: str) -> Lattice:
        """Parses cell lattice parameters from &SYSTEM block."""
        is_angstrom: bool = bool(re.search(r"ANGSTROM", content, re.IGNORECASE))
        scale_factor: float = 1.0 if is_angstrom else 0.529177210903

        # Match CELL block under &SYSTEM
        cell_match = re.search(
            r"CELL\s*\n\s*([^\n]+)", content, re.IGNORECASE
        )
        if cell_match:
            tokens: List[float] = [
                float(x) for x in cell_match.group(1).strip().split()
            ]
            if len(tokens) >= 6:
                # CPMD CELL format: a, b/a, c/a, cos(alpha), cos(beta), cos(gamma)
                a_val: float = tokens[0] * scale_factor
                b_val: float = tokens[1] * a_val
                c_val: float = tokens[2] * a_val
                alpha: float = np.arccos(tokens[3]) * 180.0 / np.pi if abs(tokens[3]) <= 1.0 else 90.0
                beta: float = np.arccos(tokens[4]) * 180.0 / np.pi if abs(tokens[4]) <= 1.0 else 90.0
                gamma: float = np.arccos(tokens[5]) * 180.0 / np.pi if abs(tokens[5]) <= 1.0 else 90.0
                return Lattice.from_parameters(a_val, b_val, c_val, alpha, beta, gamma)
            elif len(tokens) == 1:
                a_val = tokens[0] * scale_factor
                return Lattice.cubic(a_val)

        # Fallback default lattice if not explicitly parsed
        return Lattice.cubic(10.0)

    def _parse_cpmd_atoms(self, content: str) -> Tuple[List[str], List[List[float]]]:
        """Parses atomic species and cartesian coordinates from &ATOMS block."""
        is_angstrom: bool = bool(re.search(r"ANGSTROM", content, re.IGNORECASE))
        scale_factor: float = 1.0 if is_angstrom else 0.529177210903

        atoms_match = re.search(
            r"&ATOMS\s*\n(.*?)(?=\n\s*&|\Z)", content, re.DOTALL | re.IGNORECASE
        )
        if not atoms_match:
            return [], []

        atoms_block: str = atoms_match.group(1)
        lines: List[str] = atoms_block.splitlines()

        species: List[str] = []
        coords: List[List[float]] = []

        curr_symbol: str = "X"
        lines_to_read: int = 0

        for line in lines:
            stripped: str = line.strip()
            if not stripped or stripped.startswith("!"):
                continue

            if stripped.startswith("*"):
                # Pseudopotential header, e.g. *pseudo/O_MT_GIA_BLYP KLEINMAN-BYLANDER
                pseudo_name: str = stripped.split()[0]
                # Infer element symbol from pseudo header
                match_elem = re.search(r"([A-Z][a-z]?)", pseudo_name.replace("*pseudo/", "").replace("*", ""))
                curr_symbol = match_elem.group(1) if match_elem else "X"
                lines_to_read = 0
                continue

            if stripped.upper().startswith("LMAX"):
                continue

            # Number of atoms for this species
            if lines_to_read == 0 and stripped.isdigit():
                lines_to_read = int(stripped)
                continue

            # Coordinate line
            parts: List[str] = stripped.split()
            if len(parts) >= 3:
                try:
                    c_xyz: List[float] = [
                        float(parts[0]) * scale_factor,
                        float(parts[1]) * scale_factor,
                        float(parts[2]) * scale_factor,
                    ]
                    species.append(curr_symbol)
                    coords.append(c_xyz)
                    lines_to_read -= 1
                except ValueError:
                    pass

        return species, coords

    def _update_cpmd_cell(self, content: str, lattice: Lattice) -> str:
        """Updates CELL vector line under &SYSTEM block."""
        a, b, c = lattice.a, lattice.b, lattice.c
        alpha, beta, gamma = lattice.alpha, lattice.beta, lattice.gamma

        # Ensure ANGSTROM keyword exists in &SYSTEM block
        if not re.search(r"ANGSTROM", content, re.IGNORECASE):
            content = re.sub(
                r"(&SYSTEM\b)",
                r"\1\n   ANGSTROM",
                content,
                flags=re.IGNORECASE,
            )

        cos_alpha: float = float(np.cos(np.radians(alpha)))
        cos_beta: float = float(np.cos(np.radians(beta)))
        cos_gamma: float = float(np.cos(np.radians(gamma)))

        b_over_a: float = b / a if a > 0 else 1.0
        c_over_a: float = c / a if a > 0 else 1.0

        cell_line: str = (
            f"   CELL\n         {a:.10f} {b_over_a:.10f} {c_over_a:.10f}  "
            f"{cos_alpha:.10f}  {cos_beta:.10f}  {cos_gamma:.10f}"
        )

        pattern = r"CELL\s*\n\s*[^\n]+"
        if re.search(pattern, content, re.IGNORECASE):
            content = re.sub(pattern, cell_line, content, flags=re.IGNORECASE)
        else:
            content = re.sub(
                r"(&SYSTEM\b[^\n]*)",
                r"\1\n" + cell_line,
                content,
                flags=re.IGNORECASE,
            )

        return content

    def _embed_atom_map_comment(self, content: str, inverse_indices: List[int]) -> str:
        """Embeds inverse atom index mapping comment into the &INFO section or header."""
        comment_str: str = f"   ! ATOM_MAP_INVERSE: {inverse_indices}"
        if re.search(r"&INFO", content, re.IGNORECASE):
            content = re.sub(
                r"(&INFO\b[^\n]*)",
                r"\1\n" + comment_str,
                content,
                flags=re.IGNORECASE,
            )
        else:
            content = f"{comment_str}\n" + content
        return content

    def _update_cpmd_atoms(
        self, content: str, sorted_symbols: List[str], sorted_coords: List[np.ndarray]
    ) -> str:
        """Constructs and updates the &ATOMS block grouped by species."""
        # Group coordinates by chemical symbol
        grouped: Dict[str, List[np.ndarray]] = {}
        for symbol, coord in zip(sorted_symbols, sorted_coords):
            grouped.setdefault(symbol, []).append(coord)

        atoms_str: str = "&ATOMS\n"
        for symbol, coords_list in grouped.items():
            lmax_str: str = "LMAX=S" if symbol == "H" else "LMAX=P"
            atoms_str += f"*pseudo/{symbol}_BLYP.psp KLEINMAN-BYLANDER\n"
            atoms_str += f"   {lmax_str}\n"
            atoms_str += f"   {len(coords_list)}\n"
            for c in coords_list:
                atoms_str += f"  {c[0]:.10f}  {c[1]:.10f}  {c[2]:.10f}\n"
            atoms_str += "\n"
        atoms_str += "&END"

        pattern = r"&ATOMS\b.*?(?=\Z)"
        if re.search(pattern, content, re.DOTALL | re.IGNORECASE):
            content = re.sub(pattern, atoms_str, content, flags=re.DOTALL | re.IGNORECASE)
        else:
            content = content.rstrip() + "\n\n" + atoms_str

        return content

    def _get_minimal_cpmd_template(self) -> str:
        return """&INFO
   Input generated by cellify
&END

&CPMD
   OPTIMIZE GEOMETRY XYZ
   MIRROR
   FILEPATH
         ./tmp
&END

&SYSTEM
   SYMMETRY
         1
   ANGSTROM
   CELL
         10.0 1.0 1.0 0.0 0.0 0.0
   CUTOFF
         100.0
&END

&DFT
   FUNCTIONAL BLYP
&END

&ATOMS
&END
"""

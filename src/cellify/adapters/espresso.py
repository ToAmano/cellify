"""
Quantum ESPRESSO input/output adapter for cellify.
"""

import os
import re
from typing import Any, Dict, Tuple

from pymatgen.core import Structure

from cellify.adapters.base import BaseAdapter


class EspressoAdapter(BaseAdapter):
    """
    Quantum ESPRESSO input file adapter.
    Preserves calculation parameters (&CONTROL, &SYSTEM, etc.) and comment lines,
    while automatically updating nat/ntyp and replacing structure sections.
    """

    def read(self, filepath: str) -> Tuple[Structure, Dict[str, Any]]:
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Input file not found: {filepath}")

        with open(filepath, "r", encoding="utf-8") as f:
            content: str = f.read()

        # Check if it looks like a QE output log file (missing &control namelist)
        is_output: bool = not re.search(r"&control", content, re.IGNORECASE)

        # Safely parse structure using ASE espresso-in or espresso-out reader
        try:
            # pylint: disable=import-outside-toplevel
            from ase.io import read as ase_read
            from pymatgen.io.ase import AseAtomsAdaptor

            if is_output:
                atoms = ase_read(filepath, format="espresso-out", index=-1)
            else:
                atoms = ase_read(filepath, format="espresso-in")
            structure: Structure = AseAtomsAdaptor.get_structure(atoms)
        except Exception as ase_err:  # pylint: disable=broad-exception-caught
            # Fallback to try the other format
            try:
                if is_output:
                    atoms = ase_read(filepath, format="espresso-in")
                    is_output = False
                else:
                    atoms = ase_read(filepath, format="espresso-out", index=-1)
                    is_output = True
                structure = AseAtomsAdaptor.get_structure(atoms)
            except Exception as fallback_err:  # pylint: disable=broad-exception-caught
                raise ValueError(
                    f"Failed to parse structure from Quantum ESPRESSO file: {ase_err} (fallback: {fallback_err})"
                ) from fallback_err

        meta_data: Dict[str, Any] = {
            "mode": "espresso_out" if is_output else "espresso_text_replace",
            "content": content,
            "filepath": filepath,
        }
        return structure, meta_data

    def write(
        self, filepath: str, structure: Structure, meta_data: Dict[str, Any]
    ) -> None:
        if meta_data.get("mode") == "espresso_out":
            raise ValueError(
                "Cannot write a QE input file using a QE output log as a template. Please provide a template input file."
            )

        content: str = meta_data.get("content", "")
        if not content:
            # Generate a minimal default QE input template if original content is missing
            content = """&CONTROL
  calculation = 'scf'
  restart_mode = 'from_scratch'
  pseudo_dir = './'
  outdir = './'
/
&SYSTEM
  ibrav = 0
  nat = 0
  ntyp = 0
/
&ELECTRONS
/

ATOMIC_SPECIES

CELL_PARAMETERS angstrom

ATOMIC_POSITIONS crystal
"""

        # 1. Calculate new nat and ntyp
        nat_new: int = len(structure)
        ntyp_new: int = len(structure.composition.elements)

        # 2. Update calculation if specified
        calculation_override = meta_data.get("calculation")
        if calculation_override:
            content = re.sub(
                r"(\bcalculation\s*=\s*['\"]?)[a-zA-Z0-9_-]+([\'\"]?)",
                r"\g<1>" + calculation_override + r"\2",
                content,
                flags=re.IGNORECASE,
            )

        # 3. Update nat and ntyp inside namelists
        content = re.sub(
            r"(\bnat\s*=\s*)\d+", r"\g<1>" + str(nat_new), content, flags=re.IGNORECASE
        )
        content = re.sub(
            r"(\bntyp\s*=\s*)\d+",
            r"\g<1>" + str(ntyp_new),
            content,
            flags=re.IGNORECASE,
        )

        # 4. Strip old structure-related blocks from text
        cleaned_content: str = content
        struct_keywords = ["ATOMIC_SPECIES", "CELL_PARAMETERS", "ATOMIC_POSITIONS"]
        for kw in struct_keywords:
            pattern = (
                r"(?i)^\s*"
                + kw
                + r"\b.*?(?=\n\s*(?:ATOMIC_SPECIES|CELL_PARAMETERS|ATOMIC_POSITIONS|K_POINTS|KPOINTS|&[A-Za-z]+)|\Z)"
            )
            cleaned_content = re.sub(
                pattern, "", cleaned_content, flags=re.DOTALL | re.MULTILINE
            )

        # Clean extra leading/trailing whitespaces
        cleaned_content = cleaned_content.strip() + "\n\n"

        # 4. Extract existing pseudopotential information from the original file
        pseudos: Dict[str, Any] = {}
        species_match = re.search(
            r"(?i)ATOMIC_SPECIES\s*\n(.*?)(?=\n\s*(?:ATOMIC_|CELL_|K_POINTS|KPOINTS|&[A-Za-z]+)|\Z)",
            content,
            re.DOTALL,
        )
        if species_match:
            for line in species_match.group(1).strip().split("\n"):
                parts = line.split()
                if len(parts) >= 3:
                    pseudos[parts[0]] = (parts[1], parts[2])

        # 5. Reconstruct structure blocks
        # ATOMIC_SPECIES
        species_str: str = "ATOMIC_SPECIES\n"
        for el in structure.composition.elements:
            el_symbol: str = el.symbol
            mass, pseudo = pseudos.get(
                el_symbol, (str(el.atomic_mass), f"{el_symbol}.UPF")
            )
            species_str += f"  {el_symbol}  {mass}  {pseudo}\n"

        # CELL_PARAMETERS
        cell_str: str = "\nCELL_PARAMETERS angstrom\n"
        for vec in structure.lattice.matrix:
            cell_str += f"  {vec[0]:.10f}  {vec[1]:.10f}  {vec[2]:.10f}\n"

        # ATOMIC_POSITIONS
        pos_str: str = "\nATOMIC_POSITIONS crystal\n"
        for site in structure:
            pos_str += (
                f"  {site.specie.symbol}  {site.a:.10f}  {site.b:.10f}  {site.c:.10f}\n"
            )

        # 6. Save file
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(cleaned_content)
            f.write(species_str)
            f.write(cell_str)
            f.write(pos_str)

import cairosvg
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

import os
from collections import defaultdict

from openmmdl.openmmdl_analysis.core.utils import extract_ints


class FigureHighlighter:
    """
    Identifies and highlights ligand atoms involved in various types of protein-ligand interactions.

    Attributes
    ----------
    complex_pdb_file : str
        Path to the protein-ligand complex PDB file.
    ligand_no_h_pdb_file : str
        Path to the ligand PDB file without hydrogens.
    complex : mda.Universe
        MDAnalysis Universe object of the protein-ligand complex.
    ligand_no_h : mda.Universe
        MDAnalysis Universe object of the ligand without hydrogens.
    lig_noh : mda.AtomGroup
        AtomGroup of all atoms in the ligand without hydrogens.
    """

    def __init__(self, complex_pdb_file, ligand_no_h_pdb_file, ligand_name, mapping_topology_file=None):
        self.complex_pdb_file = complex_pdb_file
        self.ligand_no_h_pdb_file = ligand_no_h_pdb_file
        self.ligand_name = ligand_name

        # Universe used for drawing/selection (can be reduced/renumbered)
        self.complex = mda.Universe(complex_pdb_file)

        # Universe used for interpreting ProLIF atom codes (MUST match ProLIF topology)
        self.map_u = mda.Universe(mapping_topology_file or complex_pdb_file)

        self._map_lig_heavy = self.map_u.select_atoms(f"resname {ligand_name} and not name H*")

        self._id_to_ligidx = {int(a.id): i for i, a in enumerate(self._map_lig_heavy)}
        self._index_to_ligidx = {int(a.index): i for i, a in enumerate(self._map_lig_heavy)}
        self._index1_to_ligidx = {int(a.index) + 1: i for i, a in enumerate(self._map_lig_heavy)}

        # Name-based fallback mapping (atom names like C2, N2, O1, etc.)
        name_to_idxs = defaultdict(list)
        for i, a in enumerate(self._map_lig_heavy):
            name_to_idxs[str(a.name)].append(i)

        self._name_to_ligidx = {}
        for name, idxs in name_to_idxs.items():
            # If duplicates exist, keep the first and optionally warn
            self._name_to_ligidx[name] = idxs[0]

    def _build_code_to_ligidx(self):
        """
        Build mapping from various "codes" (Atom.index+1, Atom.id, Atom.serial if present)
        to ligand-local heavy-atom indices (0..n_heavy-1), which should match RDKit
        indices after RemoveAllHs on a ligand RDKit mol created from the same selection.
        """
        lig_heavy = self.complex.select_atoms(
            f"resname {self.ligand_name} and not name H*"
        )

        code_to_ligidx = {}
        for lig_idx, a in enumerate(lig_heavy):
            candidates = set()

            # Common ProLIF-style identifiers
            try:
                candidates.add(int(a.index) + 1)  # global 1-based index
            except Exception:
                pass
            try:
                candidates.add(int(a.id))         # topology id (often PDB serial, but not always)
            except Exception:
                pass
            # Some topologies expose serial explicitly
            if hasattr(a, "serial"):
                try:
                    candidates.add(int(a.serial))
                except Exception:
                    pass

            for c in candidates:
                code_to_ligidx[c] = lig_idx

        return code_to_ligidx

    def _lig_index_from_complex_atom_id(self, atom_code):
        try:
            return self._code_to_ligidx.get(int(atom_code))
        except Exception:
            return None

    def _tok_to_ligidxs(self, tok: str):
        t = tok.strip(",;")
        tl = t.lower()
        if tl in ("donor", "acceptor", "pi", "ni", "cation", "anion", "positive", "negative"):
            return []

        # If token contains letters, try atom-name mapping first (PLIP often does this)
        if any(c.isalpha() for c in t):
            idx = self._name_to_ligidx.get(t)
            return [] if idx is None else [idx]

        # Otherwise treat as numeric codes
        out = []
        for atom_id in extract_ints(t):
            idx = self._lig_index_from_complex_code(atom_id)
            if idx is not None:
                out.append(idx)
        return out

    def _lig_index_from_complex_code(self, code_int: int):
        # 1) Direct lookups in mapping universe ligand atoms
        if code_int in self._id_to_ligidx:
            return self._id_to_ligidx[code_int]
        if code_int in self._index_to_ligidx:
            return self._index_to_ligidx[code_int]
        if code_int in self._index1_to_ligidx:
            return self._index1_to_ligidx[code_int]

        # 2) Fallback: resolve the atom in the *mapping universe* and map by atom name
        u = self.map_u
        atom = None

        # Try as universe index (0-based and 1-based)
        for candidate in (code_int, code_int - 1):
            if 0 <= candidate < len(u.atoms):
                atom = u.atoms[candidate]
                break

        # Try as topology id (often PDB serial)
        if atom is None:
            sel = u.select_atoms(f"id {code_int}")
            if len(sel) == 1:
                atom = sel[0]

        if atom is None:
            return None

        return self._name_to_ligidx.get(atom.name)

    def split_interaction_data(self, data):
        """
        Splits the input data into multiple parts.

        Parameters
        ----------
        data : list of str
            A list containing strings with protein residue name, interacting indices and interaction type.

        Returns
        -------
        list of str
            List of separate formatted strings with separated protein residue name, interacting indices and interaction type.
        """
        split_data = []
        for item in data:
            parts = item.split("_")
            protein_partner_name = parts[0]
            numeric_codes = " ".join(parts[1:-1])
            interaction_type = parts[-1]
            split_value = f"{protein_partner_name} {numeric_codes} {interaction_type}"
            split_data.append(split_value)

        return split_data

    def highlight_numbers(self, split_data, starting_idx):
        """
        Parse split interaction strings and return per-interaction-type ligand-local
        heavy-atom indices (0..n_heavy-1) suitable for RDKit highlighting.

        IMPORTANT:
        - This version assumes self._tok_to_ligidxs(tok) returns ligand indices.
        - Therefore we MUST NOT call _lig_index_from_complex_code() again inside this method
        (that was the double-mapping bug that turned valid codes into 0..24 and then None).
        """
        highlighted_hbond_acceptor = set()
        highlighted_hbond_donor = set()
        highlighted_hydrophobic = set()
        highlighted_hbond_both = set()
        highlighted_waterbridge = set()
        highlighted_pistacking = set()
        highlighted_halogen = set()
        highlighted_ni = set()
        highlighted_pi = set()
        highlighted_pication = set()
        highlighted_metal = set()

        for item in split_data:
            parts = item.split()
            if len(parts) < 3:
                continue

            interaction_type = parts[-1].lower()
            tokens = parts[1:-1]  # exclude protein partner name and interaction_type word

            # HBond format: "<RES> <atomcodes...> <Donor|Acceptor> hbond"
            if interaction_type == "hbond":
                label = None
                for t in tokens:
                    tl = t.lower()
                    if tl in ("donor", "acceptor"):
                        label = tl
                        break

                for tok in tokens:
                    if tok.lower() in ("donor", "acceptor"):
                        continue
                    for lig_idx in self._tok_to_ligidxs(tok):
                        if label == "donor":
                            highlighted_hbond_donor.add(lig_idx)
                        elif label == "acceptor":
                            highlighted_hbond_acceptor.add(lig_idx)
                        else:
                            highlighted_hbond_both.add(lig_idx)

            elif interaction_type == "hydrophobic":
                for tok in tokens:
                    for lig_idx in self._tok_to_ligidxs(tok):
                        highlighted_hydrophobic.add(lig_idx)

            elif interaction_type == "waterbridge":
                # tokens may include "Acceptor"/"Donor"; ignore, just map indices
                for tok in tokens:
                    if tok.lower() in ("donor", "acceptor"):
                        continue
                    for lig_idx in self._tok_to_ligidxs(tok):
                        highlighted_waterbridge.add(lig_idx)

            elif interaction_type == "pistacking":
                for tok in tokens:
                    for lig_idx in self._tok_to_ligidxs(tok):
                        highlighted_pistacking.add(lig_idx)

            elif interaction_type == "halogen":
                for tok in tokens:
                    for lig_idx in self._tok_to_ligidxs(tok):
                        highlighted_halogen.add(lig_idx)

            elif interaction_type == "saltbridge":
                # ProLIF/PLIP labels vary; determine NI vs PI if present; default to PI
                sb_type = None
                for t in tokens:
                    tu = t.upper()
                    if tu in ("NI", "PI"):
                        sb_type = tu
                        break
                    if tu in ("CATION", "POSITIVE"):
                        sb_type = "PI"
                    if tu in ("ANION", "NEGATIVE"):
                        sb_type = "NI"

                target_set = highlighted_pi if sb_type != "NI" else highlighted_ni

                for tok in tokens:
                    if tok.lower() in ("cation", "anion", "positive", "negative", "pi", "ni"):
                        continue
                    for lig_idx in self._tok_to_ligidxs(tok):
                        target_set.add(lig_idx)

            elif interaction_type == "pication":
                for tok in tokens:
                    for lig_idx in self._tok_to_ligidxs(tok):
                        highlighted_pication.add(lig_idx)

            elif interaction_type == "metal":
                # safest: take first ligand index we can extract from tokens
                found = None
                found_tok = None
                for tok in tokens:
                    idxs = list(self._tok_to_ligidxs(tok))
                    if idxs:
                        found = idxs[0]
                        found_tok = tok
                        break
                if found is not None:
                    highlighted_metal.add(found)

            # else: ignore unknown interaction types silently

        # Merge donor/acceptor overlap into "both"
        overlap = highlighted_hbond_donor.intersection(highlighted_hbond_acceptor)
        if overlap:
            highlighted_hbond_donor.difference_update(overlap)
            highlighted_hbond_acceptor.difference_update(overlap)
            highlighted_hbond_both.update(overlap)

        return (
            sorted(highlighted_hbond_donor),
            sorted(highlighted_hbond_acceptor),
            sorted(highlighted_hbond_both),
            sorted(highlighted_hydrophobic),
            sorted(highlighted_waterbridge),
            sorted(highlighted_pistacking),
            sorted(highlighted_halogen),
            sorted(highlighted_ni),
            sorted(highlighted_pi),
            sorted(highlighted_pication),
            sorted(highlighted_metal),
        )


    def generate_interaction_dict(self, interaction_type, keys):
        """
        Generates a dictionary of interaction RGB color model based on the provided interaction type.

        Parameters
        ----------
        interaction_type : str
            The type of interaction (e.g., 'hydrophobic', 'hbond_donor').
        keys : list of int
            Atom indices corresponding to the given interaction type.

        Returns
        -------
        dict
            Dictionary mapping each atom index to an RGB color code tuple.
        """
        interaction_dict = {
            "hbond_acceptor": (1.0, 0.6, 0.6),  # light red / pink
            "hbond_both": (0.6, 0.0, 0.5),  # dark magenta / purple
            "hbond_donor": (0.3, 0.5, 1.0),  # light blue
            "hydrophobic": (1.0, 1.0, 0.0),  # yellow
            "waterbridge": (0.0, 1.0, 0.9),  # cyan / aqua
            "pistacking": (0.0, 0.0, 1.0),  # blue
            "halogen": (1.0, 0.0, 0.9),  # magenta / hot pink
            "ni": (1.0, 0.6, 0.0),  # orange
            "pi": (0.3, 0.9, 0.8),  # turquoise / teal
            "pication": (0.0, 0.0, 1.0),  # blue
            "metal": (1.0, 0.6, 0.0),  # orange
        }

        interaction_dict = {int(key): interaction_dict[interaction_type] for key in keys}

        return interaction_dict


class LigandImageGenerator:
    """
    Generates 2D images of the ligand structure from a protein-ligand complex with atom indices mapped.

    Attributes
    ----------
    ligand_name : str
        Name of the ligand (3 letters) in the protein-ligand complex topology.
    complex_pdb_file : str
        Path to the protein-ligand complex PDB file.
    ligand_no_h_pdb_file : str
        Path to the ligand PDB file without hydrogens.
    output_svg_filename : str
        Output filename for the generated SVG image.
    fig_type : str
        Type of image to generate. Can be "svg" or "png".
    """

    def __init__(
        self,
        ligand_name,
        complex_pdb_file,
        ligand_no_h_pdb_file,
        output_svg_filename,
        fig_type="svg",
    ):
        self.ligand_name = ligand_name
        self.complex_pdb_file = complex_pdb_file
        self.ligand_no_h_pdb_file = ligand_no_h_pdb_file
        self.output_svg_filename = output_svg_filename
        self.fig_type = fig_type

    def generate_image(self):
        """
        Generates an SVG image (or PNG) of the ligand.

        Returns
        -------
        None
            This function writes out a figure and does not return anything.

        Raises
        ------
        Exception
            If any step in the image generation pipeline fails.
        """
        try:
            # Load complex and ligand structures
            complex = mda.Universe(self.complex_pdb_file)
            ligand_no_h = mda.Universe(self.ligand_no_h_pdb_file)
            lig_noh = ligand_no_h.select_atoms("all")
            complex_lig = complex.select_atoms(f"resname {self.ligand_name}")

            # Application of RDKit Converter to obtain rdkit mol of ligand
            lig_atoms = complex_lig.convert_to("RDKIT")
            prepared_ligand = Chem.RemoveAllHs(lig_atoms)

            AllChem.Compute2DCoords(prepared_ligand)

            # Map atom indices between ligand_no_h and complex
            for atom in prepared_ligand.GetAtoms():
                atom_index = atom.GetIdx()
                for lig_atom in lig_noh:
                    lig_index = lig_atom.index
                    if atom_index == lig_index:
                        lig_atom_name = lig_atom.name
                        for comp_lig in complex_lig:
                            comp_lig_name = comp_lig.name
                            if lig_atom_name == comp_lig_name:
                                num = int(comp_lig.id)
                                atom.SetAtomMapNum(num)

            # Generate an SVG image of the ligand
            drawer = Draw.MolDraw2DSVG(5120, 3200)
            drawer.drawOptions().addStereoAnnotation = True  # Add stereo information if available
            drawer.DrawMolecule(prepared_ligand)

            # Adjust font size in the SVG output using the FontSize method
            font_size = drawer.FontSize()
            drawer.SetFontSize(font_size * 0.5)  # You can adjust the multiplier as needed

            drawer.FinishDrawing()
            svg = drawer.GetDrawingText().replace("svg:", "")

            # Save the SVG image to the specified output file
            with open(self.output_svg_filename, "w") as f:
                f.write(svg)

            # Convert to PNG if requested
            if self.fig_type == "png":
                png_filename = self.output_svg_filename.replace(".svg", ".png")
                cairosvg.svg2png(url=self.output_svg_filename, write_to=png_filename)
                print(f"PNG image saved as: {png_filename}")

        except Exception as e:
            print(f"Error: {e}")

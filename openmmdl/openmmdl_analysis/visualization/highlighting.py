import cairosvg
import MDAnalysis as mda
from rdkit import Chem
from rdkit.Chem import AllChem, Draw


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
    def __init__(self, complex_pdb_file, ligand_no_h_pdb_file):
        """
        Initialize the InteractionProcessor class.

        Args:
            complex_pdb_file (str): Path to the protein-ligand complex PDB file.
            ligand_no_h_pdb_file (str): Path to the ligand PDB file without hydrogens.
        """
        self.complex_pdb_file = complex_pdb_file
        self.ligand_no_h_pdb_file = ligand_no_h_pdb_file
        self.complex = mda.Universe(complex_pdb_file)
        self.ligand_no_h = mda.Universe(ligand_no_h_pdb_file)
        self.lig_noh = self.ligand_no_h.select_atoms("all")

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
        Extracts the data from the split_data output of the interactions and categorizes it to the correct interaction list.

        Parameters
        ----------
        split_data : list of str
            Split interaction data strings with the protein residue name, interacting indices and interaction type.
        starting_idx : list of int
            Starting indices of ligand atoms used to correctly identify atoms.

        Returns
        -------
        tuple of lists
            Lists of highlighted atom indices categorized by interaction types.
        """
        highlighted_hbond_acceptor = []
        highlighted_hbond_donor = []
        highlighted_hydrophobic = []
        highlighted_hbond_both = []
        highlighted_waterbridge = []
        highlighted_pistacking = []
        highlighted_halogen = []
        highlighted_ni = []
        highlighted_pi = []
        highlighted_pication = []
        highlighted_metal = []

        for item in split_data:
            parts = item.split()
            protein_partner_name = parts[0]
            numeric_codes = parts[1:-1]
            interaction_type = parts[-1]

            if interaction_type == "hbond":
                parts = item.split()
                protein_partner_name = parts[0]
                numeric_codes = parts[1:-2]
                type = parts[-2]
                interaction_type = parts[-1]
                for code in numeric_codes:
                    atom_index = int(code)
                    complex_id = self.complex.select_atoms(f"id {atom_index}")
                    for atom in complex_id:
                        atom_name = atom.name
                    for lig_atom in self.lig_noh:
                        if lig_atom.name == atom_name:
                            lig_real_index = lig_atom.id
                    if type == "Donor":
                        highlighted_hbond_donor.append(lig_real_index - 1)
                    elif type == "Acceptor":
                        highlighted_hbond_acceptor.append(lig_real_index - 1)

            elif interaction_type == "hydrophobic":
                for code in numeric_codes:
                    atom_index = int(code)
                    complex_id = self.complex.select_atoms(f"id {atom_index}")
                    for atom in complex_id:
                        atom_name = atom.name
                    for lig_atom in self.lig_noh:
                        if lig_atom.name == atom_name:
                            lig_real_index = lig_atom.id
                    highlighted_hydrophobic.append(lig_real_index - 1)

            elif interaction_type == "waterbridge":
                for code in numeric_codes:
                    atom_index = int(code)
                    complex_id = self.complex.select_atoms(f"id {atom_index}")
                    for atom in complex_id:
                        atom_name = atom.name
                    for lig_atom in self.lig_noh:
                        if lig_atom.name == atom_name:
                            lig_real_index = lig_atom.id
                    highlighted_waterbridge.append(lig_real_index - 1)

            elif interaction_type == "pistacking":
                split_codes = numeric_codes[0].split(",")
                for code in split_codes:
                    atom_index = int(code)
                    complex_id = self.complex.select_atoms(f"id {atom_index}")
                    for atom in complex_id:
                        atom_name = atom.name
                    for lig_atom in self.lig_noh:
                        if lig_atom.name == atom_name:
                            lig_real_index = lig_atom.id
                    highlighted_pistacking.append(lig_real_index - 1)

            elif interaction_type == "halogen":
                numeric_codes = parts[1:-2]
                for code in numeric_codes:
                    atom_index = int(code)
                    complex_id = self.complex.select_atoms(f"id {atom_index}")
                    for atom in complex_id:
                        atom_name = atom.name
                    for lig_atom in self.lig_noh:
                        if lig_atom.name == atom_name:
                            lig_real_index = lig_atom.id
                    highlighted_halogen.append(lig_real_index - 1)

            elif interaction_type == "saltbridge":
                numeric_codes = parts[1:-3]
                saltbridge_type = parts[-2]
                if saltbridge_type == "NI":
                    split_codes = numeric_codes[0].split(",")
                    for code in split_codes:
                        atom_index = int(code)
                        complex_id = self.complex.select_atoms(f"id {atom_index}")
                        for atom in complex_id:
                            atom_name = atom.name
                        for lig_atom in self.lig_noh:
                            if lig_atom.name == atom_name:
                                lig_real_index = lig_atom.id
                        highlighted_ni.append(lig_real_index - 1)
                elif saltbridge_type == "PI":
                    for code in numeric_codes:
                        atom_index = int(code)
                        complex_id = self.complex.select_atoms(f"id {atom_index}")
                        for atom in complex_id:
                            atom_name = atom.name
                        for lig_atom in self.lig_noh:
                            if lig_atom.name == atom_name:
                                lig_real_index = lig_atom.id
                        highlighted_pi.append(lig_real_index - 1)

            elif interaction_type == "pication":
                numeric_codes = parts[1:-2]
                for code in numeric_codes:
                    atom_index = int(code)
                    complex_id = self.complex.select_atoms(f"id {atom_index}")
                    for atom in complex_id:
                        atom_name = atom.name
                    for lig_atom in self.lig_noh:
                        if lig_atom.name == atom_name:
                            lig_real_index = lig_atom.id
                    highlighted_pication.append(lig_real_index - 1)

            elif interaction_type == "metal":
                ligidx = parts[1]
                atom_index = int(ligidx)
                complex_id = self.complex.select_atoms(f"id {atom_index}")
                for atom in complex_id:
                    atom_name = atom.name
                for lig_atom in self.lig_noh:
                    if lig_atom.name == atom_name:
                        lig_real_index = lig_atom.id
                highlighted_metal.append(lig_real_index - 1)

        for value in highlighted_hbond_donor[:]:
            if value in highlighted_hbond_acceptor:
                highlighted_hbond_donor.remove(value)
                highlighted_hbond_acceptor.remove(value)
                highlighted_hbond_both.append(value)

        return (
            highlighted_hbond_donor,
            highlighted_hbond_acceptor,
            highlighted_hbond_both,
            highlighted_hydrophobic,
            highlighted_waterbridge,
            highlighted_pistacking,
            highlighted_halogen,
            highlighted_ni,
            highlighted_pi,
            highlighted_pication,
            highlighted_metal,
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

        interaction_dict = {
            int(key): interaction_dict[interaction_type] for key in keys
        }
        
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
            drawer.drawOptions().addStereoAnnotation = (
                True  # Add stereo information if available
            )
            drawer.DrawMolecule(prepared_ligand)

            # Adjust font size in the SVG output using the FontSize method
            font_size = drawer.FontSize()
            drawer.SetFontSize(
                font_size * 0.5
            )  # You can adjust the multiplier as needed

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

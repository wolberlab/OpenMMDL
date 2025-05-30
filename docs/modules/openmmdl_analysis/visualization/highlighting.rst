API Documentation for highlighting
============================


.. py:class:: FigureHighlighter(complex_pdb_file, ligand_no_h_pdb_file)

   dentifies and highlights ligand atoms involved in various types of protein-ligand interactions.

    :param str complex_pdb_file: Path to the protein-ligand complex PDB file.
    :param str ligand_no_h_pdb_file: Path to the ligand PDB file without hydrogens.
    :ivar mda.Universe complex: MDAnalysis Universe object of the protein-ligand complex.
    :ivar mda.Universe ligand_no_h: MDAnalysis Universe object of the ligand without hydrogens.
    :ivar mda.AtomGroup lig_noh: AtomGroup of all atoms in the ligand without hydrogens.

    .. py:method:: split_interaction_data(data)

        Splits the input data into multiple parts.

        :param list data: A list containing strings with protein residue name, interacting indices and interaction type.
        :return: List of separate formatted strings with separated protein residue name, interacting indices and interaction type.
        :rtype: list

    .. py:method:: highlight_numbers(split_data, starting_idx)

        Extracts the data from the split_data output of the interactions and categorizes it to the correct interaction list.

        :param list split_data: Split interaction data strings with the protein residue name, interacting indices and interaction type.
        :param list starting_idx: Starting indices of ligand atoms used to correctly identify atoms.
        :return: A tuple of lists of highlighted atom indices by interaction type.
        :rtype: tuple

    .. py:method:: generate_interaction_dict(interaction_type, keys)

        Generates a dictionary of interaction RGB color model based on the provided interaction type.

        :param str interaction_type: The type of interaction (e.g., 'hydrophobic', 'hbond_donor').
        :param list keys: Atom indices corresponding to the given interaction type.
        :return: Dictionary mapping each atom index to an RGB color code tuple.
        :rtype: dict


.. py:class:: LigandImageGenerator(ligand_name, complex_pdb_file, ligand_no_h_pdb_file, output_svg_filename, fig_type='svg')

    Generates 2D images of the ligand structure from a protein-ligand complex with atom indices mapped.

    :param str ligand_name: Name of the ligand (3 letters) in the protein-ligand complex topology.
    :param str complex_pdb_file: Path to the protein-ligand complex PDB file.
    :param str ligand_no_h_pdb_file: Path to the ligand PDB file without hydrogens.
    :param str output_svg_filename: Output filename for the generated SVG image.
    :param str fig_type: Type of image to generate. Can be "svg" or "png".

    .. py:method:: generate_image()

        Generates and saves a 2D depiction of the ligand with atom indices labeled. Converts to PNG if specified.

        :return: None. This function writes out a figure and does not return anything.
        :rtype: None

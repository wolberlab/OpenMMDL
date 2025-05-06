API Documentation for FigureHighlighter and LigandImageGenerator
============================


.. py:class:: FigureHighlighter(complex_pdb_file, ligand_no_h_pdb_file)

    A class for identifying and categorizing ligand atoms involved in different types of protein-ligand interactions.

    :param str complex_pdb_file: Path to the protein-ligand complex PDB file.
    :param str ligand_no_h_pdb_file: Path to the ligand PDB file without hydrogens.

    .. py:method:: split_interaction_data(data)

        Splits interaction strings into components: residue/atom identifiers and interaction type.

        :param list data: A list of raw interaction strings.
        :return: A list of split interaction strings.
        :rtype: list

    .. py:method:: highlight_numbers(split_data, starting_idx)

        Parses interaction data to identify and categorize ligand atom indices involved in interactions.

        :param list split_data: Parsed interaction data.
        :param list starting_idx: Starting indices to help resolve atom positions.
        :return: A tuple of lists of highlighted atom indices by interaction type.
        :rtype: tuple

    .. py:method:: generate_interaction_dict(interaction_type, keys)

        Returns an RGB color mapping for a given interaction type and corresponding atoms.

        :param str interaction_type: Type of interaction (e.g., 'hydrophobic').
        :param list keys: Atom indices involved in the interaction.
        :return: Dictionary mapping atom indices to RGB tuples.
        :rtype: dict


.. py:class:: LigandImageGenerator(ligand_name, complex_pdb_file, ligand_no_h_pdb_file, output_svg_filename, fig_type='svg')

    A class for generating annotated 2D ligand structure images using RDKit, with optional export to PNG.

    :param str ligand_name: Ligand residue name in the complex.
    :param str complex_pdb_file: Path to the protein-ligand complex PDB file.
    :param str ligand_no_h_pdb_file: Path to the ligand PDB file without hydrogens.
    :param str output_svg_filename: Path to save the output SVG image.
    :param str fig_type: Type of image to save ('svg' or 'png'). Default is 'svg'.

    .. py:method:: generate_image()

        Generates and saves a 2D depiction of the ligand with atom indices labeled. Converts to PNG if specified.

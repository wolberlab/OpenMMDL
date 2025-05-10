API Documentation for bindingmodes
============================


.. py:class:: BindingModeProcesser(pdb_md, ligand, peptide, special, ligand_rings, interaction_list, threshold, total_frames)

    A class that processes protein-ligand interaction data for a given molecular dynamics (MD) simulation trajectory.
    The class performs multiple analyses on interaction data, such as filtering interactions based on frequency, 
    generating fingerprints according to the combination of interactions present in frames creating specific binding modes.

    :ivar str pdb_md: The path to the molecular dynamics (MD) PDB file for the protein-ligand complex.
    :ivar str ligand: The ligand information (e.g., name or structure).
    :ivar str or None peptide: The peptide chain ID, if applicable, or None if not considered.
    :ivar bool special: Special handling or unique identifiers related to the ligand or protein.
    :ivar list ligand_rings: List of ligand ring information used for hydrophobic interaction detection.
    :ivar pd.DataFrame interaction_list: DataFrame containing the interactions for the trajectory frames.
    :ivar float threshold: The threshold for filtering interactions based on their occurrence frequency.
    :ivar int total_frames: The total number of frames in the molecular dynamics simulation trajectory.
    :ivar dict unique_columns_rings_grouped: A dictionary containing the grouped interactions based on frames.
    :ivar pd.DataFrame interactions_all: DataFrame containing all interaction data processed.
    :ivar dict unique_data_all: Dictionary of unique data generated for all interactions across all frames.
    :ivar dict unique_data: Dictionary of unique data generated based on filtered interactions and threshold.

    .. py:method:: _gather_interactions(df)

        Process a DataFrame with the protein-ligand interaction and generate column names for each unique interaction.

        :param pd.DataFrame df: DataFrame that contains the interaction data for the whole trajectory.
        :returns: A dictionary with the keys being 'FRAME' numbers and values being dictionaries containing 
        row indices and their corresponding unique column names for interactions.
        :rtype: dict

    .. py:method:: _process_interaction_wraper(interaction_list, threshold)

        Apply filtering and interaction enumeration to an interaction DataFrame.

        :param pd.DataFrame interaction_list: Interaction data obtained from the MD simulation.
        :param float threshold: Threshold for interaction occurrence (as a fraction of total frames).
        :returns: Tuple of the processed interaction list and a dictionary of unique identifiers.
        :rtype: Tuple[pd.DataFrame, dict]

    .. py:method:: filtering_values(threshold, df)

        Apply filtering and interaction enumeration to an interaction DataFrame.

        :param float threshold: Minimum fraction of total frames an interaction must appear in.
        :param pandas.DataFrame df: The interaction list DataFrame.
        :returns: 
            - **interaction_list** (*pd.DataFrame*): Modified DataFrame including new interaction columns that contain the filtered values.
            - **unique_data** (*dict*): Dictionary containing unique filtered interaction names.
        :rtype: tuple(pd.DataFrame, dict)

    .. py:method:: unique_data_generation(filtered_values)

        Generates a dictionary of unique interactions from a filtered list of interactions.

        :param list filtered_values: List of filtered interaction identifiers.
        :returns: Dictionary mapping interaction identifiers to themselves.
        :rtype: dict

    .. py:method:: df_iteration_numbering(df, unique_data)

        Iterates over a DataFrame of interaction data and assigns binary indicators (1 or 0) to each row depending on whether the interaction matches any entry in `unique_data`.

        :param pandas.DataFrame df: DataFrame containing interaction data for all frames.
        :param dict unique_data: Dictionary of unique interaction identifiers obtained from :py:meth:`unique_data_generation`.
        :returns: Modifies the input DataFrame in-place by appending columns corresponding to recurring interactions.
        :rtype: None

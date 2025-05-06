API Documentation for BindingModeProcesser
============================


.. py:class:: BindingModeProcesser(pdb_md, ligand, peptide, special, ligand_rings, interaction_list, threshold, total_frames)

    Processes protein-ligand interaction data to extract and analyze recurring binding modes.

    :param str pdb_md: Path to the MD simulation PDB file.
    :param str ligand: Identifier for the ligand.
    :param str peptide: Identifier for the peptide chain (can be None).
    :param bool special: Special flag used in downstream logic (not documented further).
    :param list ligand_rings: List of ligand rings for identifying ring-based hydrophobic interactions.
    :param pandas.DataFrame interaction_list: DataFrame containing protein-ligand interaction data.
    :param float threshold: Threshold (percentage) used to filter recurring interactions.
    :param int total_frames: Total number of trajectory frames used in the analysis.

    .. py:method:: gather_interactions(df)

        Processes a DataFrame with protein-ligand interactions and generates unique interaction column names.

        :param pandas.DataFrame df: Interaction data for the full trajectory.
        :returns: Dictionary mapping frame numbers to index–column name mappings.
        :rtype: dict

    .. py:method:: process_interaction_wraper(interaction_list, threshold)

        Wrapper that filters interaction values, generates unique identifiers, and numbers interactions per frame.

        :param pandas.DataFrame interaction_list: The full interaction list DataFrame.
        :param float threshold: Fraction of total frames to define a recurring interaction.
        :returns: Tuple of the processed interaction list and a dictionary of unique identifiers.
        :rtype: Tuple[pandas.DataFrame, dict]

    .. py:method:: filtering_values(threshold, df)

        Filters interaction values based on occurrence frequency and appends them as columns to the input DataFrame.

        :param float threshold: Minimum fraction of total frames an interaction must appear in.
        :param pandas.DataFrame df: The interaction list DataFrame.
        :returns: A list of interaction identifiers that passed the threshold.
        :rtype: list

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

API Documentation for interactions
============================


.. py:class:: InteractionAnalyzer(pdb_md, dataframe, num_processes, lig_name, special_ligand, peptide, md_len)

    A class for analyzing protein-ligand interactions in molecular dynamics (MD) simulations. It processes the trajectory and generates interaction data.

    :param pdb_md: The MDAnalysis object representing the topology and trajectory of the simulation.
    :param dataframe: Path to a CSV file containing pre-existing interaction data. If None, the class will process the trajectory.
    :param num_processes: Number of processes to use for parallel processing of the trajectory.
    :param lig_name: Name of the ligand in the complex.
    :param special_ligand: Name of a special ligand (optional).
    :param peptide: Chain ID of the peptide (optional).
    :param md_len: The length of the MD trajectory.

    .. py:method:: __init__(self, pdb_md, dataframe, num_processes, lig_name, special_ligand, peptide, md_len)

        Initializes the InteractionAnalyzer class.

    .. py:method:: characterize_complex(self, pdb_file, binding_site_id)

        Characterizes the protein-ligand complex and returns their interaction set.

        :param pdb_file: Path to the PDB file representing the complex.
        :param binding_site_id: Identifier of the binding site to be analyzed.
        :returns: Interaction set for the specified binding site.
        :rtype: PLInteraction or None

    .. py:method:: retrieve_plip_interactions(self, pdb_file, lig_name)

        Retrieves the interactions from PLIP for a given ligand in the PDB file.

        :param pdb_file: Path to the PDB file of the complex.
        :param lig_name: Name of the ligand in the complex to be analyzed.
        :returns: Dictionary of binding sites and their interactions.
        :rtype: dict

    .. py:method:: retrieve_plip_interactions_peptide(self, pdb_file)

        Retrieves the interactions from PLIP for a peptide in the PDB file.

        :param pdb_file: Path to the PDB file of the complex.
        :returns: Dictionary of binding sites and their interactions.
        :rtype: dict

    .. py:method:: create_df_from_binding_site(self, selected_site_interactions, interaction_type="hbond")

        Creates a DataFrame from the selected binding site interactions and specified interaction type.

        :param selected_site_interactions: Pre-calculated interactions from PLIP for the selected binding site.
        :param interaction_type: The interaction type of interest (default is "hbond").
        :returns: DataFrame containing the interaction data.
        :rtype: pandas.DataFrame

    .. py:method:: change_lig_to_residue(self, file_path, new_residue_name)

        Changes the ligand residue name in the topology file to mimic an amino acid residue (for special ligands).

        :param file_path: Path to the topology file.
        :param new_residue_name: New residue name to be used in place of the ligand.

    .. py:method:: process_frame(self, frame)

        Processes a single frame of the MD trajectory and retrieves interaction data for the specified frame.

        :param frame: The index of the frame to be processed.
        :returns: DataFrame containing the interaction data for the processed frame.
        :rtype: pandas.DataFrame

    .. py:method:: process_frame_special(self, frame)

        Processes a single frame of the MD trajectory for special ligands and retrieves interaction data for the specified frame.

        :param frame: The index of the frame to be processed.
        :returns: List of DataFrames containing interaction data for the processed frame with special ligands.
        :rtype: list

    .. py:method:: process_frame_wrapper(self, args)

        Wrapper function for processing MD trajectory frames in parallel.

        :param args: Tuple containing frame index, MDAnalysis Universe, ligand name, special ligand, and peptide chain ID.
        :returns: Tuple containing the frame index and the result from :py:meth:`process_frame`.
        :rtype: tuple

    .. py:method:: fill_missing_frames(self, df)

        Fills missing frames in the DataFrame with placeholder values where no interactions occurred.

        :param df: The interaction data DataFrame.
        :returns: The DataFrame with missing frames filled.
        :rtype: pandas.DataFrame

    .. py:method:: process_trajectory(self)

        Processes the entire protein-ligand trajectory in parallel across multiple CPUs and returns the final interaction data.

        :returns: DataFrame containing interaction data for the entire trajectory.
        :rtype: pandas.DataFrame

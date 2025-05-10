API Documentation for interactions
============================


.. py:class:: InteractionAnalyzer(pdb_md, dataframe, num_processes, lig_name, special_ligand, peptide, md_len)

    Analyzes molecular interactions between a protein and a ligand/peptide 
    throughout an MD trajectory using PLIP (Protein-Ligand Interaction Profiler).

    :ivar mda.Universe pdb_md: MDAnalysis Universe object representing the topology and trajectory.
    :ivar str or None dataframe: Path to an existing interaction CSV file. If None, the trajectory will be processed anew.
    :ivar int num_processes: Number of CPU cores to use for parallel frame analysis.
    :ivar str lig_name: Residue name of the ligand in the complex.
    :ivar str special_ligand: Residue name for special ligands like metal ions (optional).
    :ivar str peptide: Chain ID of the peptide ligand (optional).
    :ivar int md_len: Number of frames in the trajectory.
    :ivar pd.DataFrame interaction_list: DataFrame storing the extracted interactions across the trajectory.

    .. py:method:: _retrieve_plip_interactions(self, pdb_file, lig_name)

        Retrieves the interactions from PLIP.

        :param str pdb_file: The path of the PDB file of the complex.
        :param str lig_name: Name of the Ligand in the complex topology that will be analyzed.
        :returns: A dictionary of the binding sites and the interactions.
        :rtype: dict

    .. py:method:: _retrieve_plip_interactions_peptide(self, pdb_file)

        Retrives the interactions from PLIP for a peptide.

        :param str pdb_file: The path of the PDB file of the complex.
        :returns: A dictionary of the binding sites and the interactions.
        :rtype: dict

    .. py:method:: _create_df_from_binding_site(self, selected_site_interactions, interaction_type="hbond")

        Creates a data frame from a binding site and interaction type.

        :param dict selected_site_interactions: Precaluclated interactions from PLIP for the selected site.
        :param str, optional interaction_type: The interaction type of interest (default set to hydrogen bond). Defaults to "hbond".
        :returns: DataFrame with information retrieved from PLIP.
        :rtype: pd.DataFrame

    .. py:method:: _change_lig_to_residue(self, file_path, new_residue_name)

        Reformats the topology file to change the ligand to a residue. This is needed for interactions with special ligands such as metal ions.

        :param str file_path: Filepath of the topology file.
        :param str new_residue_name: New residue name of the ligand now changed to mimic an amino acid residue.
        :returns: Modifies and writes out new topology file.
        :rtype: None

    .. py:method:: _process_frame(self, frame)

        Process a single frame of MD simulation.

        :param int frame: The number of the frame that will be processed.
        :returns: A dataframe conatining the interaction data for the processed frame.
        :rtype: pd.DataFrame

    .. py:method:: _process_frame_special(self, frame)

        Function extension of process_frame to process special ligands.

        :param int frame: Number of the frame that will be processed.
        :returns: List of dataframes containing the interaction data for the processed frame with the special ligand.
        :rtype: list of pd.DataFrame 

    .. py:method:: _process_frame_wrapper(self, args)

        Wrapper for the MD Trajectory procession.

        :param tuple args: Tuple containing (frame_idx: int - number of the frame to be processed).
        :returns: Tuple containing the frame index and the result of from the process_frame function.
        :rtype: tuple

    .. py:method:: :fill_missing_frames(self, df)

        Fills the frames with no interactions in the DataFrame with placeholder values.

        :param pd.DataFrame df: The input DataFrame with frames that have no interactions
        :returns: DataFrame with placeholder values in the frames with no interactions.
        :rtype: pd.DataFrame

    .. py:method:: _process_trajectory(self)

        Process protein-ligand trajectory with multiple CPUs in parallel.

        :returns: A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
        :rtype: pd.DataFrame

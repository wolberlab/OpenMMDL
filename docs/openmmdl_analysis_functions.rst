OpenMMDL Analysis functions
=============================

This page displays all the functions of **OpenMMDL Analysis**.

openmmdl_analysis.barcode_generation
------------------------------

.. py:function:: barcodegeneration(df, interaction)
    
    Generates barcodes for a given interaction.
    
    :param str protein_name: Name of the protein PDB.
    :param pandas.dataframe df: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param str interaction: name of the interaction to generate a barcode for

    :returns: binary array of wit 1 representing the interaction is present in the corresponding frame
    :rtype: numpy.array

.. py:function::  waterids_barcode_generator(df, interaction)

    Generates a barcode containing coresponding water ids for a given interaction.

    :param pandas.dataframe df: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param str interaction: name of the interaction to generate a barcode for

    :returns: returns a list of waterids for the frames where the interaction is present 0 if no interaction present
    :rtype: list


.. py:function:: plot_barcodes(barcodes, save_path)

    Generates picture of barcodes for interactions of a specific type.

    :param list barcodes: list of np arrays containing the barcodes for each interaction
    :param str save_path: name of the file to save the picture to

    :returns: None
    :rtype: None    

.. py:function:: plot_waterbridge_piechart(df_all, waterbridge_barcodes, waterbridge_interactions)

    Generates piecharts for each waterbridge interaction with the water ids of the interacting waters.

    :param pandas.dataframe df_all: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param list waterbridge_barcodes: list of np arrays containing the barcodes for each waterbridge interaction
    :param list waterbridge_interactions: list of strings containing waterbridge interactions

    :returns: None
    :rtype: None    

.. py:function:: plot_bacodes_grouped(interactions, df_all, interaction_type)

    Generates barcode figures and groups them by ligandatom, aswell as total interaction barcode for a giveen lingenatom.

    :param list interactions: list of pandas.indexes that contain the interactions to generate barcodes for
    :param pandas.dataframe df_all: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param str interaction_type: name of the interaction to generate a barcode for
    
    :returns: None 
    :rtype: None

openmmdl_analysis.binding_mode_processing
------------------------------------------

.. py:function:: gather_interactions(df, ligand_rings, peptide=None)

    Process a DataFrame with the protein-ligand interaction and generate column names for each unique interaction.

    :param pandas.dataframe df: DataFrame that contains the interaction data for the whole trajectory.
    :param list ligand_rings: List of the ligand ring information to recognize the atom numbers belonging to rings for hydrophobic interactions.
    :param str peptide: Name of the peptide chain in the protein. If None, the peptide chain is not considered.

    :returns: A dictionary with the keys being 'FRAME' numbers and values being dictionaries containing row indices and their corresponding unique column names for interactions.
    :rtype: dict

.. py:function:: remove_duplicate_values(data)
    
    Remove the duplicate values from sub-dictionaries within the input dictionary.

    :param dict data: The input dictionary containing sub-dictionaries with possible duplicate values.
    
    :returns: A dictionary without duplicate values.
    :rtype: dict

.. py:function:: combine_subdict_values(data)

    Combines the values from the individual sub-dictionaries into a single list.

    :param dict data: Dictionary with values that are sub-dictionaries.
    
    :returns: A dictionary with a single key named 'all' that contains a list of all combined values from all the sub-dictionaries.
    :rtype: dict

.. py:function:: filtering_values(threshold, frames, df, unique_columns_rings_grouped)

    Filter and append values (interactions) to a DataFrame based on occurrence counts.

    :param float threshold: A treshold value that is used for filtering of the values (interactions) based upon the occurence count.
    :param int frames: The number of frames that is used to calculate the treshold.
    :param pandas.dataframe df: DataFrame to which the filtered values (interactions) will be added.
    :param dict unique_columns_rings_grouped: Dictionary containing the grouped and unique values otained from gather_interactions.

    :returns:  A list of values, with unique values and their corresponding occurence counts.
    :rtype: list

.. py:function:: unique_data_generation(filtered_values)

    :param list filtered_values: A list of values, where the unique interactions are extracted from.

    :returns: A dictionary containing the filtered unique interactions.
    :rtype: dict

.. py:function:: df_iteration_numbering(df, unique_data, peptide=None)

    Loop through the DataFrame and assign the values 1 and 0 to the rows, depending if the corresponding interaction from unique data is present.

    :param pandas.dataframe df: DataFrame which has the interaction data for all of the frames.
    :param dict unique_data: Dictionary that contains the unique interactions obtained from unique_data_generation.
    :param str peptide: Name of the peptide chainid in the original topology. Defaults to None. If None, the peptide chain is not considered.

    :returns: None
    :rtype: None

.. py:function:: update_values(df, new, unique_data)

    Update the values in the input DataFrame based upon the frame values and an reference DataFrame.

    :param pandas.dataframe df: Input DataFrame that will be updated.
    :param pandas.dataframe new: The reference DataFrame containing values that are used to update the input DataFrame.
    :param dict unique_data: A dictionary containing keys that represent the specific unique column names that need to be updated in the input DataFrame.

    :returns: None
    :rtype: None

.. py:function:: calculate_representative_frame(traj, bmode_frame_list, lig)

    Calculates the most representative frame for a bindingmode. This is based uppon the averagwe RMSD of a frame to all other frames in the binding mode.

    :param mdanalysis.universe traj: MDAnalysis universe object containing the trajectory.
    :param list bmode_frame_list: List of frames belonging to a binding mode.
    :param str lig: Name of the ligand in the topology.

openmmdl_analysis.find_stable_waters
------------------------------------------

.. py:function:: trace_waters(topology, trajectory, output_directory)

    Trace the water molecules in a trajectory and write all which move below one Angstrom distance. To adjust the distance alter the integer

    :param str topology: Path to the topology file.
    :param str trajectory: Path to the trajectory file.
    :param str output_directory: Path to the output directory.

    :returns: DataFrame containing stable water coordinates.
    :rtype: pandas.DataFrame
    :returns: Total number of frames.
    :rtype: int

.. py:function:: perform_clustering_and_writing(stable_waters, cluster_eps, total_frames, output_directory)

    Perform DBSCAN clustering on the stable water coordinates, and write the clusters and their representatives to PDB files.

    :param pandas.DataFrame stable_waters: DataFrame containing stable water coordinates.
    :param float cluster_eps: DBSCAN clustering epsilon parameter. This is in Angstrom in this case, and defines which Water distances should be within one cluster
    :param int total_frames: Total number of frames.
    :param str output_directory: Path to the output directory.

    :returns: None
    :rtype: None

.. py:function::  write_pdb_clusters_and_representatives(clustered_waters, min_samples, output_sub_directory)

    Writes the clusters and their representatives to PDB files.

    :param pandas.dataframe clustered_waters: DataFrame containing clustered water coordinates.
    :param int min_samples: DBSCAN clustering min_samples parameter.
    :param str output_sub_directory: Path to the output subdirectory.

    :returns: None
    :rtype: None

.. py:function:: stable_waters_pipeline(topology, trajectory, water_eps, output_directory="./stableWaters")

    Function to run the pipeline to extract stable water clusters, and their representatives from a PDB & DCD file

    :param str topology: Path to the topology file.
    :param str trajectory: Path to the trajectory file.
    :param float water_eps: DBSCAN clustering epsilon parameter.
    :param str output_directory: Path to the output directory. Optional, defaults to "./stableWaters"

    :returns: None
    :rtype: None

.. py:function:: filter_and_parse_pdb(protein_pdb)

    This function reads in a PDB and returns the structure with bioparser.

    :param str protein_pdb: Path to the PDB file.

    :returns: Biopython PDB Structure object.
    :rtype: biopython.structure

.. py:function:: find_interacting_residues(structure, representative_waters, distance_threshold)

    This function maps waters (e.g. the representative waters) to interacting residues of a different PDB structure input. Use "filter_and_parse_pdb" to get the input for this function

    :param biopython.structure structure: Biopython PDB Structure object.
    :param pandas.dataframe representative_waters: DataFrame containing representative water coordinates.
    :param float distance_threshold: Threshold distance for identifying interacting residues.

    :returns:  Dictionary mapping cluster numbers to interacting residues.
    :rtype: dict

.. py:function:: read_pdb_as_dataframe(pdb_file)

    Helper function reading a PDB

    :param str pdb_file: Path to the PDB file.

    :returns: DataFrame containing PDB data.
    :rtype: pandas.dataframe

.. py:function:: analyze_protein_and_water_interaction(protein_pdb_file, representative_waters_file, cluster_eps, output_directory="./stableWaters", distance_threshold=5.0,)

    Analyse the interaction of residues to water molecules using a threshold that can be specified when calling the function

    :param str protein_pdb_file: Path to the protein PDB file without waters.
    :param str representative_waters_file: Path to the representative waters PDB file, or any PDB file containing only waters
    :param float cluster_eps: DBSCAN clustering epsilon parameter.
    :param str output_directory: Path to the output directory. Optional, defaults to "./stableWaters"
    :param float distance_threshold: Threshold distance for identifying interacting residues. Optional, defaults to 5.0

    :returns: None
    :rtype: None

openmmdl_analysis.interaction_gathering
------------------------------------------

.. py:function:: characterize_complex(pdb_file, binding_site_id)

    Characterize the protein-ligand complex and return their interaction set

    :param str pdb_file: Path to the PDB file.
    :param str binding_site_id: A string that specifies the identifier of the binding site

    :returns:  A object representing the interactions if. If Binding site is not found returns None
    :rtype: plip.pdb_complex.basic.interaction_sets

.. py:function:: retrieve_plip_interactions(pdb_file, lig_name)

    Retrieves the interactions from PLIP.

    :param str pdb_file: Path to the PDB file.
    :param str lig_name: Name of the ligand in the topology.

    :returns: A dictionary of the binding sites and the interactions.
    :rtype: dict

.. py:function:: retrieve_plip_interactions_peptide(pdb_file, peptide)

    Retrives the interactions from PLIP for a peptide.

    :param str pdb_file: Path to the PDB file.
    :param str peptide: Name of the peptide chainid in the original topology.

    :returns: A dictionary of the binding sites and the interactions.
    :rtype: dict

.. py:function:: create_df_from_binding_site(selected_site_interactions, interaction_type="hbond")

    Creates a data frame from a binding site and interaction type.

    :param dict selected_site_interactions: Precaluclated interactions from PLIP for the selected site
    :param str interaction_type: The interaction type of interest (default set to hydrogen bond). Defaults to "hbond".

    :returns: DataFrame with information retrieved from PLIP.
    :rtype: pandas.DataFrame

.. py:function:: change_lig_to_residue(file_path, old_residue_name, new_residue_name)

    Reformats the topology file to change the ligand to a residue. This is needed for interactions with special ligands such as metal ions.

    :param str file_path: Path to the topology file.
    :param str old_residue_name: Name of the ligand in the topology.
    :param str new_residue_name: New residue name of the ligand now changed to mimic an amino acid residue.

    :returns: None
    :rtype: None

.. py:function:: process_frame(frame, pdb_md, lig_name, special=None, peptide=None):
    
    Process a single frame of MD simulation.

    :param int frame: Number of frame to be processed.
    :param mdanalysis.universe pdb_md: MDAnalysis universe object containing the trajectory.
    :param str lig_name: Name of the ligand in the topology.
    :param str special: Name of the special ligand in the topology. Defaults to None.
    :param str peptide: Name of the peptide chainid in the original topology. Defaults to None.

    :returns: A dataframe conatining the interaction data for the processed frame.
    :rtype: pandas.dataframe

.. py:function:: process_frame_special(frame, pdb_md, lig_name, special=None)

    Function extension of process_frame to process special ligands.

    :param int frame: Number of the frame that will be processed.
    :param mdanalysis.universe pdb_md: MDAnalysis universe object containing the trajectory.
    :param str lig_name: Name of the ligand in the topology.
    :param str special: Name of the special ligand in the topology. Defaults to None.

    :returns: list of dataframes containing the interaction data for the processed frame with the special ligand.
    :rtype: list   

.. py:function:: process_frame_wrapper(args)

    Wrapper for the MD Trajectory procession.

    :param tuple args: Tuple containing (frame_idx: int - number of the frame to be processed,
                                        pdb_md: mda.universe - MDA Universe class representation of the topology and the trajectory of the file that is being processed,
                                        lig_name: str - Name of the ligand in the complex that will be analyzed,
                                        special_ligand: str - Name of the special ligand that will be analysed,
                                        peptide: str - Chainid of the peptide that will be analyzed)

    :returns: Tuple containing the frame index and the result of from the process_frame function.
    :rtype: tuple

.. py:function:: process_trajectory(pdb_md, dataframe, num_processes, lig_name, special_ligand, peptide)

    Process protein-ligand trajectory with multiple CPUs in parallel.

    :param mdanalysis.universe pdb_md: MDAnalysis universe object containing the trajectory.
    :param str dataframe:  Name of a CSV file as str, where the interaction data will be read from if not None.
    :param int num_processes: Number of processes to be used for the parallelization.
    :param str lig_name: Name of the ligand in the topology.
    :param str special_ligand: Name of the special ligand in the topology.
    :param str peptide: Name of the peptide chainid in the original topology.

    :returns: A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
    :rtype: pandas.dataframe

.. py:function:: fill_missing_frames(df, md_len)

    Fills the frames with no interactions in the DataFrame with placeholder values.

    :param pandas.dataframe df: The input DataFrame with frames that have no Interactions
    :param int md_len: The value that indicates the number of frames, thus allowing the function to loop through the DataFrame

    :returns: DataFrame with placeholder values in the frames with no interactions.
    :rtype: pandas.dataframe

openmmdl_analysis.markov_state_figure_generation
-------------------------------------------------

.. py:function:: min_transition_calculation(min_transition)

    Calculates a list based on the minimum transition time provided values and returns it in factors 1, 2, 5, 10.

    :param int min_transition: The minimum tranisiton time input for the generation of the factors.

    :returns: List with the minimum transition time with factors 1, 2, 5, 10.
    :rtype: list

.. py:function:: binding_site_markov_network(total_frames, min_transitions, combined_dict, font_size=12, size_node=200)

    Generate Markov Chain plots based on transition probabilities.

    :param int total_frames: The number of frames in the protein-ligand MD simulation.
    :param list min_transitions: List of transition tresholds in %. A Markov Chain plot will be generated for each of the tresholds.
    :param dict combined_dict: A dictionary with the information of the Binding Modes and their order of appearance during the simulation for all frames.
    :param int font_size: The font size for the node labels. The default value is set to 12.
    :param int size_node: The size of the nodes in the Markov Chain plot. the default value is set to 200.

    :returns: None
    :rtype: None

openmmdl_analysis.pml_writer
-----------------------------

.. py:function:: generate_pharmacophore_centers(df, interactions)
    
    Generates pharmacophore points for interactions that are points such as hydrophobic and ionic interactions

    :param pandas.dataframe df: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param list interactions: list of strings containing the interactions to generate pharmacophore points for

    :returns: Dict of interactions from which pharmacophore is generated as key and list of coordinates as value
    :rtype: dict

.. py:function:: generate_pharmacophore_vectors(df, interactions)

    Generates pharmacophore points for interactions that are vectors such as hydrogen bond donors or acceptors

    :param pandas.dataframe df: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param list interactions: list of strings containing the interactions to generate pharmacophore points for

    :returns:  Dict of interactions from which pharmacophore is generated as key and list of coordinates as value (first coords are ligand side, second are protein side)
    :rtype: dict

.. py:function:: generate_md_pharmacophore_cloudcenters(df, core_compound, output_filename, sysname, id_num=0)

    Generates pharmacophore from all interactions formed in the MD simulation.
    A feature is generated for each interaction at the center of all its ocurrences.

    :param pandas.dataframe df: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param str core_compound: Name of the ligand.
    :param str output_filename: Name of the output file.
    :param str sysname: Name of the system.
    :param int id_num: Number of the system. Defaults to 0.

    :returns: None
    :rtype: None

.. py:function:: generate_bindingmode_pharmacophore(dict_bindingmode, core_compound, sysname, outname, id_num=0)

    Generates pharmacophore from a binding mode and writes it to a .pml file

    :param dict dict_bindingmode: Dictionary containing all interactions of the bindingmode and thei coresponding ligand and protein coordinates.
    :param str core_compound: Name of the ligand.
    :param str sysname: Name of the system.
    :param str outname: Name of the output file.
    :param int id_num: Number of the system. Defaults to 0.

    :returns: None 
    :rtype: None

.. py:function:: generate_pharmacophore_centers_all_points(df, interactions)

    Generates pharmacophore points for all interactions to generate point cloud.

    :param pandas.dataframe df: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param list interactions: list of strings containing the interactions to generate pharmacophore points for.

    :returns: Dict of interactions from which pharmacophore is generated as key and list of coordinates as value
    :rtype: dict

.. py:function:: generate_point_cloud_pml(cloud_dict, sysname, outname)

    Generates pharmacophore point cloud and writes it to a .pml file

    :param dict cloud_dict: Dictionary containing all interactions of the trajectory and their corresponding ligand coordinates.
    :param str sysname: Name of the system.
    :param str outname: Name of the output file.

    :returns: None
    :rtype: None

openmmdl_analysis.preprocessing
--------------------------------

.. py:function:: increase_ring_indices(ring, lig_index)

    Increases the atom indices in a ring of the ligand obtained from the ligand to fit the atom indices present in the protein-ligand complex.

    :param str ring:  A list of atom indices belonging to a ring that need to be modified.
    :param int lig_index: An integer that is the first number of the ligand atom indices obtained from the protein-ligand, which is used to modify the ring indices

    :returns: A new list with modified atom indicies.
    :rtype: list

.. py:function:: convert_ligand_to_smiles(input_sdf, output_smi)
    
    Converts ligand structures from an SDF file to SMILES :) format

    :param str input_sdf: Path to the input SDF file.
    :param str output_smi: Path to the output SMILES file.

    :returns: None  
    :rtype: None

.. py:function:: process_pdb_file(input_pdb_filename)
    
    Process a PDB file to make it compatible with the openmmdl_analysis package.

    :param str input_pdb_filename: Path to the input PDB file.

    :returns: None
    :rtype: None

.. py:function:: extract_and_save_ligand_as_sdf(input_pdb_filename, output_filename, target_resname)
    
    Extract and save the ligand from the receptor ligand complex PDB file into a new PDB file by itself.

    :param str input_pdb_filename: Path to the input PDB file.
    :param str output_filename: Path to the output SDF file.
    :param str target_resname: Name of the ligand in the target PDB file.

    :returns: None
    :rtype: None

.. py:function:: renumber_atoms_in_residues(input_pdb_file, output_pdb_file, lig_name)
    
    Renumer the atoms of the ligand in the topology PDB file.

    :param str input_pdb_file: Path to the input PDB file.
    :param str output_pdb_file: Path to the output PDB file.
    :param str lig_name: Name of the ligand in the input PDB file.

    :returns: None
    :rtype: None

.. py:function:: replace_atom_type(data)

    Replace wrong ligand atom types in the topology PDB file.

    :param str data: Text of the initial PDB file.

    :returns: Edited text of the PDB file.
    :rtype: str

.. py:function:: process_pdb(input_file, output_file)

    Wrapper function to process a PDB file.

    :param str input_file: Path to the input PDB file.
    :param str output_file: Path to the output PDB file.

    :returns: None
    :rtype: None

.. py:function:: move_hydrogens_to_end(structure, target_residue_name)

    Moves hydrogens to the last lines of theresidue in the PDB file.

    :param biopython.structure structure: Biopython PDB Structure object.
    :param str target_residue_name: Name of the target residue in the input PDB file.

    :returns: None  
    :rtype: None

openmmdl_analysis.rdkit_figure_generation
------------------------------------------

.. py:function:: generate_ligand_image(ligand_name, complex_pdb_file, ligand_no_h_pdb_file, smiles_file, output_png_filename)
    
    Generates a PNG image of the ligand.

    :param str ligand_name: Name of the ligand in the topology.
    :param str complex_pdb_file: Path to the PDB file of the protein-ligand complex.
    :param str ligand_no_h_pdb_file: Path to the PDB file of the ligand without hydrogens.
    :param str smiles_file: Path to the SMILES file of the ligand.
    :param str output_png_filename: Path to the output PNG file.

    :returns: None
    :rtype: None

.. py:function:: split_interaction_data(data)
    
    Splits the Input which consists of the ResNr and ResType, Atom indices, interaction type in multiple parts.

    :param list data: A list of ResNr and ResType, Atom indices, interaction type that needs to be split.
    
    :returns: A new list of the interaction data that consists of three parts, being the protein_partner_name that represents the interacting protein residue, numeric codes, that represent the atom indices of the interacting atoms of the ligand and the interaction type.
    :rtype: list

.. py:function:: highlight_numbers(split_data, starting_idx)
    
    Extracts the data from the split_data output of the interactions and categorizes it to its respective list.

    :param list split_data: A list of interaction data items, where each item contains information about protein partner name, numeric codes and interaction type.
    :param list starting_idx: Starting index of the ligand atom indices used for identifying the correct atom to highlight.

    :returns: A tuple that contains list of all of the highlighted atoms of all of the interactions.
        - highlighted_hbond_donor (list of int): Atom indices for hydrogen bond donors.
        - highlighted_hbond_acceptor (list of int): Atom indices for hydrogen bond acceptors.
        - highlighted_hbond_both (list of int): Atom indices for interactions that are both donors and acceptors.
        - highlighted_hydrophobic (list of int): Atom indices for hydrophobic interactions.
        - highlighted_waterbridge (list of int): Atom indices for water-bridge interactions.
        - highlighted_pistacking (list of int): Atom indices for pi-stacking interactions.
        - highlighted_halogen (list of int): Atom indices for halogen interactions.
        - highlighted_ni (list of int): Atom indices for negative ionizable salt bridge interactions.
        - highlighted_pi (list of int): Atom indices for positive ionizable salt bridge interactions.
        - highlighted_pication (list of int): Atom indices for pi-cation interactions.
        - highlighted_metal (list of int): Atom indices for metal interactions.
    :rtype: tuple

.. py:function:: generate_interaction_dict(interaction_type, keys)
    
    Generates a dictionary of interaction RGB color model based on the provided interaction type.

    :param str interaction_type: The type of the interaction, for example 'hydrophobic'.
    :param list keys: List of the highlighted atoms that display an interaction.

    :returns: A dictionary with the interaction types are associated with their respective RGB color codes.
    :rtype: dict

.. py:function:: update_dict(target_dict, *source_dicts)

    Updates the dictionary wth the keys and values from other dictionaries.

    :param dict target_dict: The dictionary that needs to be updated with new keys and values.
    :param dict source_dicts: One or multiple dictionaries that are used to update the target dictionary with new keys and values.

    :returns: None
    :rtype: None

.. py:function:: create_and_merge_images(binding_mode, occurrence_percent, split_data, merged_image_paths)

    Create and merge images to generate a legend for binding modes.

    :param str binding_mode: Name of the binding mode.
    :param float occurrence_percent: Percentage of the binding mode occurrence.
    :param list split_data: Data of the interactions used to generate the legend.
    :param list merged_image_paths: A list with the paths to the rdkit figures.

    :returns: Paths to the merged images. 
    :rtype: list

.. py:function:: arranged_figure_generation(merged_image_paths, output_path)

    Generate an arranged figure by arranging merged images in rows and columns.

    :param list merged_image_paths: Paths of the merged images with the rdkit figure and legend.
    :param dict output_path: The paths where the arranged output should be saved.

    :returns: None
    :rtype: None

openmmdl_analysis.rmsd_calculation
-----------------------------------

.. py:function:: rmsd_for_atomgroups(prot_lig_top_file, prot_lig_traj_file, selection1, selection2=None)

    Calulate the RMSD for selected atom groups, and save the csv file and plot.

    :param str prot_lig_top_file: Path to the topology file.
    :param str prot_lig_traj_file: Path to the trajectory file.
    :param str selection1: Selection string for main atom group, also used during alignment.
    :param list selection2: Selection strings for additional atom groups. Defaults to None.

    :returns: DataFrame containing RMSD of the selected atom groups over time.
    :rtype: pandas.dataframe

.. py:function:: RMSD_dist_frames(prot_lig_top_file, prot_lig_traj_file, lig, nucleic=False)

    Calculate the RMSD between all frames in a matrix.

    :param str prot_lig_top_file: Path to the topology file.
    :param str prot_lig_traj_file: Path to the trajectory file.
    :param str lig: Name of the ligand in the topology.
    :param bool nucleic: Boolean to indicate if the receptor is a nucleic acid. Defaults to False.

    :returns: pairwise_rmsd_prot. Numpy array of RMSD values for pairwise protein structures.
    :rtype: numpy.array
    :returns: pairwise_rmsd_lig. Numpy array of RMSD values for ligand structures.
    :rtype: numpy.array

openmmdl_analysis.visualization_functions
------------------------------------------

.. py:function:: interacting_water_ids(df_all, waterbridge_interactions)

    Generates a list of all water ids that form water bridge interactions.

    :param pandas.dataframe df_all: Dataframe containing all interactions from plip analysis (typicaly df_all)
    :param list waterbridge_interactions: list of strings containing waterbridge interactions

    :returns: list of all unique water ids that form water bridge interactions
    :rtype: list

.. py:function:: save_interacting_waters_trajectory(pdb_file_path, dcd_file_path, interacting_waters, ligname, special, outputpath="./Visualization/",)

    Saves .pdb and .dcd files of the trajectory containing ligand, receptor and all interacting waters.

    :param str pdb_file_path: Path to the original PDB file.
    :param str dcd_file_path: Path to the original DCD file.
    :param list interacting_waters: List of all interacting water ids
    :param str ligname: Name of the ligand in the topology.
    :param str special: Name of the special ligand in the topology.
    :param str outputpath: Path to the output directory. Optional, defaults to "./Visualization/"

    :returns: None
    :rtype: None

.. py:function:: cloud_json_generation(df_all)
    
    Generates dict for visualization of interaction clouds. Later saved as .json file.

    :param pandas.dataframe df_all: Dataframe containing all interactions from plip analysis (typicaly df_all)

    :returns: Dict containing all interaction clouds
    :rtype: dict

.. py:function:: visualization(ligname, receptor_type="protein or nucleic", height="1000px", width="1000px")

    Generates visualization of the trajectory with the interacting waters and interaction clouds.

    :param str ligname: Name of the ligand in the topology.
    :param str receptor_type: Type of the receptor. Defaults to "protein or nucleic".
    :param str height: Height of the visualization. Defaults to "1000px".
    :param str width: Width of the visualization. Defaults to "1000px".

    :returns: Returns an nglview.widget object containing the visualization
    :rtype: nglview.widget

.. py:function:: run_visualization()

    Runs the visualization notebook in the current directory. The visualization notebook is copied from the package directory to the current directory and automaticaly started.

    :returns: None
    :rtype: None
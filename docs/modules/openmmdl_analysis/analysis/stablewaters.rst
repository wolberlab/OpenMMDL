API Documentation for stablewaters
============================

.. py:class:: StableWaters(trajectory, topology, water_eps)

    A pipeline for identifying and analyzing stable water molecules from molecular dynamics (MD) trajectories.

    This class processes MD simulations to trace water molecules that exhibit limited movement over time,
    clusters them using DBSCAN and identifies representative water positions. Additionally, it analyzes
    potential interactions between stable water clusters and protein residues.

    :param str trajectory: Path to the trajectory file.
    :param str topology: Path to the topology file.
    :param float water_eps: Epsilon parameter for DBSCAN clustering, in Angstrom.
    :ivar mda.Universe u: Universe object created from the topology and trajectory files.
    :ivar float water_eps: Epsilon parameter for DBSCAN clustering, in Angstrom.


    .. py:method:: _trace_waters(output_directory)

        Trace the water molecules in a trajectory and write all which move below one Angstrom distance.
        To adjust the distance alter the integer.

        :param str output_directory: Directory where output files will be saved.
        :returns:
            - **stable_waters** (*pd.DataFrame*): DataFrame containing stable water coordinates.
            - **total_frames** (*int*): Total number of frames.
        :rtype: Tuple[pd.DataFrame, int]


    .. py:method:: _perform_clustering_and_writing(stable_waters, cluster_eps, total_frames, output_directory)

        Perform DBSCAN clustering on the stable water coordinates, and write the clusters and their representatives to PDB files.

        :param pd.DataFrame stable_waters: DataFrame containing stable water coordinates.
        :param float cluster_eps: DBSCAN clustering epsilon parameter. This is in Angstrom in this case, and defines which Water distances should be within one cluster.
        :param int total_frames: Total number of frames.
        :param str output_directory: Directory where output files will be saved.
        :returns: None. Writes out the clusters into their representatives PDB files.
        :rtype: None

    .. py:method:: _write_pdb_clusters_and_representatives(clustered_waters, min_samples, output_sub_directory)

        Writes the clusters and their representatives to PDB files.

        :param pd.DataFrame clustered_waters: DataFrame containing clustered water coordinates.
        :param int min_samples: Minimum number of samples for DBSCAN clustering.
        :param str output_sub_directory: Subdirectory where output PDB files will be saved.
        :returns: None. Writes clusters out as PDB files.
        :rtype: None

    .. py:method:: stable_waters_pipeline(output_directory="./stableWaters")

        Function to run the pipeline to extract stable water clusters, and their representatives from a PDB & DCD file.

        :param str, optional output_directory: Directory where output files will be saved. Default is "./stableWaters".
        :returns: None. This function does not return anything and saves the files.
        :rtype: None


    .. py:method:: _find_interacting_residues(structure, representative_waters, distance_threshold)

        This function maps waters (e.g. the representative waters) to interacting residues of a different PDB structure input.
        Use "filter_and_parse_pdb" to get the input for this function.

        :param Bio.PDB.Structure.Structure: Biopython PDB structure object.
        :param pd.DataFrame representative_waters: DataFrame containing representative water coordinates.
        :param float distance_threshold: Threshold distance for identifying interacting residues.
        :returns: Dictionary mapping cluster numbers to interacting residues.
        :rtype: dict


    .. py:method:: analyze_protein_and_water_interaction(protein_pdb_file, representative_waters_file, cluster_eps, output_directory="./stableWaters", distance_threshold=5.0)

        Analyse the interaction of residues to water molecules using a threshold that can be specified when calling the function.

        :param str protein_pdb_file: Path to the protein PDB file without waters.
        :param str representative_waters_file: Path to the representative waters PDB file, or any PDB file containing only waters.
        :param float cluster_eps: DBSCAN clustering epsilon parameter.
        :param str, optional output_directory: Directory where output files will be saved. Default is "./stableWaters".
        :param float, optional distance_threshold: Threshold distance for identifying interacting residues. Default is 5.0 (Angstrom).
        :returns: None. This function does not return anything and saves the data in a Dataframe.
        :rtype: None

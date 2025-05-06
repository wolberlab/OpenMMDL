API Documentation for stablewaters
============================

.. py:class:: StableWaters(trajectory, topology, water_eps)

    Class to analyze stable water molecules in molecular dynamics trajectories and their interactions with proteins.

    :param str trajectory: Path to the trajectory file (e.g. DCD).
    :param str topology: Path to the topology file (e.g. PDB).
    :param float water_eps: DBSCAN clustering epsilon (in Ångström).


    .. py:method:: trace_waters(output_directory)

        Identifies water molecules that move less than 1 Å between consecutive frames and saves their coordinates.

        :param str output_directory: Directory where the results will be saved.
        :returns: A DataFrame with stable water coordinates and the total number of frames.
        :rtype: Tuple[pandas.DataFrame, int]


    .. py:method:: perform_clustering_and_writing(stable_waters, cluster_eps, total_frames, output_directory)

        Performs DBSCAN clustering on stable waters and writes cluster files and representative waters to PDB.

        :param pandas.DataFrame stable_waters: DataFrame containing stable water coordinates.
        :param float cluster_eps: Clustering epsilon distance for DBSCAN.
        :param int total_frames: Total number of trajectory frames.
        :param str output_directory: Output path for clustered files.


    .. py:method:: write_pdb_clusters_and_representatives(clustered_waters, min_samples, output_sub_directory)

        Writes each DBSCAN cluster and its average (representative) water position to separate PDB files.

        :param pandas.DataFrame clustered_waters: Clustered water coordinates.
        :param int min_samples: Minimum samples parameter used in DBSCAN.
        :param str output_sub_directory: Directory to save the PDB files.


    .. py:method:: stable_waters_pipeline(output_directory="./stableWaters")

        Runs the full pipeline: stable water identification, clustering, and representative water export.

        :param str output_directory: Directory for output files. Default is "./stableWaters".


    .. py:method:: analyze_protein_and_water_interaction(protein_pdb_file, representative_waters_file, cluster_eps, output_directory="./stableWaters", distance_threshold=5.0)

        Analyzes interactions between protein residues and representative water molecules based on a distance threshold.

        :param str protein_pdb_file: Path to the protein PDB file (without water).
        :param str representative_waters_file: PDB file containing only water molecules (usually representative waters).
        :param float cluster_eps: Clustering epsilon value.
        :param str output_directory: Base directory for outputs. Default is "./stableWaters".
        :param float distance_threshold: Distance threshold in Å for interaction detection. Default is 5.0 Å.


    .. py:method:: filter_and_parse_pdb(protein_pdb)

        Parses a protein PDB file, filtering out water and heteroatoms.

        :param str protein_pdb: Path to a PDB file.
        :returns: Parsed Biopython structure.
        :rtype: Bio.PDB.Structure.Structure


    .. py:method:: find_interacting_residues(structure, representative_waters, distance_threshold)

        Finds residues within a distance threshold of representative waters.

        :param Bio.PDB.Structure.Structure structure: Biopython structure object.
        :param pandas.DataFrame representative_waters: DataFrame of water coordinates.
        :param float distance_threshold: Distance in Å for considering an interaction.
        :returns: Dictionary mapping cluster number to interacting residue tuples.
        :rtype: dict


    .. py:method:: read_pdb_as_dataframe(pdb_file)

        Reads water atom coordinates from a PDB file and returns a DataFrame.

        :param str pdb_file: Path to a PDB file.
        :returns: DataFrame with coordinates.
        :rtype: pandas.DataFrame

API Documentation for RMSDAnalyzer
============================


.. py:class:: RMSDAnalyzer(prot_lig_top_file, prot_lig_traj_file)

    Class to analyze RMSD (Root Mean Square Deviation) of protein-ligand systems using MDAnalysis.

    :param str prot_lig_top_file: Path to the topology file (e.g., PDB).
    :param str prot_lig_traj_file: Path to the trajectory file (e.g., DCD).

    .. py:method:: rmsd_for_atomgroups(fig_type, selection1, selection2=None)

        Calculates RMSD over time for specified atom groups and saves results as CSV and a plot.

        :param str fig_type: Format to save the figure (e.g., "png", "jpg").
        :param str selection1: Atom selection string for alignment and RMSD calculation.
        :param list selection2: Optional. Additional selection strings for other atom groups.
        :returns: DataFrame containing RMSD values over trajectory frames.
        :rtype: pandas.DataFrame

    .. py:method:: rmsd_dist_frames(fig_type, lig, nucleic=False)

        Computes pairwise RMSD matrix for protein/nucleic acid and ligand over all frames, and plots heatmaps.

        :param str fig_type: Format to save the figure (e.g., "png", "jpg").
        :param str lig: Residue name of the ligand.
        :param bool nucleic: Whether the receptor contains nucleic acids. Default is False.
        :returns: Pairwise RMSD matrices for protein/nucleic acid and ligand.
        :rtype: Tuple[numpy.ndarray, numpy.ndarray]

    .. py:method:: calc_rmsd_2frames(ref, frame)

        Calculates RMSD between a reference and a single trajectory frame using JIT-accelerated computation.

        :param numpy.ndarray ref: Reference atomic coordinates.
        :param numpy.ndarray frame: Atomic coordinates of a frame.
        :returns: RMSD value.
        :rtype: float

    .. py:method:: calculate_distance_matrix(selection)

        Calculates a full pairwise RMSD distance matrix for a selected group of atoms across all trajectory frames.

        :param str selection: Atom selection string.
        :returns: Symmetric RMSD distance matrix.
        :rtype: numpy.ndarray

    .. py:method:: calculate_representative_frame(bmode_frames, DM)

        Determines the most representative frame for a group of frames by average RMSD.

        :param list bmode_frames: List of frame indices for a binding mode.
        :param numpy.ndarray DM: Pairwise RMSD distance matrix.
        :returns: Frame index of the most representative structure.
        :rtype: int


.. py:function:: calc_rmsd_2frames_jit(ref, frame)

    JIT-accelerated function to compute RMSD between two sets of atomic coordinates.

    :param numpy.ndarray ref: Reference atomic coordinates (N x 3 array).
    :param numpy.ndarray frame: Atomic coordinates of the comparison frame (N x 3 array).
    :returns: Root Mean Square Deviation between the two coordinate sets.
    :rtype: float

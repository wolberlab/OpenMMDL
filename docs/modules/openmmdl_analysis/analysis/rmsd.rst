API Documentation for RMSDAnalyzer
==================================

.. py:class:: RMSDAnalyzer(top_file, traj_file)

    A class for Root-Mean-Square Deviation (RMSD) analysis in molecular dynamics simulations.
    Provides methods to compute RMSD over time, between frames, and to determine representative structures.

    :param str top_file: Path to the topology file.
    :param str traj_file: Path to the trajectory file.
    
    :ivar mda.Universe universe: MDAnalysis Universe initialized with the input files.

    .. py:method:: rmsd_for_atomgroups(fig_type, selection1, selection2=None)

        Calculate RMSD over time for selected atom groups and save the results as CSV and a plot.

        :param str fig_type: Format to save the figure (e.g., "png", "jpg").
        :param str selection1: Atom selection string for alignment and RMSD calculation.
        :param list selection2: (Optional) Additional selection strings for other atom groups. Defaults to None.
        :returns: DataFrame containing RMSD values over trajectory frames.
        :rtype: pandas.DataFrame

    .. py:method:: rmsd_dist_frames(fig_type, lig, nucleic=False)

        Computes pairwise RMSD matrix for protein/nucleic acid and ligand over all frames, and plots heatmaps.

        :param str fig_type: Format to save the figure (e.g., "png", "jpg").
        :param str lig: Residue name of the ligand (used in atom selection).
        :param bool nucleic: Whether the receptor contains nucleic acids. Defaults to False.
        :returns: Tuple containing RMSD matrices for the protein/nucleic acid and the ligand.
        :rtype: Tuple[numpy.ndarray, numpy.ndarray]

    .. py:method:: calc_rmsd_2frames(ref, frame)

        Wrapper for JIT-accelerated RMSD calculation between two sets of atomic coordinates.

        :param numpy.ndarray ref: Reference atomic coordinates.
        :param numpy.ndarray frame: Atomic coordinates of the comparison frame.
        :returns: RMSD value.
        :rtype: float

    .. py:method:: calculate_distance_matrix(selection)

        Computes a full pairwise RMSD distance matrix for a selected group of atoms across all trajectory frames.

        :param str selection: Atom selection string (e.g., "protein", "resname LIG").
        :returns: Symmetric RMSD distance matrix.
        :rtype: numpy.ndarray

    .. py:method:: calculate_representative_frame(bmode_frames, DM)

        Identifies the most representative frame from a list of frames based on average RMSD to others.

        :param list bmode_frames: List of frame indices belonging to a binding mode.
        :param numpy.ndarray DM: Pairwise RMSD distance matrix.
        :returns: Index of the most representative frame.
        :rtype: int


.. py:function:: calc_rmsd_2frames_jit(ref, frame)

    JIT-accelerated function to compute RMSD between two sets of atomic coordinates.

    :param numpy.ndarray ref: Reference atomic coordinates (N x 3 array).
    :param numpy.ndarray frame: Atomic coordinates of the comparison frame (N x 3 array).
    :returns: Root Mean Square Deviation between the two coordinate sets.
    :rtype: float

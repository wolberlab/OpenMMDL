API Documentation for rmsd
==================================

.. py:class:: RMSDAnalyzer(top_file, traj_file)

    A class responsible for the Root-Mean-Square Deviation (RMSD) analysis throughout the molecular dynamics simulation. 
    The class provides functionalities to calculate RMSD over time, compute pairwise RMSD between trajectory frames 
    and identify the representative frames within clusters of the binding modes.

    :param str top_file: Path to the topology file.
    :param str traj_file: Path to the trajectory file.
    
    :ivar mda.Universe universe: MDAnalysis Universe object initialized with the provided topology and trajectory files.

    .. py:method:: rmsd_for_atomgroups(fig_type, selection1, selection2=None)

        Calculate the RMSD for selected atom groups, and save the CSV file and plot.

        :param str fig_type: Type of the figure to save (e.g., 'png', 'jpg').
        :param str selection1: Selection string for main atom group, also used during alignment.
        :param list of str, optional selection2: Selection strings for additional atom groups. Defaults to None.
        :returns: DataFrame containing RMSD of the selected atom groups over time.
        :rtype: pd.DataFrame

    .. py:method:: rmsd_dist_frames(fig_type, lig, nucleic=False)

        Calculate the RMSD between all frames in a matrix.

        :param str fig_type: Type of the figure to save (e.g., 'png', 'jpg').
        :param str lig: Ligand name saved in the above PDB file. Selection string for the MDAnalysis AtomGroup to be investigated, also used during alignment.
        :param bool, optional nucleic: Bool indicating if the receptor to be analyzed contains nucleic acids. Defaults to False.
        :returns:
            - **pairwise_rmsd_prot** (*np.ndarray*): Numpy array of RMSD values for pairwise protein structures.
            - **pairwise_rmsd_lig** (*np.ndarray*): Numpy array of RMSD values for ligand structures.
        :rtype: Tuple[numpy.ndarray, numpy.ndarray]

    .. py:method:: _calc_rmsd_2frames(ref, frame)

        Calculates the RMSD between a reference and a target frame.
    
        This method serves as a wrapper for the `calc_rmsd_2frames_jit` function, 
        which dpes the actual RMSD calculation between two sets of coordinates.

        :param np.ndarray ref: Numpy array representing the reference atom positions, shape (N, 3).
        :param np.ndarray frame: Numpy array representing the atom positions of the target frame, shape (N, 3).
        :returns: The RMSD value between the reference and target frame.
        :rtype: float

    .. py:method:: calculate_distance_matrix(selection)

        Calculates the pairwise RMSD-based distance matrix for all trajectory frames 
        for the selected atom selection.

        :param str selection: Selection string for the atoms (e.g., 'protein', 'resname LIG') used to compute the RMSD between frames.
        :returns: Numpy array containing RMSD values between all pairs of frames.
        :rtype: np.ndarray

    .. py:method:: calculate_representative_frame(bmode_frames, DM)

        Calculates the most representative frame for a bindingmode. 
        This is based uppon the averagwe RMSD of a frame to all other frames in the binding mode.

        :param list of int bmode_frames: List of frames belonging to a binding mode.
        :param np.ndarray DM: Distance matrix of trajectory.
        :returns: Number of the most representative frame.
        :rtype: int


.. py:function:: calc_rmsd_2frames_jit(ref, frame)

    Calculates the RMSD between two frames of atomic coordinates.

    :param np.ndarray ref: Numpy array containing the reference atomic coordinates.
    :param np.ndarray frame: Numpy array containing the atomic coordinates of the target frame.
    :returns: The RMSD value between the reference and target frame.
    :rtype: float

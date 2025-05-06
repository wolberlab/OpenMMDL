API Documentation for TrajectorySaver
============================


.. py:class:: TrajectorySaver(pdb_md, ligname, special, nucleic)

    Class to handle saving of trajectory frames and interacting water molecules from MDAnalysis simulations.

    :param mda.Universe pdb_md: MDAnalysis Universe object containing the trajectory.
    :param str ligname: Residue name of the ligand in the PDB file.
    :param str special: Residue name of a special molecule (e.g., HEM).
    :param bool nucleic: Indicates if the receptor contains nucleic acids.

    .. py:method:: save_interacting_waters_trajectory(interacting_waters, outputpath)

        Saves a trajectory (.pdb and .dcd) including protein/nucleic acids, ligand, special residues, and selected water molecules.

        :param list interacting_waters: List of water residue IDs to include.
        :param str outputpath: Path to directory where output files will be saved (default: "./Visualization/").

    .. py:method:: save_frame(frame, outpath, selection=False)

        Saves a single frame of the trajectory to a file.

        :param int frame: Index of the frame to save.
        :param str outpath: Output path for the saved frame.
        :param str selection: Optional atom selection string. If not provided, all atoms are included.

API Documentation for TrajectorySaver
============================


.. py:class:: TrajectorySaver(pdb_md, ligname, special, nucleic)

    A class to handle the saving of molecular dynamics trajectory data for ligand–receptor complexes,
    with or without interacting water molecules.

    :param mda.Universe pdb_md: An MDAnalysis Universe object containing the trajectory and topology.
    :param str ligname: The residue name (3 letters) of the ligand in the PDB file.
    :param str special: The residue name (3 letters) of any special ligand or cofactor (e.g., HEM) in the PDB file.
    :param bool nucleic: Flag indicating if the receptor is nucleic acid (DNA/RNA). Used for proper atom selection.

    .. py:method:: save_interacting_waters_trajectory(interacting_waters, outputpath)

        Saves .pdb and .dcd files of the trajectory containing ligand, receptor and all interacting waters.

        :param list interacting_waters: List of all interacting water IDs.
        :param str outputpath: Path to directory where output files will be saved (default: "./Visualization/").
        :returns: None. This function writes output directly to a new PDB or DCD file and does not return anything.
        :rtype: None

    .. py:method:: save_frame(frame, outpath, selection=False)

        Saves a single frame of the trajectory to a file.

        :param int frame: Index of the frame to save.
        :param str outpath: Output path for the saved frame.
        :param str selection: Optional atom selection string. If not provided, all atoms are included.

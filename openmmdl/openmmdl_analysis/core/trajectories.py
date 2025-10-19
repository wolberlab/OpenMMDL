import MDAnalysis as mda


class TrajectorySaver:
    """
    A class to handle the saving of molecular dynamics trajectory data for ligand–receptor complexes, with or without interacting water molecules.

    Attributes
    ----------
    pdb_md : mda.Universe
        An MDAnalysis Universe object containing the trajectory and topology.
    ligname : str
        The residue name (3 letters) of the ligand in the PDB file.
    special : str
        The residue name (3 letters) of any special ligand or cofactor (e.g., HEM) in the PDB file.
    nucleic : bool
        Flag indicating if the receptor is nucleic acid (DNA/RNA). Used for proper atom selection.
    """

    def __init__(self, pdb_md, ligname, special, nucleic):
        self.pdb_md = pdb_md
        self.ligname = ligname
        self.special = special
        self.nucleic = nucleic

    def save_interacting_waters_trajectory(self, interacting_waters, outputpath):
        """
        Saves .pdb and .dcd files of the trajectory containing ligand, receptor and all interacting waters.

        Parameters
        ----------
        interacting_waters : list of int
            List of all interacting water IDs.
        outputpath : str
            Filepath to output new PDB and DCD files.

        Returns
        -------
        None
            This function writes output directly to a new PDB or DCD file and does not return anything.
        """
        water_atoms = self.pdb_md.select_atoms(
            f"protein or nucleic or resname {self.ligname} or resname {self.special}"
        )

        for water in interacting_waters:
            add_water_atoms = self.pdb_md.select_atoms(f"resname HOH and resid {water}")
            water_atoms = water_atoms + add_water_atoms

        water_atoms.write(f"{outputpath}interacting_waters.pdb")

        with mda.Writer(f"{outputpath}interacting_waters.dcd", water_atoms.n_atoms) as W:
            for ts in self.pdb_md.trajectory:
                W.write(water_atoms)

    def save_frame(self, frame, outpath, selection=False):
        """
        Saves a single frame of the trajectory.

        Parameters
        ----------
        frame : int
            Number of the frame to save.
        outpath : str
            Path to save the frame to.
        selection : str or bool, optional
            A MDAnalysis selection string to specify which atoms to save. Defaults to False.

        Returns
        -------
        None
            This function writes output directly to a new PDB file and does not return anything.
        """
        self.pdb_md.trajectory[frame]
        if selection:
            frame_atomgroup = self.pdb_md.atoms[selection]
        else:
            frame_atomgroup = self.pdb_md.atoms
        frame_atomgroup.write(outpath)

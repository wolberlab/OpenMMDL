import MDAnalysis as mda


class TrajectorySaver:
    def __init__(self, pdb_md, ligname, special, nucleic):
        """Initializes the TrajectorySaver with a mda.Universe object, ligand name, special residue name and if receptor is nucleic.

        Args:
            pdb_md (mda.Universe): MDAnalysis Universe object containing the trajectory.
            ligname (str): name of the ligand in the pdb file.
            special (str): name of the special residue/ligand in the pdb file (e.g. HEM).
            nucleic (bool): True if receptor is nucleic, False otherwise.
        """
        self.pdb_md = pdb_md
        self.ligname = ligname
        self.special = special
        self.nucleic = nucleic

    def save_interacting_waters_trajectory(self, interacting_waters, outputpath):
        """Saves .pdb and .dcd files of the trajectory containing ligand, receptor and all interacting waters.

        Args:
            interacting_waters (list): list of all interacting water ids.
            outputpath (str, optional): filepath to output new pdb and dcd files. Defaults to './Visualization/'.
        """
        water_atoms = self.pdb_md.select_atoms(
            f"protein or nucleic or resname {self.ligname} or resname {self.special}"
        )

        for water in interacting_waters:
            add_water_atoms = self.pdb_md.select_atoms(f"resname HOH and resid {water}")
            water_atoms = water_atoms + add_water_atoms

        water_atoms.write(f"{outputpath}interacting_waters.pdb")

        with mda.Writer(
            f"{outputpath}interacting_waters.dcd", water_atoms.n_atoms
        ) as W:
            for ts in self.pdb_md.trajectory:
                W.write(water_atoms)

    def save_frame(self, frame, outpath, selection=False):
        """Saves a single frame of the trajectory.

        Args:
            frame (int): Number of the frame to save.
            outpath (str): Path to save the frame to.
            selection (str, optional): MDAnalysis selection string. Defaults to False.
        """
        self.pdb_md.trajectory[frame]
        if selection:
            frame_atomgroup = self.pdb_md.atoms[selection]
        else:
            frame_atomgroup = self.pdb_md.atoms
        frame_atomgroup.write(outpath)

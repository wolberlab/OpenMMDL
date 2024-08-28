import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.distances import dist


class MdtrajConverter:
    """
    A class to handle the conversion of molecular dynamics trajectories into centered trajectories in a desired format.
    
    The MdtrajConverter class is designed to center molecules and apply periodic boundary conditions for each frame
    of the trajectory, then save the output in the desired format(s).

    Attributes:
        pdb_file (str): The name of the PDB file used for the topology.
        mdtraj_output (str): The selected format(s) for the output of the topology and trajectory.
    """

    def __init__(self, pdb_file, mdtraj_output):
        """
        Initialize the MdtrajConverter with the provided PDB file and output format.

        Args:
            pdb_file (str): Name of the PDB file.
            mdtraj_output (str): The selected format that will be used as an output of the topology and trajectory.
        """
        self.pdb_file = pdb_file
        self.mdtraj_output = mdtraj_output

    def mdtraj_conversion(self):
        """Recenter and apply periodic boundary conditions to the molecules in each frame of the trajectory."""
        mdtraj_frames = md.load_dcd("trajectory.dcd", top=self.pdb_file)
        mdtraj_frames.image_molecules()

        if "dcd" in self.mdtraj_output:
            mdtraj_frames.save_dcd("centered_old_coordinates.dcd")
        if "xtc" in self.mdtraj_output:
            mdtraj_frames.save_xtc("centered_old_coordinates.xtc")

        mdtraj_first_frame = mdtraj_frames[0:1]

        if "pdb" in self.mdtraj_output:
            mdtraj_first_frame.save_pdb("centered_old_coordinates_top.pdb")
        if "gro" in self.mdtraj_output:
            mdtraj_first_frame.save_gro("centered_old_coordinates_top.gro")


class MDAnalysisConversion:
    """
    A class to handle the conversion of molecular dynamics trajectories using MDAnalysis.
    
    The MDAnalysisConversion class is designed to translate the trajectory so that all frames coincide with the center of geometry,
    and output the results in the specified formats.

    Attributes:
        post_mdtraj_pdb_file (str): Name of the post-MDtraj PDB file.
        post_mdtraj_dcd_file (str): Name of the post-MDtraj DCD file.
        mda_output (str): Selection of output formats.
        output_selection (str): Selection of topologies with specific atom selections that will be created.
        ligand_name (str, optional): Ligand name saved in the PDB file.
        special_ligname (str, optional): Special residue name saved in the PDB file.
    """
    def __init__(
        self,
        post_mdtraj_pdb_file,
        post_mdtraj_dcd_file,
        mda_output,
        output_selection,
        ligand_name=None,
        special_ligname=None,
    ):
        """
        Initialize the MDAnalysisConversion class with given parameters.

        Args:
            post_mdtraj_pdb_file (str): Name of the post-MDtraj PDB file.
            post_mdtraj_dcd_file (dcd file): Name of the post-MDtraj DCD File.
            ligand_name (str): Ligand name saved in the PDB file.
            special_ligname (str): Special residue name saved in the PDB file.
            mda_output (str): Selection of output formats.
            output_selection (str): Selection of topologies with specific atom selections that will be created.
        """
        self.post_mdtraj_pdb_file = post_mdtraj_pdb_file
        self.post_mdtraj_dcd_file = post_mdtraj_dcd_file
        self.mda_output = mda_output
        self.output_selection = output_selection
        self.ligand_name = ligand_name
        self.special_ligname = special_ligname

    def mdanalysis_conversion(self):
        """
        Translate the trajectory so that all frames coincide with its center of geometry and output the results in specified formats.
        """
        topology_trajectory = mda.Universe(
            self.post_mdtraj_pdb_file, self.post_mdtraj_dcd_file
        )
        topology_trajectory_all_atoms = topology_trajectory.select_atoms("all")
        # Translate the trajectory so that all frames coincide with its center of geometry
        topology_trajectory_all_atoms.atoms.translate(
            topology_trajectory_all_atoms.center_of_mass()
        )
        topology_trajectory_protein_ligand = topology_trajectory.select_atoms(
            f"protein or resname {self.ligand_name} or resname {self.special_ligname}"
        )

        if "pdb" in self.mda_output:
            if self.output_selection != "mda_prot_lig":
                with mda.Writer(
                    f"centered_traj_unaligned.dcd",
                    topology_trajectory_all_atoms.n_atoms,
                ) as w:
                    for ts in topology_trajectory.trajectory[1:]:
                        w.write(topology_trajectory_all_atoms)
                topology_trajectory_all_atoms.write(f"centered_top.pdb")
                topology_ref_all_pdb = mda.Universe(f"centered_top.pdb")
                mobile_all_pdb = mda.Universe(
                    f"centered_top.pdb", f"centered_traj_unaligned.dcd"
                )
                alignment_all_pdb = align.AlignTraj(
                    mobile_all_pdb,
                    topology_ref_all_pdb,
                    select="protein and name CA",
                    weights="mass",
                    filename=f"centered_traj.dcd",
                )
                alignment_all_pdb.run()

            if self.output_selection != "mda_all":
                with mda.Writer(
                    f"prot_lig_traj_unaligned.dcd",
                    topology_trajectory_protein_ligand.n_atoms,
                ) as w:
                    for ts in topology_trajectory.trajectory[1:]:
                        w.write(topology_trajectory_protein_ligand)
                topology_trajectory_protein_ligand.write(f"prot_lig_top.pdb")
                topology_ref_prot_lig_pdb = mda.Universe(f"prot_lig_top.pdb")
                mobile_prot_lig_pdb = mda.Universe(
                    f"prot_lig_top.pdb", f"prot_lig_traj_unaligned.dcd"
                )
                alignment_prot_lig_pdb = align.AlignTraj(
                    mobile_prot_lig_pdb,
                    topology_ref_prot_lig_pdb,
                    select="protein and name CA",
                    weights="mass",
                    filename=f"prot_lig_traj.dcd",
                )
                alignment_prot_lig_pdb.run()

        if "gro" in self.mda_output:
            if self.output_selection != "mda_prot_lig":
                with mda.Writer(
                    f"centered_traj_unaligned.xtc",
                    topology_trajectory_all_atoms.n_atoms,
                ) as w:
                    for ts in topology_trajectory.trajectory[1:]:
                        w.write(topology_trajectory_all_atoms)
                topology_trajectory_all_atoms.write(f"centered_top.gro")
                topology_ref_all_gro = mda.Universe(f"centered_top.gro")
                mobile_all_gro = mda.Universe(
                    f"centered_top.gro", f"centered_traj_unaligned.xtc"
                )
                alignment_all_gro = align.AlignTraj(
                    mobile_all_gro,
                    topology_ref_all_gro,
                    select="protein and name CA",
                    weights="mass",
                    filename=f"centered_traj.xtc",
                )
                alignment_all_gro.run()

            if self.output_selection != "mda_all":
                with mda.Writer(
                    f"prot_lig_traj_unaligned.xtc",
                    topology_trajectory_protein_ligand.n_atoms,
                ) as w:
                    for ts in topology_trajectory.trajectory[1:]:
                        w.write(topology_trajectory_protein_ligand)
                topology_trajectory_protein_ligand.write(f"prot_lig_top.gro")
                topology_ref_prot_lig_gro = mda.Universe(f"prot_lig_top.gro")
                mobile_prot_lig_gro = mda.Universe(
                    f"prot_lig_top.gro", f"prot_lig_traj_unaligned.xtc"
                )
                alignment_prot_lig_gro = align.AlignTraj(
                    mobile_prot_lig_gro,
                    topology_ref_prot_lig_gro,
                    select="protein and name CA",
                    weights="mass",
                    filename=f"prot_lig_traj.xtc",
                )
                alignment_prot_lig_gro.run()


class MDPostProcessingHandler:
    def __init__(self, config_parser):
        """
        Initialize the MDPostProcessingHandler with the necessary parameters.

        Parameters:
        - config_parser (ConfigParser): The configuration parser instance.
        - session (dict): The session dictionary containing post-processing configuration.
        - protein (str): The protein identifier.
        - nmLigName (str, optional): Name of the normal ligand (optional).
        - spLigName (str, optional): Name of the special ligand (optional).
        """
        self.config_parser = config_parser
        self.postprocessing = config_parser.postprocessing
        self.input_filetype = config_parser.input_file_type
        self.protein = config_parser.protein
        self.prmtop_file = config_parser.prmtop_file
        self.old_output = config_parser.old_output
        self.ligand = config_parser.ligand
        self.ligand_name = config_parser.nmLig
        self.special_ligname = config_parser.spLigName
        self.mda_output = config_parser.mda_output
        self.mda_selection = config_parser.mda_selection

    def process(self):
        """
        Perform the MD post-processing based on the configuration.
        """
        if self.postprocessing == "True":
            print("yy")
            self._handle_mdtraj_conversion()
            self._handle_mdanalysis_conversion()

    def _handle_mdtraj_conversion(self):
        """
        Handle the MDTraj conversion process based on the file type.
        """
        if self.input_filetype == "pdb":
            print("do it")
            converter = MdtrajConverter(f'Equilibration_{self.protein}', self.old_output)
            converter.mdtraj_conversion()
        elif self.input_file_type == "amber":
            converter = MdtrajConverter(self.prmtop_file, self.old_output)
            converter.mdtraj_conversion()

    def _handle_mdanalysis_conversion(self):
        """
        Handle the MDAnalysis conversion process based on the file type and session parameters.
        """
        if self.input_filetype == "pdb":
            print("dup dup do it")
            self._perform_mdanalysis_conversion(
                'centered_old_coordinates_top.pdb',
                'centered_old_coordinates.dcd',
                ligand_name='UNK' if self.ligand else None
            )
        elif self.input_filetype == "amber":
            ligand_name = self.ligand_name
            special_ligname = self.special_ligname
            self._perform_mdanalysis_conversion(
                'centered_old_coordinates_top.pdb',
                'centered_old_coordinates.dcd',
                ligand_name=ligand_name,
                special_ligname=special_ligname
            )

    def _perform_mdanalysis_conversion(self, pdb_file, dcd_file, ligand_name=None, special_ligname=None):
        """
        Perform the MDAnalysis conversion.

        Args:
        - pdb_file (str): The PDB file for the post-MDtraj processing.
        - dcd_file (str): The DCD file for the post-MDtraj processing.
        - ligand_name (str, optional): The ligand name in the PDB file.
        - special_ligname (str, optional): The special residue name in the PDB file.
        """
        converter = MDAnalysisConversion(
            pdb_file,
            dcd_file,
            self.mda_output,
            self.mda_selection,
            ligand_name=ligand_name,
            special_ligname=special_ligname,
        )
        converter.mdanalysis_conversion()

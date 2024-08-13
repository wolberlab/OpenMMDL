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

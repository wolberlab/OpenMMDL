import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import align


def mdtraj_conversion(pdb_file, mdtraj_output):
    """Recenter and apply periodic boundary conditions to the molecules in each frame of the trajectory, and save the centered trajectory and its first frame.

    Args:
        pdb_file (str): Name of the PDB file. This PDB file stores the extracted frames from the MD trajectory.
        mdtraj_output (str): The selected format that will be used as an output of the topology and trajectory.

    Returns:
        None
    """
    mdtraj_frames = md.load_dcd("trajectory.dcd", top=pdb_file)
    mdtraj_frames.image_molecules()
    if "dcd" in mdtraj_output:
        mdtraj_frames.save_dcd("centered_old_coordinates.dcd")
    if "xtc" in mdtraj_output:
        mdtraj_frames.save_xtc("centered_old_coordinates.xtc")
    mdtraj_first_frame = mdtraj_frames[0:1]
    if "pdb" in mdtraj_output:
        mdtraj_first_frame.save_pdb("centered_old_coordinates_top.pdb")
    if "gro" in mdtraj_output:
        mdtraj_first_frame.save_gro("centered_old_coordinates_top.gro")


def MDanalysis_conversion(
    post_mdtraj_pdb_file,
    post_mdtraj_dcd_file,
    mda_output,
    output_selection,
    ligand_name=None,
    special_ligname=None,
):
    """Translate the trajectory so that all frames coincide with its center of geometry.

    Args:
        post_mdtraj_pdb_file (str): Name of the post-MDtraj PDB file.
        post_mdtraj_dcd_file (dcd file): Name of the post-MDtraj DCD File.
        ligand_name (str): Ligand name saved in the PDB file.
        special_ligname (str): Special residue name saved in the PDB file.
        mda_output (str): Selection of output formats.
        output_selection (str): Selection of topologies with specific atom selections that will be created.

    Returns:
        None
    """
    topology_trajectory = mda.Universe(post_mdtraj_pdb_file, post_mdtraj_dcd_file)
    topology_trajectory_all_atoms = topology_trajectory.select_atoms("all")
    # Translate the trajectory so that all frames coincide with its center of geometry
    topology_trajectory_all_atoms.atoms.translate(topology_trajectory_all_atoms.center_of_mass())
    topology_trajectory_protein_ligand = topology_trajectory.select_atoms(
        f"protein or resname {ligand_name} or resname {special_ligname}"
    )

    if "pdb" in mda_output:
        if output_selection != "mda_prot_lig":
            with mda.Writer("centered_traj_unaligned.dcd", topology_trajectory_all_atoms.n_atoms) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_all_atoms)
            topology_trajectory_all_atoms.write("centered_top.pdb")
            topology_ref_all_pdb = mda.Universe("centered_top.pdb")
            mobile_all_pdb = mda.Universe("centered_top.pdb", "centered_traj_unaligned.dcd")
            alignment_all_pdb = align.AlignTraj(
                mobile_all_pdb,
                topology_ref_all_pdb,
                select="protein and name CA",
                weights="mass",
                filename="centered_traj.dcd",
            )
            alignment_all_pdb.run()

        if output_selection != "mda_all":
            with mda.Writer(
                "prot_lig_traj_unaligned.dcd",
                topology_trajectory_protein_ligand.n_atoms,
            ) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_protein_ligand)
            topology_trajectory_protein_ligand.write("prot_lig_top.pdb")
            topology_ref_prot_lig_pdb = mda.Universe("prot_lig_top.pdb")
            mobile_prot_lig_pdb = mda.Universe("prot_lig_top.pdb", "prot_lig_traj_unaligned.dcd")
            alignment_prot_lig_pdb = align.AlignTraj(
                mobile_prot_lig_pdb,
                topology_ref_prot_lig_pdb,
                select="protein and name CA",
                weights="mass",
                filename="prot_lig_traj.dcd",
            )
            alignment_prot_lig_pdb.run()

    if "gro" in mda_output:
        if output_selection != "mda_prot_lig":
            with mda.Writer("centered_traj_unaligned.xtc", topology_trajectory_all_atoms.n_atoms) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_all_atoms)
            topology_trajectory_all_atoms.write("centered_top.gro")
            topology_ref_all_gro = mda.Universe("centered_top.gro")
            mobile_all_gro = mda.Universe("centered_top.gro", "centered_traj_unaligned.xtc")
            alignment_all_gro = align.AlignTraj(
                mobile_all_gro,
                topology_ref_all_gro,
                select="protein and name CA",
                weights="mass",
                filename="centered_traj.xtc",
            )
            alignment_all_gro.run()

        if output_selection != "mda_all":
            with mda.Writer(
                "prot_lig_traj_unaligned.xtc",
                topology_trajectory_protein_ligand.n_atoms,
            ) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_protein_ligand)
            topology_trajectory_protein_ligand.write("prot_lig_top.gro")
            topology_ref_prot_lig_gro = mda.Universe("prot_lig_top.gro")
            mobile_prot_lig_gro = mda.Universe("prot_lig_top.gro", "prot_lig_traj_unaligned.xtc")
            alignment_prot_lig_gro = align.AlignTraj(
                mobile_prot_lig_gro,
                topology_ref_prot_lig_gro,
                select="protein and name CA",
                weights="mass",
                filename="prot_lig_traj.xtc",
            )
            alignment_prot_lig_gro.run()

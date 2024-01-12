import mdtraj as md
import MDAnalysis as mda

from MDAnalysis.analysis import align
from MDAnalysis.analysis.distances import dist

######################################################################
# prepare the topology and trajectory file for further analysis
def mdtraj_conversion(pdb_file, mdtraj_output):
    """
    Recenter and apply periodic boundary conditions to the molecules in each frame of the trajectory, and save the centered trajectory and its first frame.

    Parameters
    ----------
    pdb_file: str
        name of the pdb file. This pdb file stores the extracted frames from the MD trajectory.
    mdtraj_output : str
        The selected format that will be used as an output of the topology and trajectory

    Returns
    -------
    None
    """    
    mdtraj_frames = md.load_dcd("trajectory.dcd", top=pdb_file)
    mdtraj_frames.image_molecules()
    if "dcd" in mdtraj_output:
        mdtraj_frames.save_dcd(f'centered_old_coordinates.dcd')
    if "xtc" in mdtraj_output:
        mdtraj_frames.save_xtc(f'centered_old_coordinates.xtc')
    mdtraj_first_frame = mdtraj_frames[0:1]
    if "pdb" in mdtraj_output:
        mdtraj_first_frame.save_pdb(f'centered_old_coordinates_top.pdb')
    if "gro" in mdtraj_output:
        mdtraj_first_frame.save_gro(f'centered_old_coordinates_top.gro')
    
def MDanalysis_conversion(post_mdtraj_pdb_file, post_mdtraj_dcd_file, mda_output, output_selection, ligand_name=None, special_ligname=None):
    """
    translate the trajectory so that all frames coincide with its center of geometry.

    Parameters
    ----------
    post_mdtraj_pdb_file: str
        Name of the post-MDtraj PDB file.
    post_mdtraj_dcd_file: dcd file
        Name of the post-MDtraj DCD File.
    ligand_name: str
        ligand name saved in PDB file.
    special_ligname: str
        special residue name saved in PDB file.
    mda_output: str
        Selection of Output formats.
    output_selection: str
        Selection of Topologies with specific atom selections that will be created

    Returns
    -------
    None
    """    
    topology_trajectory = mda.Universe(post_mdtraj_pdb_file, post_mdtraj_dcd_file)
    topology_trajectory_all_atoms = topology_trajectory.select_atoms("all")
    # translate the trajectoy so that all frames coincide with its center of geometry
    topology_trajectory_all_atoms.atoms.translate(topology_trajectory_all_atoms.center_of_mass())
    topology_trajectory_protein_ligand = topology_trajectory.select_atoms(f'protein or resname {ligand_name} or resname {special_ligname}')
    
    
    if "pdb" in mda_output:
        if output_selection != "mda_prot_lig":
            with mda.Writer(f'centered_traj_unaligned.dcd', topology_trajectory_all_atoms.n_atoms) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_all_atoms)
            topology_trajectory_all_atoms.write(f'centered_top.pdb')
            topology_ref_all_pdb = mda.Universe(f'centered_top.pdb')
            mobile_all_pdb = mda.Universe(f'centered_top.pdb', f'centered_traj_unaligned.dcd')
            alignment_all_pdb = align.AlignTraj(mobile_all_pdb, topology_ref_all_pdb, select="protein and name CA", weights="mass", filename=f'centered_traj.dcd')
            alignment_all_pdb.run()

        if output_selection != "mda_all":
            with mda.Writer(f'prot_lig_traj_unaligned.dcd', topology_trajectory_protein_ligand.n_atoms) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_protein_ligand)
            topology_trajectory_protein_ligand.write(f'prot_lig_top.pdb')
            topology_ref_prot_lig_pdb = mda.Universe(f'prot_lig_top.pdb')
            mobile_prot_lig_pdb = mda.Universe(f'prot_lig_top.pdb', f'prot_lig_traj_unaligned.dcd')
            alignment_prot_lig_pdb = align.AlignTraj(mobile_prot_lig_pdb, topology_ref_prot_lig_pdb, select="protein and name CA", weights="mass", filename=f'prot_lig_traj.dcd')
            alignment_prot_lig_pdb.run()

    if "gro" in mda_output:
        if output_selection != "mda_prot_lig":
            with mda.Writer(f'centered_traj_unaligned.xtc', topology_trajectory_all_atoms.n_atoms) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_all_atoms)
            topology_trajectory_all_atoms.write(f'centered_top.gro')
            topology_ref_all_gro = mda.Universe(f'centered_top.gro')
            mobile_all_gro = mda.Universe(f'centered_top.gro', f'centered_traj_unaligned.xtc')
            alignment_all_gro = align.AlignTraj(mobile_all_gro, topology_ref_all_gro, select="protein and name CA", weights="mass", filename=f'centered_traj.xtc')
            alignment_all_gro.run()

        if output_selection != "mda_all":
            with mda.Writer(f'prot_lig_traj_unaligned.xtc', topology_trajectory_protein_ligand.n_atoms) as w:
                for ts in topology_trajectory.trajectory[1:]:
                    w.write(topology_trajectory_protein_ligand)
            topology_trajectory_protein_ligand.write(f'prot_lig_top.gro')
            topology_ref_prot_lig_gro = mda.Universe(f'prot_lig_top.gro')
            mobile_prot_lig_gro = mda.Universe(f'prot_lig_top.gro', f'prot_lig_traj_unaligned.xtc')
            alignment_prot_lig_gro = align.AlignTraj(mobile_prot_lig_gro, topology_ref_prot_lig_gro, select="protein and name CA", weights="mass", filename=f'prot_lig_traj.xtc')
            alignment_prot_lig_gro.run()

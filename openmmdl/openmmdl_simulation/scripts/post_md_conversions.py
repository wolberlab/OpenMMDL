import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import mdtraj as md
import MDAnalysis as mda
import MDAnalysis.transformations as trans

from MDAnalysis.analysis import rms, diffusionmap, align
from MDAnalysis.analysis.distances import dist
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

######################################################################
# prepare the topology and trajectory file for further analysis
def mdtraj_conversion(pdb_file):
    """
    Recenter and apply periodic boundary conditions to the molecules in each frame of the trajectory, and save the centered trajectory and its first frame.

    Parameters
    ----------
    pdb_file: str
        name of the pdb file. This pdb file stores the extracted frames from the MD trajectory.

    Returns
    -------
    None
    """    
    mdtraj_frames = md.load_dcd("trajectory.dcd", top=pdb_file)
    mdtraj_frames.image_molecules()
    mdtraj_frames.save_dcd(f'centered_old_coordinates.dcd')
    mdtraj_frames.save_xtc(f'centered_old_coordinates.xtc')
    mdtraj_first_frame = mdtraj_frames[0:1]
    mdtraj_first_frame.save_pdb(f'centered_old_coordinates_top.pdb')
    mdtraj_first_frame.save_gro(f'centered_old_coordinates_top.gro')
    
def MDanalysis_conversion(post_mdtraj_pdb_file, post_mdtraj_dcd_file, ligand_name):
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

    Returns
    -------
    None
    """    
    topology_trajectory = mda.Universe(post_mdtraj_pdb_file, post_mdtraj_dcd_file)
    topology_trajectory_all_atoms = topology_trajectory.select_atoms("all")
    # translate the trajectoy so that all frames coincide with its center of geometry
    topology_trajectory_all_atoms.atoms.translate(topology_trajectory_all_atoms.center_of_mass())
    topology_trajectory_protein_ligand = topology_trajectory.select_atoms(f'protein or resname {ligand_name}')

    with mda.Writer(f'centered_traj.dcd', topology_trajectory_all_atoms.n_atoms) as w:
        for ts in topology_trajectory.trajectory[1:]:
            w.write(topology_trajectory_all_atoms)
    topology_trajectory_all_atoms.write(f'centered_top.pdb')

    with mda.Writer(f'centered_traj.xtc', topology_trajectory_all_atoms.n_atoms) as w:
        for ts in topology_trajectory.trajectory[1:]:
            w.write(topology_trajectory_all_atoms)
    topology_trajectory_all_atoms.write(f'centered_top.gro')
    
    with mda.Writer(f'prot_lig_traj.dcd', topology_trajectory_protein_ligand.n_atoms) as w:
        for ts in topology_trajectory.trajectory[1:]:
            w.write(topology_trajectory_protein_ligand)
    topology_trajectory_protein_ligand.write(f'prot_lig_top.pdb')

    with mda.Writer(f'prot_lig_traj.xtc', topology_trajectory_protein_ligand.n_atoms) as w:
        for ts in topology_trajectory.trajectory[1:]:
            w.write(topology_trajectory_protein_ligand)
    topology_trajectory_protein_ligand.write(f'prot_lig_top.gro')

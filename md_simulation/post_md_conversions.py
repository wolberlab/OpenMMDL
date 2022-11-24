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
        name of the pdb file, i.e. "output...pdb". This pdb file actually stores the whole extracted frames from the MD trajectory"

    Returns
    -------
    'centered_old_coord.dcd': dcd file
        dcd file containing the recentered coordinates
    'centered_old_coord.pdb': pdb file
        pdb file, the first frame of recentered MD trajectories
    """    
    mdtraj_frames = md.load(pdb_file)
    mdtraj_frames.image_molecules()
    mdtraj_frames.save_dcd(f'centered_old_coordinates.dcd')
    mdtraj_first_frame = mdtraj_frames[0:1]
    mdtraj_first_frame.save_pdb(f'centered_old_coorinates.pdb')
    
def MDanalysis_conversion(post_mdtraj_pdb_file, post_mdtraj_dcd_file, ligand_name):
    """
    translate the trajectoy so that all frames conincdie its center of geomertry

    Parameters
    ----------
    post_mdtraj_pdb_file: pdb file
        name of the pdb file, i.e. 'centered_old_coord.pdb'
    post_mdtraj_dcd_file: dcd file
        name of the dcd file, i.e. 'centered_old_coord.dcd'
    ligand_name: str
        ligand name saved in above pdb file.

    Returns
    -------
    'centered_traj.dcd': dcd file
        dcd file containing all particles, including protein, ligand, solvent, ions.
    'centered_top.pdb': pdb file
        pdb file containing all particiles. This is used to stroe the topology information
    'prot_lig_traj.dcd': dcd file
        dcd file containing only the complex, i.e, protein and ligand.
    'prot_lig_top.pdb'
        pdb file containing only the complex. This is used to stroe the topology information

    """    
    topology_trajectory = mda.Universe(post_mdtraj_pdb_file, post_mdtraj_dcd_file)
    topology_trajectory_all_atoms = topology_trajectory.select_atoms("all")
    # translate the trajectoy so that all frames conincdie its center of geomertry
    topology_trajectory_all_atoms.atoms.translate(topology_trajectory_all_atoms.center_of_mass())
    topology_trajectory_protein_ligand = topology_trajectory.select_atoms(f'protein or resname {ligand_name}')

    with mda.Writer(f'centered_traj.dcd', topology_trajectory_all_atoms.n_atoms) as w:
        for ts in topology_trajectory.trajectory[1:]:
            w.write(topology_trajectory_all_atoms)
    topology_trajectory_all_atoms.write(f'centered_top.pdb')
    
    with mda.Writer(f'prot_lig_traj.dcd', topology_trajectory_protein_ligand.n_atoms) as w:
        for ts in topology_trajectory.trajectory[1:]:
            w.write(topology_trajectory_protein_ligand)
    topology_trajectory_protein_ligand.write(f'prot_lig_top.pdb')


# RMSD calculation
def rmsd_for_atomgroups(prot_lig_top_file, prot_lig_traj_file, selection1, selection2=None):
    """Calulate the RMSD for selected atom groups, and save the csv file and plot.

    Parameters
    ----------
    prot_lig_top_file: pdb file
        name of the pdb file, i.e. 'prot_lig_top.pdb'
    prot_lig_traj_file:
        name of the dcd file, i.e. 'prot_lig_traj.dcd'
    selection1: str
        Selection string for main atom group, also used during alignment.
    selection2: list of str, optional
        Selection strings for additional atom groups.

    Returns
    -------
    rmsd_df: pandas.core.frame.DataFrame
        DataFrame containing RMSD of the selected atom groups over time.
    
    """ 
    universe=mda.Universe(prot_lig_top_file, prot_lig_traj_file)    
    universe.trajectory[0]
    ref = universe
    rmsd_analysis = rms.RMSD(universe, ref, select=selection1, groupselections=selection2)
    rmsd_analysis.run()
    columns = [selection1, *selection2] if selection2 else [selection1]
    rmsd_df = pd.DataFrame(np.round(rmsd_analysis.rmsd[:, 2:], 2), columns=columns)
    rmsd_df.index.name = "frame"

    rmsd_df.to_csv('RMSD_over_time.csv', sep=' ')

    rmsd_df.plot(title="RMSD of protein and ligand")
    plt.ylabel("RMSD (Å)")
    plt.savefig('RMSD_over_time.png')

    return rmsd_df

def RMSD_dist_frames(prot_lig_top_file, prot_lig_traj_file, lig):
    """Calculate the RMSD between all frames in a matrix.

    Parameters
    ----------
    prot_lig_top_file: pdb file
        name of the pdb file, i.e. 'prot_lig_top.pdb'.
    prot_lig_traj_file:
        name of the dcd file, i.e. 'prot_lig_traj.dcd'.
    lig: str
        ligand name saved in above pdb file. Selection string for the atomgroup to be investigated, also used during alignment.

    Returns
    -------
    pairwise_rmsd_prot: np.ndarray
        Numpy array of RMSD values for pairwise protein structures.
    pairwise_rmsd_lig: np.ndarray
        Numpy array of RMSD values for ligand structures.

    """
    universe=mda.Universe(prot_lig_top_file, prot_lig_traj_file)
    pairwise_rmsd_prot = diffusionmap.DistanceMatrix(universe, select="protein").run().dist_matrix
    pairwise_rmsd_lig = diffusionmap.DistanceMatrix(universe, f"resname {lig}").run().dist_matrix
    print(type(pairwise_rmsd_lig))

    max_dist = max(np.amax(pairwise_rmsd_lig), np.amax(pairwise_rmsd_prot))
    
    fig, ax = plt.subplots(1,2)
    fig.suptitle("RMSD between the frames")

    # protein image
    img1 = ax[0].imshow(pairwise_rmsd_prot, cmap="viridis", vmin=0, vmax=max_dist)
    ax[0].title.set_text("protein")
    ax[0].set_xlabel("frames")
    ax[0].set_ylabel("frames")
    
    # ligand image
    img2 = ax[1].imshow(pairwise_rmsd_lig, cmap="viridis", vmin=0, vmax=max_dist)
    ax[1].title.set_text("ligand")
    ax[1].set_xlabel("frames")

    fig.colorbar(img1, ax=ax, orientation="horizontal", fraction=0.1, label="RMSD (Å)")

    plt.savefig('RMSD_between_the_frames.png')
    return pairwise_rmsd_prot, pairwise_rmsd_lig

# Interaction analysis
def atomic_distance(
    prot_lig_top_file, 
    prot_lig_traj_file, 
    prot_resid, 
    prot_atom_name, 
    lig_name, 
    lig_id, 
    lig_atom_name):
    """ 
    calculate the distance between the interested pair atoms, and plot the distance vaule along the trajectory.

    Parameters
    ----------
    prot_lig_top_file: pdb file
        name of the pdb file, i.e. 'prot_lig_top.pdb'.
    prot_lig_traj_file:
        name of the dcd file, i.e. 'prot_lig_traj.dcd'.

    Returns
    -------

    """
    # select the atomgroup of protein and ligand separately.
    universe=mda.Universe(prot_lig_top_file, prot_lig_traj_file)
    atomgroup_prot = universe.select_atoms(
        f"resid {prot_resid} and name {prot_atom_name}"
    )
    atomgroup_lig = universe.select_atoms(f"resname {ligand_name} and name {lig_atom_name}")
    
    # calculate the distance along the trajectory
    distances = []
    for frame in universe.trajectory:
        distance = dist(atomgroup_prot, atomgroup_lig)
        distances.append(distance[2][0])
    
    # plot
    plt.plot(distances)
    plt.ylabel("distance (Å)")
    plt.xlabel("frame")
    plt.title(f"Atomic distance between prot-{prot_resid}-{prot_atom_name} and lig-{lig_id}-{lig_atom_name}")
    plt.savefig(f'Atomic_distance_between_prot-{prot_resid}-{prot_atom_name}_and_lig-{lig_id}-{lig_atom_name}.png')

    return


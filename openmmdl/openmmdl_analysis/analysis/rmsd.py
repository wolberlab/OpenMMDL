import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
import matplotlib.pyplot as plt
from numba import jit
from tqdm import tqdm
from MDAnalysis.analysis import rms, diffusionmap


class RMSDAnalyzer:
    """
    A class responsible for the Root-Mean-Square Deviation (RMSD) analysis throughout the molecular dynamics simulation. 
    The class provides functionalities to calculate RMSD over time, compute pairwise RMSD between trajectory frames 
    and identify the representative frames within clusters of the binding modes.

    Parameters
    ----------
    top_file : str
        Path to the topology file.
    traj_file : str
        Path to the trajectory file.

    Attributes
    ----------
    universe : mda.Universe
        MDAnalysis Universe object initialized with the provided topology and trajectory files.
    """
    def __init__(self, top_file, traj_file):
        self.universe = mda.Universe(top_file, traj_file)

    def rmsd_for_atomgroups(self, fig_type, selection1, selection2=None):
        """
        Calculate the RMSD for selected atom groups, and save the CSV file and plot.

        Parameters
        ----------
        fig_type : str
            Type of the figure to save (e.g., 'png', 'jpg').
        selection1 : str 
            Selection string for main atom group, also used during alignment.
        selection2 : list of str, optional 
            Selection strings for additional atom groups. Defaults to None.

        Returns
        -------
        pd.DataFrame
            DataFrame containing RMSD of the selected atom groups over time.
        """
        self.universe.trajectory[0]
        ref = self.universe
        rmsd_analysis = rms.RMSD(
            self.universe, ref, select=selection1, groupselections=selection2
        )
        rmsd_analysis.run()
        columns = [selection1, *selection2] if selection2 else [selection1]
        rmsd_df = pd.DataFrame(np.round(rmsd_analysis.rmsd[:, 2:], 2), columns=columns)
        rmsd_df.index.name = "frame"

        # Create the directory if it doesn't exist
        output_directory = "./RMSD/"
        os.makedirs(output_directory, exist_ok=True)

        # Save the RMSD values to a CSV file in the created directory
        rmsd_df.to_csv("./RMSD/RMSD_over_time.csv", sep=" ")

        # Plot and save the RMSD over time as a PNG file
        rmsd_df.plot(title="RMSD of protein and ligand")
        plt.ylabel("RMSD (Å)")
        plt.savefig(f"./RMSD/RMSD_over_time.{fig_type}")

        return rmsd_df

    def rmsd_dist_frames(self, fig_type, lig, nucleic=False):
        """
        Calculate the RMSD between all frames in a matrix.

        Parameters
        ----------
        fig_type : str 
            Type of the figure to save (e.g., 'png', 'jpg').
        lig : str 
            Ligand name saved in the above PDB file. Selection string for the MDAnalysis AtomGroup to be investigated, also used during alignment.
        nucleic : bool, optional 
            Bool indicating if the receptor to be analyzed contains nucleic acids. Defaults to False.

        Returns
        -------
        pairwise_rmsd_prot: np.ndarray
            Numpy array of RMSD values for pairwise protein structures.
        pairwise_rmsd_lig: np.ndarray 
            Numpy array of RMSD values for ligand structures.
        """
        if nucleic:
            pairwise_rmsd_prot = (
                diffusionmap.DistanceMatrix(self.universe, select="nucleic")
                .run()
                .dist_matrix
            )
        else:
            pairwise_rmsd_prot = (
                diffusionmap.DistanceMatrix(self.universe, select="protein")
                .run()
                .dist_matrix
            )
        pairwise_rmsd_lig = (
            diffusionmap.DistanceMatrix(self.universe, f"resname {lig}")
            .run()
            .dist_matrix
        )

        max_dist = max(np.amax(pairwise_rmsd_lig), np.amax(pairwise_rmsd_prot))

        fig, ax = plt.subplots(1, 2)
        fig.suptitle("RMSD between the frames")

        # protein image
        img1 = ax[0].imshow(pairwise_rmsd_prot, cmap="viridis", vmin=0, vmax=max_dist)
        if nucleic:
            ax[0].title.set_text("nucleic")
        else:
            ax[0].title.set_text("protein")
        ax[0].set_xlabel("frames")
        ax[0].set_ylabel("frames")

        # ligand image
        img2 = ax[1].imshow(pairwise_rmsd_lig, cmap="viridis", vmin=0, vmax=max_dist)
        ax[1].title.set_text("ligand")
        ax[1].set_xlabel("frames")

        fig.colorbar(
            img1, ax=ax, orientation="horizontal", fraction=0.1, label="RMSD (Å)"
        )

        plt.savefig(f"./RMSD/RMSD_between_the_frames.{fig_type}")
        
        return pairwise_rmsd_prot, pairwise_rmsd_lig

    def calculate_distance_matrix(self, selection):
        """
        Calculates the pairwise RMSD-based distance matrix for all trajectory frames 
        for the selected atom selection.
    
        Parameters
        ----------
        selection : str 
            Selection string for the atoms (e.g., 'protein', 'resname LIG') 
            used to compute the RMSD between frames.
    
        Returns
        -------
        np.ndarray
            Numpy array containing RMSD values between all pairs of frames.
        """
        distances = np.zeros(
            (len(self.universe.trajectory), len(self.universe.trajectory))
        )
        # calculate distance matrix
        for i in tqdm(
            range(len(self.universe.trajectory)),
            desc="\033[1mCalculating distance matrix:\033[0m",
        ):
            self.universe.trajectory[i]
            frame_i = self.universe.select_atoms(selection).positions
            # distances[i] = md.rmsd(traj_aligned, traj_aligned, frame=i)
            for j in range(i + 1, len(self.universe.trajectory)):
                self.universe.trajectory[j]
                frame_j = self.universe.select_atoms(selection).positions
                rmsd = self._calc_rmsd_2frames(frame_i, frame_j)
                distances[i][j] = rmsd
                distances[j][i] = rmsd
                
        return distances

    def calculate_representative_frame(self, bmode_frames, DM):
        """
        Calculates the most representative frame for a bindingmode. 
        This is based uppon the averagwe RMSD of a frame to all other frames in the binding mode.

        Parameters
        ----------
        bmode_frames : list of int 
            List of frames belonging to a binding mode.
        DM : np.ndarray 
            Distance matrix of trajectory.

        Returns
        -------
        int
            Number of the most representative frame.
        """
        frames = bmode_frames
        mean_rmsd_per_frame = {}
        # first loop  : first frame
        for frame_i in frames:
            mean_rmsd_per_frame[frame_i] = 0
            # we will add the rmsd between theses 2 frames and then calcul the
            # mean
            for frame_j in frames:
                # We don't want to calcul the same frame.
                if not frame_j == frame_i:
                    # we add to the corresponding value in the list of all rmsd
                    # the RMSD betwween frame_i and frame_j
                    mean_rmsd_per_frame[frame_i] += DM[frame_i - 1, frame_j - 1]
            # mean calculation
            mean_rmsd_per_frame[frame_i] /= len(frames)

            # Representative frame = frame with lower RMSD between all other
            # frame of the cluster
            repre = min(mean_rmsd_per_frame, key=mean_rmsd_per_frame.get)

        return repre

    def _calc_rmsd_2frames(self, ref, frame):
        """
        Calculates the RMSD between a reference and a target frame.
    
        This method serves as a wrapper for the `calc_rmsd_2frames_jit` function, 
        which dpes the actual RMSD calculation between two sets of coordinates.
    
        Parameters
        ----------
        ref : np.ndarray 
            Numpy array representing the reference atom positions, shape (N, 3).
        frame : np.ndarray 
            Numpy array representing the atom positions of the target frame, shape (N, 3).
    
        Returns
        -------
        float 
            The RMSD value between the reference and target frame.
        """
        return calc_rmsd_2frames_jit(ref, frame)


@jit(nopython=True, parallel=True, nogil=True)
def calc_rmsd_2frames_jit(ref, frame):
    """
    Calculates the RMSD between two frames of atomic coordinates.

    Parameters
    ----------
    ref : np.ndarray
        Numpy array containing the reference atomic coordinates.
    frame : np.ndarray 
        Numpy array containing the atomic coordinates of the target frame.

    Returns
    -------
    float
        The RMSD value between the reference and target frame.
    """
    dist = np.zeros(len(frame))
    for atom in range(len(frame)):
        dist[atom] = (
            (ref[atom][0] - frame[atom][0]) ** 2
            + (ref[atom][1] - frame[atom][1]) ** 2
            + (ref[atom][2] - frame[atom][2]) ** 2
        )

    return np.sqrt(dist.mean())

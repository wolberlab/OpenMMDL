import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

import MDAnalysis as mda

from MDAnalysis.analysis import rms, diffusionmap
from MDAnalysis.analysis.distances import dist


class RMSDAnalyzer:
    def __init__(self, prot_lig_top_file, prot_lig_traj_file):
        self.prot_lig_top_file = prot_lig_top_file
        self.prot_lig_traj_file = prot_lig_traj_file
        self.universe = mda.Universe(prot_lig_top_file, prot_lig_traj_file)
        
    def rmsd_for_atomgroups(self, fig_type, selection1, selection2=None):
        """Calculate the RMSD for selected atom groups, and save the csv file and plot.

        Args:
            fig_type (str): Type of the figure to save (e.g., 'png', 'jpg').
            selection1 (str): Selection string for main atom group, also used during alignment.
            selection2 (list, optional): Selection strings for additional atom groups. Defaults to None.

        Returns:
            pandas dataframe: rmsd_df. DataFrame containing RMSD of the selected atom groups over time.
        """
        self.universe.trajectory[0]
        ref = self.universe
        rmsd_analysis = rms.RMSD(self.universe, ref, select=selection1, groupselections=selection2)
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
        """Calculate the RMSD between all frames in a matrix.

        Args:
            fig_type (str): Type of the figure to save (e.g., 'png', 'jpg').
            lig (str): ligand name saved in the above pdb file. Selection string for the atomgroup to be investigated, also used during alignment.
            nucleic (bool, optional): Bool indicating if the receptor to be analyzed contains nucleic acids. Defaults to False.

        Returns:
            np.array: pairwise_rmsd_prot. Numpy array of RMSD values for pairwise protein structures.
            np.array: pairwise_rmsd_lig. Numpy array of RMSD values for ligand structures.
        """
        if nucleic:
            pairwise_rmsd_prot = diffusionmap.DistanceMatrix(self.universe, select="nucleic").run().dist_matrix
        else:
            pairwise_rmsd_prot = diffusionmap.DistanceMatrix(self.universe, select="protein").run().dist_matrix
        pairwise_rmsd_lig = diffusionmap.DistanceMatrix(self.universe, f"resname {lig}").run().dist_matrix

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

        fig.colorbar(img1, ax=ax, orientation="horizontal", fraction=0.1, label="RMSD (Å)")

        plt.savefig(f"./RMSD/RMSD_between_the_frames.{fig_type}")
        return pairwise_rmsd_prot, pairwise_rmsd_lig

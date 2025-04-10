import os
import pytest
import contextlib
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md

from openmmdl.openmmdl_analysis.rmsd_calculation import (
    rmsd_for_atomgroups,
    RMSD_dist_frames,
)

test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/0_unk_hoh.pdb"
trajectory_file = f"{test_data_directory}/all_50.dcd"
fig_type = "png"
selection1 = "protein"
selection2 = ("resname UNK", "")
ligand_name = "UNK"


def test_rmsd_for_atomgroups():

    # Call the function
    rmsd_df = rmsd_for_atomgroups(
        topology_file, trajectory_file, fig_type, selection1, selection2
    )

    # Check if the output DataFrame has the correct structure
    assert isinstance(rmsd_df, pd.DataFrame)
    assert rmsd_df.index.name == "frame"

    # Define file paths
    csv_path = os.path.join("RMSD", "RMSD_over_time.csv")
    plot_path = os.path.join("RMSD", "RMSD_over_time.png")

    print("Checking CSV file:", csv_path)
    # Check if the CSV file exists
    assert os.path.exists(csv_path), f"CSV file does not exist at {csv_path}"

    print("Checking plot file:", plot_path)
    # Check if the plot file exists
    assert os.path.exists(plot_path), f"Plot file does not exist at {plot_path}"

    # Cleanup created files after the test
    os.remove(csv_path)
    os.remove(plot_path)


def test_rmsd_dist_frames():

    # Call the function
    pairwise_rmsd_prot, pairwise_rmsd_lig = RMSD_dist_frames(
        topology_file, trajectory_file, fig_type, ligand_name
    )

    # Check if the function returns numpy arrays for pairwise RMSD
    assert isinstance(pairwise_rmsd_prot, np.ndarray)
    assert isinstance(pairwise_rmsd_lig, np.ndarray)

    # Define file paths
    plot_path = "./RMSD/RMSD_between_the_frames.png"

    print("Checking plot file:", plot_path)
    # Check if the plot file exists
    assert os.path.exists(plot_path), f"Plot file does not exist at {plot_path}"

    # Cleanup created files after the test
    with contextlib.suppress(FileNotFoundError):
        os.remove(plot_path)

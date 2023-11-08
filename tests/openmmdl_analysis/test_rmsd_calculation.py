import os
import pytest
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md

from openmmdl.openmmdl_analysis.rmsd_calculation import rmsd_for_atomgroups, RMSD_dist_frames

test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/0_unk_hoh.pdb"
trajectory_file = f"{test_data_directory}/all_50.dcd"
selection1 = "protein"
selection2 = ["resname UNK"]
ligand_name = "UNK"

def test_rmsd_for_atomgroups():

    # Call the function
    rmsd_df = rmsd_for_atomgroups(topology_file, trajectory_file, selection1, selection2)

    # Check if the output DataFrame has the correct structure
    assert isinstance(rmsd_df, pd.DataFrame)
    assert rmsd_df.index.name == "frame"
    
    # Check if the CSV file exists
    assert os.path.exists("RMSD_over_time.csv")
    
    # Check if the plot file exists
    assert os.path.exists("RMSD_over_time.png")
    
    # Cleanup created files after the test
    os.remove("RMSD_over_time.csv")
    os.remove("RMSD_over_time.png")

def test_rmsd_dist_frames():

    # Call the function
    pairwise_rmsd_prot, pairwise_rmsd_lig = RMSD_dist_frames(topology_file, trajectory_file, ligand_name)

    # Check if the function returns numpy arrays for pairwise RMSD
    assert isinstance(pairwise_rmsd_prot, np.ndarray)
    assert isinstance(pairwise_rmsd_lig, np.ndarray)

    # Check if the plot file exists
    assert os.path.exists("RMSD_between_the_frames.png")

    # Cleanup created files after the test
    os.remove("RMSD_between_the_frames.png")

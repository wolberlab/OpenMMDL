import pytest
import pandas as pd
import os
import shutil
from Bio.PDB.Structure import Structure
from pathlib import Path
from unittest.mock import patch, mock_open
from openmmdl.openmmdl_analysis.find_stable_waters import (perform_clustering_and_writing, process_trajectory_and_cluster, 
                              filter_and_parse_pdb, find_interacting_residues, 
                              read_pdb_as_dataframe, analyze_protein_and_water_interaction)

# Fixtures and mock data setup

test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/metal_top.pdb"
trajectory_file = f"{test_data_directory}/metal_traj_25.dcd"
output_dirs = []

def test_process_trajectory_and_cluster():
    water_eps_values = [0.5,1.0,2.0]  # Example epsilon values

    for water_eps in water_eps_values:
        output_directory = f"./test_output"
        process_trajectory_and_cluster(topology_file, trajectory_file, water_eps, output_directory)

        strEps = str(water_eps).replace(".", "")
        output_directory = f"./test_output_clusterEps_{strEps}"
        # Check if the expected output directory is created
        
        assert os.path.isdir(output_directory), f"Directory {output_directory} was not created"
        output_dirs.append(output_directory)

        # Check if stable_waters.csv is created
        print("Check-CSV")
        csv_file = os.path.join(output_directory, "stable_waters.csv")
        assert os.path.isfile(csv_file)

        # Load and verify the data in stable_waters.csv
        print("Check-CSV-content")
        stable_waters_df = pd.read_csv(csv_file)
        assert not stable_waters_df.empty
        assert set(stable_waters_df.columns) == {"Frame", "Residue", "Oxygen_X", "Oxygen_Y", "Oxygen_Z"}

    # Cleanup: remove created directories and files
    for dir in output_dirs:
        shutil.rmtree(dir)




def test_filter_and_parse_pdb():
    # Call the function with the sample PDB file
    structure = filter_and_parse_pdb(topology_file)

    # Check if the returned object is a Structure
    assert isinstance(structure, Structure), "The returned object is not a Structure"




#### FOR THE NEXT ONE PLEASE DON'T KEEP THE NESTED FUNCTIONS NESTED!
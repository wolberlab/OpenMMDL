import os
import pytest
import shutil
import tempfile
from pathlib import Path
import pandas as pd
import numpy as np
import mdtraj as md
import MDAnalysis as mda
import unittest
from unittest.mock import Mock, patch
from plip.structure.preparation import PDBComplex, LigandFinder, Mol, PLInteraction

from openmmdl.openmmdl_analysis.interaction_gathering import *


test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/complex.pdb"
frame_file = f"{test_data_directory}/processing_frame_1.pdb"
topology_metal = f"{test_data_directory}/metal_top.pdb"
trajetory_metal = f"{test_data_directory}/metal_traj_25.dcd"

binding_site_id = "UNK:X:0"
lig_name = "UNK"
peptide = "X"

# Test the function
def test_characterize_complex():
    # Call the function
    interaction_set = characterize_complex(topology_file, binding_site_id)

    # Check if the function returns a PLInteraction object
    assert isinstance(interaction_set, PLInteraction)

def test_retrieve_plip_interactions():
    # Call the function
    interactions = retrieve_plip_interactions(topology_file, lig_name)

    # Check if the function returns a dictionary
    assert isinstance(interactions, dict)

def test_retrieve_plip_interactions_peptide():
    # Call the function
    interactions = retrieve_plip_interactions_peptide(topology_file, peptide)

    # Check if the function returns a dictionary
    assert isinstance(interactions, dict)

# Define test data
sample_interactions = {
    "hydrophobic": [["Column1", "Column2"], [1, 2], [3, 4]],
    "hbond": [["ColumnA", "ColumnB"], ["A", "B"], ["C", "D"]],
}

def test_create_df_from_binding_site():
    # Test with valid interaction type
    df = create_df_from_binding_site(sample_interactions, interaction_type="hydrophobic")
    assert isinstance(df, pd.DataFrame)
    assert df.shape == (2, 2)
    assert list(df.columns) == ["Column1", "Column2"]

    # Test with default interaction type
    df_default = create_df_from_binding_site(sample_interactions)
    assert isinstance(df_default, pd.DataFrame)
    assert df_default.shape == (2, 2)
    assert list(df_default.columns) == ["ColumnA", "ColumnB"]

    # Test with an invalid interaction type (should default to 'hbond')
    df_invalid = create_df_from_binding_site(sample_interactions, interaction_type="invalid_type")
    assert isinstance(df_invalid, pd.DataFrame)
    assert df_invalid.shape == (2, 2)
    assert list(df_invalid.columns) == ["ColumnA", "ColumnB"]


@pytest.fixture
def input_pdb_filename(tmp_path):
    input_pdb_filename = tmp_path / "input.pdb"
    
    # Create a mock PDB file with 10 atoms
    input_pdb_content = """ATOM      1  N   UNK A 454      43.493  48.319  35.835  1.00  0.00      A    N  
ATOM      2  N1  UNK A 454      44.740  47.862  35.697  1.00  0.00      A    N  
ATOM      3  C14 UNK A 454      44.608  46.866  34.829  1.00  0.00      A    C  
ATOM      4  N2  UNK A 454      43.265  46.644  34.450  1.00  0.00      A    N  
ATOM      5  C7  UNK A 454      42.607  47.556  35.077  1.00  0.00      A    C  
ATOM      6  H5  UNK A 454      41.542  47.701  34.954  1.00  0.00      A    H  
ATOM      7  H10 UNK A 454      45.308  46.132  34.453  1.00  0.00      A    H  
ATOM      8  C   UNK A 454      43.168  49.513  36.656  1.00  0.00      A    C  
ATOM      9  C2  UNK A 454      42.743  50.705  35.818  1.00  0.00      A    C  
ATOM     10  C4  UNK A 454      43.545  51.052  34.671  1.00  0.00      A    C"""

    input_pdb_filename.write_text(input_pdb_content)
    return input_pdb_filename

def test_change_lig_to_residue():
    topology_file = f"{test_data_directory}/complex.pdb"
    shutil.copy(str(topology_file), '.')
    topology_file = "complex.pdb"

    # Change ligand to residue
    change_lig_to_residue(str(topology_file), 'UNK', 'NEW')

    # Read the output PDB file and check if residues are modified
    with open(topology_file, 'r') as output_file:
        modified_lines = output_file.readlines()
        assert any('NEW' in line for line in modified_lines)
        assert all('UNK' not in line for line in modified_lines)


def test_process_frame_with_sample_data():
    # Define a sample frame number
    frame_number = 1

    destination_file = "processing_frame_1.pdb"

    shutil.copy(frame_file, destination_file)

    # Load the sample PDB file into an MDAnalysis Universe
    sample_universe = mda.Universe(topology_file)

    # Call the process_frame function with the sample data
    result = process_frame(frame_number, sample_universe, lig_name)

    # Define the expected columns you want to check
    expected_columns = ["FRAME", "INTERACTION"]  # Add the specific columns you want to validate

    # Check if the result is a Pandas DataFrame
    assert isinstance(result, pd.DataFrame)

    # Check if all expected columns are present in the result
    for column in expected_columns:
        assert column in result.columns

def test_process_trajectory():
    topology_file = f"{test_data_directory}/0_unk_hoh.pdb"
    trajectory_file = f"{test_data_directory}/all_50.dcd"
    pdb_md = mda.Universe(topology_file,trajectory_file)
    dataframe = None
    num_processes = 2
    lig_name = "UNK"

    interaction_list = pd.DataFrame(columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "DIST", "LIGCARBONIDX", "PROTCARBONIDX", "LIGCOO", "PROTCOO"])

    interaction_list = process_trajectory(pdb_md, dataframe, num_processes, lig_name, special_ligand=None, peptide=None)

    assert interaction_list is not None
    assert len(interaction_list) > 10
    
def test_process_frame_special_with_files():
    test_data_directory = "openmmdl/tests/data/in"  # Replace with the actual path to your test data directory
    topology_metal = f"{test_data_directory}/metal_top.pdb"
    trajetory_metal = f"{test_data_directory}/metal_traj_25.dcd"

    # Load PDB and DCD files using mdanalysis.Universe
    import MDAnalysis as mda
    u = mda.Universe(topology_metal, trajetory_metal)

    lig_name = "UNK"  # Replace with the actual ligand name
    special = "HEM"  # Replace with the actual special residue name
    frame = 0

    result = process_frame_special(frame, u, lig_name, special)

    assert isinstance(result, list)
    assert all(isinstance(df, pd.DataFrame) for df in result)

    # Add more specific assertions based on the expected behavior of the function
    # For example, check if the columns in the DataFrame are as expected, or if certain conditions hold

    # Clean up any temporary files created during the test
    for frame in range(len(u.trajectory)):
        temp_file = f'processing_frame_{frame}.pdb'
        if os.path.exists(temp_file):
            os.remove(temp_file)

def test_process_frame_wrapper():

    test_data_directory = "openmmdl/tests/data/in"  # Replace with the actual path to your test data directory
    topology_metal = f"{test_data_directory}/metal_top.pdb"
    trajetory_metal = f"{test_data_directory}/metal_traj_25.dcd"
    ligand_special = f"{test_data_directory}/ligand_special.pdb"
    shutil.copy(str(topology_metal), '.')
    shutil.copy(str(trajetory_metal), '.')
    shutil.copy(str(ligand_special), '.')
    topology_metal = "metal_top.pdb"
    trajetory_metal = "metal_traj_25.dcd"

    # Load PDB and DCD files using MDAnalysis
    pdb_md = mda.Universe(topology_metal, trajetory_metal)
    lig_name = "UNK"  # Replace with the actual ligand name
    special_ligand = "HEM"  # Replace with the actual special ligand name
    peptide = None  # Replace with the actual peptide name
    frame_idx = 2

    args = (frame_idx, pdb_md, lig_name, special_ligand, peptide)
    result = process_frame_wrapper(args)

    # Perform assertions based on the expected behavior of the process_frame_special function
    assert isinstance(result, tuple)
    assert len(result) == 2
    assert isinstance(result[0], int)


def test_fill_missing_frames():
    # Create a sample DataFrame with missing frames
    data = {'FRAME': [1, 2, 4, 5],
            'Value1': ['A', 'B', 'C', 'D']}
    df = pd.DataFrame(data)
    md_len = 6


    # Call the fill_missing_frames function
    filled_df = fill_missing_frames(df, md_len)  # md_len = 6, should include frames 1 to 5

    # Assert that all frame numbers from 1 to 5 are present in the 'FRAME' column
    assert all(filled_df['FRAME'] == [1, 2, 3, 4, 5])

    # Assert that missing frames have 'Value1' column set to "skip"
    assert all(filled_df.loc[filled_df['FRAME'] == 3, 'Value1'] == 'skip')

    assert isinstance(filled_df, pd.DataFrame)
if __name__ == "__main":
    pytest.main()

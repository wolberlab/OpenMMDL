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

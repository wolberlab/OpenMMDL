import numpy as np
import pandas as pd
import re
import os
from pathlib import Path
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import pytest
from openmmdl.openmmdl_analysis.pml_writer import *

#pml_writer tests
@pytest.fixture
def sample_dataframe_generate_pharmacophore_centers():
    data = {
        'Hydrophobic': [1, 1, 0, 1, 0],
        'Ionic': [0, 1, 0, 0, 1],
        'LIGCOO': ["(1.0, 2.0, 3.0)", "(2.0, 3.0, 4.0)", "(3.0, 4.0, 5.0)", "(4.0, 5.0, 6.0)", "(5.0, 6.0, 7.0)"]
    }
    df = pd.DataFrame(data)
    return df

@pytest.fixture
def sample_interactions_generate_pharmacophore_centers():
    return ['Hydrophobic', 'Ionic']

def test_generate_pharmacophore_centers(sample_dataframe_generate_pharmacophore_centers, sample_interactions_generate_pharmacophore_centers):
    result = generate_pharmacophore_centers(sample_dataframe_generate_pharmacophore_centers, sample_interactions_generate_pharmacophore_centers)
    
    expected_pharmacophore = {
        'Hydrophobic': [2.333, 3.333, 4.333],
        'Ionic': [3.5, 4.5, 5.5]
    }
    
    assert result == expected_pharmacophore


@pytest.fixture
def sample_dataframe_generate_pharmacophore_vectors():
    # Create a sample dataframe for testing
    data = {
        'HBDonors': [1, 0, 1, 0, 1],
        'HBAcceptors': [0, 1, 0, 1, 0],
        'LIGCOO': [
            "(1.0, 2.0, 3.0)",
            "(2.0, 3.0, 4.0)",
            "(3.0, 4.0, 5.0)",
            "(4.0, 5.0, 6.0)",
            "(5.0, 6.0, 7.0)"
        ],
        'PROTCOO': [
            "(0.5, 1.5, 2.5)",
            "(1.5, 2.5, 3.5)",
            "(2.5, 3.5, 4.5)",
            "(3.5, 4.5, 5.5)",
            "(4.5, 5.5, 6.5)"
        ]
    }
    df = pd.DataFrame(data)
    return df

@pytest.fixture
def sample_interactions_generate_pharmacophore_vectors():
    return ['HBDonors', 'HBAcceptors']

def test_generate_pharmacophore_vectors(sample_dataframe_generate_pharmacophore_vectors, sample_interactions_generate_pharmacophore_vectors):
    result = generate_pharmacophore_vectors(sample_dataframe_generate_pharmacophore_vectors, sample_interactions_generate_pharmacophore_vectors)
    
    expected_pharmacophore = {
        'HBDonors': [
            [3.0, 4.0, 5.0],
            [2.5, 3.5, 4.5]
        ],
        'HBAcceptors': [
            [3.0, 4.0, 5.0],
            [2.5, 3.5, 4.5]
        ]
    }

    assert result == expected_pharmacophore

def test_generate_md_pharmacophore_cloudcenters(tmp_path):
    # Sample data for the DataFrame
    data = {
        'Acceptor_hbond_1': [1, 0, 1, 0, 1],
        'Donor_hbond_1': [0, 1, 0, 1, 0],
        'pistacking_1': [1, 0, 0, 1, 1],
        'hydrophobic_1': [0, 1, 0, 1, 0],
        'PI_saltbridge_1': [1, 0, 1, 0, 1],
        'NI_saltbridge_1': [0, 1, 0, 1, 0],
        'LIGCOO': ['(1.0, 2.0, 3.0)', '(2.0, 3.0, 4.0)', '(3.0, 4.0, 5.0)', '(4.0, 5.0, 6.0)', '(5.0, 6.0, 7.0)'],
        'PROTCOO': ['(7.0, 6.0, 5.0)', '(6.0, 5.0, 4.0)', '(5.0, 4.0, 3.0)', '(4.0, 3.0, 2.0)', '(3.0, 2.0, 1.0)'],
    }

    df = pd.DataFrame(data)

    # Output file paths
    output_filename = tmp_path / "test_output.pml"

    # Call the function
    generate_md_pharmacophore_cloudcenters(df, 'core_compound', output_filename, 'system_name', id_num=0)

    # Check if the output file is created
    assert os.path.isfile(output_filename), f"File {output_filename} not found."

    # Check if the generated XML is valid
    try:
        ET.parse(output_filename)
    except ET.ParseError:
        pytest.fail(f"Invalid XML in {output_filename}")


def test_generate_pharmacophore_centers_all_points():
    # Sample data for the DataFrame
    data = {
        'interaction1': [1, 0, 1, 0, 1],
        'interaction2': [0, 1, 0, 1, 0],
        'LIGCOO': ['(1.0, 2.0, 3.0)', '(2.0, 3.0, 4.0)', '(3.0, 4.0, 5.0)', '(4.0, 5.0, 6.0)', '(5.0, 6.0, 7.0)'],
    }

    df = pd.DataFrame(data)

    # Sample interactions
    interactions = ['interaction1', 'interaction2']

    # Call the function
    pharmacophore = generate_pharmacophore_centers_all_points(df, interactions)

    # Check if the generated pharmacophore has the expected structure
    assert isinstance(pharmacophore, dict), "Pharmacophore should be a dictionary."
    
    for interaction in interactions:
        assert interaction in pharmacophore, f"{interaction} not found in the generated pharmacophore."

        points = pharmacophore[interaction]
        assert isinstance(points, list), f"Pharmacophore points for {interaction} should be a list."
        
        # Check if the points have the expected structure
        for point in points:
            assert isinstance(point, list) and len(point) == 3, "Each point should be a list of three coordinates."

def test_generate_point_cloud_pml(tmp_path):
    # Sample data for the cloud_dict
    cloud_dict = {
        'feature1': {
            'interaction1': [(1.0, 2.0, 3.0), (1.5, 2.5, 3.5), (2.0, 3.0, 4.0)],
            'interaction2': [(2.0, 3.0, 4.0), (2.5, 3.5, 4.5), (3.0, 4.0, 5.0)],
        },
        'feature2': {
            'interaction3': [(3.0, 4.0, 5.0), (3.5, 4.5, 5.5), (4.0, 5.0, 6.0)],
        },
    }

    # Output file paths
    outname = tmp_path / "test_output"
    outname_pml = f"{outname}.pml"

    # Call the function
    generate_point_cloud_pml(cloud_dict, "system_name", outname)

    # Check if the output file is created
    assert os.path.isfile(outname_pml), f"File {outname_pml} not found."

    # Check if the generated XML is valid
    try:
        ET.parse(outname_pml)
    except ET.ParseError:
        pytest.fail(f"Invalid XML in {outname_pml}")
        

def test_generate_bindingmode_pharmacophore(tmp_path):
    # Prepare inputs
    dict_bindingmode = {
        "Acceptor_hbond": {
            "PROTCOO": [[1, 2, 3]],
            "LIGCOO": [[4, 5, 6]]
        }
    }
    core_compound = "ligand"
    sysname = "system"
    outname = "test_output"
    id_num = 0

    # Change the current working directory to tmp_path
    os.chdir(tmp_path)
    os.mkdir("Binding_Modes_Markov_States")

    # Call the function
    generate_bindingmode_pharmacophore(dict_bindingmode, core_compound, sysname, outname, id_num)

    # Check if the output file is created
    outname_pml = f"./Binding_Modes_Markov_States/{outname}.pml"
    assert os.path.isfile(outname_pml), f"File {outname_pml} not found."

    # Check if the generated XML is valid
    try:
        ET.parse(outname_pml)
    except ET.ParseError:
        pytest.fail(f"Invalid XML in {outname_pml}")
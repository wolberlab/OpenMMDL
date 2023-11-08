import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
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
    
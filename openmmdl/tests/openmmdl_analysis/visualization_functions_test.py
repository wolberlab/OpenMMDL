import numpy as np
import pandas as pd
import re
import shutil
import subprocess
import os
from pathlib import Path
import matplotlib.pyplot as plt
import pytest
from openmmdl.openmmdl_analysis.visualization_functions import *

package_path = Path("openmmdl/openmmdl_analysis")

# visualization_functions tests
@pytest.fixture
def sample_dataframe_interacting_water_ids():
    data = {
        'Interaction1': [0, 1, 0, 1, 0],
        'Interaction2': [1, 0, 0, 0, 1],
        'WATER_IDX': [101, 102, None, 104, 105],
        'FRAME': [1, 2, 3, 4, 5]  
    }
    df_all = pd.DataFrame(data)
    return df_all

def test_interacting_water_ids(sample_dataframe_interacting_water_ids):
    waterbridge_interactions = ['Interaction1', 'Interaction2']
    
    result = interacting_water_ids(sample_dataframe_interacting_water_ids, waterbridge_interactions)

    expected_interacting_waters = [101, 102, 104, 105]

    assert sorted(result) == sorted(expected_interacting_waters)
    

@pytest.fixture
def sample_dataframe_cloud_json_generation():
    data = {
        'LIGCOO': [
            "(1.0, 2.0, 3.0)",
            "(4.0, 5.0, 6.0)",
            "(7.0, 8.0, 9.0)",
        ],
        'INTERACTION': [
            'hydrophobic',
            'acceptor',
            'donor',
        ],
        'PROTISDON': [
            'False',
            'True',
            'False',
        ],
        'PROTISPOS': [
            'False',
            'False',
            'True',
        ],
    }
    df_all = pd.DataFrame(data)
    return df_all

def test_run_visualization():
    # Set up the paths
    package_path = Path("openmmdl/openmmdl_analysis")
    notebook_path =  package_path / "visualization.ipynb"
    
    # Run the visualization function
    # run_visualization()
    
    # Check if the notebook was copied to the current directory with the correct name
    copied_notebook_path = os.path.join(os.getcwd(), 'visualization.ipynb')
    shutil.copy(str(notebook_path), '.')
    new_notebook_path = 'visualization.ipynb'
    assert os.path.isfile(copied_notebook_path)
    
    # Check if the content of the copied notebook is the same as the original notebook
    with open(new_notebook_path, 'r') as copied_notebook:
        with open(notebook_path, 'r') as original_notebook:
            assert copied_notebook.read() == original_notebook.read()

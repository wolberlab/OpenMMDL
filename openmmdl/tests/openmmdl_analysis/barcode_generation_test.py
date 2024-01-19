import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
import pytest
from openmmdl.openmmdl_analysis.barcode_generation import *

# Barcode generation tests
@pytest.fixture
def sample_dataframe_barcode_generation():
    data = {
        'FRAME': [1, 1, 2, 2, 3],
        'Interaction1': [1, 0, 1, 0, 0],
        'Interaction2': [0, 0, 0, 1, 1],
        'WATER_IDX': [101, 102, 103, 104, 105],
    }
    return pd.DataFrame(data)

def test_barcodegeneration(sample_dataframe_barcode_generation):
    interaction = 'Interaction1'
    barcode = barcodegeneration(sample_dataframe_barcode_generation, interaction)
    
    assert isinstance(barcode, np.ndarray)
    
    expected_barcode = np.array([1, 1, 0])
    assert np.array_equal(barcode, expected_barcode)
    
def test_waterids_barcode_generator(sample_dataframe_barcode_generation):
    interaction = 'Interaction2'
    waterid_barcode = waterids_barcode_generator(sample_dataframe_barcode_generation, interaction)
    
    # Test if the output is a list
    assert isinstance(waterid_barcode, list)
    
    # Test the expected waterid barcode for the sample dataframe and interaction
    expected_waterid_barcode = [0, 104, 105]
    assert waterid_barcode == expected_waterid_barcode

def test_plot_barcodes():
    # create barcode data
    working_directory = os.getcwd()
    # Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory before:", files_in_working_directory)

    # Test case 1: No barcode
    plot_barcodes({}, "no_barcodes.png")
    assert not os.path.isfile("no_barcodes.png")

    # Test case 2: Single barcode
    barcode_data = {'166ARGA_4220,4221_Carboxylate_NI_saltbridge': np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])}
    plot_barcodes(barcode_data, "single_barcode.png")
    single_barcode = "single_barcode.png"
    assert single_barcode is not None
    
    barcodes = {
        "Barcode 1": np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]),
        "Barcode 2": np.array([0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]),
        # Include more barcodes as needed
    }
    plot_barcodes(barcodes, "multiple_barcodes.png")

    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory after:", files_in_working_directory)
    save_path = "multiple_barcodes.png"
    
    assert save_path is not None


    def test_plot_waterbridge_piechart(tmp_path):
        # Prepare inputs
        df_all = pd.DataFrame({
            'interaction1': [1, 0, 1, 0],
            'interaction2': [0, 1, 0, 1],
            'waterid': [1, 2, 1, 2]
        })
        waterbridge_barcodes = [np.array([1, 0, 1, 0]), np.array([0, 1, 0, 1])]
        waterbridge_interactions = ['interaction1', 'interaction2']

        # Change the current working directory to tmp_path
        os.chdir(tmp_path)

        # Use os.makedirs instead of os.mkdir
        os.makedirs("./Barcodes/Waterbridge_Piecharts/", exist_ok=True)

        # Call the function
        plot_waterbridge_piechart(df_all, waterbridge_barcodes, waterbridge_interactions)

        # Check if the output files are created
        for interaction in waterbridge_interactions:
            outname_png = f"./Barcodes/Waterbridge_Piecharts/{interaction}.png"
            assert os.path.isfile(outname_png), f"File {outname_png} not found."
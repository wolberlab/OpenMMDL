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

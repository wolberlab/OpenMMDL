import numpy as np
import pandas as pd
import rdkit
import MDAnalysis as mda
import re
import os
import matplotlib.pyplot as plt
import pytest

from openmmdl.openmmdl_analysis.binding_mode_processing import *


# binding_mode_processing tests
@pytest.fixture
def sample_dataframe_bindingmode_processing():
    data = {
        'FRAME': {0: 1, 1: 2, 2: 3, 3: 2},
        'Prot_partner': {0: 'A', 1: 'B', 2: 'C', 3: 'A'},
        'INTERACTION': {0: 'hydrophobic', 1: 'hbond', 2: 'saltbridge', 3: 'hydrophobic'},
        'LIGCARBONIDX': {0: 101, 1: 102, 2: 103, 3: 102},
        'DONORIDX': {0: 201, 1: 202, 2: 203, 3: 202},
        'ACCEPTORIDX': {0: 301, 1: 302, 2: 303, 3: 302},
        'PROTISDON': {0: True, 1: False, 2: True, 3: False},
        'LIG_IDX_LIST': {0: [1, 2], 1: [3, 4], 2: [5, 6], 3: [3, 4]},
        'LIG_GROUP': {0: 'Group1', 1: 'Group2', 2: 'Group3', 3: 'Group1'},
        'PROTISPOS': {0: True, 1: False, 2: True, 3: True},
        'DON_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONORTYPE': {0: 0, 1: 0, 2: 0, 3: 0},
        'ACCEPTOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'LOCATION': {0: 0, 1: 0, 2: 0, 3: 0},
        'METAL_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'METAL_TYPE': {0: 0, 1: 0, 2: 0, 3: 0},
        'RESTYPE_LIG': {0: 0, 1: 0, 2: 0, 3: 0},
        'TARGET_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'COORDINATION': {0: 0, 1: 0, 2: 0, 3: 0}
    }

    # Add 'halogen' and 'hbond' data to the existing DataFrame
    data['FRAME'][4] = 4  # Add a new 'FRAME' value
    data['Prot_partner'][4] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][4] = 'halogen'  # Add 'halogen' interaction
    data['DON_IDX'][4] = 501  # DON_IDX for 'halogen'
    data['DONORTYPE'][4] = 'F'  # Halogen type
    data['ACCEPTOR_IDX'][4] = 0
    data['DONOR_IDX'][4] = 0
    data['LIG_IDX_LIST'][4] = 0
    data['LIG_GROUP'][4] = 0  # LIG_GROUP for 'pication
    data['RESTYPE_LIG'][4] = 0 
    data['TARGET_IDX'][4] = 0

    data['FRAME'][5] = 5  # Add a new 'FRAME' value
    data['Prot_partner'][5] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][5] = 'hbond'  # Add 'hbond' interaction
    data['ACCEPTORIDX'][5] = 301  # ACCEPTORIDX for 'hbond'
    data['DON_IDX'][5] = 0  # DON_IDX
    data['DONORTYPE'][5] = 0  # DON_IDX
    data['PROTISDON'][5] = True  # PROTISDON is True for 'hbond'
    data['ACCEPTOR_IDX'][5] = 0
    data['LIG_IDX_LIST'][5] = 0
    data['DONOR_IDX'][5] = 0
    data['LIG_GROUP'][5] = 0  # LIG_GROUP for 'pication
    data['RESTYPE_LIG'][5] = 0 
    data['TARGET_IDX'][5] = 0

    # Add 'waterbridge' cases where PROTISDON is both True and False
    data['FRAME'][6] = 6  # Add a new 'FRAME' value
    data['Prot_partner'][6] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][6] = 'waterbridge'  # Add 'waterbridge' interaction
    data['ACCEPTOR_IDX'][6] = 401  # ACCEPTOR_IDX for 'waterbridge'
    data['DON_IDX'][6] = 0  # DON_IDX
    data['DONORTYPE'][6] = 0  # DON_IDX
    data['DONOR_IDX'][6] = 0
    data['LIG_IDX_LIST'][6] = 0
    data['PROTISDON'][6] = True  # PROTISDON is True for 'waterbridge'
    data['LIG_GROUP'][6] = 0  # LIG_GROUP for 'pication
    data['RESTYPE_LIG'][6] = 0
    data['TARGET_IDX'][6] = 0

    data['FRAME'][7] = 7  # Add a new 'FRAME' value
    data['Prot_partner'][7] = 'B'  # Add a new 'Prot_partner' value
    data['INTERACTION'][7] = 'waterbridge'  # Add 'waterbridge' interaction
    data['DONOR_IDX'][7] = 501  # DONOR_IDX for 'waterbridge'
    data['DON_IDX'][7] = 0  # DON_IDX
    data['DONORTYPE'][7] = 0  # DON_IDX
    data['PROTISDON'][7] = False  # PROTISDON is False for 'waterbridge'
    data['ACCEPTOR_IDX'][7] = 0
    data['LIG_IDX_LIST'][7] = 0 # LIG_IDX_LIST for 'pication'
    data['LIG_GROUP'][7] = 0  # LIG_GROUP for 'pication
    data['RESTYPE_LIG'][7] = 0 
    data['TARGET_IDX'][7] = 0

    # Add 'pistacking' case
    data['FRAME'][8] = 8  # Add a new 'FRAME' value
    data['Prot_partner'][8] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][8] = 'pistacking'  # Add 'pistacking' interaction
    data['LIG_IDX_LIST'][8] = [7, 8]  # LIG_IDX_LIST for 'pistacking'
    data['LIG_GROUP'][8] = 0  # LIG_GROUP for 'pication
    data['ACCEPTOR_IDX'][8] = 0
    data['DON_IDX'][8] = 0  # DON_IDX
    data['DONOR_IDX'][8] = 0 
    data['PROTISDON'][8] = False
    data['DONORTYPE'][8] = 0  # DON_IDX
    data['RESTYPE_LIG'][8] = 0 
    data['TARGET_IDX'][8] = 0

    # Add 'pication' case
    data['FRAME'][9] = 9  # Add a new 'FRAME' value
    data['Prot_partner'][9] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][9] = 'pication'  # Add 'pication' interaction
    data['LIG_IDX_LIST'][9] = [9, 10]  # LIG_IDX_LIST for 'pication'
    data['LIG_GROUP'][9] = 'Group4'  # LIG_GROUP for 'pication'
    data['ACCEPTOR_IDX'][9] = 0
    data['DON_IDX'][9] = 0  # DON_IDX
    data['PROTISDON'][9] = False
    data['DONOR_IDX'][9] = 0 
    data['DONORTYPE'][9] = 0  # DON_IDX
    data['RESTYPE_LIG'][9] = 0 
    data['TARGET_IDX'][9] = 0
    
    # Add 'metal' interaction case
    data['FRAME'][10] = 10  # Add a new 'FRAME' value
    data['Prot_partner'][10] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][10] = 'metal'  # Add 'metal' interaction
    data['METAL_IDX'][10] = 401  # METAL_IDX for 'metal'
    data['METAL_TYPE'][10] = 'Fe'  # Metal type
    data['LOCATION'][10] = 'site1'  # Location
    data['ACCEPTOR_IDX'][10] = 0
    data['DONOR_IDX'][10] = 0
    data['RESTYPE_LIG'][10] = 'A' 
    data['TARGET_IDX'][10] = 401
    data['COORDINATION'][10] = "site1"

    data['FRAME'][11] = 11  # Add a new 'FRAME' value
    data['Prot_partner'][11] = 'A'  # Add a new 'Prot_partner' value
    data['INTERACTION'][11] = 'saltbridge'  # Add 'saltbridge' interaction
    data['LIG_IDX_LIST'][11] = [7, 8]  # Ligand index list for 'saltbridge PI'
    data['LIG_GROUP'][11] = 'Group4'  # Ligand group for 'saltbridge PI'
    data['PROTISPOS'][11] = False  # PROTISPOS is False for 'saltbridge PI'
    data['RESTYPE_LIG'][11] = 0
    data['TARGET_IDX'][11] = 0

    # Add 'hydrophobic' case where 'ring_found' is False
    data['FRAME'][12] = 12  # Add a new 'FRAME' value
    data['Prot_partner'][12] = 'C'  # Add a new 'Prot_partner' value
    data['INTERACTION'][12] = 'hydrophobic'  # Add 'hydrophobic' interaction
    data['LIGCARBONIDX'][12] = 104  # LIGCARBONIDX for 'hydrophobic' (not in any ring)
    data['RESTYPE_LIG'][12] = 0 
    data['TARGET_IDX'][12] = 0

    return pd.DataFrame(data)


 
@pytest.fixture
def sample_dataframe_bindingmode_processing_with_peptides():
    data = {
        'FRAME': {0: 1, 1: 2, 2: 3, 3: 2},
        'Prot_partner': {0: '62VAL', 1: 'SER144', 2: 'GLU321', 3: 'ILE432'},
        'INTERACTION': {0: 'hydrophobic', 1: 'hbond', 2: 'saltbridge', 3: 'hydrophobic'},
        'LIGCARBONIDX': {0: 101, 1: 102, 2: 103, 3: 102},
        'DONORIDX': {0: 201, 1: 202, 2: 203, 3: 202},
        'ACCEPTORIDX': {0: 301, 1: 302, 2: 303, 3: 302},
        'PROTISDON': {0: True, 1: False, 2: True, 3: False},
        'LIG_IDX_LIST': {0: [1, 2], 1: [3, 4], 2: [5, 6], 3: [3, 4]},
        'LIG_GROUP': {0: 'Group1', 1: 'Group2', 2: 'Group3', 3: 'Group1'},
        'PROTISPOS': {0: True, 1: False, 2: True, 3: True},
        'DON_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONORTYPE': {0: 0, 1: 0, 2: 0, 3: 0},
        'ACCEPTOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'DONOR_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'LOCATION': {0: 0, 1: 0, 2: 0, 3: 0},
        'METAL_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'METAL_TYPE': {0: 0, 1: 0, 2: 0, 3: 0},
        'RESTYPE_LIG': {0: 'ILE', 1: 'TYR', 2: 'ARG', 3: 'VAL'},
        'TARGET_IDX': {0: 0, 1: 0, 2: 0, 3: 0},
        'COORDINATION': {0: 0, 1: 0, 2: 0, 3: 0},
        'RESNR_LIG' : {0: '101', 1: '202', 2: '155', 3: '102'}
    }
     
    # Additional data for peptide interactions
    data['FRAME'][4] = 4
    data['Prot_partner'][4] = 'LEU248'
    data['INTERACTION'][4] = 'halogen'
    data['DON_IDX'][4] = 501
    data['DONORTYPE'][4] = 'F'
    data['ACCEPTOR_IDX'][4] = 0
    data['DONOR_IDX'][4] = 0
    data['LIG_IDX_LIST'][4] = 0
    data['LIG_GROUP'][4] = 0
    data['RESTYPE_LIG'][4] = 'ILE'
    data['TARGET_IDX'][4] = 0
    data['RESNR_LIG'][4] = '501'
 
    data['FRAME'][5] = 5
    data['Prot_partner'][5] = 'SER300'
    data['INTERACTION'][5] = 'hbond'
    data['ACCEPTORIDX'][5] = 301
    data['DON_IDX'][5] = 0
    data['DONORTYPE'][5] = 0
    data['PROTISDON'][5] = True
    data['ACCEPTOR_IDX'][5] = 0
    data['LIG_IDX_LIST'][5] = 0
    data['DONOR_IDX'][5] = 0
    data['LIG_GROUP'][5] = 0
    data['RESTYPE_LIG'][5] = 'HIS'
    data['TARGET_IDX'][5] = 0
    data['RESNR_LIG'][5] = '301'
 
    data['FRAME'][6] = 6
    data['Prot_partner'][6] = 'TYR343'
    data['INTERACTION'][6] = 'waterbridge'
    data['ACCEPTOR_IDX'][6] = 401
    data['DON_IDX'][6] = 0
    data['DONORTYPE'][6] = 0
    data['DONOR_IDX'][6] = 0
    data['LIG_IDX_LIST'][6] = 0
    data['PROTISDON'][6] = True
    data['LIG_GROUP'][6] = 0
    data['RESTYPE_LIG'][6] = 'SER'
    data['TARGET_IDX'][6] = 0
    data['RESNR_LIG'][6] = '455'

    data['FRAME'][7] = 7
    data['Prot_partner'][7] = 'ILE178'
    data['INTERACTION'][7] = 'waterbridge'
    data['DONOR_IDX'][7] = 501
    data['DON_IDX'][7] = 0
    data['DONORTYPE'][7] = 0
    data['PROTISDON'][7] = False
    data['ACCEPTOR_IDX'][7] = 0
    data['LIG_IDX_LIST'][7] = 0
    data['LIG_GROUP'][7] = 0
    data['RESTYPE_LIG'][7] = 'TYR'
    data['TARGET_IDX'][7] = 0
    data['RESNR_LIG'][7] = '467'

    data['FRAME'][8] = 8
    data['Prot_partner'][8] = 'PHE344'
    data['INTERACTION'][8] = 'pistacking'
    data['LIG_IDX_LIST'][8] = [7, 8]
    data['LIG_GROUP'][8] = 0
    data['ACCEPTOR_IDX'][8] = 0
    data['DON_IDX'][8] = 0
    data['DONOR_IDX'][8] = 0
    data['PROTISDON'][8] = False
    data['DONORTYPE'][8] = 0
    data['RESTYPE_LIG'][8] = 'PHE'
    data['TARGET_IDX'][8] = 0
    data['RESNR_LIG'][8] = '398'
 
    data['FRAME'][9] = 9
    data['Prot_partner'][9] = 'PHE754'
    data['INTERACTION'][9] = 'pication'
    data['LIG_IDX_LIST'][9] = [9, 10]
    data['LIG_GROUP'][9] = 'B'
    data['ACCEPTOR_IDX'][9] = 0
    data['DON_IDX'][9] = 0
    data['PROTISDON'][9] = False
    data['DONOR_IDX'][9] = 0
    data['DONORTYPE'][9] = 0
    data['RESTYPE_LIG'][9] = 'ARG'
    data['TARGET_IDX'][9] = 0
    data['RESNR_LIG'][9] = '245'
     
    data['FRAME'][10] = 10
    data['Prot_partner'][10] = 'LYS567'
    data['INTERACTION'][10] = 'pication'
    data['LIG_IDX_LIST'][10] = [9, 10]
    data['LIG_GROUP'][10] = 'B'
    data['RESTYPE_LIG'][10] = 'PHE'
    data['TARGET_IDX'][10] = 501
    data['RESNR_LIG'][10] = '228'
 
    data['FRAME'][11] = 11
    data['Prot_partner'][11] = 'LYS567'
    data['INTERACTION'][11] = 'saltbridge'
    data['LIG_IDX_LIST'][11] = [7, 8]
    data['LIG_GROUP'][11] = 'Group4'
    data['PROTISPOS'][11] = False
    data['RESTYPE_LIG'][11] = 'GLU'
    data['TARGET_IDX'][11] = 0
    data['RESNR_LIG'][11] = '423'
 
    data['FRAME'][12] = 12
    data['Prot_partner'][12] = 'HEM144'
    data['INTERACTION'][12] = 'metal'
    data['METAL_IDX'][12] = 401  # METAL_IDX for 'metal'
    data['METAL_TYPE'][12] = 'Fe'  # Metal type
    data['LOCATION'][12] = 'site1'  # Location
    data['ACCEPTOR_IDX'][12] = 0
    data['DON_IDX'][12] = 0
    data['RESTYPE_LIG'][12] = 'HIS'  # Assuming 'A' as the RESTYPE_LIG for the metal interaction
    data['TARGET_IDX'][12] = 401
    data['COORDINATION'][12] = "site1"
    data['RESNR_LIG'][12] = '256'
 
    return pd.DataFrame(data)


def test_gather_interactions(sample_dataframe_bindingmode_processing):
    df = sample_dataframe_bindingmode_processing
    ligand_rings = [[101], [102], [103]]  # Define sample ligand rings for testing

    result = gather_interactions(df, ligand_rings)

    # Assert that the result is a dictionary
    
    assert isinstance(result, dict)

    # Check specific values in the generated dictionary for known interactions based on the updated fixture
    expected_result = {
    1: {0: 'A_101_hydrophobic'},
    2: {1: 'B_202_Donor_hbond', 3: 'A_102_hydrophobic'},
    3: {2: 'C_[5, 6]_Group3_NI_saltbridge'},
    4: {4: 'A_501_F_halogen'},
    5: {5: 'A_301_Acceptor_hbond'},
    6: {6: 'A_401_Acceptor_waterbridge'},
    7: {7: 'B_501_Donor_waterbridge'},
    8: {8: 'A_[7, 8]_pistacking'},
    9: {9: 'A_[9_ 10]_Group4_pication'},
    10: {10: 'A_401_Fe_site1_metal'},
    11: {11: 'A_[7, 8]_Group4_PI_saltbridge'},
    12: {12: 'C_104_hydrophobic'}
}
    # Check if the actual result matches the expected result
    assert result == expected_result

def test_gather_interactions_with_peptides(sample_dataframe_bindingmode_processing_with_peptides):
    df = sample_dataframe_bindingmode_processing_with_peptides
    ligand_rings = [[101], [102], [103]]  # Define sample ligand rings for testing

    result = gather_interactions(df, ligand_rings, peptide=True)
 
    # Assert that the result is a dictionary
    assert isinstance(result, dict)
 
    # Check specific values in the generated dictionary for known interactions based on the updated fixture
    expected_result = {
        1: {0: '62VAL_101ILE_hydrophobic'},
        2: {1: 'SER144_202TYR_Donor_hbond', 3: 'ILE432_102VAL_hydrophobic'},
        3: {2: 'GLU321_155ARG_ARG_NI_saltbridge'},
        4: {4: 'LEU248_501ILE_F_halogen'},
        5: {5: 'SER300_301HIS_Acceptor_hbond'},
        6: {6: 'TYR343_455SER_Acceptor_waterbridge'},
        7: {7: 'ILE178_467TYR_Donor_waterbridge'},
        8: {8: 'PHE344_398PHE_pistacking'},
        9: {9: 'PHE754_245ARG_ARG_pication'},
        10: {10: 'LYS567_228PHE_PHE_pication'},
        11: {11: 'LYS567_423GLU_GLU_PI_saltbridge'},
        12: {12: 'HIS_256HIS_Fe_site1_metal'},
    }
 
    # Check if the actual result matches the expected result
    assert result == expected_result



@pytest.fixture
def test_remove_duplicates_data():
    input_data = {
        'a': {'x': 1, 'y': 2, 'z': 1},
        'b': {'p': 3, 'q': 3, 'r': 4}
    }
    expected_output = {
        'a': {'x': 1, 'y': 2},
        'b': {'p': 3, 'r': 4}
    }
    return input_data, expected_output

def test_unique_data_generation():
    # Test case 1: Check if the function returns an empty dictionary for an empty list
    result = unique_data_generation([])
    assert result == {}

    # Test case 2: Check if the function correctly identifies and stores unique values
    input_list = [1, 2, 2, 3, 3, 4, 5, 5]
    result = unique_data_generation(input_list)
    expected_result = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}
    assert result == expected_result

    # Test case 3: Check if the function handles strings
    input_list = ["apple", "banana", "apple", "cherry"]
    result = unique_data_generation(input_list)
    expected_result = {"apple": "apple", "banana": "banana", "cherry": "cherry"}
    assert result == expected_result

    
# Define a test case that uses the fixture
def test_remove_duplicate_values(test_remove_duplicates_data):
    input_data, expected_output = test_remove_duplicates_data
    assert remove_duplicate_values(input_data) == expected_output

def test_combine_subdict_values():
    # Test case 1: Empty input dictionary
    data = {}
    result = combine_subdict_values(data)
    assert result == {'all': []}

    # Test case 2: Input dictionary with sub-dictionaries
    data = {
        'dict1': {'a': 1, 'b': 2},
        'dict2': {'c': 3, 'd': 4},
        'dict3': {'e': 5, 'f': 6},
    }
    result = combine_subdict_values(data)
    assert result == {'all': [1, 2, 3, 4, 5, 6]}

    # Test case 3: Input dictionary with empty sub-dictionaries
    data = {
        'dict1': {},
        'dict2': {},
    }
    result = combine_subdict_values(data)
    assert result == {'all': []}

    # Test case 4: Input dictionary with sub-dictionaries containing various data types
    data = {
        'dict1': {'a': 1, 'b': 'text', 'c': [1, 2, 3]},
        'dict2': {'d': None, 'e': 5.5},
    }
    result = combine_subdict_values(data)
    assert result == {'all': [1, 'text', [1, 2, 3], None, 5.5]}

# Define a sample DataFrame for testing
sample_data = {
    'A': [1, 2, 3, 4, 5],
    'B': [2, 3, 4, 5, 6],
    'C': [3, 4, 5, 6, 7],
    'D': [4, 5, 6, 7, 8],
}
sample_df = pd.DataFrame(sample_data)

# Define the provided 'unique_columns_rings_grouped' data for testing
unique_columns_rings_grouped = {
    1: {0: 'A_101_hydrophobic'},
    2: {1: 'B_202_Donor_hbond', 3: 'A_102_hydrophobic'},
    3: {2: 'C_[5, 6]_Group3_NI_saltbridge'},
    4: {4: 'A_501_F_halogen'},
    5: {5: 'A_301_Acceptor_hbond'},
    6: {6: 'A_401_Acceptor_waterbridge'},
    7: {7: 'B_501_Donor_waterbridge'},
    8: {8: 'A_[7, 8]_pistacking'},
    9: {9: 'A_[9_ 10]_Group4_pication'},
    10: {10: 'A_401_Fe_site1_metal'},
    11: {11: 'A_[7, 8]_Group4_PI_saltbridge'},
    12: {12: 'C_104_hydrophobic'}
}

def test_filtering_values_with_provided_data():
    # Test case 1: Check if the function returns a list of filtered values
    threshold = 0.2  # 20% threshold
    frames = 1000  # Some arbitrary number of frames
    df = pd.DataFrame()  # Create an empty DataFrame for testing
    result = filtering_values(threshold, frames, df, unique_columns_rings_grouped)

    assert isinstance(result, list)
    

def test_df_iteration_numbering():
    # Sample DataFrame for testing
    data = {
        'Unnamed: 0': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
        'RESNR': [98, 63, 162, 161, 166, 165, 125, 166, 211, 227, 223, 165, 100, 59, 98, 207, 164],
        'RESTYPE': ['PHE', 'ARG', 'ALA', 'PHE', 'ARG', 'ASP', 'TYR', 'ARG', 'PHE', 'LEU', 'THR', 'ASP', 'ASP', 'ARG', 'PHE', 'PHE', 'LYS'],
        'RESCHAIN': ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'],
        'RESNR_LIG': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'RESTYPE_LIG': ['UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK', 'UNK'],
        'RESCHAIN_LIG': ['X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X'],
        'DIST': [3.46, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 3.36, 3.61, 3.84, 3.62, 3.72, 3.62, 3.99, 3.65, 3.70, 5.16],
        'LIGCARBONIDX': [4196.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4206.0, 4207.0, 4207.0, 4215.0, 4217.0, 4217.0, 4194.0, 4208.0, 0.0],
        '162ALAA_4214_4215_4216_4217_4218_4213_hydrophobic': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '166ARGA_4220,4221_Carboxylate_NI_saltbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '98PHEA_4194_hydrophobic': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '63ARGA_4201_Acceptor_waterbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '164LYSA_4213_4214_4215_4216_4217_4218_Aromatic_pication': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '166ARGA_4220_Acceptor_hbond': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '98PHEA_4225_Donor_hbond': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '207PHEA_4213,4214,4215,4216,4217,4218_pistacking': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '100ASPA_4005_Donor_waterbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        '59ARGA_4222_Acceptor_waterbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    }

    df = pd.DataFrame(data)

    interactions = [
    'hbond',
    'waterbridge',
    'hbond',
    'hbond',
    'hbond',
    'hbond',
    'hbond',
    'saltbridge',
    'hydrophobic',
    'hydrophobic',
    'hydrophobic',
    'hydrophobic',
    'waterbridge',
    'waterbridge',
    'hydrophobic',
    'pistacking',
    'pication'
]
    df['INTERACTION'] = interactions

    
    # Define the values for the "PROTISDON" column
    protisdon_values = [False, True, True, True, True, True, True, 0, 0, 0, 0, 0, False, True, 0, 0, 0]

    # Update the "PROTISDON" column in the DataFrame
    df['PROTISDON'] = protisdon_values

    # Define the values for the "Prot_partner" column
    prot_partner_values = ['98PHEA', '63ARGA', '162ALAA', '161PHEA', '166ARGA', '165ASPA', '125TYRA', '166ARGA', '211PHEA', '227LEUA', '223THRA', '165ASPA', '100ASPA', '59ARGA', '98PHEA', '207PHEA', '164LYSA']

    # Update the "Prot_partner" column in the DataFrame
    df['Prot_partner'] = prot_partner_values

    # Define the values for the "ACCEPTORIDX" column
    acceptoridx_values = [0.0, 0.0, 4221.0, 4221.0, 4220.0, 4220.0, 4192.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Update the "ACCEPTORIDX" column in the DataFrame
    df['ACCEPTORIDX'] = acceptoridx_values

    # Define the values for the "DONORIDX" column
    donoridx_values = [4225.0, 0.0, 2417.0, 2397.0, 2468.0, 2456.0, 1828.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    # Update the "ACCEPTORIDX" column in the DataFrame
    df['DONORIDX'] = donoridx_values

    # Define the values for the "ACCEPTOR_IDX" column
    acceptor_idx_values = [0.0, 4201.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4222.0, 0.0, 0.0, 0.0]

    # Add the "ACCEPTOR_IDX" column to the DataFrame
    df['ACCEPTOR_IDX'] = acceptor_idx_values

    # Define the values for the "DONOR_IDX" column
    donor_idx_values = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4005.0, 0.0, 0.0, 0.0, 0.0]

    # Add the "ACCEPTOR_IDX" column to the DataFrame
    df['DONOR_IDX'] = donor_idx_values
    
    # Define the values for the "LIG_IDX_LIST" column
    lig_idx_list_values = [0, 0, 0, 0, 0, 0, 0, "4220,4221", 0, 0, 0, 0, 0, 0, 0, "4213,4214,4215,4216,4217,4218", "4213,4214,4215,4216,4217,4218"]

    # Add the "LIG_IDX_LIST" column to the DataFrame
    df['LIG_IDX_LIST'] = lig_idx_list_values

    # Define the values for the "LIG_GROUP" column
    lig_group_values = [0, 0, 0, 0, 0, 0, 0, "Carboxylate", 0, 0, 0, 0, 0, 0, 0, "Aromatic", "Aromatic"]

    # Add the "LIG_GROUP" column to the DataFrame
    df['LIG_GROUP'] = lig_group_values
    
    # Updated unique_data dictionary
    unique_data = {
        '63ARGA_4201_Acceptor_waterbridge': '63ARGA_4201_Acceptor_waterbridge',
        '166ARGA_4220_Acceptor_hbond': '166ARGA_4220_Acceptor_hbond',
        '166ARGA_4220,4221_Carboxylate_NI_saltbridge': '166ARGA_4220,4221_Carboxylate_NI_saltbridge',
        '162ALAA_4214_4215_4216_4217_4218_4213_hydrophobic': '162ALAA_4214_4215_4216_4217_4218_4213_hydrophobic',
        '98PHEA_4194_hydrophobic': '98PHEA_4194_hydrophobic',
        '98PHEA_4225_Donor_hbond': '98PHEA_4225_Donor_hbond',
        '164LYSA_4213_4214_4215_4216_4217_4218_Aromatic_pication': '164LYSA_4213_4214_4215_4216_4217_4218_Aromatic_pication',
        '207PHEA_4213,4214,4215,4216,4217,4218_pistacking': '207PHEA_4213,4214,4215,4216,4217,4218_pistacking',
        '59ARGA_4222_Acceptor_waterbridge': '59ARGA_4222_Acceptor_waterbridge',
        '100ASPA_4005_Donor_waterbridge': '100ASPA_4005_Donor_waterbridge'
    }


    # Call the function with the sample DataFrame and unique_data
    df_iteration_numbering(df, unique_data)

    expected_98PHEA_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
    assert (df['98PHEA_4194_hydrophobic'] == expected_98PHEA_values).all()

    expected_166ARGA_4220_Acceptor_hbond_values = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (df['166ARGA_4220_Acceptor_hbond'] == expected_166ARGA_4220_Acceptor_hbond_values).all()

    expected_63ARGA_4201_Acceptor_waterbridge_values = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (df['63ARGA_4201_Acceptor_waterbridge'] == expected_63ARGA_4201_Acceptor_waterbridge_values).all()

    expected_98PHEA_4225_Donor_hbond_values = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (df['98PHEA_4225_Donor_hbond'] == expected_98PHEA_4225_Donor_hbond_values).all()

    expected_164LYSA_4213_4214_4215_4216_4217_4218_Aromatic_pication_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    assert (df['164LYSA_4213_4214_4215_4216_4217_4218_Aromatic_pication'] == expected_164LYSA_4213_4214_4215_4216_4217_4218_Aromatic_pication_values).all()

    expected_207PHEA_4213_4214_4215_4216_4217_4218_pistacking_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    assert (df['207PHEA_4213,4214,4215,4216,4217,4218_pistacking'] == expected_207PHEA_4213_4214_4215_4216_4217_4218_pistacking_values).all()

    expected_59ARGA_4222_Acceptor_waterbridge_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    assert (df['59ARGA_4222_Acceptor_waterbridge'] == expected_59ARGA_4222_Acceptor_waterbridge_values).all()

    expected_100ASPA_4005_Donor_waterbridge_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    assert (df['100ASPA_4005_Donor_waterbridge'] == expected_100ASPA_4005_Donor_waterbridge_values).all()


@pytest.fixture
def sample_dataframe_it_peptides():
    # Create a sample DataFrame for testing
    data = {
        'Unnamed: 0': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
        'Prot_partner': ['62VAL', 'SER144', 'GLU321', 'ILE432', 'LEU248', 'SER300', 'TYR343', 'ILE178', 'PHE344', 'PHE754', 'LYS567', 'LYS567', 'HIS'],
        'LIGCARBONIDX': [101, 202, 155, 102, 501, 301, 467, 467, 398, 245, 228, 423, 256],
        'INTERACTION': ['hydrophobic', 'hbond', 'saltbridge', 'hydrophobic', 'halogen', 'hbond', 'waterbridge', 'waterbridge', 'pistacking', 'pication', 'pication', 'saltbridge', 'metal'],
        'PROTISDON': [None, True, None, None, None, True, True, False, None, None, None, False, None],
        'ACCEPTORIDX': [None, 202, None, None, None, 301, None, None, None, None, None, None, None],
        'RESNR_LIG': [101, 202, 155, 102, 501, 301, 455, 467, 398, 245, 228, 423, 256],
        'DONORIDX': [None, None, None, None, None, None, None, None, None, None, None, None, None],
        'DONORTYPE': [None, None, None, None, 'F', None, None, None, None, None, None, None, None],
        'LIG_IDX_LIST': [None, None, None, None, None, None, None, None, None, None, None, None, None],
        'RESTYPE_LIG': ['VAL', 'TYR', 'ARG', 'VAL', 'ILE', 'HIS', 'SER', 'TYR', 'PHE', 'ARG', 'PHE', 'GLU', 'HIS'],
        'TARGET_IDX': [None, None, None, None, None, None, None, None, None, None, None, None, None],
        'METAL_TYPE': [None, None, None, None, None, None, None, None, None, None, None, None, 'Fe'],
        '62VAL_101ILE_hydrophobic': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'SER144_202TYR_Donor_hbond': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'GLU321_155ARG_ARG_NI_saltbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'ILE432_102VAL_hydrophobic': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'LEU248_501ILE_F_halogen': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'SER300_301HIS_Acceptor_hbond': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'TYR343_455SER_Acceptor_waterbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'ILE178_467TYR_Donor_waterbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'PHE344_398PHE_pistacking': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'PHE754_245ARG_ARG_pication': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'LYS567_228PHE_PHE_pication': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'LYS567_423GLU_GLU_PI_saltbridge': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        'HIS_256HIS_Fe_site1_metal': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    }
    
    df = pd.DataFrame(data)
    return df

def test_df_iteration_numbering_with_peptide(sample_dataframe_it_peptides):
    unique_data = {
        0: '62VAL_101ILE_hydrophobic',
        1: 'SER144_202TYR_Donor_hbond',
        2: 'GLU321_155ARG_ARG_NI_saltbridge',
        3: 'ILE432_102VAL_hydrophobic',
        4: 'LEU248_501ILE_F_halogen',
        5: 'SER300_301HIS_Acceptor_hbond',
        6: 'TYR343_455SER_Acceptor_waterbridge',
        7: 'ILE178_467TYR_Donor_waterbridge',
        8: 'PHE344_398PHE_pistacking',
        9: 'PHE754_245ARG_ARG_pication',
        10: 'LYS567_228PHE_PHE_pication',
        11: 'LYS567_423GLU_GLU_PI_saltbridge',
        12: 'HIS_256HIS_Fe_site1_metal',
    }

    df_iteration_numbering(sample_dataframe_it_peptides, unique_data, peptide=True)

    # Assertions similar to the provided ones
    expected_101ILE_hydrophobic_values = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['62VAL_101ILE_hydrophobic'] == expected_101ILE_hydrophobic_values).all()

    expected_202TYR_Donor_hbond_values = [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['SER144_202TYR_Donor_hbond'] == expected_202TYR_Donor_hbond_values).all()

    expected_155ARG_ARG_NI_saltbridge_values = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['GLU321_155ARG_ARG_NI_saltbridge'] == expected_155ARG_ARG_NI_saltbridge_values).all()

    expected_102VAL_hydrophobic_values = [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['ILE432_102VAL_hydrophobic'] == expected_102VAL_hydrophobic_values).all()

    expected_501ILE_halogen_values = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['LEU248_501ILE_F_halogen'] == expected_501ILE_halogen_values).all()

    expected_301HIS_Acceptor_hbond_values = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['SER300_301HIS_Acceptor_hbond'] == expected_301HIS_Acceptor_hbond_values).all()

    expected_455SER_Acceptor_waterbridge_values = [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['TYR343_455SER_Acceptor_waterbridge'] == expected_455SER_Acceptor_waterbridge_values).all()

    expected_467TYR_Donor_waterbridge_values = [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['ILE178_467TYR_Donor_waterbridge'] == expected_467TYR_Donor_waterbridge_values).all()

    expected_398PHE_pistacking_values = [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0]
    assert (sample_dataframe_it_peptides['PHE344_398PHE_pistacking'] == expected_398PHE_pistacking_values).all()

    expected_245ARG_ARG_pication_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]
    assert (sample_dataframe_it_peptides['PHE754_245ARG_ARG_pication'] == expected_245ARG_ARG_pication_values).all()

    expected_228PHE_PHE_pication_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
    assert (sample_dataframe_it_peptides['LYS567_228PHE_PHE_pication'] == expected_228PHE_PHE_pication_values).all()

    expected_423GLU_GLU_PI_saltbridge_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    assert (sample_dataframe_it_peptides['LYS567_423GLU_GLU_PI_saltbridge'] == expected_423GLU_GLU_PI_saltbridge_values).all()

    expected_256HIS_Fe_site1_metal_values = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    assert (sample_dataframe_it_peptides['HIS_256HIS_Fe_site1_metal'] == expected_256HIS_Fe_site1_metal_values).all()

@pytest.fixture
def sample_data():
    # Create sample data for testing
    df = pd.DataFrame({
        'FRAME': [1, 2, 3],
        'Column1': [10, 20, 30],
        'Column2': [40, 50, 60]
    })

    new_df = pd.DataFrame({
        'FRAME': [1, 2, 3],
        'Column1': [100, 200, 300],
        'Column2': [400, 500, 600]
    })

    unique_data = {'Column1': 'Column1', 'Column2': 'Column2'}

    return df, new_df, unique_data

def test_update_values(sample_data):
    # Arrange
    df, new_df, unique_data = sample_data

    # Set 'FRAME' as the index for new_df
    new_df = new_df.set_index('FRAME')

    # Act
    update_values(df, new_df, unique_data)

    # Assert
    expected_df = pd.DataFrame({
        'FRAME': [1, 2, 3],
        'Column1': [100, 200, 300],
        'Column2': [400, 500, 600]
    })

    # Check if the specific values are updated
    assert df[['Column1', 'Column2']].equals(expected_df[['Column1', 'Column2']])

    # Check if the specific values are updated
    assert df[['Column1', 'Column2']].equals(expected_df[['Column1', 'Column2']])

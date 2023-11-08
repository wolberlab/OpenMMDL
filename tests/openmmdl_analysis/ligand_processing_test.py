import pytest
import openmmdl
from pathlib import Path
from openmmdl.openmmdl_analysis.ligand_processing import increase_ring_indices, convert_ligand_to_smiles
import os


test_data_directory = Path("openmmdl/tests/data/in")
TEST_LIGAND_FILE = f"{test_data_directory}/CVV.sdf"
TEST_OUTPUT_FILE = "CVV.smi"

def test_increase_ring_indices():
    # Test case 1: Check if ring indices are correctly increased
    ring = [1, 2, 3]
    lig_index = 10
    result = increase_ring_indices(ring, lig_index)
    assert result == [11, 12, 13]

    # Test case 2: Check with a different lig_index
    ring = [3, 4, 5]
    lig_index = 20
    result = increase_ring_indices(ring, lig_index)
    assert result == [23, 24, 25]


def test_convert_ligand_to_smiles():
    # Convert the ligand structure to SMILES in the same directory as the input SDF file
    convert_ligand_to_smiles(TEST_LIGAND_FILE, TEST_OUTPUT_FILE)

    # Verify that the output SMILES file was created in the same directory as the input file
    assert os.path.exists(TEST_OUTPUT_FILE)

    # Optionally, you can also read and validate the content of the output SMILES file
    with open(TEST_OUTPUT_FILE, "r") as smi_file:
        smiles_lines = smi_file.readlines()
        assert len(smiles_lines) > 0  # Check that there are SMILES representations



if __name__ == '__main__':
    pytest.main()

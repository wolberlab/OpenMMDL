import os
import shutil
import pytest
from pathlib import Path
from unittest.mock import mock_open, patch
from openmmdl.openmmdl_simulation.scripts.cleaning_procedures import cleanup, create_directory_if_not_exists, post_md_file_movement, copy_file, create_directory_if_not_exists, organize_files

@pytest.fixture
def test_protein_name():
    return "test_protein"

@pytest.fixture
def test_directory_path():
    return "test_directory"

def test_cleanup(test_protein_name):
    # Create a dummy file to be removed
    with open(f'output_{test_protein_name}', 'w') as dummy_file:
        dummy_file.write("Dummy content")

    # Call the cleanup function
    cleanup(test_protein_name)

    # Check if the file has been removed
    assert not os.path.exists(f'output_{test_protein_name}')

def test_create_directory_if_not_exists(test_directory_path):
    # Create a test directory
    create_directory_if_not_exists(test_directory_path)

    # Check if the directory exists
    assert os.path.exists(test_directory_path)

    # Call the function again, it should not raise an error
    create_directory_if_not_exists(test_directory_path)

    # Cleanup: Remove the test directory
    shutil.rmtree(test_directory_path)
    assert not os.path.exists(test_directory_path)


@patch("os.path.exists")
@patch("shutil.copy")
def test_copy_file(mock_copy, mock_exists):
    
    src = "source_file.txt"
    dest = "destination_directory"

    # Mock the os.path.exists to return True, indicating the source file exists
    mock_exists.return_value = True

    # Call the copy_file function
    copy_file(src, dest)

    # Check that os.path.exists was called with the source file
    mock_exists.assert_called_with(src)

    # Check that shutil.copy was called with the source file and destination directory
    mock_copy.assert_called_with(src, dest)

# Mock the os.path.exists and os.rename functions
@patch("os.path.exists")
@patch("os.rename")
def test_organize_files(mock_rename, mock_exists):
    source = ["file1.txt", "file2.txt", "file3.txt"]
    destination = "destination_directory"

    # Mock os.path.exists to return True for all source files
    mock_exists.side_effect = [True] * len(source)

    # Call the organize_files function
    organize_files(source, destination)

    # Print the calls made to os.rename
    for call in mock_rename.call_args_list:
        print(call)

def test_post_md_file_movement():
    # Get the absolute path to the test data directory
    test_data_directory = Path("openmmdl/tests/data/in")

    # Define the full path to the input files
    ligand = test_data_directory / 'CVV.sdf'
    protein_name = test_data_directory / '6b73.pdb'
    prmtop = test_data_directory / '6b73.prmtop'
    inpcrd = test_data_directory / '6b73.inpcrd'
    
    # Assert that the input files exist before moving
    assert os.path.exists(ligand)
    assert os.path.exists(protein_name)
    assert os.path.exists(prmtop)
    assert os.path.exists(inpcrd)

    # Call the post_md_file_movement function
    post_md_file_movement(protein_name, prmtop, inpcrd, ligand)
    
    # Check if the files have been organized and moved to the correct directories
    input_files_dir = Path("Input_Files")

    assert os.path.exists(input_files_dir)
    assert os.path.exists(input_files_dir / "6b73.pdb")
    assert os.path.exists(input_files_dir / "6b73.prmtop")
    assert os.path.exists(input_files_dir / "6b73.inpcrd")
    assert os.path.exists(input_files_dir / "CVV.sdf")

# Run the tests
if __name__ == "__main__":
    pytest.main()

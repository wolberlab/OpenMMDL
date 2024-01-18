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

#def test_post_md_file_movement():
#    # Get the absolute path to the test data directory
#    test_data_directory = Path("openmmdl/tests/data/in")
#
#    # Define the full path to the input files
#    ligand = test_data_directory / 'CVV.sdf'
#    protein_name = test_data_directory / '6b73.pdb'
#    prmtop = test_data_directory / '6b73.prmtop'
#    inpcrd = test_data_directory / '6b73.inpcrd'
#    
#    # Assert that the input files exist before moving
#    assert os.path.exists(ligand)
#    assert os.path.exists(protein_name)
#    assert os.path.exists(prmtop)
#    assert os.path.exists(inpcrd)
#
#    # Call the post_md_file_movement function
#    post_md_file_movement(protein_name, prmtop, inpcrd, ligand)
#    
#    # Check if the files have been organized and moved to the correct directories
#    input_files_dir = Path("Input_Files")
#
#    assert os.path.exists(input_files_dir)
#    assert os.path.exists(input_files_dir / "6b73.pdb")
#    assert os.path.exists(input_files_dir / "6b73.prmtop")
#    assert os.path.exists(input_files_dir / "6b73.inpcrd")
#    assert os.path.exists(input_files_dir / "CVV.sdf")

def test_post_md_file_movement():
    # Get the absolute path to the test data directory
    test_data_directory = Path("openmmdl/tests/data/in")

    # Define the full path to the input files
    ligand = test_data_directory / 'CVV.sdf'
    protein_name = test_data_directory / '6b73.pdb'
    prmtop = test_data_directory / '6b73.prmtop'
    inpcrd = test_data_directory / '6b73.inpcrd'
    protein_no_solvent = test_data_directory / 'prepared_no_solvent_6b73.pdb'
    protein_solvent = test_data_directory / 'solvent_padding_6b73.pdb'
    protein_equilibration = test_data_directory / 'Equilibration_6b73.pdb'
    protein_minimization = test_data_directory / 'Energyminimization_6b73.pdb'
    output_pdb = test_data_directory / 'output_6b73.pdb'
    mdtraj_top = test_data_directory / 'centered_old_coordinates_top.pdb'
    prot_lig_top = test_data_directory / 'prot_lig_top.pdb'
    checkpoint = test_data_directory / 'checkpoint.chk'
    checkpoint_10x = test_data_directory / '10x_checkpoint.chk'
    
    # Assert that the input files exist before moving
    assert os.path.exists(ligand)
    assert os.path.exists(protein_name)
    assert os.path.exists(prmtop)
    assert os.path.exists(inpcrd)
    assert os.path.exists(protein_no_solvent)

    shutil.copy(str(protein_no_solvent), '.')
    shutil.copy(str(protein_solvent), '.')
    shutil.copy(str(protein_equilibration), '.')
    shutil.copy(str(protein_minimization), '.')
    shutil.copy(str(output_pdb), '.')
    shutil.copy(str(mdtraj_top), '.')
    shutil.copy(str(prot_lig_top), '.')
    shutil.copy(str(checkpoint), '.')
    shutil.copy(str(checkpoint_10x), '.')
    shutil.copy(str(protein_name), '.')
    protein_name = '6b73.pdb'

    # Call the post_md_file_movement function
    post_md_file_movement(str(protein_name), str(prmtop), str(inpcrd), [str(ligand)])
    
    # Check if the files have been organized and moved to the correct directories
    input_files_dir = Path("Input_Files")
    md_files_dir = Path("MD_Files")
    md_postprocessing_dir = Path("MD_Postprocessing")
    final_output_dir = Path("Final_Output")
    checkpoints_dir = Path("Checkpoints")

    assert os.path.exists(input_files_dir)
    assert os.path.exists(md_files_dir / "Pre_MD")
    assert os.path.exists(md_files_dir / "Pre_MD" / "prepared_no_solvent_6b73.pdb")
    assert os.path.exists(md_files_dir / "Pre_MD" / "solvent_padding_6b73.pdb")
    assert os.path.exists(md_files_dir / "Minimization_Equilibration" / "Equilibration_6b73.pdb")
    assert os.path.exists(md_files_dir / "Minimization_Equilibration" / "Energyminimization_6b73.pdb")
    assert os.path.exists(md_files_dir / "MD_Output" / "output_6b73.pdb")
    assert os.path.exists(md_postprocessing_dir / "centered_old_coordinates_top.pdb")
    assert os.path.exists(final_output_dir / "Prot_Lig" / "prot_lig_top.pdb")
    assert os.path.exists(checkpoints_dir / "checkpoint.chk")
    assert os.path.exists(checkpoints_dir / "10x_checkpoint.chk")

# Run the tests
if __name__ == "__main__":
    pytest.main()

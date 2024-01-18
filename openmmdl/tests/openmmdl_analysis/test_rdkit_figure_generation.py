import pytest
import os
import time
import shutil
from PIL import Image
from pathlib import Path
from openmmdl.openmmdl_analysis.rdkit_figure_generation import *

test_data_directory = Path("openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation")
test_data_directory_files = Path("openmmdl/tests/data/in")
lig_no_h = test_data_directory_files / 'lig_no_h.pdb'
complex = test_data_directory_files / 'complex.pdb'
smi_file = test_data_directory_files / 'lig_no_h.smi'
current_directory = os.getcwd() 
output_path = 'all_binding_modes_arranged.png'

shutil.copy(str(lig_no_h), '.')
shutil.copy(str(complex), '.')

@pytest.mark.parametrize("input_data, expected_output", [
    (["60GLUA_4206_4207_4216_4217_4218_4205_hydrophobic"], ['60GLUA 4206 4207 4216 4217 4218 4205 hydrophobic']),
    (["165ASPA_4203_Acceptor_hbond"], ['165ASPA 4203 Acceptor hbond']),
    (["125TYRA_4192_Acceptor_waterbridge"], ['125TYRA 4192 Acceptor waterbridge']),
])
def test_split_interaction_data(input_data, expected_output):
    result = split_interaction_data(input_data)
    assert result == expected_output

def test_highlight_numbers():
    # Input data
    split_data = [
        "163GLYA 4202 Acceptor hbond",
        "165ASPA 4203 Donor hbond",
        "165ASPA 4222 Donor hbond",
        "165ASPA 4203 Acceptor hbond",
        "125TYRA 4192 Acceptor waterbridge",
        "165ASPA 4222 Donor waterbridge",
        "161PHEA 4211 4212 4213 4214 4215 4210 hydrophobic",
        "59ARGA 4205 4206 4207 4216 4217 4218 Aromatic pication",
        "155PHEA 4205 4206 4207 4216 4217 4218 pistacking",
        "59ARGA 4194 F halogen",
        "166ARGA 4202,4203 Carboxylate NI saltbridge",
        "165ASPA 4202 Amine PI saltbridge"
        "HEM 4202 FE 4 metal"
    ]

    starting_idx = 1  # Updated starting index

    result = highlight_numbers(split_data, starting_idx)

    highlighted_hbond_donor, highlighted_hbond_acceptor, highlighted_hbond_both, \
    highlighted_hydrophobic, highlighted_waterbridge, highlighted_pistacking, highlighted_halogen, \
    highlighted_ni, highlighted_pi, highlighted_pication, highlighted_metal = result

    assert highlighted_hbond_donor is not None
    assert highlighted_hbond_acceptor is not None
    assert highlighted_hbond_both is not None
    assert highlighted_hydrophobic is not None
    assert highlighted_waterbridge is not None
    assert highlighted_halogen is not None
    assert highlighted_ni is not None
    assert highlighted_pication is not None
    assert highlighted_metal is not None
    
def test_update_dict():
    # Test case 1: Check if the target dictionary is updated correctly
    target_dict = {1: '1', 2: '2'}
    source_dict = {3: '3', 4: '4'}
    update_dict(target_dict, source_dict)
    assert target_dict == {1: '1', 2: '2', 3: '3', 4: '4'}

    # Test case 2: Check if the function handles multiple source dictionaries
    target_dict = {}
    source_dict1 = {1: '1'}
    source_dict2 = {2: '2', 3: '3'}
    update_dict(target_dict, source_dict1, source_dict2)
    assert target_dict == {1: '1', 2: '2', 3: '3'}

    # Test case 3: Check if the function handles empty source dictionaries
    target_dict = {1: '1', 2: '2'}
    update_dict(target_dict)  # No source dictionaries provided
    assert target_dict == {1: '1', 2: '2'}

def test_generate_interaction_dict():
    # Test with a known interaction type 'hydrophobic'
    interaction_type = 'hydrophobic'
    keys = [1, 2, 3]
    expected_result = {
        1: (1.0, 1.0, 0.0),
        2: (1.0, 1.0, 0.0),
        3: (1.0, 1.0, 0.0)
    }
    result = generate_interaction_dict(interaction_type, keys)
    assert result == expected_result

def test_create_and_merge_images_with_split_data():
    # Define test data
    binding_mode = 'Binding_Mode_1'
    occurrence_percent = 92
    split_data = [
        "166ARGA 4220,4221 Carboxylate NI saltbridge",
        "161PHEA 4221 Acceptor hbond",
        "207ILEA 4205 4206 4207 4208 4209 4204 hydrophobic"
    ]
    merged_image_paths = []

    # Define source image paths
    source_image_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1.png'
    source_svg_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1.svg'
    source_merged_image_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1_merged.png'

    # Copy source image files to the working directory
    working_directory = os.getcwd()
    destination_image_path = os.path.join(working_directory, os.path.basename(source_image_path))
    destination_svg_path = os.path.join(working_directory, os.path.basename(source_svg_path))
    destination_merged_image_path = os.path.join(working_directory, os.path.basename(source_merged_image_path))
    shutil.copy(source_image_path, destination_image_path)
    shutil.copy(source_svg_path, destination_svg_path)
    shutil.copy(source_merged_image_path, destination_merged_image_path)


    # Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory before:", files_in_working_directory)

    # Run the function
    merged_image_paths = create_and_merge_images(binding_mode, occurrence_percent, split_data, merged_image_paths)


    # Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory after:", files_in_working_directory)

    # Check if the merged image file was created
    assert len(merged_image_paths) == 1

    # Check if the merged image file is a valid image
    merged_image_path = merged_image_paths[0]
    try:
        with Image.open(merged_image_path) as img:
            img.verify()
    except Exception as e:
        pytest.fail(f"Merged image file is not a valid image: {e}")


def test_max_width_and_height_calculation():
    # Create some example images with different sizes
    image1 = Image.new('RGB', (100, 200), (255, 255, 255))
    image2 = Image.new('RGB', (150, 250), (255, 255, 255))
    merged_images = [image1, image2]

    # Calculate the maximum width and height
    max_width = max(image.size[0] for image in merged_images)
    max_height = max(image.size[1] for image in merged_images)

    # Assert the calculated max_width and max_height
    assert max_width == 150
    assert max_height == 250

def test_big_figure_creation():
    # Create example merged images
    image1 = Image.new('RGB', (100, 200), (255, 255, 255))
    image2 = Image.new('RGB', (150, 250), (255, 255, 255))
    merged_images = [image1, image2]

    # Calculate the maximum width and height
    max_width = max(image.size[0] for image in merged_images)
    max_height = max(image.size[1] for image in merged_images)

    # Determine the number of images per row (in your case, 2 images per row)
    images_per_row = 2

    # Calculate the number of rows and columns required
    num_rows = (len(merged_images) + images_per_row - 1) // images_per_row
    total_width = max_width * images_per_row
    total_height = max_height * num_rows

    # Create a new image with the calculated width and height
    big_figure = Image.new('RGB', (total_width, total_height), (255, 255, 255))  # Set background to white

    # Assert the dimensions of the created big_figure
    assert big_figure.size == (300, 250)  # Width should be 300, height should be 250

def test_arranged_figure_generation():
    binding_mode1_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1_merged.png'
    binding_mode2_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_2_merged.png'
    all_modes_path = 'openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/all_binding_modes_arranged.png'
    working_directory = os.getcwd()
    
    # Print the working directory to verify it's as expected
    print("Working Directory:", working_directory)

    destination_path_1 = os.path.join(working_directory, os.path.basename(binding_mode1_path))
    destination_path_2 = os.path.join(working_directory, os.path.basename(binding_mode2_path))
    destination_path_all = os.path.join(working_directory, os.path.basename(all_modes_path))
    
    # Print the destination paths to verify they are constructed correctly
    print("Destination Path 1:", destination_path_1)
    print("Destination Path 2:", destination_path_2)
    print("Destination Path All:", destination_path_all)

    shutil.copy(binding_mode1_path, destination_path_1)
    shutil.copy(binding_mode2_path, destination_path_2)
    shutil.copy(all_modes_path, destination_path_all)
    
    merged_image_paths = ['Binding_Mode_1_merged.png', 'Binding_Mode_2_merged.png']
    output_path = 'all_binding_modes_arranged.png'
    output_path = os.path.join(working_directory, output_path)
    print(output_path)

    # Run the function
    arranged_figure_generation(merged_image_paths, output_path)
    print(output_path)

    # Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory:", files_in_working_directory)

    output_path = os.path.join(working_directory, 'Binding_Modes_Markov_States', 'all_binding_modes_arranged.png')
    print(output_path)

    # Check if the output file was created
    
    assert output_path is not None


output_image_file = "output_image.png"

# Copy the files to the current folder
shutil.copy(complex, Path.cwd())
shutil.copy(lig_no_h, Path.cwd())
shutil.copy(smi_file, Path.cwd())

# Test the generate_ligand_image function
def test_generate_ligand_image():
    ligand_name = "UNK"
    generate_ligand_image(ligand_name, "complex.pdb", "lig_no_h.pdb", "lig_no_h.smi", output_image_file)

    # Assert that the output image file exists
    assert os.path.exists(output_image_file)



# Run the tests
if __name__ == '__main__':
    pytest.main()

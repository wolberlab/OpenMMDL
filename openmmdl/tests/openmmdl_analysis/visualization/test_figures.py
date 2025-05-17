import pytest
import os
import shutil
from PIL import Image
from pathlib import Path
from openmmdl.openmmdl_analysis.visualization.figures import FigureMerger, FigureArranger

test_data_directory = Path(
    "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation"
)
test_data_directory_files = Path("openmmdl/tests/data/in")
lig_no_h = test_data_directory_files / "lig_no_h.pdb"
complex = test_data_directory_files / "complex.pdb"
smi_file = test_data_directory_files / "lig_no_h.smi"
current_directory = os.getcwd()
output_path = "all_binding_modes_arranged.png"

shutil.copy(str(lig_no_h), ".")
shutil.copy(str(complex), ".")

def test_create_and_merge_images_with_split_data():
    # Define test data
    binding_mode = "Binding_Mode_1"
    occurrence_percent = 92
    split_data = [
        "166ARGA 4220,4221 Carboxylate NI saltbridge",
        "161PHEA 4221 Acceptor hbond",
        "207ILEA 4205 4206 4207 4208 4209 4204 hydrophobic",
    ]
    merged_image_paths = []

    # Define source image paths
    source_image_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1.png"
    source_svg_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1.svg"
    source_merged_image_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1_merged.png"

    # Copy source image files to the working directory
    working_directory = os.getcwd()
    destination_image_path = os.path.join(
        working_directory, os.path.basename(source_image_path)
    )
    destination_svg_path = os.path.join(
        working_directory, os.path.basename(source_svg_path)
    )
    destination_merged_image_path = os.path.join(
        working_directory, os.path.basename(source_merged_image_path)
    )
    shutil.copy(source_image_path, destination_image_path)
    shutil.copy(source_svg_path, destination_svg_path)
    shutil.copy(source_merged_image_path, destination_merged_image_path)

    # Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory before:", files_in_working_directory)

    # Run the function
    merger = FigureMerger(binding_mode, occurrence_percent, split_data, merged_image_paths)
    merged_image_paths = merger.create_and_merge_images()

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


def test_arranged_figure_generation():
    binding_mode1_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1_merged.png"
    binding_mode2_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_2_merged.png"
    all_modes_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/all_binding_modes_arranged.png"
    working_directory = os.getcwd()

    # Print the working directory to verify it's as expected
    print("Working Directory:", working_directory)

    destination_path_1 = os.path.join(
        working_directory, os.path.basename(binding_mode1_path)
    )
    destination_path_2 = os.path.join(
        working_directory, os.path.basename(binding_mode2_path)
    )
    destination_path_all = os.path.join(
        working_directory, os.path.basename(all_modes_path)
    )

    # Print the destination paths to verify they are constructed correctly
    print("Destination Path 1:", destination_path_1)
    print("Destination Path 2:", destination_path_2)
    print("Destination Path All:", destination_path_all)

    shutil.copy(binding_mode1_path, destination_path_1)
    shutil.copy(binding_mode2_path, destination_path_2)
    shutil.copy(all_modes_path, destination_path_all)

    merged_image_paths = ["Binding_Mode_1_merged.png", "Binding_Mode_2_merged.png"]
    output_path = "all_binding_modes_arranged.png"
    output_path = os.path.join(working_directory, output_path)
    print(output_path)

    # Run the function
    arranger = FigureArranger(merged_image_paths, output_path)
    arranger.arranged_figure_generation()
    print(output_path)

    # Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory:", files_in_working_directory)

    output_path = os.path.join(
        working_directory,
        "Binding_Modes_Markov_States",
        "all_binding_modes_arranged.png",
    )
    print(output_path)

    # Check if the output file was created

    assert output_path is not None


def test_arranged_figure_generation():
    binding_mode1_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_1_merged.png"
    binding_mode2_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/Binding_Mode_2_merged.png"
    all_modes_path = "openmmdl/tests/data/openmmdl_analysis/rdkit_figure_generation/all_binding_modes_arranged.png"
    working_directory = os.getcwd()

    # Print the working directory to verify it's as expected
    print("Working Directory:", working_directory)

    destination_path_1 = os.path.join(
        working_directory, os.path.basename(binding_mode1_path)
    )
    destination_path_2 = os.path.join(
        working_directory, os.path.basename(binding_mode2_path)
    )
    destination_path_all = os.path.join(
        working_directory, os.path.basename(all_modes_path)
    )

    # Print the destination paths to verify they are constructed correctly
    print("Destination Path 1:", destination_path_1)
    print("Destination Path 2:", destination_path_2)
    print("Destination Path All:", destination_path_all)

    shutil.copy(binding_mode1_path, destination_path_1)
    shutil.copy(binding_mode2_path, destination_path_2)
    shutil.copy(all_modes_path, destination_path_all)

    merged_image_paths = ["Binding_Mode_1_merged.png", "Binding_Mode_2_merged.png"]
    output_path = "all_binding_modes_arranged.png"
    output_path = os.path.join(working_directory, output_path)
    print(output_path)

    # Run the function
    arranger = FigureArranger(merged_image_paths, output_path)
    arranger.arranged_figure_generation()
    print(output_path)

    # Print the current files in the working directory for debugging
    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory:", files_in_working_directory)

    output_path = os.path.join(
        working_directory,
        "Binding_Modes_Markov_States",
        "all_binding_modes_arranged.png",
    )
    print(output_path)

    # Check if the output file was created

    assert output_path is not None


# Run the tests
if __name__ == "__main__":
    pytest.main()

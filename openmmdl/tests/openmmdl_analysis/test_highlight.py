import pytest
import os
import time
import shutil
from PIL import Image
from pathlib import Path
from openmmdl.openmmdl_analysis.highlight import FigureHighlighter, LigandImageGenerator

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


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        (
            ["60GLUA_4206_4207_4216_4217_4218_4205_hydrophobic"],
            ["60GLUA 4206 4207 4216 4217 4218 4205 hydrophobic"],
        ),
        (["165ASPA_4203_Acceptor_hbond"], ["165ASPA 4203 Acceptor hbond"]),
        (["125TYRA_4192_Acceptor_waterbridge"], ["125TYRA 4192 Acceptor waterbridge"]),
    ],
)

def test_split_interaction_data(input_data, expected_output):
    interaction = FigureHighlighter(complex, lig_no_h)
    result = interaction.split_interaction_data(input_data)
    assert result == expected_output


def test_update_dict():
    interaction = FigureHighlighter(complex, lig_no_h)
    # Test case 1: Check if the target dictionary is updated correctly
    target_dict = {1: "1", 2: "2"}
    source_dict = {3: "3", 4: "4"}
    interaction.update_dict(target_dict, source_dict)
    assert target_dict == {1: "1", 2: "2", 3: "3", 4: "4"}

    # Test case 2: Check if the function handles multiple source dictionaries
    target_dict = {}
    source_dict1 = {1: "1"}
    source_dict2 = {2: "2", 3: "3"}
    interaction.update_dict(target_dict, source_dict1, source_dict2)
    assert target_dict == {1: "1", 2: "2", 3: "3"}

    # Test case 3: Check if the function handles empty source dictionaries
    target_dict = {1: "1", 2: "2"}
    interaction.update_dict(target_dict)  # No source dictionaries provided
    assert target_dict == {1: "1", 2: "2"}


def test_generate_interaction_dict():
    # Test with a known interaction type 'hydrophobic'
    interaction_type = "hydrophobic"
    keys = [1, 2, 3]
    expected_result = {1: (1.0, 1.0, 0.0), 2: (1.0, 1.0, 0.0), 3: (1.0, 1.0, 0.0)}
    interaction = FigureHighlighter(complex, lig_no_h)
    result = interaction.generate_interaction_dict(interaction_type, keys)
    assert result == expected_result


# Test the generate_ligand_image function
def test_generate_ligand_image():
    ligand_name = "UNK"
    shutil.copy(complex, Path.cwd())
    shutil.copy(lig_no_h, Path.cwd())
    shutil.copy(smi_file, Path.cwd())
    output_image_file = "output_image.png"
    image = LigandImageGenerator(ligand_name, "complex.pdb", "lig_no_h.pdb", output_image_file)
    image.generate_image()

    # Assert that the output image file exists
    assert os.path.exists(output_image_file)


# Run the tests
if __name__ == "__main__":
    pytest.main()

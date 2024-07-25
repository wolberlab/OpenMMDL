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
        "FRAME": [1, 1, 2, 2, 3],
        "Interaction1": [1, 0, 1, 0, 0],
        "Interaction2": [0, 0, 0, 1, 1],
        "WATER_IDX": [101, 102, 103, 104, 105],
    }
    return pd.DataFrame(data)


def test_barcodegeneration(sample_dataframe_barcode_generation):
    interaction = "Interaction1"
    barcode = barcodegeneration(sample_dataframe_barcode_generation, interaction)

    assert isinstance(barcode, np.ndarray)

    expected_barcode = np.array([1, 1, 0])
    assert np.array_equal(barcode, expected_barcode)


def test_waterids_barcode_generator(sample_dataframe_barcode_generation):
    interaction = "Interaction2"
    waterid_barcode = waterids_barcode_generator(
        sample_dataframe_barcode_generation, interaction
    )

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
    barcode_data = {
        "166ARGA_4220,4221_Carboxylate_NI_saltbridge": np.array(
            [
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
            ]
        )
    }
    plot_barcodes(barcode_data, "single_barcode.png")
    single_barcode = "single_barcode.png"
    assert single_barcode is not None

    barcodes = {
        "Barcode 1": np.array(
            [
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
            ]
        ),
        "Barcode 2": np.array(
            [
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
                0,
                1,
            ]
        ),
        # Include more barcodes as needed
    }
    plot_barcodes(barcodes, "multiple_barcodes.png")

    files_in_working_directory = os.listdir(working_directory)
    print("Files in Working Directory after:", files_in_working_directory)
    save_path = "multiple_barcodes.png"

    assert save_path is not None


def test_plot_waterbridge_piechart(tmp_path):
    # Prepare inputs
    df_all = pd.DataFrame(
        {
            "interaction1": [1, 0, 1, 0],
            "interaction2": [0, 1, 0, 1],
            "WATER_IDX": [1, 2, 1, 2],  # changed 'waterid' to 'WATER_IDX'
            "FRAME": [0, 1, 2, 3],  # added 'FRAME' column
        }
    )
    waterbridge_barcodes = [np.array([1, 0, 1, 0]), np.array([0, 1, 0, 1])]
    waterbridge_interactions = ["interaction1", "interaction2"]
    fig_type = "png"

    # Change the current working directory to tmp_path

    # Use os.makedirs
    os.makedirs(f"{tmp_path}/Barcodes/Waterbridge_Piecharts/", exist_ok=True)

    # Call the function
    plot_waterbridge_piechart(
        df_all, waterbridge_barcodes, waterbridge_interactions, fig_type
    )

    # Check if the output files are created
    for interaction in waterbridge_interactions:
        outname_png = f"./Barcodes/Waterbridge_Piecharts/{interaction}.png"
        assert os.path.isfile(outname_png), f"File {outname_png} not found."

        # Additional assertions for content or specific properties of the generated files
        img = plt.imread(outname_png)
        assert img is not None, f"Unable to read image file {outname_png}."

        # Retrieve the percentage directly from the Axes object
        percentage_text = plt.gca().texts[0].get_text()
        assert percentage_text is not None, "Percentage text is None."

        # Retrieve the title directly from the Axes object
        title_text = plt.gca().get_title()
        assert title_text is not None, "Title text is None."

        # You can add more assertions based on your specific requirements
        # For example, check if the file size is greater than zero, etc.
        assert os.path.getsize(outname_png) > 0, f"File {outname_png} is empty."


def test_plot_bacodes_grouped(tmp_path):
    # Create a mock dataframe with all necessary columns
    df_all = pd.DataFrame(
        {
            "column1": [1, 2, 3],
            "column2": ["a", "b", "c"],
            "FRAME": [0, 1, 2],
            "atom1_atom2_interaction": [1, 0, 1],
            "atom3_atom4_interaction": [0, 1, 1],
        }
    )

    # Define interactions and interaction_type
    interactions = ["atom1_atom2_interaction", "atom3_atom4_interaction"]
    interaction_type = "interaction"
    fig_type = "png"

    working_directory = os.getcwd()
    plot_barcodes_grouped(interactions, df_all, interaction_type, fig_type)
    # Check if the output files were created
    assert os.path.exists(
        os.path.join(
            working_directory,
            "Barcodes",
            "atom2",
            f"atom2_{interaction_type}_barcodes.png",
        )
    )
    assert os.path.exists(
        os.path.join(
            working_directory,
            "Barcodes",
            "atom4",
            f"atom4_{interaction_type}_barcodes.png",
        )
    )
    assert os.path.exists(
        os.path.join(
            working_directory, "Barcodes", f"{interaction_type}_interactions.png"
        )
    )

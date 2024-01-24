import MDAnalysis as mda
import pytest
import pandas as pd
import re
import os
import shutil
from Bio.PDB.Structure import Structure
from Bio.PDB import PDBParser
from pathlib import Path
from unittest.mock import patch, mock_open
from openmmdl.openmmdl_analysis.find_stable_waters import (
    perform_clustering_and_writing,
    stable_waters_pipeline,
    trace_waters,
    filter_and_parse_pdb,
    write_pdb_clusters_and_representatives,
    find_interacting_residues,
    read_pdb_as_dataframe,
    analyze_protein_and_water_interaction,
)

# Fixtures and mock data setup

test_data_directory = Path("openmmdl/tests/data/in")
topology_file = f"{test_data_directory}/metal_top.pdb"
trajectory_file = f"{test_data_directory}/metal_traj_25.dcd"
csv_file_path = f"{test_data_directory}/stable_waters.csv"
repwat_file_path = f"{test_data_directory}/representative_waters.pdb"
output_dirs = []


def test_stable_waters_pipeline():
    water_eps_values = [0.5, 1.0, 2.0]  # Example epsilon values

    for water_eps in water_eps_values:
        output_directory = f"./test_output"
        stable_waters_pipeline(
            topology_file, trajectory_file, water_eps, output_directory
        )

        strEps = str(water_eps).replace(".", "")
        output_directory = f"./test_output_clusterEps_{strEps}"
        # Check if the expected output directory is created
        assert os.path.isdir(
            output_directory
        ), f"Directory {output_directory} was not created"
        output_dirs.append(output_directory)

        # Check if stable_waters.csv is created
        csv_file = os.path.join(output_directory, "stable_waters.csv")
        assert os.path.isfile(csv_file)

        # Load and verify the data in stable_waters.csv
        stable_waters_df = pd.read_csv(csv_file)
        assert not stable_waters_df.empty
        assert set(stable_waters_df.columns) == {
            "Frame",
            "Residue",
            "Oxygen_X",
            "Oxygen_Y",
            "Oxygen_Z",
        }

    # Cleanup: remove created directories and files
    for dir in output_dirs:
        shutil.rmtree(dir)


def test_perform_clustering():
    # Load the stable_waters data from the CSV file
    stable_waters_df = pd.read_csv(csv_file_path)

    # Define test parameters
    cluster_eps = 2

    u = mda.Universe(topology_file, trajectory_file)
    # Get the total number of frames for the progress bar
    total_frames = len(u.trajectory)
    output_directory = "./test_output_clustering"

    # Run the function
    perform_clustering_and_writing(
        stable_waters_df, cluster_eps, total_frames, output_directory
    )

    # Define the regular expression pattern for matching the line
    pattern = re.compile(
        r"ATOM\s+\d+\s+O\s+WAT\s+A\s+\d+\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+1\.00\s+0\.00\s+O"
    )

    # Assert subdirectory creation and file creation
    percentage_values = [25, 50, 75, 90, 99]
    for percent in percentage_values:
        min_samples = int((percent / 100) * total_frames)
        sub_directory = os.path.join(output_directory, f"clusterSize{min_samples}")

        assert os.path.isdir(sub_directory), f"Subdirectory for {percent}% not created"

        # Assuming the names of files created, adjust as necessary
        expected_files = [
            "cluster_0.pdb",
            "cluster_1.pdb",
        ]  # Replace with actual expected file names
        for file_name in expected_files:
            file_path = os.path.join(sub_directory, file_name)
            assert os.path.isfile(
                file_path
            ), f"File {file_name} was not created in {sub_directory}"

        # Check the contents of the files
        for file_name in expected_files:
            file_path = os.path.join(sub_directory, file_name)
            with open(file_path, "r") as file:
                # Read file and search for the pattern
                if not any(pattern.match(line) for line in file):
                    assert (
                        False
                    ), f"File {file_name} does not contain the required line format"

    # Cleanup
    shutil.rmtree(output_directory)


def test_write_pdb_clusters_and_representatives():
    # Mock data setup
    data = {
        "Oxygen_X": [1.0, 2.0, 3.0],
        "Oxygen_Y": [4.0, 5.0, 6.0],
        "Oxygen_Z": [7.0, 8.0, 9.0],
        "Cluster_Label": [0, 0, 1],
    }
    clustered_waters = pd.DataFrame(data)
    min_samples = 2
    output_sub_directory = "test_write_representatives"

    if os.path.exists(output_sub_directory):
        shutil.rmtree(output_sub_directory)
    os.makedirs(output_sub_directory, exist_ok=True)

    # Run the function
    write_pdb_clusters_and_representatives(
        clustered_waters, min_samples, output_sub_directory
    )

    # Assert file creation
    unique_labels = clustered_waters["Cluster_Label"].unique()
    for label in unique_labels:
        filename = os.path.join(output_sub_directory, f"cluster_{label}.pdb")
        assert os.path.isfile(filename), f"File {filename} not created"

    # Assert representative_waters.pdb creation and contents
    rep_file = os.path.join(output_sub_directory, "representative_waters.pdb")
    assert os.path.isfile(rep_file), "representative_waters.pdb not created"

    # Cleanup
    shutil.rmtree(output_sub_directory)


def test_filter_and_parse_pdb():
    # Call the function with the sample PDB file
    structure = filter_and_parse_pdb(topology_file)

    # Check if the returned object is a Structure
    assert isinstance(structure, Structure), "The returned object is not a Structure"


def test_find_interacting_residues():
    representative_waters_file = test_data_directory / "representative_waters.pdb"
    distance_threshold = 2.0  # Example threshold

    # Parse structure.pdb
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(topology_file))

    # Read representative_waters.pdb into a DataFrame
    waters_data = []
    with open(representative_waters_file, "r") as file:
        for line in file:
            if line.startswith("ATOM"):
                parts = line.split()
                x, y, z = map(float, parts[5:8])
                waters_data.append([x, y, z])
    representative_waters = pd.DataFrame(
        waters_data, columns=["Oxygen_X", "Oxygen_Y", "Oxygen_Z"]
    )

    # Run find_interacting_residues
    interacting_residues = find_interacting_residues(
        structure, representative_waters, distance_threshold
    )

    # Assert the results
    assert isinstance(interacting_residues, dict)
    # Example: Check if a specific water molecule interacts with any residues. We now fromthe test data that water 17 should interact.
    assert 17 in interacting_residues


def test_read_pdb_as_dataframe():
    # Mock PDB file content
    mock_pdb_content = (
        "ATOM      1  O   WAT A   1      26.091  60.495  24.828  1.00  0.00           O\n"
        "ATOM      2  O   WAT A   2      30.000  50.000  40.000  1.00  0.00           O\n"
    )

    # Expected data
    expected_data = [[26.091, 60.495, 24.828], [30.000, 50.000, 40.000]]
    expected_df = pd.DataFrame(
        expected_data, columns=["Oxygen_X", "Oxygen_Y", "Oxygen_Z"]
    )

    # Mock open function
    with patch("builtins.open", mock_open(read_data=mock_pdb_content)):
        # Call the function
        result_df = read_pdb_as_dataframe("dummy_path.pdb")

    # Assert DataFrame content
    pd.testing.assert_frame_equal(result_df, expected_df)


def test_analyze_protein_and_water_interaction():
    # Paths to the real PDB files

    protein_pdb_file = topology_file
    representative_waters_file = (
        "representative_waters.pdb"  # Assuming this is the correct name
    )

    # Setup output directory
    cluster_eps = 1.0  # Example value, adjust as needed
    strEps = str(cluster_eps).replace(".", "")
    output_directory = Path("testprotwatint/output_clusterEps_" + strEps)
    if output_directory.exists():
        shutil.rmtree(output_directory)
    os.makedirs(output_directory, exist_ok=True)

    # Create subdirectories and copy representative_waters.pdb into each
    mock_subdirectories = ["subdir1", "subdir2"]
    for subdir in mock_subdirectories:
        sub_path = output_directory / subdir
        os.makedirs(sub_path, exist_ok=True)
        shutil.copy(test_data_directory / representative_waters_file, sub_path)

    test_output_directory = Path("testprotwatint/output")
    os.makedirs(test_output_directory, exist_ok=True)
    # Run the function
    analyze_protein_and_water_interaction(
        str(protein_pdb_file),
        representative_waters_file,
        cluster_eps,
        str(test_output_directory),
        distance_threshold=5.0,
    )

    # Assert file creation in each subdirectory
    for subdir in mock_subdirectories:
        result_file = output_directory / subdir / "interacting_residues.csv"
        assert result_file.is_file(), f"File {result_file} not created"

    # Cleanup
    shutil.rmtree(output_directory)
    shutil.rmtree(test_output_directory)

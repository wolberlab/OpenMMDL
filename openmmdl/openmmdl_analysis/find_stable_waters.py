import MDAnalysis as mda
import os
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from tqdm import tqdm
from io import StringIO
from Bio.PDB import PDBParser


def trace_waters(topology, trajectory, output_directory):
    """trace the water molecules in a trajectory and write all which move below one Angstrom distance. To adjust the distance alter the integer
    Args:
        topology (pdb_file_name): Path to the topology file.
        trajectory (dcd_file_name): Path to the trajectory file.
        output_directory (str): Directory where output files will be saved.

    Returns:
        pd.DataFrame: DataFrame containing stable water coordinates.
        int: Total number of frames.
    """
    u = mda.Universe(topology, trajectory)
    # Get the total number of frames for the progress bar
    total_frames = len(u.trajectory)
    # Create an empty DataFrame to store stable water coordinates
    stable_waters = pd.DataFrame(
        columns=["Frame", "Residue", "Oxygen_X", "Oxygen_Y", "Oxygen_Z"]
    )

    # Initialize variables for the previous frame's data
    prev_frame_coords = {}

    # Iterate through frames with tqdm for the progress bar
    for ts in tqdm(
        u.trajectory, total=total_frames, desc="Processing frames for the wateranalysis"
    ):
        frame_num = ts.frame
        frame_coords = {}

        # Iterate through oxygen atoms of the specified water type
        # for atom in u.select_atoms(f"resname {water_type} and name O"):
        for atom in u.select_atoms(f"resname HOH and name O"):
            frame_coords[atom.index] = (
                atom.position[0],
                atom.position[1],
                atom.position[2],
            )

        # Check if it's not the first frame
        if frame_num > 0:
            stable_coords = []

            # Iterate through the oxygen atoms in the current frame
            for atom_index, coords in frame_coords.items():
                prev_coords = prev_frame_coords.get(
                    atom_index, coords
                )  # Get previous coordinates or use current if not found

                # Calculate the distance between current and previous coordinates
                distance = np.linalg.norm(np.array(coords) - np.array(prev_coords))

                # If the distance is less than 1 Angstrom, consider it a stable water
                if distance < 1:
                    stable_coords.append(
                        (frame_num, atom_index, coords[0], coords[1], coords[2])
                    )

            # Append stable water coordinates to the stable_waters DataFrame
            if stable_coords:  # Check if stable_coords is not empty
                stable_waters = pd.concat(
                    [
                        stable_waters,
                        pd.DataFrame(
                            stable_coords,
                            columns=[
                                "Frame",
                                "Residue",
                                "Oxygen_X",
                                "Oxygen_Y",
                                "Oxygen_Z",
                            ],
                        ),
                    ]
                )

        # Update the previous frame's coordinates
        prev_frame_coords = frame_coords

    stable_waters.to_csv(
        os.path.join(output_directory, "stable_waters.csv"), index=False
    )
    return stable_waters, total_frames


def perform_clustering_and_writing(
    stable_waters, cluster_eps, total_frames, output_directory
):
    """
    Args:
        stable_waters (pd.DataFrame): DataFrame containing stable water coordinates.
        cluster_eps (float): DBSCAN clustering epsilon parameter. This is in Angstrom in this case, and defines which Water distances should be within one cluster
        total_frames (int): Total number of frames.
        output_directory (str): Directory where output files will be saved.

    Returns:
        None, it writes files.
    """
    # Feature extraction: XYZ coordinates
    X = stable_waters[["Oxygen_X", "Oxygen_Y", "Oxygen_Z"]]

    # List of percentages to iterate over
    percentage_values = [25, 50, 75, 90, 99]

    for percent in percentage_values:
        min_percent = percent / 100
        min_samples = int(min_percent * total_frames)
        dbscan = DBSCAN(eps=cluster_eps, min_samples=min_samples)
        labels = dbscan.fit_predict(X)

        clustered_waters = stable_waters.copy()
        clustered_waters["Cluster_Label"] = labels
        print(clustered_waters["Cluster_Label"])
        clustered_waters = clustered_waters[clustered_waters["Cluster_Label"] != -1]

        output_sub_directory = os.path.join(
            output_directory, f"clusterSize{min_samples}"
        )
        os.makedirs(output_sub_directory, exist_ok=True)
        print("cluster_eps:")
        print(cluster_eps)
        write_pdb_clusters_and_representatives(
            clustered_waters, min_samples, output_sub_directory
        )


def write_pdb_clusters_and_representatives(
    clustered_waters, min_samples, output_sub_directory
):
    """
    Args:
        clustered_waters (pd.DataFrame): DataFrame containing clustered water coordinates.
        min_samples (int): Minimum number of samples for DBSCAN clustering.
        output_sub_directory (str): Subdirectory where output PDB files will be saved.

    Returns:
        None, it will output PDB files.
    """
    atom_counter = 1
    pdb_file_counter = 1
    print("minsamples:")
    print(min_samples)
    os.makedirs(output_sub_directory, exist_ok=True)
    with pd.option_context("display.max_rows", None):  # Temporarily set display options
        for label, cluster in clustered_waters.groupby("Cluster_Label"):
            pdb_lines = []
            for _, row in cluster.iterrows():
                x, y, z = row["Oxygen_X"], row["Oxygen_Y"], row["Oxygen_Z"]
                atom_counter = 1 if atom_counter > 9999 else atom_counter
                pdb_line = f"ATOM{atom_counter:6}  O   WAT A{atom_counter:4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
                pdb_lines.append(pdb_line)
                atom_counter += 1

            # Write the current cluster to a new PDB file
            output_filename = os.path.join(output_sub_directory, f"cluster_{label}.pdb")
            with open(output_filename, "w") as pdb_file:
                pdb_file.write("".join(pdb_lines))
                print(f"Cluster {label} written")

            pdb_file_counter += 1

        # Write representative water molecules to a PDB file
        representative_waters = clustered_waters.groupby("Cluster_Label").mean()
        representative_waters.reset_index(inplace=True)
        representative_filename = os.path.join(
            output_sub_directory, "representative_waters.pdb"
        )
        with open(representative_filename, "w") as pdb_file:
            for index, row in representative_waters.iterrows():
                x, y, z = row["Oxygen_X"], row["Oxygen_Y"], row["Oxygen_Z"]
                pdb_line = f"ATOM{index + 1:6}  O   WAT A{index + 1:4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
                pdb_file.write(pdb_line)


# Example usage
# stable_waters_pipeline("topology_file", "trajectory_file", 0.5)
def stable_waters_pipeline(
    topology, trajectory, water_eps, output_directory="./stableWaters"
):
    """Function to run the pipeline to extract stable water clusters, and their representatives from a PDB & DCD file
    Args:
        topology (PDB_file_name): Path to the topology file.
        trajectory (DCD_file_name): Path to the trajectory file.
        water_eps (float): DBSCAN clustering epsilon parameter.
        output_directory (str, optional): Directory where output files will be saved. Default is "./stableWaters".

    Returns:
        None, it starts the pipeline which will create output files.
    """
    # Load the PDB and DCD files
    output_directory += "_clusterEps_"
    strEps = str(water_eps).replace(".", "")
    output_directory += strEps
    os.makedirs(output_directory, exist_ok=True)
    # Create a stable waters list by calling the process_trajectory_and_cluster function
    stable_waters, total_frames = trace_waters(topology, trajectory, output_directory)
    # Now call perform_clustering_and_writing with the returned values
    perform_clustering_and_writing(
        stable_waters, water_eps, total_frames, output_directory
    )


def filter_and_parse_pdb(protein_pdb):
    """This function reads in a PDB and returns the structure with bioparser.
    Args:
        protein_pdb (PDB_file_path): Path to a protein PDB file.

    Returns:
        Structure: PDB structure object.
    """
    with open(protein_pdb, "r") as pdb_file:
        lines = [
            line
            for line in pdb_file
            if (
                line.startswith("ATOM")
                and line[17:20].strip() not in ["HOH", "WAT", "T4P", "T3P"]
                and line[22:26]
                .strip()
                .isdigit()  # Exclude lines with non-numeric sequence identifiers
            )
        ]

    # Convert the list of lines to a string buffer
    pdb_string = "".join(lines)
    pdb_buffer = StringIO(pdb_string)

    # Now parse the filtered lines
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_buffer)

    return structure


def find_interacting_residues(structure, representative_waters, distance_threshold):
    """This function maps waters (e.g. the representative waters) to interacting residues of a different PDB structure input. Use "filter_and_parse_pdb" to get the input for this function
    Args:
        structure (Structure): PDB structure object.
        representative_waters (pd.DataFrame): DataFrame containing representative water coordinates.
        distance_threshold (float): Threshold distance for identifying interacting residues.

    Returns:
        dict: Dictionary mapping cluster numbers to interacting residues.
    """
    interacting_residues = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                # Check if the residue is a protein residue (not a heteroatom or water molecule)
                if (
                    residue.id[0] == " "
                    and residue.id[2] == " "
                    and residue.resname not in ["HOH", "WAT"]
                ):
                    for wat_index, wat_row in representative_waters.iterrows():
                        wat_coords = np.array(
                            [
                                wat_row["Oxygen_X"],
                                wat_row["Oxygen_Y"],
                                wat_row["Oxygen_Z"],
                            ]
                        )
                        residue_coords = np.array(
                            [atom.get_coord() for atom in residue.get_atoms()][0]
                        )

                        distance = np.linalg.norm(wat_coords - residue_coords)
                        if distance < distance_threshold:
                            key = wat_index  # Assuming wat_index is the number of the water cluster
                            if key not in interacting_residues:
                                interacting_residues[key] = []
                            interacting_residues[key].append((chain.id, residue.id[1]))

    return interacting_residues


def read_pdb_as_dataframe(pdb_file):
    """Helper function reading a PDB
    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        pd.DataFrame: DataFrame containing PDB data.
    """
    lines = []
    with open(pdb_file, "r") as f:
        lines = f.readlines()

    # Extract relevant information from PDB file lines
    data = []
    for line in lines:
        if line.startswith("ATOM"):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            data.append([x, y, z])

    # Create a DataFrame
    columns = ["Oxygen_X", "Oxygen_Y", "Oxygen_Z"]
    representative_waters = pd.DataFrame(data, columns=columns)

    return representative_waters


# Analyse protein and water interaction, get the residues and the corresponding weater molecules that interact.
def analyze_protein_and_water_interaction(
    protein_pdb_file,
    representative_waters_file,
    cluster_eps,
    output_directory="./stableWaters",
    distance_threshold=5.0,
):
    """Analyse the interaction of residues to water molecules using a threshold that can be specified when calling the function
    Args:
        protein_pdb_file (str): Path to the protein PDB file without waters.
        representative_waters_file (str): Path to the representative waters PDB file, or any PDB file containing only waters
        cluster_eps (float): DBSCAN clustering epsilon parameter.
        output_directory (str, optional): Directory where output files will be saved. Default is "./stableWaters".
        distance_threshold (float, optional): Threshold distance for identifying interacting residues. Default is 5.0 (Angstrom).

    Returns:
        None, it will write a csv file.
    """
    output_directory += "_clusterEps_"
    strEps = str(cluster_eps).replace(".", "")
    output_directory += strEps

    # Iterate over subdirectories
    for subdirectory in os.listdir(output_directory):
        subdirectory_path = os.path.join(output_directory, subdirectory)
        if os.path.isdir(subdirectory_path):
            # Perform operations within each subdirectory
            representative_waters = read_pdb_as_dataframe(
                os.path.join(subdirectory_path, representative_waters_file)
            )
            filtered_structure = filter_and_parse_pdb(protein_pdb_file)
            interacting_residues = find_interacting_residues(
                filtered_structure, representative_waters, distance_threshold
            )
            result_df = pd.DataFrame(
                interacting_residues.items(),
                columns=["Cluster_Number", "Interacting_Residues"],
            )
            # Save result to each subdirectory
            result_df.to_csv(
                os.path.join(subdirectory_path, "interacting_residues.csv"), index=False
            )
            print(f"Exported interacting_residues.csv in {subdirectory_path}")

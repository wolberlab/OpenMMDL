import MDAnalysis as mda
import os
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from tqdm import tqdm
from io import StringIO
from Bio.PDB import PDBParser, Structure
from typing import Tuple, Dict, List, Optional


class StableWaters:
    def __init__(self, trajectory: str, topology: str, water_eps: float) -> None:
        self.trajectory = trajectory
        self.topology = topology
        self.u = mda.Universe(self.topology, self.trajectory)
        self.water_eps = water_eps

    def trace_waters(self, output_directory: str) -> Tuple[pd.DataFrame, int]:
        """Trace the water molecules in a trajectory and write all which move below one Angstrom distance. To adjust the distance alter the integer
        Args:
            output_directory (str): Directory where output files will be saved.

        Returns:
            Tuple[pd.DataFrame, int]: DataFrame containing stable water coordinates and total number of frames.
        """
        # Get the total number of frames for the progress bar
        total_frames = len(self.u.trajectory)
        
        # Create an empty DataFrame to store stable water coordinates
        stable_waters = pd.DataFrame(
            columns=["Frame", "Residue", "Oxygen_X", "Oxygen_Y", "Oxygen_Z"]
        )

        # Initialize variables for the previous frame's data
        prev_frame_coords = {}

        # Iterate through frames with tqdm for the progress bar
        for ts in tqdm(
            self.u.trajectory,
            total=total_frames,
            desc="Processing frames for the water analysis",
        ):
            frame_num = ts.frame
            frame_coords = {}

            # Iterate through oxygen atoms of the specified water type
            for atom in self.u.select_atoms("resname HOH and name O"):
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
                if stable_coords:
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
        self,
        stable_waters: pd.DataFrame,
        cluster_eps: float,
        total_frames: int,
        output_directory: str,
    ) -> None:
        """
        Perform DBSCAN clustering on the stable water coordinates, and write the clusters and their representatives to PDB files.

        Args:
            stable_waters (pd.DataFrame): DataFrame containing stable water coordinates.
            cluster_eps (float): DBSCAN clustering epsilon parameter. This is in Angstrom in this case
                                 and defines which Water distances should be within one cluster.
            total_frames (int): Total number of frames.
            output_directory (str): Directory where output files will be saved.
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
            clustered_waters = clustered_waters[clustered_waters["Cluster_Label"] != -1]

            output_sub_directory = os.path.join(
                output_directory, f"clusterSize{min_samples}"
            )
            os.makedirs(output_sub_directory, exist_ok=True)
            print("cluster_eps:")
            print(cluster_eps)
            self.write_pdb_clusters_and_representatives(
                clustered_waters, min_samples, output_sub_directory
            )

    def write_pdb_clusters_and_representatives(
        self,
        clustered_waters: pd.DataFrame,
        min_samples: int,
        output_sub_directory: str,
    ) -> None:
        """
        Writes the clusters and their representatives to PDB files.

        Args:
            clustered_waters (pd.DataFrame): DataFrame containing clustered water coordinates.
            min_samples (int): Minimum number of samples for DBSCAN clustering.
            output_sub_directory (str): Subdirectory where output PDB files will be saved.
        """
        atom_counter = 1
        pdb_file_counter = 1
        print("minsamples:")
        print(min_samples)
        os.makedirs(output_sub_directory, exist_ok=True)
        with pd.option_context("display.max_rows", None):
            for label, cluster in clustered_waters.groupby("Cluster_Label"):
                pdb_lines = []
                for _, row in cluster.iterrows():
                    x, y, z = row["Oxygen_X"], row["Oxygen_Y"], row["Oxygen_Z"]
                    atom_counter = 1 if atom_counter > 9999 else atom_counter
                    pdb_line = f"ATOM{atom_counter:6}  O   WAT A{atom_counter:4}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n"
                    pdb_lines.append(pdb_line)
                    atom_counter += 1

                # Write the current cluster to a new PDB file
                output_filename = os.path.join(
                    output_sub_directory, f"cluster_{label}.pdb"
                )
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

    def stable_waters_pipeline(self, output_directory: str = "./stableWaters") -> None:
        """Function to run the pipeline to extract stable water clusters, and their representatives from a PDB & DCD file
        Args:
            output_directory (str, optional): Directory where output files will be saved. Default is "./stableWaters".
        """
        # Load the PDB and DCD files
        output_directory += "_clusterEps_"
        strEps = str(self.water_eps).replace(".", "")
        output_directory += strEps
        os.makedirs(output_directory, exist_ok=True)
        
        # Create a stable waters list by calling the process_trajectory_and_cluster function
        stable_waters, total_frames = self.trace_waters(output_directory)
        # Now call perform_clustering_and_writing with the returned values
        self.perform_clustering_and_writing(
            stable_waters, self.water_eps, total_frames, output_directory
        )

    def filter_and_parse_pdb(self, protein_pdb: str) -> Structure:
        """Reads in a PDB and returns the structure with bioparser.
        Args:
            protein_pdb (str): Path to a protein PDB file.

        Returns:
            Structure: Biopython PDB structure object.
        """
        with open(protein_pdb, "r") as pdb_file:
            lines = [
                line
                for line in pdb_file
                if (
                    line.startswith("ATOM")
                    and line[17:20].strip() not in ["HOH", "WAT", "T4P", "T3P"]
                    and line[22:26].strip().isdigit()
                )
            ]

        # Convert the list of lines to a string buffer
        pdb_string = "".join(lines)
        pdb_buffer = StringIO(pdb_string)

        # Now parse the filtered lines
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_buffer)

        return structure

    def find_interacting_residues(
        self,
        structure: Structure,
        representative_waters: pd.DataFrame,
        distance_threshold: float,
    ) -> Dict[int, List[Tuple[str, int]]]:
        """Maps waters (e.g. the representative waters) to interacting residues of a different PDB structure input. 
           Use "filter_and_parse_pdb" to get the input for this function.
        Args:
            structure (Structure): Biopython PDB structure object.
            representative_waters (pd.DataFrame): DataFrame containing representative water coordinates.
            distance_threshold (float): Threshold distance for identifying interacting residues.

        Returns:
            Dict[int, List[Tuple[str, int]]]: Dictionary mapping cluster numbers to interacting residues.
        """
        interacting_residues = {}

        for model in structure:
            for chain in model:
                # Check if the residue is a protein residue (not a heteroatom or water molecule)
                for residue in chain:
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
                                key = wat_index # Assuming wat_index is the number of the water cluster
                                if key not in interacting_residues:
                                    interacting_residues[key] = []
                                interacting_residues[key].append(
                                    (chain.id, residue.id[1])
                                )

        return interacting_residues

    def read_pdb_as_dataframe(self, pdb_file: str) -> pd.DataFrame:
        """Helper function reading a PDB
        Args:
            pdb_file (str): Path to the PDB file.

        Returns:
            pd.DataFrame: DataFrame containing PDB data.
        """
        lines = []
        with open(pdb_file, "r") as f:
            lines = f.readlines()

        data = []
        for line in lines:
            if line.startswith("ATOM"):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                data.append([x, y, z])

        columns = ["Oxygen_X", "Oxygen_Y", "Oxygen_Z"]
        representative_waters = pd.DataFrame(data, columns=columns)

        return representative_waters

    def analyze_protein_and_water_interaction(
        self,
        protein_pdb_file: str,
        representative_waters_file: str,
        cluster_eps: float,
        output_directory: str = "./stableWaters",
        distance_threshold: float = 5.0,
    ) -> None:
        """Analyze the interaction of residues to water molecules using a threshold that can be specified when calling the function
        Args:
            protein_pdb_file (str): Path to the protein PDB file without waters.
            representative_waters_file (str): Path to the representative waters PDB file, or any PDB file containing only waters.
            cluster_eps (float): DBSCAN clustering epsilon parameter.
            output_directory (str, optional): Directory where output files will be saved. Default is "./stableWaters".
            distance_threshold (float, optional): Threshold distance for identifying interacting residues. Default is 5.0 (Angstrom).
        """
        output_directory += "_clusterEps_"
        strEps = str(cluster_eps).replace(".", "")
        output_directory += strEps

        # Iterate over subdirectories
        for subdirectory in os.listdir(output_directory):
            subdirectory_path = os.path.join(output_directory, subdirectory)
            if os.path.isdir(subdirectory_path):
                # Perform operations within each subdirectory
                representative_waters = self.read_pdb_as_dataframe(
                    os.path.join(subdirectory_path, representative_waters_file)
                )
                filtered_structure = self.filter_and_parse_pdb(protein_pdb_file)
                interacting_residues = self.find_interacting_residues(
                    filtered_structure, representative_waters, distance_threshold
                )
                result_df = pd.DataFrame(
                    interacting_residues.items(),
                    columns=["Cluster_Number", "Interacting_Residues"],
                )
                # Save result to each subdirectory
                result_df.to_csv(
                    os.path.join(subdirectory_path, "interacting_residues.csv"),
                    index=False,
                )
                print(f"Exported interacting_residues.csv in {subdirectory_path}")

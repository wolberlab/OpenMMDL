import os
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from plip.basic import config
from plip.structure.preparation import PDBComplex, PLInteraction
from plip.exchange.report import BindingSiteReport
from multiprocessing import Pool

config.KEEPMOD = True


class InteractionAnalyzer:
    """
    Analyzes molecular interactions between a protein and a ligand/peptide 
    throughout an MD trajectory using PLIP (Protein-Ligand Interaction Profiler).

    Attributes
    ----------
    pdb_md : mda.Universe
        MDAnalysis Universe object representing the topology and trajectory.
    dataframe : str or None
        Path to an existing interaction CSV file. If None, the trajectory will be processed anew.
    num_processes : int
        Number of CPU cores to use for parallel frame analysis.
    lig_name : str
        Residue name of the ligand in the complex.
    special : str
        Residue name for special ligands like metal ions (optional).
    peptide : str
        Chain ID of the peptide ligand (optional).
    md_len : int
        Number of frames in the trajectory.
    ineraction_list : pd.DataFrame
        DataFrame storing the extracted interactions across the trajectory.
    """
    def __init__(
        self,
        pdb_md,
        dataframe,
        num_processes,
        lig_name,
        special_ligand,
        peptide,
        md_len,
    ):
        self.pdb_md = pdb_md
        self.dataframe = dataframe
        self.num_processes = num_processes
        self.lig_name = lig_name
        self.special = special_ligand
        self.peptide = peptide
        self.md_len = md_len
        self.ineraction_list = self._process_trajectory()

    def _retrieve_plip_interactions(self, pdb_file, lig_name):
        """
        Retrieves the interactions from PLIP.

        Parameters
        ----------
        pdb_file : str 
            The path of the PDB file of the complex.
        lig_name : str 
            Name of the Ligand in the complex topology that will be analyzed.

        Returns
        -------
        dict
            A dictionary of the binding sites and the interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)  # load the pdb file
        for ligand in protlig.ligands:
            if str(ligand.longname) == lig_name:
                protlig.characterize_complex(
                    ligand
                )  # find ligands and analyze interactions
        sites = {}
        # loop over binding sites
        for key, site in sorted(protlig.interaction_sets.items()):
            binding_site = BindingSiteReport(site)  # collect data about interactions
            # tuples of *_features and *_info will be converted to pandas DataFrame
            keys = (
                "hydrophobic",
                "hbond",
                "waterbridge",
                "saltbridge",
                "pistacking",
                "pication",
                "halogen",
                "metal",
            )
            # interactions is a dictionary which contains relevant information for each
            # of the possible interactions: hydrophobic, hbond, etc. in the considered
            # binding site.
            interactions = {
                k: [getattr(binding_site, k + "_features")]
                + getattr(binding_site, k + "_info")
                for k in keys
            }
            sites[key] = interactions

        return sites

    def _retrieve_plip_interactions_peptide(self, pdb_file):
        """
        Retrives the interactions from PLIP for a peptide.

        Parameters
        ----------
        pdb_file : str 
            The path of the PDB file of the complex.

        Returns
        -------
        dict 
            A dictionary of the binding sites and the interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)  # load the pdb file
        protlig.characterize_complex(
            protlig.ligands[-1]
        )  # find ligands and analyze interactions
        sites = {}
        # loop over binding sites
        for key, site in sorted(protlig.interaction_sets.items()):
            binding_site = BindingSiteReport(site)  # collect data about interactions
            # tuples of *_features and *_info will be converted to pandas DataFrame
            keys = (
                "hydrophobic",
                "hbond",
                "waterbridge",
                "saltbridge",
                "pistacking",
                "pication",
                "halogen",
                "metal",
            )
            # interactions is a dictionary which contains relevant information for each
            # of the possible interactions: hydrophobic, hbond, etc. in the considered
            # binding site.
            interactions = {
                k: [getattr(binding_site, k + "_features")]
                + getattr(binding_site, k + "_info")
                for k in keys
            }
            sites[key] = interactions

        return sites

    def _create_df_from_binding_site(
        self, selected_site_interactions, interaction_type="hbond"
    ):
        """
        Creates a data frame from a binding site and interaction type.

        Parameters
        ----------
        selected_site_interactions : dict 
            Precaluclated interactions from PLIP for the selected site.
        interaction_type : str, optional
            The interaction type of interest (default set to hydrogen bond). Defaults to "hbond".

        Returns
        -------
        pd.DataFrame 
            DataFrame with information retrieved from PLIP.
        """
        # check if interaction type is valid:
        valid_types = [
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        ]

        if interaction_type not in valid_types:
            print(
                "\033[1m!!! Wrong interaction type specified. Hbond is chosen by default!!!\033[0m\n"
            )
            interaction_type = "hbond"

        df = pd.DataFrame.from_records(
            # data is stored AFTER the column names
            selected_site_interactions[interaction_type][1:],
            # column names are always the first element
            columns=selected_site_interactions[interaction_type][0],
        )
        return df

    def _change_lig_to_residue(self, file_path, new_residue_name):
        """
        Reformats the topology file to change the ligand to a residue. This is needed for interactions with special ligands such as metal ions.

        Parameters
        ----------
        file_path : str
            Filepath of the topology file.
        new_residue_name : str 
            New residue name of the ligand now changed to mimic an amino acid residue.

        Returns
        -------
        None
        """
        with open(file_path, "r") as file:
            lines = file.readlines()

        with open(file_path, "w") as file:
            for line in lines:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    # Assuming the standard PDB format for simplicity
                    # You may need to adapt this part based on your specific PDB file
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()

                    # Check if the residue name matches the one to be changed
                    if residue_name == self.lig_name:
                        # Change the residue name to the new one
                        modified_line = line[:17] + new_residue_name + line[20:]
                        file.write(modified_line)
                    else:
                        file.write(line)
                else:
                    file.write(line)

    def _process_frame(self, frame):
        """
        Process a single frame of MD simulation.

        Parameters
        ----------
        frame : int
            The number of the frame that will be processed.

        Returns
        -------
        pd.DataFrame 
            A dataframe conatining the interaction data for the processed frame.
        """
        atoms_selected = self.pdb_md.select_atoms(
            f"protein or nucleic or resname {self.lig_name} or (resname HOH and around 10 resname {self.lig_name}) or resname {self.special}"
        )
        for num in self.pdb_md.trajectory[(frame) : (frame + 1)]:
            atoms_selected.write(f"processing_frame_{frame}.pdb")
        if self.peptide is None:
            interactions_by_site = self._retrieve_plip_interactions(
                f"processing_frame_{frame}.pdb", self.lig_name
            )
            index_of_selected_site = -1
            selected_site = list(interactions_by_site.keys())[index_of_selected_site]

            interaction_types = [
                "hydrophobic",
                "hbond",
                "waterbridge",
                "saltbridge",
                "pistacking",
                "pication",
                "halogen",
                "metal",
            ]

            interaction_list = pd.DataFrame()
            for interaction_type in interaction_types:
                tmp_interaction = self._create_df_from_binding_site(
                    interactions_by_site[selected_site],
                    interaction_type=interaction_type,
                )
                tmp_interaction["FRAME"] = int(frame)
                tmp_interaction["INTERACTION"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])
            if os.path.exists(f"processing_frame_{frame}.pdb"):
                os.remove(f"processing_frame_{frame}.pdb")

            if self.special is not None:
                combi_lig_special = mda.Universe("ligand_special.pdb")
                complex = mda.Universe("complex.pdb")
                complex_all = complex.select_atoms("all")
                result = self._process_frame_special(frame)
                results_df = pd.concat(result, ignore_index=True)
                results_df = results_df[results_df["LOCATION"] == "protein.sidechain"]
                results_df["RESTYPE"] = results_df["RESTYPE"].replace(
                    ["HIS", "SER", "CYS"], self.lig_name
                )
                results_df["LOCATION"] = results_df["LOCATION"].replace(
                    "protein.sidechain", "ligand"
                )
                updated_target_idx = []

                for index, row in results_df.iterrows():
                    ligand_special_int_nr = int(row["TARGET_IDX"])
                    ligand_special_int_nr_atom = combi_lig_special.select_atoms(
                        f"id {ligand_special_int_nr}"
                    )
                    for atom in ligand_special_int_nr_atom:
                        atom_name = atom.name
                        # Adjust atom_name based on the specified conditions
                        if atom_name in ["N", "C", "O", "S"]:
                            atom_name = f"{atom_name}1"
                        else:
                            # Assuming the format is a single letter followed by a number
                            base_name, atom_number = atom_name[:-1], int(atom_name[-1])
                            new_atom_number = atom_number + 1
                            atom_name = f"{base_name}{new_atom_number}"
                        for complex_atom in complex_all:
                            complex_atom_name = complex_atom.name
                            if atom_name == complex_atom_name:
                                true_number = complex_atom.id
                                break  # Exit the loop once a match is found
                        updated_target_idx.append(true_number)

                # Update 'TARGET_IDX' in interaction_list
                results_df["TARGET_IDX"] = updated_target_idx
                interaction_list["TARGET_IDX"] = interaction_list["TARGET_IDX"]

                # Concatenate the updated results_df to interaction_list
                interaction_list = pd.concat([interaction_list, results_df])
        if self.peptide is not None:
            interactions_by_site = self._retrieve_plip_interactions_peptide(
                f"processing_frame_{frame}.pdb"
            )
            index_of_selected_site = -1
            selected_site = list(interactions_by_site.keys())[index_of_selected_site]

            interaction_types = [
                "hydrophobic",
                "hbond",
                "waterbridge",
                "saltbridge",
                "pistacking",
                "pication",
                "halogen",
                "metal",
            ]

            interaction_list = pd.DataFrame()
            for interaction_type in interaction_types:
                tmp_interaction = self._create_df_from_binding_site(
                    interactions_by_site[selected_site],
                    interaction_type=interaction_type,
                )
                tmp_interaction["FRAME"] = int(frame)
                tmp_interaction["INTERACTION"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])
            if os.path.exists(f"processing_frame_{frame}.pdb"):
                os.remove(f"processing_frame_{frame}.pdb")

        return interaction_list

    def _process_frame_special(self, frame):
        """
        Function extension of process_frame to process special ligands.

        Parameters
        ----------
        frame : int 
            Number of the frame that will be processed.

        Returns
        -------
        list of pd.DataFrame 
            List of dataframes containing the interaction data for the processed frame with the special ligand.
        """
        res_renaming = ["HIS", "SER", "CYS"]
        interaction_dfs = []
        for res in res_renaming:
            self.pdb_md.trajectory[frame]
            atoms_selected = self.pdb_md.select_atoms(
                f"resname {self.lig_name} or resname {self.special}"
            )
            atoms_selected.write(f"processing_frame_{frame}.pdb")
            self._change_lig_to_residue(f"processing_frame_{frame}.pdb", res)
            interactions_by_site = self._retrieve_plip_interactions(
                f"processing_frame_{frame}.pdb", self.special
            )
            index_of_selected_site = -1
            selected_site = list(interactions_by_site.keys())[index_of_selected_site]
            interaction_types = ["metal"]
            interaction_list = pd.DataFrame()
            for interaction_type in interaction_types:
                tmp_interaction = self._create_df_from_binding_site(
                    interactions_by_site[selected_site],
                    interaction_type=interaction_type,
                )
                tmp_interaction["FRAME"] = int(frame)
                tmp_interaction["INTERACTION"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])
            interaction_dfs.append(interaction_list)
            os.remove(f"processing_frame_{frame}.pdb")
        return interaction_dfs

    def _process_frame_wrapper(self, args):
        """
        Wrapper for the MD Trajectory procession.

        Parameters
        ----------
        args : tuple
            Tuple containing (frame_idx: int - number of the frame to be processed)

        Returns
        -------
        tuple
            Tuple containing the frame index and the result of from the process_frame function.
        """
        frame_idx, pdb_md, lig_name, special_ligand, peptide = args

        return frame_idx, self._process_frame(frame_idx)

    def _fill_missing_frames(self, df):
        """
        Fills the frames with no interactions in the DataFrame with placeholder values.

        Parameters
        ----------
        df : pd.DataFrame 
            The input DataFrame with frames that have no interactions.

        Returns
        -------
        pd.DataFrame 
            DataFrame with placeholder values in the frames with no interactions.
        """

        # Create a set containing all unique values in the 'FRAME' column
        existing_frames = set(df["FRAME"])

        # Create a list to store new rows for missing numbers
        missing_rows = []

        # Iterate through numbers from 0 to md_len
        for frame_number in range(1, self.md_len):
            if frame_number not in existing_frames:
                # Create a new row with 'FRAME' set to the missing number and other columns set to "skip"
                missing_row = {"FRAME": frame_number}
                for col in df.columns:
                    if col != "FRAME":
                        missing_row[col] = "skip"
                missing_rows.append(missing_row)

        # Concatenate the missing rows with the original DataFrame
        df = pd.concat([df, pd.DataFrame(missing_rows)], ignore_index=True)

        # Sort the DataFrame by the 'FRAME' column
        df.sort_values(by="FRAME", inplace=True)

        return df

    def _process_trajectory(self):
        """
        Process protein-ligand trajectory with multiple CPUs in parallel.

        Returns
        -------
        pd.DataFrame 
            A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
        """
        if self.dataframe is None:
            print("\033[1mProcessing protein-ligand trajectory\033[0m")
            print(f"\033[1mUsing {self.num_processes} CPUs\033[0m")

            with Pool(processes=self.num_processes) as pool:
                frame_args = [
                    (i, self.pdb_md, self.lig_name, self.special, self.peptide)
                    for i in range(1, self.md_len)
                ]

                # Initialize the progress bar with the total number of frames
                pbar = tqdm(
                    total=self.md_len - 1,
                    ascii=True,
                    desc="\033[1mAnalyzing frames\033[0m",
                )

                results = []
                for result in pool.imap(self._process_frame_wrapper, frame_args):
                    results.append(result)
                    pbar.update(1)  # Update the progress manually

            # Close the progress bar
            pbar.close()

            # Extract the results and sort them by frame index
            results.sort(key=lambda x: x[0])
            interaction_lists = [result[1] for result in results]

            interaction_list = pd.concat(interaction_lists)

            interaction_list.to_csv("interactions_gathered.csv")

        elif self.dataframe is not None:
            print(f"\033[1mGathering data from {self.dataframe}\033[0m")
            interaction_tmp = pd.read_csv(self.dataframe)
            interaction_list = interaction_tmp.drop(interaction_tmp.columns[0], axis=1)

        interaction_list["Prot_partner"] = (
            interaction_list["RESNR"].astype(str)
            + interaction_list["RESTYPE"]
            + interaction_list["RESCHAIN"]
        )

        interaction_list = self._fill_missing_frames(
            interaction_list,
        )

        print("\033[1mProtein-ligand trajectory processed\033[0m")

        return interaction_list

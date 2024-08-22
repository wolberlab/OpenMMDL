import os
import itertools
import pandas as pd
import numpy as np
from MDAnalysis.analysis import rms, diffusionmap
from MDAnalysis.analysis.distances import dist
from tqdm import tqdm
from pathlib import Path
from numba import jit
from typing import Dict, List, Union, Optional, Tuple


class BindingModeProcesser:
    def __init__(
        self,
        pdb_md,
        ligand: str,
        peptide: Optional[str],
        special: str,
        ligand_rings: List[List[int]],
        interaction_list: pd.DataFrame,
        threshold: float,
    ) -> None:
        self.pdb_md = pdb_md
        self.ligand = ligand
        self.peptide = peptide
        self.special = special
        self.ligand_rings = ligand_rings
        self.unique_columns_rings_grouped = self.gather_interactions(interaction_list)
        self.interaction_list, self.unique_data = self.process_interaction_wrapper(
            interaction_list, threshold / 100
        )
        self.interactions_all, self.unique_data_all = self.process_interaction_wrapper(
            interaction_list.copy(), 0.00001
        )

    def process_interaction_wrapper(
        self, interaction_list: pd.DataFrame, threshold: float
    ) -> Tuple[pd.DataFrame, Dict[str, str]]:
        filtered_values = self.filtering_values(threshold, interaction_list)
        interaction_list.fillna(0, inplace=True)
        unique_data = self.unique_data_generation(filtered_values)
        self.df_iteration_numbering(interaction_list, unique_data)
        return interaction_list, unique_data

    def gather_interactions(self, df: pd.DataFrame) -> Dict[int, Dict[int, str]]:
        """Process a DataFrame with the protein-ligand interaction and generate column names for each unique interaction.

        Args:
            df (pandas.DataFrame): DataFrame that contains the interaction data for the whole trajectory.

        Returns:
            dict: A dictionary with the keys being 'FRAME' numbers and values being dictionaries containing row indices and their corresponding unique column names for interactions.
        """
        unique_columns_rings = {}
        unique_columns_rings_grouped = {}

        if self.peptide is None:
            for index, row in df.iterrows():
                interaction = row["INTERACTION"]
                prot_partner = row["Prot_partner"]

                if interaction == "hydrophobic":
                    ligcarbonidx = row["LIGCARBONIDX"]
                    ring_found = False
                    for ligand_ring in self.ligand_rings:
                        if ligcarbonidx in ligand_ring:
                            numbers_as_strings = [str(idx) for idx in ligand_ring]
                            name_with_numbers = "_".join(numbers_as_strings)
                            col_name = (
                                f"{prot_partner}_{name_with_numbers}_{interaction}"
                            )
                            ring_found = True
                            break
                    if not ring_found:
                        col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
                elif interaction == "hbond":
                    if row["PROTISDON"]:
                        ligcarbonidx = int(row["ACCEPTORIDX"])
                        type_ = "Acceptor"
                    else:
                        ligcarbonidx = int(row["DONORIDX"])
                        type_ = "Donor"
                    col_name = f"{prot_partner}_{ligcarbonidx}_{type_}_{interaction}"
                elif interaction == "halogen":
                    ligcarbonidx = int(row["DON_IDX"])
                    halogen = row["DONORTYPE"]
                    col_name = f"{prot_partner}_{ligcarbonidx}_{halogen}_{interaction}"
                elif interaction == "waterbridge":
                    if row["PROTISDON"]:
                        ligcarbonidx = int(row["ACCEPTOR_IDX"])
                        type_ = "Acceptor"
                    else:
                        ligcarbonidx = int(row["DONOR_IDX"])
                        type_ = "Donor"
                    col_name = f"{prot_partner}_{ligcarbonidx}_{type_}_{interaction}"
                elif interaction == "pistacking":
                    ligcarbonidx = row["LIG_IDX_LIST"]
                    col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
                elif interaction == "pication":
                    ligidx = row["LIG_IDX_LIST"]
                    ligtype = row["LIG_GROUP"]
                    col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
                    col_name = col_name.replace(",", "_")
                elif interaction == "saltbridge":
                    ligidx = row["LIG_IDX_LIST"]
                    lig_group = row["LIG_GROUP"]
                    if row["PROTISPOS"]:
                        type_ = "NI"
                    else:
                        type_ = "PI"
                    col_name = (
                        f"{prot_partner}_{ligidx}_{lig_group}_{type_}_{interaction}"
                    )
                elif interaction == "metal":
                    special_ligand = row["RESTYPE_LIG"]
                    ligcarbonidx = int(row["TARGET_IDX"])
                    metal_type = row["METAL_TYPE"]
                    coordination = row["COORDINATION"]
                    col_name = f"{special_ligand}_{ligcarbonidx}_{metal_type}_{coordination}_{interaction}"

                frame_value = row["FRAME"]
                if frame_value not in unique_columns_rings_grouped:
                    unique_columns_rings_grouped[frame_value] = {}
                if interaction != "skip":
                    unique_columns_rings_grouped[frame_value][index] = col_name
                    unique_columns_rings[index] = col_name

        else:  # If peptide is not None
            for index, row in df.iterrows():
                interaction = row["INTERACTION"]
                prot_partner = row["Prot_partner"]
                peptide_partner = f"{row['RESNR_LIG']}{row['RESTYPE_LIG']}"

                if interaction == "hydrophobic":
                    col_name = f"{prot_partner}_{peptide_partner}_{interaction}"
                elif interaction == "hbond":
                    if row["PROTISDON"]:
                        type_ = "Acceptor"
                    else:
                        type_ = "Donor"
                    col_name = f"{prot_partner}_{peptide_partner}_{type_}_{interaction}"
                elif interaction == "halogen":
                    halogen = row["DONORTYPE"]
                    col_name = (
                        f"{prot_partner}_{peptide_partner}_{halogen}_{interaction}"
                    )
                elif interaction == "waterbridge":
                    if row["PROTISDON"]:
                        type_ = "Acceptor"
                    else:
                        type_ = "Donor"
                    col_name = f"{prot_partner}_{peptide_partner}_{type_}_{interaction}"
                elif interaction == "pistacking":
                    col_name = f"{prot_partner}_{peptide_partner}_{interaction}"
                elif interaction == "pication":
                    ligidx = f"{row['RESNR_LIG']}{row['RESTYPE_LIG']}"
                    ligtype = row["RESTYPE_LIG"]
                    col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
                    col_name = col_name.replace(",", "_")
                elif interaction == "saltbridge":
                    ligidx = f"{row['RESNR_LIG']}{row['RESTYPE_LIG']}"
                    lig_group = row["RESTYPE_LIG"]
                    if row["PROTISPOS"]:
                        type_ = "NI"
                    else:
                        type_ = "PI"
                    col_name = (
                        f"{prot_partner}_{ligidx}_{lig_group}_{type_}_{interaction}"
                    )
                elif interaction == "metal":
                    special_ligand = row["RESTYPE_LIG"]
                    ligcarbonidx = f"{row['RESNR_LIG']}{row['RESTYPE_LIG']}"
                    metal_type = row["METAL_TYPE"]
                    coordination = row["COORDINATION"]
                    col_name = f"{special_ligand}_{ligcarbonidx}_{metal_type}_{coordination}_{interaction}"

                frame_value = row["FRAME"]
                if frame_value not in unique_columns_rings_grouped:
                    unique_columns_rings_grouped[frame_value] = {}
                if interaction != "skip":
                    unique_columns_rings_grouped[frame_value][index] = col_name
                    unique_columns_rings[index] = col_name

        print("\033[1minteraction partners generated\033[0m")
        return unique_columns_rings_grouped

    def remove_duplicate_values(
        self, data: Dict[int, Dict[int, str]]
    ) -> Dict[int, Dict[int, str]]:
        """Remove the duplicate values from sub-dictionaries within the input dictionary.

        Args:
            data (dict): The input dictionary containing sub-dictionaries with possible duplicate values.

        Returns:
            dict: A dictionary without duplicate values.
        """
        unique_data = {}

        for key, sub_dict in data.items():
            unique_sub_dict = {}
            seen_values = set()

            for sub_key, value in sub_dict.items():
                if value not in seen_values:
                    unique_sub_dict[sub_key] = value
                    seen_values.add(value)

            if unique_sub_dict:
                unique_data[key] = unique_sub_dict

        return unique_data

    def combine_subdict_values(
        self, data: Dict[int, Dict[int, str]]
    ) -> Dict[str, List[str]]:
        """Combines the values from the individual sub-dictionaries into a single list.

        Args:
            data (dict): Dictionary with values that are sub-dictionaries.

        Returns:
            dict: A dictionary with a single key named 'all' that contains a list of all combined values from all the sub-dictionaries.
        """
        combined_data = {"all": []}
        for sub_dict in data.values():
            combined_data["all"].extend(sub_dict.values())

        return combined_data

    def filtering_values(self, threshold: float, df: pd.DataFrame) -> List[str]:
        """Filter and append values (interactions) to a DataFrame based on occurrence counts.

        Args:
            threshold (float): A threshold value used for filtering the values (interactions) based on occurrence count.
            df (pandas.DataFrame): DataFrame to which the filtered values (interactions) will be added.

        Returns:
            list: A list of values, with unique values and their corresponding occurrence counts.
        """
        frames = len(self.pdb_md.trajectory) - 1
        unique_data = self.remove_duplicate_values(self.unique_columns_rings_grouped)
        unique_colums_rings_all = self.combine_subdict_values(unique_data)
        all_values = unique_colums_rings_all["all"]

        occurrences = {}
        for value in all_values:
            occurrences[value] = occurrences.get(value, 0) + 1

        threshold = threshold * frames
        filtered_values = [
            value for value, count in occurrences.items() if count >= threshold
        ]

        for value in filtered_values:
            df[value] = None

        return filtered_values

    def unique_data_generation(self, filtered_values: List[str]) -> Dict[str, str]:
        """Generate a dictionary containing the unique interactions from a list of filtered values obtained by filtering_values.

        Args:
            filtered_values (list): A list of values, where the unique interactions are extracted from.

        Returns:
            dict: A dictionary containing the filtered unique interactions.
        """
        unique_data = {}
        for element in filtered_values:
            if element not in unique_data:
                unique_data[element] = element

        return unique_data

    def df_iteration_numbering(
        self, df: pd.DataFrame, unique_data: Dict[str, str]
    ) -> None:
        """Loop through the DataFrame and assign the values 1 and 0 to the rows, depending if the corresponding interaction from unique data is present.

        Args:
            df (pd.DataFrame): DataFrame which has the interaction data for all of the frames.
            unique_data (Dict[str, str]): Dictionary that contains the unique interactions obtained from unique_data_generation.
        """
        if self.peptide is None:
            for index, row in df.iterrows():
                if row["INTERACTION"] == "hydrophobic":
                    for col in unique_data.values():
                        if "hydrophobic" in col:
                            ligcarbonidx_check = str(int(row["LIGCARBONIDX"]))
                            if ligcarbonidx_check in col:
                                parts = col.split("_")
                                prot_partner = parts[0]
                                interaction = parts[-1]
                                condition = (row["Prot_partner"] == prot_partner) & (
                                    row["INTERACTION"] == interaction
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "hbond":
                    if row["PROTISDON"] == True:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, _type, interaction = (
                                    col.split("_")
                                )
                                ligcarbonidx = int(ligcarbonidx)
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["ACCEPTORIDX"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                    elif row["PROTISDON"] == False:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, _type, interaction = (
                                    col.split("_")
                                )
                                ligcarbonidx = int(ligcarbonidx)
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["DONORIDX"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "halogen":
                    for col in unique_data.values():
                        if "halogen" in col:
                            prot_partner, ligcarbonidx, halogen, interaction = (
                                col.split("_")
                            )
                            ligcarbonidx = int(ligcarbonidx)
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (int(row["DON_IDX"]) == ligcarbonidx)
                                & (row["DONORTYPE"] == halogen)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pistacking":
                    for col in unique_data.values():
                        if "pistacking" in col:
                            prot_partner, ligcarbonidx, interaction = col.split("_")
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["LIG_IDX_LIST"] == ligcarbonidx)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "waterbridge":
                    for col in unique_data.values():
                        if "waterbridge" in col:
                            if row["PROTISDON"]:
                                prot_partner, ligcarbonidx, _type, interaction = (
                                    col.split("_")
                                )
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["ACCEPTOR_IDX"]) == int(ligcarbonidx))
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                            else:
                                prot_partner, ligcarbonidx, _type, interaction = (
                                    col.split("_")
                                )
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["DONOR_IDX"]) == int(ligcarbonidx))
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pication":
                    for col in unique_data.values():
                        if "pication" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1:-2]
                            ligidx = ",".join(ligidx)
                            ligtype = parts[-2]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["LIG_IDX_LIST"] == ligidx)
                                & (row["LIG_GROUP"] == ligtype)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "saltbridge":
                    for col in unique_data.values():
                        if "saltbridge" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1:-3]
                            ligidx = ",".join(ligidx)
                            lig_group = parts[-3]
                            _type = parts[-2]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["LIG_IDX_LIST"] == ligidx)
                                & (row["LIG_GROUP"] == lig_group)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "metal":
                    for col in unique_data.values():
                        if "metal" in col:
                            parts = col.split("_")
                            special_ligand = parts[0]
                            ligidx = int(parts[1])
                            metal_type = parts[2]
                            interaction = parts[-1]
                            condition = (
                                (row["RESTYPE_LIG"] == special_ligand)
                                & (int(row["TARGET_IDX"]) == ligidx)
                                & (row["METAL_TYPE"] == metal_type)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0

        else:  # If self.peptide is not None
            for index, row in df.iterrows():
                if row["INTERACTION"] == "hydrophobic":
                    for col in unique_data.values():
                        if "hydrophobic" in col:
                            ligcarbonidx_check = str(int(row["RESNR_LIG"]))
                            if ligcarbonidx_check in col:
                                parts = col.split("_")
                                prot_partner = parts[0]
                                interaction = parts[-1]
                                condition = (row["Prot_partner"] == prot_partner) & (
                                    row["INTERACTION"] == interaction
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "hbond":
                    if row["PROTISDON"]:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, _type, interaction = (
                                    col.split("_")
                                )
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (
                                        (str(row["RESNR_LIG"]) + row["RESTYPE_LIG"])
                                        == ligcarbonidx
                                    )
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                    else:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, _type, interaction = (
                                    col.split("_")
                                )
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (
                                        ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]))
                                        == ligcarbonidx
                                    )
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "halogen":
                    for col in unique_data.values():
                        if "halogen" in col:
                            prot_partner, ligcarbonidx, halogen, interaction = (
                                col.split("_")
                            )
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (
                                    (str(row["RESNR_LIG"]) + row["RESTYPE_LIG"])
                                    == ligcarbonidx
                                )
                                & (row["DONORTYPE"] == halogen)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pistacking":
                    for col in unique_data.values():
                        if "pistacking" in col:
                            prot_partner, ligcarbonidx, interaction = col.split("_")
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (
                                    (str(row["RESNR_LIG"]) + row["RESTYPE_LIG"])
                                    == ligcarbonidx
                                )
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "waterbridge":
                    for col in unique_data.values():
                        if "waterbridge" in col:
                            prot_partner, ligcarbonidx, _type, interaction = col.split(
                                "_"
                            )
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (
                                    (str(row["RESNR_LIG"]) + row["RESTYPE_LIG"])
                                    == ligcarbonidx
                                )
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pication":
                    for col in unique_data.values():
                        if "pication" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1:-2]
                            ligidx = ",".join(ligidx)
                            ligtype = parts[-2]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["RESNR_LIG"] == ligidx)
                                & (row["RESTYPE_LIG"] == ligtype)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "saltbridge":
                    for col in unique_data.values():
                        if "saltbridge" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1:-3]
                            ligidx = ",".join(ligidx)
                            lig_group = parts[-3]
                            _type = parts[-2]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["RESNR_LIG"] == ligidx)
                                & (row["RESTYPE_LIG"] == lig_group)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "metal":
                    for col in unique_data.values():
                        if "metal" in col:
                            parts = col.split("_")
                            special_ligand = parts[0]
                            ligidx = int(parts[1])
                            metal_type = parts[2]
                            interaction = parts[-1]
                            condition = (
                                (row["RESTYPE_LIG"] == special_ligand)
                                & (int(row["TARGET_IDX"]) == ligidx)
                                & (row["METAL_TYPE"] == metal_type)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0

    def update_values(
        self, df: pd.DataFrame, new: pd.DataFrame, unique_data: Dict[str, str]
    ) -> None:
        """Update the values in the input DataFrame based upon the frame values and a reference DataFrame.

        Args:
        df (pandas.DataFrame): Input DataFrame that will be updated.
        new (pandas.DataFrame): The reference DataFrame containing values that are used to update the input DataFrame.
        unique_data (Dict[str, str]): A dictionary containing keys that represent the specific unique column names that need to be updated in the input DataFrame.
        """
        for idx, row in df.iterrows():
            frame_value = row["FRAME"]
            values_to_update = new.loc[frame_value, list(unique_data.values())]
            df.loc[idx, list(unique_data.values())] = values_to_update

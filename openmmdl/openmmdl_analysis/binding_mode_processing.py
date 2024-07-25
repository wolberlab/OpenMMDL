import os
import itertools
import pandas as pd
import numpy as np
from MDAnalysis.analysis import rms, diffusionmap
from MDAnalysis.analysis.distances import dist
from tqdm import tqdm
from pathlib import Path
from numba import jit


def gather_interactions(df, ligand_rings, peptide=None):
    """Process a DataFrame with the protein-ligand interaction and generate column names for each unique interaction.

    Args:
        df (pandas dataframe): DataFrame that contains the interaction data for the whole trajectory.
        ligand_rings (list): A list of the ligand ring information to recognize the atom numbers belonging to rings for hydrophobic interactions.
        peptide (str, optional): chainid of the peptide in the topology. Defaults to None.

    Returns:
        dict: A dictionary with the keys being 'FRAME' numbers and values being dictionaries containing row indices and their corresponding unique column names for interactions.
    """
    unique_columns_rings = {}
    unique_columns_rings_grouped = {}

    # Iterate through the rows of the DataFrame
    if peptide is None:
        for index, row in df.iterrows():
            # Check if the 'INTERACTION' is 'hydrophobic'
            if row["INTERACTION"] == "hydrophobic":
                # Get the values from the current row
                prot_partner = row["Prot_partner"]
                ligcarbonidx = row["LIGCARBONIDX"]
                interaction = row["INTERACTION"]
                ring_found = False
                # Concatenate the values to form a unique column name
                for ligand_ring in ligand_rings:
                    if ligcarbonidx in ligand_ring:
                        numbers_as_strings = [
                            str(ligcarbonidx) for ligcarbonidx in ligand_ring
                        ]
                        # Create the name with numbers separated by underscores
                        name_with_numbers = "_".join(numbers_as_strings)
                        col_name = f"{prot_partner}_{name_with_numbers}_{interaction}"
                        ring_found = True
                        break
                if not ring_found:
                    ligcarbonidx = int(row["LIGCARBONIDX"])
                    col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
            elif row["INTERACTION"] == "hbond":
                if row["PROTISDON"] == True:
                    prot_partner = row["Prot_partner"]
                    ligcarbonidx = int(row["ACCEPTORIDX"])
                    interaction = row["INTERACTION"]
                    type = "Acceptor"
                elif row["PROTISDON"] == False:
                    prot_partner = row["Prot_partner"]
                    ligcarbonidx = int(row["DONORIDX"])
                    interaction = row["INTERACTION"]
                    type = "Donor"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
            elif row["INTERACTION"] == "halogen":
                prot_partner = row["Prot_partner"]
                ligcarbonidx = int(row["DON_IDX"])
                halogen = row["DONORTYPE"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligcarbonidx}_{halogen}_{interaction}"
            elif row["INTERACTION"] == "waterbridge":
                if row["PROTISDON"] == True:
                    prot_partner = row["Prot_partner"]
                    ligcarbonidx = int(row["ACCEPTOR_IDX"])
                    interaction = row["INTERACTION"]
                    type = "Acceptor"
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
                elif row["PROTISDON"] == False:
                    prot_partner = row["Prot_partner"]
                    ligcarbonidx = int(row["DONOR_IDX"])
                    interaction = row["INTERACTION"]
                    type = "Donor"
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
            elif row["INTERACTION"] == "pistacking":
                prot_partner = row["Prot_partner"]
                ligcarbonidx = row["LIG_IDX_LIST"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
            elif row["INTERACTION"] == "pication":
                prot_partner = row["Prot_partner"]
                ligidx = row["LIG_IDX_LIST"]
                ligtype = row["LIG_GROUP"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
                col_name = col_name.replace(",", "_")
            elif row["INTERACTION"] == "saltbridge":
                prot_partner = row["Prot_partner"]
                ligidx = row["LIG_IDX_LIST"]
                lig_group = row["LIG_GROUP"]
                interaction = row["INTERACTION"]
                if row["PROTISPOS"] == True:
                    type = "NI"
                    # Concatenate the values to form a unique column name
                    col_name = (
                        f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                    )
                elif row["PROTISPOS"] == False:
                    type = "PI"
                    col_name = (
                        f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                    )
            elif row["INTERACTION"] == "metal":
                special_ligand = row["RESTYPE_LIG"]
                ligcarbonidx = int(row["TARGET_IDX"])
                metal_type = row["METAL_TYPE"]
                coordination = row["COORDINATION"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{special_ligand}_{ligcarbonidx}_{metal_type}_{coordination}_{interaction}"
            frame_value = row["FRAME"]
            if frame_value not in unique_columns_rings_grouped:
                unique_columns_rings_grouped[frame_value] = {}
            if row["INTERACTION"] != "skip":
                unique_columns_rings_grouped[frame_value][index] = col_name
                # Add the column name and its value to the dictionary
                unique_columns_rings[index] = col_name
    if peptide is not None:
        for index, row in df.iterrows():
            # Check if the 'INTERACTION' is 'hydrophobic'
            if row["INTERACTION"] == "hydrophobic":
                # Get the values from the current row
                prot_partner = row["Prot_partner"]
                peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                interaction = row["INTERACTION"]
                ring_found = False
                col_name = f"{prot_partner}_{peptide_partner}_{interaction}"
            elif row["INTERACTION"] == "hbond":
                if row["PROTISDON"] == True:
                    prot_partner = row["Prot_partner"]
                    peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    type = "Acceptor"
                elif row["PROTISDON"] == False:
                    prot_partner = row["Prot_partner"]
                    peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    type = "Donor"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{peptide_partner}_{type}_{interaction}"
            elif row["INTERACTION"] == "halogen":
                prot_partner = row["Prot_partner"]
                peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                halogen = row["DONORTYPE"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{peptide_partner}_{halogen}_{interaction}"
            elif row["INTERACTION"] == "waterbridge":
                if row["PROTISDON"] == True:
                    prot_partner = row["Prot_partner"]
                    peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    type = "Acceptor"
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{peptide_partner}_{type}_{interaction}"
                elif row["PROTISDON"] == False:
                    prot_partner = row["Prot_partner"]
                    peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    type = "Donor"
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{peptide_partner}_{type}_{interaction}"
            elif row["INTERACTION"] == "pistacking":
                prot_partner = row["Prot_partner"]
                peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{peptide_partner}_{interaction}"
            elif row["INTERACTION"] == "pication":
                prot_partner = row["Prot_partner"]
                ligidx = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                ligtype = row["RESTYPE_LIG"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
                col_name = col_name.replace(",", "_")
            elif row["INTERACTION"] == "saltbridge":
                prot_partner = row["Prot_partner"]
                ligidx = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                lig_group = row["RESTYPE_LIG"]
                interaction = row["INTERACTION"]
                if row["PROTISPOS"] == True:
                    type = "NI"
                    # Concatenate the values to form a unique column name
                    col_name = (
                        f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                    )
                elif row["PROTISPOS"] == False:
                    type = "PI"
                    col_name = (
                        f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                    )
            elif row["INTERACTION"] == "metal":
                special_ligand = row["RESTYPE_LIG"]
                ligcarbonidx = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                metal_type = row["METAL_TYPE"]
                coordination = row["COORDINATION"]
                interaction = row["INTERACTION"]
                # Concatenate the values to form a unique column name
                col_name = f"{special_ligand}_{ligcarbonidx}_{metal_type}_{coordination}_{interaction}"
            frame_value = row["FRAME"]
            if frame_value not in unique_columns_rings_grouped:
                unique_columns_rings_grouped[frame_value] = {}
            if row["INTERACTION"] != "skip":
                unique_columns_rings_grouped[frame_value][index] = col_name
                # Add the column name and its value to the dictionary
                unique_columns_rings[index] = col_name
    print("\033[1minteraction partners generated\033[0m")

    return unique_columns_rings_grouped


def remove_duplicate_values(data):
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


def combine_subdict_values(data):
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


def filtering_values(threshold, frames, df, unique_columns_rings_grouped):
    """Filter and append values (interactions) to a DataFrame based on occurrence counts.

    Args:
        threshold (float): A treshold value that is used for filtering of the values (interactions) based upon the occurence count.
        frames (int): The number of frames that is used to calculate the treshold.
        df (pandas dataframe): DataFrame to which the filtered values (interactions) will be added.
        unique_columns_rings_grouped (dict): Dictionary containing the grouped and unique values otained from gather_interactions.

    Returns:
        list: A list of values, with unique values and their corresponding occurence counts.
    """
    # Call the function to remove duplicate keys
    unique_data = remove_duplicate_values(unique_columns_rings_grouped)

    # Call the function to combine sub-dictionary values
    unique_colums_rings_all = combine_subdict_values(unique_data)

    # Flatten the list of values
    all_values = unique_colums_rings_all["all"]

    # Count the occurrences of each value
    occurrences = {}
    for value in all_values:
        if value in occurrences:
            occurrences[value] += 1
        else:
            occurrences[value] = 1

    # Calculate the threshold (20% of 1000)
    threshold = threshold * frames

    # Filter out values that appear in less than 20% of 1000
    filtered_values = [
        value for value, count in occurrences.items() if count >= threshold
    ]

    # Append the filtered values as columns to the DataFrame
    for value in filtered_values:
        df[value] = None

    return filtered_values


def unique_data_generation(filtered_values):
    """Generate a dictionary conataing the unique interactions from a list of filtered values obtained by filtering_values.

    Args:
        filtered_values (list): A list of values, where the unique interactions are extracted from.

    Returns:
        dict: A dictionary containing the filtered unique interactions.
    """
    # Create a new dictionary to store unique values
    unique_data = {}

    # Loop through the elements of the data_list
    for element in filtered_values:
        # Check if the element is not already in the unique_data dictionary keys
        if element not in unique_data:
            # If not, add the element as both key and value
            unique_data[element] = element

    return unique_data


def df_iteration_numbering(df, unique_data, peptide=None):
    """Loop through the DataFrame and assign the values 1 and 0 to the rows, depending if the corresponding interaction from unique data is present.

    Args:
        df (pandas dataframe): DataFrame which has the interaction data for all of the frames.
        unique_data (dict): Dictionary that contains the unique interactions obtained from unique_data_generation.
        peptide (str, optional): name of the peptide chainid in the original topology. Defaults to None.
    """
    if peptide is None:
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
                        else:
                            continue
            elif row["INTERACTION"] == "hbond":
                if row["PROTISDON"] == True:
                    for col in unique_data.values():
                        if "hbond" in col:
                            prot_partner, ligcarbonidx, type, interaction = col.split(
                                "_"
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
                            prot_partner, ligcarbonidx, type, interaction = col.split(
                                "_"
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
                        prot_partner, ligcarbonidx, halogen, interaction = col.split(
                            "_"
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
                        if row["PROTISDON"] == True:
                            prot_partner, ligcarbonidx, type, interaction = col.split(
                                "_"
                            )
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (int(row["ACCEPTOR_IDX"]) == int(ligcarbonidx))
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                        elif row["PROTISDON"] == False:
                            prot_partner, ligcarbonidx, type, interaction = col.split(
                                "_"
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
                        type = parts[-2]
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
    if peptide is not None:
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
                        else:
                            continue
            elif row["INTERACTION"] == "hbond":
                if row["PROTISDON"] == True:
                    for col in unique_data.values():
                        if "hbond" in col:
                            prot_partner, ligcarbonidx, type, interaction = col.split(
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
                elif row["PROTISDON"] == False:
                    for col in unique_data.values():
                        if "hbond" in col:
                            prot_partner, ligcarbonidx, type, interaction = col.split(
                                "_"
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
                        prot_partner, ligcarbonidx, halogen, interaction = col.split(
                            "_"
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
                        if row["PROTISDON"] == True:
                            prot_partner, ligcarbonidx, type, interaction = col.split(
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
                        elif row["PROTISDON"] == False:
                            prot_partner, ligcarbonidx, type, interaction = col.split(
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
                        ligidx = parts[1]
                        ligtype = parts[-2]
                        interaction = parts[-1]
                        condition = (
                            (row["Prot_partner"] == prot_partner)
                            & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligidx)
                            & (row["INTERACTION"] == interaction)
                        )
                        df.at[index, col] = 1 if condition else 0
            elif row["INTERACTION"] == "saltbridge":
                for col in unique_data.values():
                    if "saltbridge" in col:
                        parts = col.split("_")
                        prot_partner = parts[0]
                        ligidx = parts[1]
                        lig_group = parts[-3]
                        type = parts[-2]
                        interaction = parts[-1]
                        condition = (
                            (row["Prot_partner"] == prot_partner)
                            & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligidx)
                            & (row["INTERACTION"] == interaction)
                        )
                        df.at[index, col] = 1 if condition else 0
            elif row["INTERACTION"] == "metal":
                for col in unique_data.values():
                    if "metal" in col:
                        parts = col.split("_")
                        special_ligand = parts[0]
                        ligidx = parts[1]
                        metal_type = parts[2]
                        interaction = parts[-1]
                        condition = (
                            (row["RESTYPE_LIG"] == special_ligand)
                            & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligidx)
                            & (row["METAL_TYPE"] == metal_type)
                            & (row["INTERACTION"] == interaction)
                        )
                        df.at[index, col] = 1 if condition else 0


def update_values(df, new, unique_data):
    """Update the values in the input DataFrame based upon the frame values and an reference DataFrame.

    Args:
        df (pandas dataframe): Input DataFrame that will be updated.
        new (pandas dataframe): The reference DataFrame containing values that are used to update the input DataFrame.
        unique_data (dict): A dictionary containing keys that represent the specific unique column names that need to be updated in the input DataFrame.
    """
    for idx, row in df.iterrows():
        frame_value = row["FRAME"]
        values_to_update = new.loc[frame_value, list(unique_data.values())]
        df.loc[idx, list(unique_data.values())] = values_to_update


@jit
def calc_rmsd_2frames(ref, frame):
    """
    RMSD calculation between a reference and a frame.
    """
    dist = np.zeros(len(frame))
    for atom in range(len(frame)):
        dist[atom] = (
            (ref[atom][0] - frame[atom][0]) ** 2
            + (ref[atom][1] - frame[atom][1]) ** 2
            + (ref[atom][2] - frame[atom][2]) ** 2
        )

    return np.sqrt(dist.mean())


def calculate_distance_matrix(pdb_md, selection):
    distances = np.zeros((len(pdb_md.trajectory), len(pdb_md.trajectory)))
    # calculate distance matrix
    for i in tqdm(range(len(pdb_md.trajectory))):
        pdb_md.trajectory[i]
        frame_i = pdb_md.select_atoms(selection).positions
        # distances[i] = md.rmsd(traj_aligned, traj_aligned, frame=i)
        for j in range(i + 1, len(pdb_md.trajectory)):
            pdb_md.trajectory[j]
            frame_j = pdb_md.select_atoms(selection).positions
            rmsd = calc_rmsd_2frames(frame_i, frame_j)
            distances[i][j] = rmsd
            distances[j][i] = rmsd
    return distances


def calculate_representative_frame(bmode_frames, DM):
    """Calculates the most representative frame for a bindingmode. This is based uppon the averagwe RMSD of a frame to all other frames in the binding mode.

    Args:
        traj (mda universe): MDAnalysis univers containing the trajectory.
        bmode_frame_list (list): List of frames belonging to a binding mode.
        lig (str): Name of the ligand in the topology.

    Returns:
        int: Number of the most representative frame.
    """
    frames = bmode_frames
    mean_rmsd_per_frame = {}
    # first loop  : first frame
    for frame_i in frames:
        mean_rmsd_per_frame[frame_i] = 0
        # we will add the rmsd between theses 2 frames and then calcul the
        # mean
        for frame_j in frames:
            # We don't want to calcul the same frame.
            if not frame_j == frame_i:
                # we add to the corresponding value in the list of all rmsd
                # the RMSD betwween frame_i and frame_j
                mean_rmsd_per_frame[frame_i] += DM[frame_i - 1, frame_j - 1]
        # mean calculation
        mean_rmsd_per_frame[frame_i] /= len(frames)

        # Representative frame = frame with lower RMSD between all other
        # frame of the cluster
        repre = min(mean_rmsd_per_frame, key=mean_rmsd_per_frame.get)

    return repre

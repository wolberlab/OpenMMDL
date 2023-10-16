# Import libraries
from pathlib import Path
import pandas as pd
import itertools
import os

def gather_interactions(df, ligand_rings):
    """
    Process a DataFrame with the protein-ligand interaction and generate column names for each unique interaction.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame that contains the interaction data for the whole trajectory.
    ligand_rings : list
        A list of the ligand ring information to recognize the atom numbers belonging to rings for hydrophobic interactions.

    Returns
    -------
    dict :
        A dictionary with the keys being 'FRAME' numbers and values being dictionaries containing row indices and their corresponding unique column names for interactions.
    """
    unique_columns_rings = {}
    unique_columns_rings_grouped = {}

    # Iterate through the rows of the DataFrame
    for index, row in df.iterrows():
        # Check if the 'INTERACTION' is 'hydrophobic'
        if row['INTERACTION'] == 'hydrophobic':
            # Get the values from the current row
            prot_partner = row['Prot_partner']
            ligcarbonidx = row['LIGCARBONIDX']
            interaction = row['INTERACTION']
            ring_found = False
            # Concatenate the values to form a unique column name
            for ligand_ring in ligand_rings:
                if ligcarbonidx in ligand_ring:
                    numbers_as_strings = [str(ligcarbonidx) for ligcarbonidx in ligand_ring]
                    # Create the name with numbers separated by underscores
                    name_with_numbers = '_'.join(numbers_as_strings)
                    col_name = f"{prot_partner}_{name_with_numbers}_{interaction}"
                    ring_found = True
                    break
            if not ring_found:
                ligcarbonidx = int(row['LIGCARBONIDX'])
                col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
        elif row['INTERACTION'] == 'hbond':
            if row['PROTISDON'] == True:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['ACCEPTORIDX'])
                interaction = row['INTERACTION']
                type = "Acceptor"
            elif row['PROTISDON'] == False:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['DONORIDX'])
                interaction = row['INTERACTION']
                type = "Donor"
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
        elif row['INTERACTION'] == 'halogen':
            prot_partner = row['Prot_partner']
            ligcarbonidx = int(row['DON_IDX'])
            halogen = row['DONORTYPE']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{halogen}_{interaction}"
        elif row['INTERACTION'] == 'waterbridge':
            if row['PROTISDON'] == True:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['ACCEPTOR_IDX'])
                interaction = row['INTERACTION']
                type = "Acceptor"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
            elif row['PROTISDON'] == False:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['DONOR_IDX'])
                interaction = row['INTERACTION']
                type = "Donor"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
        elif row['INTERACTION'] == 'pistacking':
            prot_partner = row['Prot_partner']
            ligcarbonidx = row['LIG_IDX_LIST']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
        elif row['INTERACTION'] == 'pication':
            prot_partner = row['Prot_partner']
            ligidx = row['LIG_IDX_LIST']
            ligtype = row['LIG_GROUP']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
            col_name = col_name.replace(',', '_')
        elif row['INTERACTION'] == 'saltbridge':
            prot_partner = row['Prot_partner']
            ligidx = row['LIG_IDX_LIST']
            lig_group = row['LIG_GROUP']
            interaction = row['INTERACTION']
            if row['PROTISPOS'] == True:
                type = "NI"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
            elif row['PROTISPOS'] == False:
                type = "PI"
                col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
        elif row['INTERACTION'] == 'metal':
            prot_partner = row['Prot_partner']
            ligcarbonidx = row['METAL_IDX']
            metal_type = row['METAL_TYPE']
            location = row['LOCATION']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{metal_type}_{location}_{interaction}"
        frame_value = row['FRAME']
        if frame_value not in unique_columns_rings_grouped:
            unique_columns_rings_grouped[frame_value] = {}
        unique_columns_rings_grouped[frame_value][index] = col_name
            
        # Add the column name and its value to the dictionary
        unique_columns_rings[index] = col_name
    print("\033[1minteraction partners generated\033[0m")

    return unique_columns_rings_grouped


def remove_duplicate_values(data):
    """
    Remove the duplicate values from sub-dictionaries within the input dictionary.

    Parameters
    ----------
    data : dict
        The input dictionary containing sub-dictionaries with possible duplicate values.

    Returns
    -------
    dict :
        A dictionary without duplicate values.
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
    """
    Combines the values from the individual sub-dictionaries into a single list.

    Parameters
    ----------
    data : dict
        Dictionary with values that are sub-dictionaries.
        

    Returns
    -------
    dict :
        A dictionary with a single key named 'all' that contains a list of all combined values from all the sub-dictionaries.
    """
    combined_data = {'all': []}
    for sub_dict in data.values():
        combined_data['all'].extend(sub_dict.values())

    return combined_data

def filtering_values(threshold, frames, df, unique_columns_rings_grouped):
    """
    Filter and append values (interactions) to a DataFrame based on occurrence counts.

    Parameters
    ----------
    treshold : flaot
        A treshold value that is used for filtering of the values (interactions) based upon the occurence count.
    frames : int
        The number of frames that is used to calculate the treshold.
    df : pandas.DataFrame
        DataFrame to which the filtered values (interactions) will be added.
    unique_columns_rings_grouped : dict
        Dictionary containing the grouped and unique values otained from gather_interactions.

    Returns
    -------
    dict :
        A dictionary with a single key named 'all' that contains a list of all combined values from all the sub-dictionaries.
    """
    # Call the function to remove duplicate keys
    unique_data = remove_duplicate_values(unique_columns_rings_grouped)

    # Call the function to combine sub-dictionary values
    unique_colums_rings_all = combine_subdict_values(unique_data)

    # Flatten the list of values
    all_values = unique_colums_rings_all['all']

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
    filtered_values = [value for value, count in occurrences.items() if count >= threshold]

    # Append the filtered values as columns to the DataFrame
    for value in filtered_values:
        df[value] = None
    
    return filtered_values


def unique_data_generation(filtered_values):
    """
    Generate a dictionary conataing the unique interactions from a list of filtered values obtained by filtering_values.

    Parameters
    ----------
    filtered_values : list
       A list of values, where the unique interactions are extracted from.
    
    Returns
    -------
    dict :
        A dictionary containing the filtered unique interactions.
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


def df_iteration_numbering(df,unique_data):
    """
    Loop through the DataFrame and assign the values 1 and 0 to the rows, depending if the corresponding interaction from unique data is present.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame which has the interaction data for all of the frames.
    unique_data : dict
        Dictionary that contains the unique interactions obtained from unique_data_generation.

    Returns
    -------
    None
    """
    for index, row in df.iterrows():
        if row['INTERACTION'] == "hydrophobic":
            for col in unique_data.values():
                if "hydrophobic" in col:
                    ligcarbonidx_check = str(int(row['LIGCARBONIDX']))
                    if ligcarbonidx_check in col:
                        parts = col.split('_')
                        prot_partner = parts[0]
                        interaction = parts[-1]
                        condition = (row['Prot_partner'] == prot_partner) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
                    else:
                        continue
        elif row['INTERACTION'] == "hbond":
            if row['PROTISDON'] == True:
                for col in unique_data.values():
                    if "hbond" in col:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        ligcarbonidx = int(ligcarbonidx)
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['ACCEPTORIDX']) == ligcarbonidx) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
            elif row['PROTISDON'] == False:
                for col in unique_data.values():
                    if "hbond" in col:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        ligcarbonidx = int(ligcarbonidx)
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['DONORIDX']) == ligcarbonidx) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "halogen":
            for col in unique_data.values():
                if "halogen" in col:
                    prot_partner, ligcarbonidx, halogen, interaction = col.split('_')
                    ligcarbonidx = int(ligcarbonidx)
                    condition = (row['Prot_partner'] == prot_partner) & (int(row['DON_IDX']) == ligcarbonidx) & (row['DONORTYPE'] == halogen) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "pistacking":
            for col in unique_data.values():
                if "pistacking" in col:
                    prot_partner, ligcarbonidx, interaction = col.split('_')
                    condition = (row['Prot_partner'] == prot_partner) & (row['LIG_IDX_LIST'] == ligcarbonidx) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "waterbridge":
            for col in unique_data.values():
                if "waterbridge" in col:
                    if row['PROTISDON'] == True:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['ACCEPTOR_IDX']) == int(ligcarbonidx)) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
                    elif row['PROTISDON'] == False:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['DONOR_IDX']) == int(ligcarbonidx)) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "pication":
            for col in unique_data.values():
                if "pication" in col:
                    parts = col.split('_')
                    prot_partner = parts[0]
                    ligidx = parts[1:-2]
                    ligidx = ','.join(ligidx)
                    ligtype = parts[-2]
                    interaction = parts[-1]
                    condition = (row['Prot_partner'] == prot_partner) & (row['LIG_IDX_LIST'] == ligidx) & (row['LIG_GROUP'] == ligtype) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0

        elif row['INTERACTION'] == "saltbridge":
            for col in unique_data.values():
                if "saltbridge" in col:
                    parts = col.split('_')
                    prot_partner = parts[0]
                    ligidx = parts[1:-3]
                    ligidx = ','.join(ligidx)
                    lig_group = parts[-3]
                    type = parts[-2]
                    interaction = parts[-1]
                    condition = (row['Prot_partner'] == prot_partner) & (row['LIG_IDX_LIST'] == ligidx) & (row['LIG_GROUP'] == lig_group) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0


def update_values(df, new, unique_data):
    """
    Update the values in the input DataFrame based upon the frame values and an reference DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame that will be updated.
    new_df : pandas.DataFrame
        The reference DataFrame containing values that are used to update the input DataFrame.
    unique_data : dict
        A dictionary containing keys that represent the specific unique column names that need to be updated in the input DataFrame.

    Returns
    -------
    None
    """
    for idx, row in df.iterrows():
        frame_value = row['FRAME']
        values_to_update = new.loc[frame_value, list(unique_data.values())]
        df.loc[idx, list(unique_data.values())] = values_to_update



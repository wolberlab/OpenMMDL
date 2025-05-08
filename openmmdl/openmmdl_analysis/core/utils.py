import pandas as pd
from io import StringIO
from Bio.PDB import PDBParser


def update_dict(target_dict, *source_dicts):
    """
    Updates the dictionary with the keys and values from other dictionaries.
    Only new keys (not already present in the target) are added.

    Parameters
    ----------
    target_dict : dict 
        The dictionary that needs to be updated with new keys and values.
    source_dicts : dict 
        One or multiple dictionaries that are used to update the target dictionary with new keys and values.

    Returns
    -------
    None
    """
    for source_dict in source_dicts:
        for key, value in source_dict.items():
            int_key = int(key)
            if int_key not in target_dict:
                target_dict[int_key] = value


def combine_subdict_values(data):
    """
    Combines the values from the individual sub-dictionaries into a single list.

    Parameters
    ----------
    data : dict
        Dictionary with values that are sub-dictionaries.

    Returns
    -------
    dict
        A dictionary with a single key named 'all' that contains a list of all combined values from all the sub-dictionaries.
    """
    combined_data = {"all": []}
    for sub_dict in data.values():
        combined_data["all"].extend(sub_dict.values())

    return combined_data


def update_values(df, new, unique_data, row_name):
    """
    Update the values in the input DataFrame based upon the frame values and an reference DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame that will be updated.
    new : pd.DataFrame 
        The reference DataFrame containing values that are used to update the input DataFrame.
    unique_data : dict 
        A dictionary containing keys that represent the specific unique column names that need to be updated in the input DataFrame.
    row_name : str
        The name of the column in the DataFrame used to index into new DataFrame.

    Returns
    -------
    None
        This function updates a dataframe and does not return anything.
    """
    for idx, row in df.iterrows():
        frame_value = row[row_name]
        values_to_update = new.loc[frame_value, list(unique_data.values())]
        df.loc[idx, list(unique_data.values())] = values_to_update


def remove_duplicate_values(data):
    """
    Remove the duplicate values from sub-dictionaries within the input dictionary.

    Parameters
    ----------
    data : dict
        The input dictionary containing sub-dictionaries with possible duplicate values.

    Returns
    -------
    dict
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

def read_pdb_as_dataframe(pdb_file):
    """
    Helper function reading a PDB

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing PDB data of the x, y, z coordinates of atoms.

    Notes
    -----
    This function extracts only lines starting with 'ATOM' and parses the
    x, y, z coordinates based on selected fields in the PDB format.
    Assumes coordinates are located at columns 31–54.
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

def filter_and_parse_pdb(protein_pdb):
    """
    This function reads in a PDB and returns the structure with bioparser.
    
    Parameters
    ----------
    protein_pdb : str
        Path to a protein PDB file.

    Returns
    -------
    Bio.PDB.Structure.Structure
        Parsed PDB structure object containing protein atoms.

    Notes
    -----
    The function:
    - Includes only lines starting with 'ATOM'.
    - Excludes water molecules (residue names 'HOH', 'WAT') and terminal phosphates ('T4P', 'T3P').
    - Skips lines with non-numeric residue sequence identifiers.
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

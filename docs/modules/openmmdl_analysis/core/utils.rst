API Documentation for utils
============================

.. py:function:: update_dict(target_dict, *source_dicts)

    Updates the dictionary with the keys and values from other dictionaries.
    Only new keys (not already present in the target) are added.

    :param dict target_dict: The dictionary that needs to be updated with new keys and values.
    :param dict source_dicts: One or multiple dictionaries that are used to update the target dictionary with new keys and values.
    :returns: None
    :rtype: None


.. py:function:: combine_subdict_values(data)

    Combines the values from the individual sub-dictionaries into a single list.

    :param dict data: Dictionary with values that are sub-dictionaries.
    :returns: A dictionary with a single key named 'all' that contains a list of all combined values from all the sub-dictionaries.
    :rtype: dict


.. py:function:: update_values(df, new, unique_data, row_name)

    Update the values in the input DataFrame based upon the frame values and an reference DataFrame.

    :param pd.DataFrame df: Input DataFrame that will be updated.
    :param pd.DataFrame new: The reference DataFrame containing values that are used to update the input DataFrame.
    :param dict unique_data: A dictionary containing keys that represent the specific unique column names that need to be updated in the input DataFrame.
    :param str row_name: The name of the column in the DataFrame used to index into new DataFrame.
    :returns: None. This function updates a dataframe and does not return anything.
    :rtype: None


.. py:function:: remove_duplicate_values(data)

    Remove the duplicate values from sub-dictionaries within the input dictionary.

    :param dict data: The input dictionary containing sub-dictionaries with possible duplicate values.
    :returns: A dictionary without duplicate values.
    :rtype: dict


.. py:function:: read_pdb_as_dataframe(pdb_file)

    Helper function reading a PDB

    :param str pdb_file: Path to the PDB file.
    :returns: DataFrame containing PDB data of the x, y, z coordinates of atoms.
    :rtype: pd.DataFrame

    .. note::
       This function extracts only lines starting with 'ATOM' and parses the
       x, y, z coordinates based on selected fields in the PDB format.
       Assumes coordinates are located at columns 31–54.


.. py:function:: filter_and_parse_pdb(protein_pdb)

    This function reads in a PDB and returns the structure with bioparser.

    :param str protein_pdb: Path to a protein PDB file.
    :returns: Parsed PDB structure object containing protein atoms.
    :rtype: Bio.PDB.Structure.Structure

    .. note::
       The function:
       - Includes only lines starting with 'ATOM'.
       - Excludes water molecules (residue names 'HOH', 'WAT') and terminal phosphates ('T4P', 'T3P').
       - Skips lines with non-numeric residue sequence identifiers.

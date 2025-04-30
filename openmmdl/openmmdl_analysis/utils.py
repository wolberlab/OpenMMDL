

def update_dict(target_dict, *source_dicts):
    """Updates the dictionary with the keys and values from other dictionaries.

    Args:
        target_dict (dict): The dictionary that needs to be updated with new keys and values.
        source_dicts (dict): One or multiple dictionaries that are used to update the target dictionary with new keys and values.
    """
    for source_dict in source_dicts:
        for key, value in source_dict.items():
            int_key = int(key)
            if int_key not in target_dict:
                target_dict[int_key] = value


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


def update_values(self, df, new, unique_data, row_name):
    """Update the values in the input DataFrame based upon the frame values and an reference DataFrame.

    Args:
        df (pandas dataframe): Input DataFrame that will be updated.
        new (pandas dataframe): The reference DataFrame containing values that are used to update the input DataFrame.
        unique_data (dict): A dictionary containing keys that represent the specific unique column names that need to be updated in the input DataFrame.
    """
    for idx, row in df.iterrows():
        frame_value = row[row_name]
        values_to_update = new.loc[frame_value, list(unique_data.values())]
        df.loc[idx, list(unique_data.values())] = values_to_update


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

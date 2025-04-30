

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

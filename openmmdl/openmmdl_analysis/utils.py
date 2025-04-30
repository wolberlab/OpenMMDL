

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

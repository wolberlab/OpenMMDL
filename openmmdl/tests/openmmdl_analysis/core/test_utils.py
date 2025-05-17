import pytest
from openmmdl.openmmdl_analysis.core.utils import update_dict, update_values


def test_update_dict():
    # Test case 1: Check if the target dictionary is updated correctly
    target_dict = {1: "1", 2: "2"}
    source_dict = {3: "3", 4: "4"}
    update_dict(target_dict, source_dict)
    assert target_dict == {1: "1", 2: "2", 3: "3", 4: "4"}

    # Test case 2: Check if the function handles multiple source dictionaries
    target_dict = {}
    source_dict1 = {1: "1"}
    source_dict2 = {2: "2", 3: "3"}
    update_dict(target_dict, source_dict1, source_dict2)
    assert target_dict == {1: "1", 2: "2", 3: "3"}

    # Test case 3: Check if the function handles empty source dictionaries
    target_dict = {1: "1", 2: "2"}
    update_dict(target_dict)  # No source dictionaries provided
    assert target_dict == {1: "1", 2: "2"}

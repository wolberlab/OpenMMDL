import pytest
from openmmdl.openmmdl_analysis.core.utils import update_dict, update_values, combine_subdict_values, extract_ints, coord_str, remove_duplicate_values


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

def test_combine_subdict_values_flattens_subdict_values_in_order():
    data = {1: {"a": 1, "b": 2}, 2: {"c": 3}}
    assert combine_subdict_values(data) == {"all": [1, 2, 3]}


def test_extract_ints_basic_and_leading_zeros():
    assert extract_ints("A12_B034") == [12, 34]

def test_remove_duplicate_values_keeps_first_key_for_each_value_and_drops_empty_subdicts():
    data = {
        1: {"a": 1, "b": 1, "c": 2},   # duplicate 1 -> keep "a"
        2: {"d": 3, "e": 3},          # duplicate 3 -> keep "d"
        3: {"x": 4, "y": 5},          # no duplicates
        4: {"p": 9, "q": 9, "r": 9},  # becomes {"p": 9}
    }

    assert remove_duplicate_values(data) == {
        1: {"a": 1, "c": 2},
        2: {"d": 3},
        3: {"x": 4, "y": 5},
        4: {"p": 9},
    }
def test_extract_ints_no_digits_and_non_string():
    assert extract_ints("no_digits") == []
    assert extract_ints(105) == [105]


def test_coord_str_formats_and_skip():
    assert coord_str([1, 2.34567, -0.1]) == "(1.000, 2.346, -0.100)"
    assert coord_str(None) == "skip"

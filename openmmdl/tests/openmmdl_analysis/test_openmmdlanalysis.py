import argparse

import pytest

from openmmdl.openmmdl_analysis.openmmdlanalysis import parse_bool_flag


@pytest.mark.parametrize(
    "value, expected",
    [
        (True, True),
        (False, False),
    ],
)
def test_parse_bool_flag_returns_bool_inputs_unchanged(value, expected):
    assert parse_bool_flag(value) is expected


@pytest.mark.parametrize(
    "value",
    [
        "1",
        1,
        "true",
        "TRUE",
        " t ",
        "yes",
        "Y",
        " on ",
    ],
)
def test_parse_bool_flag_accepts_truthy_values(value):
    assert parse_bool_flag(value) is True


@pytest.mark.parametrize(
    "value",
    [
        "0",
        0,
        "false",
        "FALSE",
        " f ",
        "no",
        "N",
        " off ",
    ],
)
def test_parse_bool_flag_accepts_falsy_values(value):
    assert parse_bool_flag(value) is False


@pytest.mark.parametrize(
    "value",
    ["", "maybe", "2", "enable", None],
)
def test_parse_bool_flag_rejects_invalid_values(value):
    with pytest.raises(argparse.ArgumentTypeError) as exc_info:
        parse_bool_flag(value)

    assert str(exc_info.value) == (
        "Boolean value expected. Use one of: true/false, yes/no, y/n, 1/0, on/off."
    )

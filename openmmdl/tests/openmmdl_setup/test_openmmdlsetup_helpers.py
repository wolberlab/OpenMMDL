from openmmdl.openmmdl_setup.openmmdlsetup import (
    app,
    configureDefaultOptions,
    _normalize_resname,
    _resnames_are_unique,
)


def test_normalize_resname_accepts_three_alnum_uppercase():
    assert _normalize_resname("unk", "UNK") == "UNK"
    assert _normalize_resname("l01", "UNK") == "L01"
    assert _normalize_resname("Ab9", "UNK") == "AB9"


def test_normalize_resname_falls_back_when_too_short():
    assert _normalize_resname("ab", "UNK") == "UNK"
    assert _normalize_resname("", "UNK") == "UNK"
    assert _normalize_resname(None, "UNK") == "UNK"


def test_normalize_resname_strips_invalid_chars_and_truncates():
    assert _normalize_resname("a_b", "UNK") == "UNK"
    assert _normalize_resname("abcd", "UNK") == "ABC"
    assert _normalize_resname("l-01", "UNK") == "L01"


def test_resnames_are_unique_true_for_unique_names():
    assert _resnames_are_unique(["UNK", "L01", "L02"]) is True


def test_resnames_are_unique_false_for_duplicates():
    assert _resnames_are_unique(["UNK", "UNK", "L02"]) is False
    assert _resnames_are_unique(["L01", "L02", "L01"]) is False

def test_configure_default_options_sets_postprocessing_defaults():
    with app.test_request_context("/"):
        session["fileType"] = "pdb"
        session["waterModel"] = "explicit"

        configureDefaultOptions()

        assert session["mda_selection"] == "mda_prot_lig_all"
        assert session["remove_postprocessing_intermediates"] == "False"

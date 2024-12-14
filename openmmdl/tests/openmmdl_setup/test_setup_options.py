import pytest
from typing import List, Dict
from flask import session
from werkzeug.datastructures import ImmutableMultiDict
from openmmdl.openmmdl_setup.setup_options import SetupOptionsConfigurator, SessionDict, RequestSessionManager

@pytest.fixture
def default_session() -> SessionDict:
    """
    Provides a default mock session with minimal required data for testing.
    """
    return {
        "fileType": "pdb",
        "waterModel": "explicit",
    }

@pytest.fixture
def configurator(default_session: SessionDict) -> SetupOptionsConfigurator:
    """
    Provides a SetupOptionsConfigurator instance initialized with the default session.
    """
    return SetupOptionsConfigurator(session=default_session)

def test_configure_default_options(configurator: SetupOptionsConfigurator, default_session: SessionDict):
    """
    Test the `configure_default_options` method to ensure default values are properly set.
    """
    configurator.configure_default_options()

    assert default_session["restart_checkpoint"] is False
    assert default_session["mdtraj_output"] == "mdtraj_pdb_dcd"
    assert default_session["mda_output"] == "mda_pdb_dcd"
    assert default_session["analysis_selection"] == "analysis_all"
    assert default_session["binding_mode"] == "40"
    assert default_session["ensemble"] == "npt"  # since waterModel is 'explicit'
    assert default_session["platform"] == "CUDA"
    assert default_session["precision"] == "mixed"
    assert default_session["cutoff"] == "1.0"  # as waterModel is not 'implicit'
    assert default_session["hmr"] is True
    assert default_session["writeDCD"] is True
    assert default_session["dataFields"] == ["step", "speed", "progress", "potentialEnergy", "temperature"]
    assert default_session["writeCheckpoint"] is True

def test_configure_default_options_with_implicit_water(configurator: SetupOptionsConfigurator, default_session: SessionDict):
    """
    Test `configure_default_options` when the water model is set to implicit.
    """
    default_session["waterModel"] = "implicit"
    configurator.configure_default_options()

    assert default_session["ensemble"] == "nvt"  # should switch to nvt due to implicit water
    assert default_session["cutoff"] == "2.0"  # cutoff should change
    assert default_session["nonbondedMethod"] == "CutoffNonPeriodic"

def test_configure_default_amber_options(configurator: SetupOptionsConfigurator, default_session: SessionDict):
    """
    Test the `configureDefaultAmberOptions` method to ensure Amber options are set correctly.
    """
    configurator.configureDefaultAmberOptions()

    assert default_session["lig_ff"] == "gaff2"
    assert default_session["charge_value"] == "0"
    assert default_session["charge_method"] == "bcc"
    assert default_session["prot_ff"] == "ff19SB"
    assert default_session["dna_ff"] == "OL15"
    assert default_session["rna_ff"] == "OL3"
    assert default_session["carbo_ff"] == "GLYCAM_06j"
    assert default_session["addType"] == "addWater"
    assert default_session["boxType"] == "cube"
    assert default_session["lipid_tp"] == "POPC"
    assert default_session["dist2Border"] == "15"
    assert default_session["water_ff"] == "opc"
    assert default_session["pos_ion"] == "Na+"
    assert default_session["neg_ion"] == "Cl-"
    assert default_session["ionConc"] == "0.15"


@pytest.fixture
def app_context():
    """
    Provides a Flask test request context for session management.
    """
    from flask import Flask

    app = Flask(__name__)
    app.secret_key = 'test_secret_key'  # Required to use sessions
    with app.test_request_context():
        yield

@pytest.fixture
def default_form() -> ImmutableMultiDict:
    """
    Provides default mock form data for testing.
    """
    return ImmutableMultiDict({
        "rcpType": "protein",
        "prot_ff": "ff14SB",
        "other_prot_ff_input": "custom_prot_ff",
        "dna_ff": "OL15",
        "other_dna_ff_input": "custom_dna_ff",
        "rna_ff": "OL3",
        "other_rna_ff_input": "custom_rna_ff",
        "carbo_ff": "GLYCAM_06j",
        "addType": "addWater",
        "boxType": "geometry",
        "geomPadding": "10",
        "ionicstrength": "0.15",
        "positiveion": "Na+",
        "negativeion": "Cl-",
        "forcefield": "amber14",
        "ml_forcefield": "openff",
        "waterModel": "tip3p",
        "smallMoleculeForceField": "gaff",
        "ligandMinimization": "Yes",
        "ligandSanitization": "Yes",
        "writeDCD": "True",
        "writeData": "True",
        "writeCheckpoint": "True",
        "dataFields": "step,speed,temperature",
        "hmr": "True",
    })

@pytest.fixture
def request_manager(default_form) -> RequestSessionManager:
    """
    Provides a RequestSessionManager instance initialized with the default form data.
    """
    return RequestSessionManager(form=default_form)

def test_set_amber_options_rcp_session(request_manager: RequestSessionManager, app_context):
    """
    Test the `setAmberOptions_rcp_session` method.
    """
    request_manager.setAmberOptions_rcp_session()

    assert session["rcpType"] == "protein"
    assert session["prot_ff"] == "ff14SB"
    assert session["other_prot_ff_input"] == "custom_prot_ff"
    assert session["dna_ff"] == "OL15"
    assert session["other_dna_ff_input"] == "custom_dna_ff"
    assert session["rna_ff"] == "OL3"
    assert session["other_rna_ff_input"] == "custom_rna_ff"
    assert session["carbo_ff"] == "GLYCAM_06j"

def test_set_amber_options_water_membrane_session(request_manager: RequestSessionManager, app_context):
    """
    Test the `setAmberOptions_water_membrane_session` method.
    """
    request_manager.setAmberOptions_water_membrane_session()

    assert session["addType"] == "addWater"
    assert session["boxType"] == "geometry"
    assert session["dist"] == ""
    assert session["lipid_tp"] == ""
    assert session["other_lipid_tp_input"] == ""
    assert session["lipid_ratio"] == ""
    assert session["lipid_ff"] == ""
    assert session["dist2Border"] == ""
    assert session["padDist"] == ""

def test_simulationoptions_add_general_settings(request_manager: RequestSessionManager, app_context):
    """
    Test the `simulationoptions_add_general_settings` method to ensure general simulation settings are correctly added to the session.
    """
    request_manager.simulationoptions_add_general_settings()

    assert session["forcefield"] == "amber14"
    assert session["ml_forcefield"] == "openff"
    assert session["waterModel"] == "tip3p"
    assert session["smallMoleculeForceField"] == "gaff"
    assert session["ligandMinimization"] == "Yes"
    assert session["ligandSanitization"] == "Yes"
    assert session["writeDCD"] is True
    assert session["writeData"] is True
    assert session["writeCheckpoint"] is True
    assert session["dataFields"] == ["step,speed,temperature"]
    assert session["hmr"] is True

def test_configure_files_add_forcefield_ligand_settings(request_manager: RequestSessionManager, app_context):
    """
    Test the `configureFiles_add_forcefield_ligand_settings` method to ensure forcefield and ligand settings are added to the session.
    """
    request_manager.configureFiles_add_forcefield_ligand_settings()

    assert session["forcefield"] == "amber14"
    assert session["ml_forcefield"] == "openff"
    assert session["waterModel"] == "tip3p"
    assert session["smallMoleculeForceField"] == "gaff"
    assert session["ligandMinimization"] == "Yes"
    assert session["ligandSanitization"] == "Yes"

def test_parse_float(request_manager: RequestSessionManager):
    """
    Test the `_parse_float` helper function.
    """
    assert request_manager._parse_float("10.5") == 10.5
    assert request_manager._parse_float(None) is None
    assert request_manager._parse_float("invalid") is None

"""Bespoke ligand parameters in the generated scripts of both setup paths."""

import ast
import io

import pytest
from flask import session

from openmmdl.openmmdl_setup import openmmdlsetup as M


@pytest.fixture
def pdb_state(monkeypatch):
    uploads = {
        "file": [(io.BytesIO(b"END\n"), "protein.pdb")],
        "sdfFile": [(io.BytesIO(b"ligand\n"), "ligand.sdf")],
        "ligandOffxml": [(io.BytesIO(b"<SMIRNOFF/>"), "custom.offxml")],
    }
    monkeypatch.setattr(M, "uploadedFiles", uploads)
    with M.app.test_request_context():
        session.update(
            fileType="pdb", pdbType="pdb", sdfFile="ligand.sdf", sdfResname="UNK",
            companionFiles=[], companionResnames=[], waterModel="TIP3P", forcefield="AMBER14",
            smallMoleculeForceField="gaff", smallMoleculeForceFieldVersion="gaff-2.11",
            ligandMinimization="False", ligandSanitization="True", solvent=True, add_membrane=False,
            water_padding=True, water_padding_distance=1.0, water_boxShape="cube",
            water_ionicstrength=0.0, water_positive="Na+", water_negative="Cl-",
        )
        M.configureDefaultOptions()
        yield uploads


def test_pdb_script_without_bespoke_parameters_is_unchanged(pdb_state):
    session["ligandParameters"] = {"mode": "standard"}
    script = M.createScript()
    ast.parse(script)
    assert "ligand_parameters" not in script and "ligand_offxml" not in script


def test_pdb_script_imports_the_uploaded_force_field(pdb_state):
    session["ligandParameters"] = {"mode": "import", "offxml": "custom.offxml"}
    script = M.createScript()
    ast.parse(script)
    assert "from openmmdl.ligand_parameters.core import FitOptions, prepare_parameters" in script
    assert "ligand_offxml = prepare_parameters(ligand, 'ligand_parameters', offxml='custom.offxml')" in script
    assert script.count("rdkit_mol=prepared_ligands, ligand_offxml=ligand_offxml)") == 2


def test_pdb_script_fits_with_the_selected_options(pdb_state):
    options = {"executable": "/env/bin/openff-bespoke", "qc_program": "xtb", "qc_method": "gfn2xtb",
               "qc_basis": "none", "workers": 1, "cores": 4, "memory_gb": 8}
    session["ligandParameters"] = {"mode": "fit", "options": options}
    script = M.createScript()
    ast.parse(script)
    assert "ligand_offxml = prepare_parameters(ligand, 'ligand_parameters', options=FitOptions(**%r))" % options in script


def test_pdb_form_reports_invalid_input_on_the_page(monkeypatch):
    monkeypatch.setattr(M, "uploadedFiles", {})
    with M.app.test_client() as client:
        with client.session_transaction() as state:
            state["fileType"] = "pdb"
        result = client.post(
            "/configureFiles",
            data={
                "file": (io.BytesIO(b"END\n"), "protein.pdb"),
                "sdfFile": (io.BytesIO(b"not a molecule\n"), "ligand.sdf"),
                "forcefield": "AMBER14", "smallMoleculeMode": "single", "ligandParameterMode": "import",
            },
        )
        assert result.status_code == 200
        assert b"alert-danger" in result.data and b"Could not read the ligand" in result.data


def _amber_state(uploads, ligand_parameters):
    uploads["protFile"] = uploads.pop("file")
    uploads["nmLigFile"] = uploads.pop("sdfFile")
    M.configureDefaultAmberOptions()
    session.update(
        fileType="amber", has_files="no", rcpType="protRcp", nmLig=True, spLig=False,
        addType="addWater", boxType="cube", water_ff="tip3p", ligandParameters=ligand_parameters,
    )


def test_amber_script_prepares_bespoke_ligand_instead_of_antechamber(pdb_state):
    _amber_state(pdb_state, {"mode": "import", "offxml": "custom.offxml"})
    script = M.createAmberBashScript()
    assert "openmmdl ligand import --ligand ligand.sdf --output bespoke_ligand --target amber --offxml custom.offxml || exit 1" in script
    assert "antechamber -" not in script and "parmchk2 -" not in script
    assert script.count("loadoff bespoke_ligand/ligand.lib") == 2
    assert "nmLig = loadpdb bespoke_ligand/ligand.pdb" in script
    assert script.index("openmmdl ligand import") < script.index("pdb4amber")
    assert script.rstrip().endswith("--ligand-dir bespoke_ligand || exit 1")
    simulation = M.createScript()
    ast.parse(simulation)
    assert "ligand_name = 'BSP'" in simulation


def test_amber_script_keeps_gaff_workflow_for_standard_ligand(pdb_state):
    _amber_state(pdb_state, {"mode": "standard"})
    script = M.createAmberBashScript()
    assert "antechamber" in script and "bespoke_ligand" not in script


def test_amber_script_shows_error_until_input_is_fixed(pdb_state):
    _amber_state(pdb_state, {"mode": "error", "error": "Upload the SMIRNOFF force field (.offxml) to import"})
    assert M.createAmberBashScript() == "# Upload the SMIRNOFF force field (.offxml) to import\n"

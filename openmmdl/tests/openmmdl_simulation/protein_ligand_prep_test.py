import pytest
import os
import rdkit
from rdkit import Chem
import simtk.openmm.app as app
from simtk.openmm.app import PDBFile, Modeller
from simtk.openmm import unit
from simtk.openmm import Vec3
import mdtraj as md
import numpy as np
import simtk
from pathlib import Path
import pdbfixer
from openmm.app import PDBFile
from pdbfixer import PDBFixer
import simtk.openmm.app as app


from simtk.openmm.app import (
    PDBFile,
    Modeller,
    PDBReporter,
    StateDataReporter,
    DCDReporter,
    CheckpointReporter,
)
from simtk.openmm import unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator
from simtk.openmm import Vec3
import simtk.openmm as mm


from openmmdl.openmmdl_simulation.scripts.forcefield_water import (
    ff_selection,
    water_forcefield_selection,
    water_model_selection,
    generate_forcefield,
    generate_transitional_forcefield,
)
from openmmdl.openmmdl_simulation.scripts.protein_ligand_prep import (
    prepare_ligand,
    rdkit_to_openmm,
    merge_protein_and_ligand,
    water_padding_solvent_builder,
    water_absolute_solvent_builder,
    write_ligand_with_partial_charges,
    membrane_builder,
    water_conversion,
)
from openmmdl.openmmdl_simulation.scripts.post_md_conversions import (
    mdtraj_conversion,
    MDanalysis_conversion,
)


protein = "6b73.pdb"
ligand = "CVV.sdf"
ligand_name = "UNK"
minimization = False
sanitization = False
ff = "AMBER14"
water = "SPC/E"
add_membrane = False
Water_Box = "Buffer"
water_padding_distance = 1.0
water_boxShape = "cube"
water_ionicstrength = 0.15
water_positive_ion = "Na+"
water_negative_ion = "Cl-"

water_box_x = 6.873
water_box_y = 7.0
water_box_z = 9.132

# Print current working directory
print("Current working directory:", os.getcwd())

# Assuming that 'test_data_directory' is properly defined in your test setup
test_data_directory = "openmmdl/tests/data/in"


test_data_directory = Path("openmmdl/tests/data/in")


# Define the full path to the input SDF file
TEST_LIGAND_FILE = f"{test_data_directory}/CVV.sdf"
TEST_MOL_FILE = f"{test_data_directory}/CVV.mol"
TEST_MOL2_FILE = f"{test_data_directory}/CVV.mol2"
TEST_PROTEIN = f"{test_data_directory}/6b73.pdb"

protein_pdb = pdbfixer.PDBFixer(str(TEST_PROTEIN))


ligand_prepared = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=minimization)
omm_ligand = rdkit_to_openmm(ligand_prepared, ligand_name)
forcefield_selected = ff_selection(ff)
water_selected = water_forcefield_selection(
    water=water, forcefield_selection=ff_selection(ff)
)
model_water = water_model_selection(water=water, forcefield_selection=ff_selection(ff))
forcefield = generate_forcefield(
    protein_ff=forcefield_selected,
    solvent_ff=water_selected,
    add_membrane=add_membrane,
    rdkit_mol=ligand_prepared,
    smallMoleculeForceField="gaff",
    smallMoleculeForceFieldVersion="gaff-2.1"
)
complex_topology, complex_positions = merge_protein_and_ligand(protein_pdb, omm_ligand)
modeller = app.Modeller(complex_topology, complex_positions)


@pytest.fixture
def ligand_charge_test_system():
    system = forcefield.createSystem(complex_topology)
    return complex_topology, system, complex_positions

# Test the prepare_ligand function
def test_prepare_ligand():
    # Test the function with the sample ligand file.
    rdkit_mol_sdf = prepare_ligand(TEST_LIGAND_FILE, minimize_molecule=False)
    rdkit_mol_mol2_2 = prepare_ligand(TEST_MOL2_FILE, minimize_molecule=True)
    rdkit_mol_mol = prepare_ligand(TEST_MOL_FILE, minimize_molecule=False)
    rdkit_mol_mol2 = prepare_ligand(TEST_MOL2_FILE, minimize_molecule=False)

    # Add your assertions here to check if the preparation worked as expected
    assert rdkit_mol_sdf is not None  # Check if the result is not None
    assert rdkit_mol_mol2_2 is not None  # Check if the result is not None
    assert rdkit_mol_mol is not None  # Check if the result is not None
    assert rdkit_mol_mol2 is not None  # Check if the result is not None


def test_rdkit_to_openmm():
    omm_ligand = rdkit_to_openmm(ligand_prepared, ligand_name)
    assert isinstance(omm_ligand, simtk.openmm.app.Modeller)


def test_merge_protein_and_ligand():
    complex_topology, complex_positions = merge_protein_and_ligand(
        protein_pdb, omm_ligand
    )
    assert complex_topology is not None
    assert complex_positions is not None

def test_write_ligand_with_partial_charges(tmp_path):
    system = forcefield.createSystem(complex_topology)

    output_file = tmp_path / f"{ligand_name}_pc.mol2"

    output = write_ligand_with_partial_charges(
        complex_topology,
        system,
        complex_positions,
        ligand_name=ligand_name,
        output_file=str(output_file),
    )

    assert output == str(output_file)
    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_write_ligand_with_partial_charges_without_name():
    system = forcefield.createSystem(complex_topology)

    output = write_ligand_with_partial_charges(
        complex_topology,
        system,
        complex_positions,
        ligand_name=None,
    )

    assert output is None

def test_write_ligand_with_partial_charges_single_returns_string(
    ligand_charge_test_system, tmp_path
):
    topology, system, positions = ligand_charge_test_system
    output_file = tmp_path / "UNK_pc.mol2"

    output = write_ligand_with_partial_charges(
        topology,
        system,
        positions,
        ligand_name="UNK",
        ligand_names=["UNK"],
        output_file=str(output_file),
    )

    assert output == str(output_file)

def test_write_ligand_with_partial_charges_without_name_returns_none(
    ligand_charge_test_system,
):
    topology, system, positions = ligand_charge_test_system

    output = write_ligand_with_partial_charges(
        topology,
        system,
        positions,
        ligand_name=None,
        ligand_names=[],
    )

    assert output is None

def test_water_padding_solvent_builder():
    protein_buffer_solved = water_padding_solvent_builder(
        model_water,
        forcefield,
        water_padding_distance,
        protein_pdb,
        modeller,
        water_positive_ion,
        water_negative_ion,
        water_ionicstrength,
        protein,
    )
    assert protein_buffer_solved is not None


def test_water_absolute_solvent_builder():
    test_data_directory = Path("openmmdl/tests/data/in")
    TEST_PROTEIN = f"{test_data_directory}/6b73.pdb"
    protein_pdb = pdbfixer.PDBFixer(str(TEST_PROTEIN))
    protein_absolute_solved = water_absolute_solvent_builder(
        model_water,
        forcefield,
        water_box_x,
        water_box_y,
        water_box_z,
        protein_pdb,
        modeller,
        water_positive_ion,
        water_negative_ion,
        water_ionicstrength,
        protein,
    )
    assert protein_absolute_solved is not None

def test_water_padding_solvent_builder_opc():
    opc_ff = "AMBER19"
    opc_water = "OPC"
    opc_forcefield_selected = ff_selection(opc_ff)
    opc_water_selected = water_forcefield_selection(
        water=opc_water,
        forcefield_selection=opc_forcefield_selected,
    )
    opc_model_water = water_model_selection(
        water=opc_water,
        forcefield_selection=opc_forcefield_selected,
    )
    opc_forcefield = generate_forcefield(
        protein_ff=opc_forcefield_selected,
        solvent_ff=opc_water_selected,
        add_membrane=False,
        rdkit_mol=ligand_prepared,
        smallMoleculeForceField="gaff",
        smallMoleculeForceFieldVersion="gaff-2.1",
    )
    opc_modeller = app.Modeller(complex_topology, complex_positions)

    protein_buffer_solved = water_padding_solvent_builder(
        opc_model_water,
        opc_forcefield,
        water_padding_distance,
        protein_pdb,
        opc_modeller,
        water_positive_ion,
        water_negative_ion,
        water_ionicstrength,
        protein,
    )

    assert protein_buffer_solved is not None


def test_water_padding_solvent_builder_opc3():
    opc3_ff = "AMBER19"
    opc3_water = "OPC3"
    opc3_forcefield_selected = ff_selection(opc3_ff)
    opc3_water_selected = water_forcefield_selection(
        water=opc3_water,
        forcefield_selection=opc3_forcefield_selected,
    )
    opc3_model_water = water_model_selection(
        water=opc3_water,
        forcefield_selection=opc3_forcefield_selected,
    )
    opc3_forcefield = generate_forcefield(
        protein_ff=opc3_forcefield_selected,
        solvent_ff=opc3_water_selected,
        add_membrane=False,
        rdkit_mol=ligand_prepared,
        smallMoleculeForceField="gaff",
        smallMoleculeForceFieldVersion="gaff-2.1",
    )
    opc3_modeller = app.Modeller(complex_topology, complex_positions)

    protein_buffer_solved = water_padding_solvent_builder(
        opc3_model_water,
        opc3_forcefield,
        water_padding_distance,
        protein_pdb,
        opc3_modeller,
        water_positive_ion,
        water_negative_ion,
        water_ionicstrength,
        protein,
    )

    assert protein_buffer_solved is not None


def test_water_absolute_solvent_builder_opc():
    opc_ff = "AMBER19"
    opc_water = "OPC"
    opc_forcefield_selected = ff_selection(opc_ff)
    opc_water_selected = water_forcefield_selection(
        water=opc_water,
        forcefield_selection=opc_forcefield_selected,
    )
    opc_model_water = water_model_selection(
        water=opc_water,
        forcefield_selection=opc_forcefield_selected,
    )
    opc_forcefield = generate_forcefield(
        protein_ff=opc_forcefield_selected,
        solvent_ff=opc_water_selected,
        add_membrane=False,
        rdkit_mol=ligand_prepared,
        smallMoleculeForceField="gaff",
        smallMoleculeForceFieldVersion="gaff-2.1",
    )
    opc_modeller = app.Modeller(complex_topology, complex_positions)

    protein_absolute_solved = water_absolute_solvent_builder(
        opc_model_water,
        opc_forcefield,
        water_box_x,
        water_box_y,
        water_box_z,
        protein_pdb,
        opc_modeller,
        water_positive_ion,
        water_negative_ion,
        water_ionicstrength,
        protein,
    )

    assert protein_absolute_solved is not None


def test_water_absolute_solvent_builder_opc3():
    opc3_ff = "AMBER19"
    opc3_water = "OPC3"
    opc3_forcefield_selected = ff_selection(opc3_ff)
    opc3_water_selected = water_forcefield_selection(
        water=opc3_water,
        forcefield_selection=opc3_forcefield_selected,
    )
    opc3_model_water = water_model_selection(
        water=opc3_water,
        forcefield_selection=opc3_forcefield_selected,
    )
    opc3_forcefield = generate_forcefield(
        protein_ff=opc3_forcefield_selected,
        solvent_ff=opc3_water_selected,
        add_membrane=False,
        rdkit_mol=ligand_prepared,
        smallMoleculeForceField="gaff",
        smallMoleculeForceFieldVersion="gaff-2.1",
    )
    opc3_modeller = app.Modeller(complex_topology, complex_positions)

    protein_absolute_solved = water_absolute_solvent_builder(
        opc3_model_water,
        opc3_forcefield,
        water_box_x,
        water_box_y,
        water_box_z,
        protein_pdb,
        opc3_modeller,
        water_positive_ion,
        water_negative_ion,
        water_ionicstrength,
        protein,
    )

    assert protein_absolute_solved is not None

class DummyForceField:
    def __init__(self, name):
        self.name = name


class DummyProteinPDB:
    topology = object()
    positions = object()


class DummyModeller:
    topology = object()
    positions = object()

    def __init__(self):
        self.calls = []

    def addMembrane(self, forcefield, **kwargs):
        self.calls.append(("addMembrane", forcefield.name))

    def convertWater(self, model_water):
        self.calls.append(("convertWater", model_water))

    def addExtraParticles(self, forcefield):
        self.calls.append(("addExtraParticles", forcefield.name))


def _run_membrane_builder(monkeypatch, tmp_path, model_water):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(membrane_builder.__globals__["PDBFile"], "writeFile", lambda *args, **kwargs: None)

    modeller = DummyModeller()
    forcefield = DummyForceField("final")
    transitional_forcefield = DummyForceField("transitional")

    returned = membrane_builder(
        ff="AMBER99SB",
        model_water=model_water,
        forcefield=forcefield,
        transitional_forcefield=transitional_forcefield,
        protein_pdb=DummyProteinPDB(),
        modeller=modeller,
        membrane_lipid_type="POPC",
        membrane_padding=1.0,
        membrane_positive_ion="Na+",
        membrane_negative_ion="Cl-",
        membrane_ionicstrength=0.15,
        protein_name="protein.pdb",
    )

    assert returned is modeller
    return modeller.calls


def test_membrane_builder_converts_tip4pew_after_transitional_build(monkeypatch, tmp_path):
    calls = _run_membrane_builder(monkeypatch, tmp_path, "tip4pew")

    assert calls == [
        ("addMembrane", "transitional"),
        ("convertWater", "tip4pew"),
    ]


def test_membrane_builder_adds_tip5p_extra_particles_after_transitional_build(monkeypatch, tmp_path):
    calls = _run_membrane_builder(monkeypatch, tmp_path, "tip5p")

    assert calls == [
        ("addMembrane", "transitional"),
        ("addExtraParticles", "final"),
    ]


def test_membrane_builder_adds_opc_extra_particles_after_transitional_build(monkeypatch, tmp_path):
    calls = _run_membrane_builder(monkeypatch, tmp_path, "opc")

    assert calls == [
        ("addMembrane", "transitional"),
        ("addExtraParticles", "final"),
    ]


if __name__ == "__main__":
    pytest.main()

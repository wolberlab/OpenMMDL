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


if __name__ == "__main__":
    pytest.main()

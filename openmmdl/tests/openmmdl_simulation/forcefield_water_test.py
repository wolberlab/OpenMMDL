import pytest
import simtk.openmm.app as app
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from openmmdl.openmmdl_simulation.scripts.forcefield_water import (
    ff_selection,
    water_forcefield_selection,
    water_model_selection,
    generate_forcefield,
    generate_transitional_forcefield,
)

# Replace 'your_module' with the actual name of the module containing your functions.


@pytest.fixture
def sample_rdkit_molecule():
    """
    A sample RDKit molecule for testing.
    """
    from rdkit import Chem

    mol = Chem.MolFromSmiles("CCO")
    return mol


def test_ff_selection():
    assert ff_selection("AMBER14") == "amber14-all.xml"
    assert ff_selection("AMBER99SB") == "amber99sb.xml"
    assert ff_selection("AMBER99SB-ILDN") == "amber99sbildn.xml"
    assert ff_selection("AMBER03") == "amber03.xml"
    assert ff_selection("AMBER10") == "amber10.xml"
    assert ff_selection("CHARMM36") == "charmm36.xml"
    assert ff_selection("NonexistentFF") is None


def test_water_forcefield_selection():
    # Test cases for 'amber14-all.xml' force field
    assert water_forcefield_selection("TIP3P", "amber14-all.xml") == "amber14/tip3p.xml"
    assert (
        water_forcefield_selection("TIP3P-FB", "amber14-all.xml")
        == "amber14/tip3pfb.xml"
    )
    assert water_forcefield_selection("SPC/E", "amber14-all.xml") == "amber14/spce.xml"
    assert (
        water_forcefield_selection("TIP4P-Ew", "amber14-all.xml")
        == "amber14/tip4pew.xml"
    )
    assert (
        water_forcefield_selection("TIP4P-FB", "amber14-all.xml")
        == "amber14/tip4pfb.xml"
    )
    assert water_forcefield_selection("TIP5P", "amber14-all.xml") is None
    assert water_forcefield_selection("NonexistentWater", "amber14-all.xml") is None
    assert water_forcefield_selection("TIP3P", "NonexistentFF") is None

    # Test cases for 'charmm36.xml' force field
    assert (
        water_forcefield_selection("CHARMM default", "charmm36.xml")
        == "charmm36/water.xml"
    )
    assert (
        water_forcefield_selection("TIP3P-PME-B", "charmm36.xml")
        == "charmm36/tip3p-pme-b.xml"
    )
    assert (
        water_forcefield_selection("TIP3P-PME-F", "charmm36.xml")
        == "charmm36/tip3p-pme-f.xml"
    )
    assert water_forcefield_selection("SPC/E", "charmm36.xml") == "charmm36/spce.xml"
    assert (
        water_forcefield_selection("TIP4P-Ew", "charmm36.xml") == "charmm36/tip4pew.xml"
    )
    assert (
        water_forcefield_selection("TIP4P-2005", "charmm36.xml")
        == "charmm36/tip4p2005.xml"
    )
    assert water_forcefield_selection("TIP5P", "charmm36.xml") == "charmm36/tip5p.xml"
    assert (
        water_forcefield_selection("TIP5P-Ew", "charmm36.xml") == "charmm36/tip5pew.xml"
    )
    assert water_forcefield_selection("NonexistentWater", "charmm36.xml") is None
    assert water_forcefield_selection("NonexistentFF", "charmm36.xml") is None


def test_water_model_selection():
    assert water_model_selection("TIP3P", "amber99sb.xml") == "tip3p"
    assert water_model_selection("TIP3P", "amber99sbildn.xml") == "tip3p"
    assert water_model_selection("TIP3P", "amber03.xml") == "tip3p"
    assert water_model_selection("TIP3P", "amber10.xml") == "tip3p"

    assert water_model_selection("SPC/E", "amber99sb.xml") == "spce"
    assert water_model_selection("SPC/E", "amber99sbildn.xml") == "spce"
    assert water_model_selection("SPC/E", "amber03.xml") == "spce"
    assert water_model_selection("SPC/E", "amber10.xml") == "spce"

    assert water_model_selection("TIP4P-Ew", "amber99sb.xml") == "tip4pew"
    assert water_model_selection("TIP4P-Ew", "amber99sbildn.xml") == "tip4pew"
    assert water_model_selection("TIP4P-Ew", "amber03.xml") == "tip4pew"
    assert water_model_selection("TIP4P-Ew", "amber10.xml") == "tip4pew"

    assert water_model_selection("TIP4P-FB", "amber99sb.xml") == "tip4pfb"
    assert water_model_selection("TIP4P-FB", "amber99sbildn.xml") == "tip4pfb"
    assert water_model_selection("TIP4P-FB", "amber03.xml") == "tip4pfb"
    assert water_model_selection("TIP4P-FB", "amber10.xml") == "tip4pfb"

    assert water_model_selection("TIP5P", "amber99sb.xml") is None
    assert water_model_selection("TIP5P", "amber99sbildn.xml") is None
    assert water_model_selection("TIP5P", "amber03.xml") is None
    assert water_model_selection("TIP5P", "amber10.xml") is None
    assert (
        water_model_selection("TIP5P", "amber14-all.xml") is None
    )  # Missing in the initial version

    assert water_model_selection("TIP3P", "amber14-all.xml") == "tip3p"

    assert water_model_selection("CHARMM default", "charmm36.xml") == "charmm"
    assert water_model_selection("TIP3P-PME-B", "charmm36.xml") == "charmm"
    assert water_model_selection("TIP3P-PME-F", "charmm36.xml") == "charmm"
    assert water_model_selection("SPC/E", "charmm36.xml") == "charmm"
    assert water_model_selection("TIP4P-Ew", "charmm36.xml") == "tip4pew"
    assert water_model_selection("TIP4P-2005", "charmm36.xml") == "tip4pew"
    assert water_model_selection("TIP5P", "charmm36.xml") == "tip5p"
    assert water_model_selection("TIP5P-Ew", "charmm36.xml") == "tip5p"

    assert water_model_selection("TIP3P", "NonexistentFF") is None


def test_generate_forcefield_with_membrane(sample_rdkit_molecule):
    forcefield = generate_forcefield(
        "amber14-all.xml", "amber14/tip3p.xml", True, "Gaff", sample_rdkit_molecule
    )
    assert isinstance(forcefield, app.ForceField)
    # Add additional assertions specific to the case with a membrane


def test_generate_forcefield_without_membrane(sample_rdkit_molecule):
    forcefield = generate_forcefield(
        "amber14-all.xml", "amber14/tip3p.xml", False, "Gaff", sample_rdkit_molecule
    )
    assert isinstance(forcefield, app.ForceField)
    # Add additional assertions specific to the case without a membrane


def test_generate_forcefield_with_old_amber_forcefield(sample_rdkit_molecule):
    forcefield = generate_forcefield(
        "amber99sb.xml", "amber14/tip3p.xml", True, "Gaff", sample_rdkit_molecule
    )
    assert isinstance(forcefield, app.ForceField)
    # Add additional assertions specific to the case with an old Amber forcefield


def test_generate_forcefield_without_small_molecule():
    forcefield = generate_forcefield("amber14-all.xml", "amber14/tip3p.xml", False)
    assert isinstance(forcefield, app.ForceField)
    # Add additional assertions specific to the case without a small molecule


def test_generate_forcefield_membrane_logic(sample_rdkit_molecule):
    forcefield_1 = generate_forcefield(
        "amber10.xml", "tip3p.xml", True, "Gaff", sample_rdkit_molecule
    )
    forcefield_2 = generate_forcefield(
        "amber14-all.xml", "amber14/tip3p.xml", True, "Gaff", sample_rdkit_molecule
    )
    forcefield_3 = generate_forcefield(
        "amber14-all.xml", "amber14/tip3p.xml", False, "Gaff", sample_rdkit_molecule
    )
    forcefield_4 = generate_forcefield(
        "amber03.xml", "tip3p.xml", False, "Gaff", sample_rdkit_molecule
    )

    assert isinstance(forcefield_1, app.ForceField)
    assert isinstance(forcefield_2, app.ForceField)
    assert isinstance(forcefield_3, app.ForceField)
    assert isinstance(forcefield_4, app.ForceField)

    # Additional tests for different force field combinations
    forcefield_5 = generate_forcefield(
        "amber14-all.xml", "tip3p.xml", True, "Gaff", sample_rdkit_molecule
    )
    forcefield_6 = generate_forcefield(
        "amber03.xml", "amber14/tip3p.xml", False, "Gaff", sample_rdkit_molecule
    )

    assert isinstance(forcefield_5, app.ForceField)
    assert isinstance(forcefield_6, app.ForceField)

    # Additional tests for membrane flag logic
    forcefield_7 = generate_forcefield(
        "amber10.xml", "tip3p.xml", True, "Gaff", sample_rdkit_molecule
    )
    forcefield_8 = generate_forcefield(
        "amber14-all.xml", "tip3p.xml", False, "Gaff", sample_rdkit_molecule
    )

    assert isinstance(forcefield_7, app.ForceField)
    assert isinstance(forcefield_8, app.ForceField)


def test_generate_transitional_forcefield(sample_rdkit_molecule):
    transitional_forcefield = generate_transitional_forcefield(
        "amber14-all.xml", "tip3p.xml", True, "Gaff", sample_rdkit_molecule
    )
    assert isinstance(transitional_forcefield, app.ForceField)

    # Additional tests for different force field combinations
    transitional_forcefield_2 = generate_transitional_forcefield(
        "amber03.xml", "amber14/tip3p.xml", False, "Gaff", sample_rdkit_molecule
    )
    assert isinstance(transitional_forcefield_2, app.ForceField)

    # Additional tests for membrane flag logic
    transitional_forcefield_3 = generate_transitional_forcefield(
        "amber14-all.xml", "tip3p.xml", False, "Gaff", sample_rdkit_molecule
    )
    assert isinstance(transitional_forcefield_3, app.ForceField)

    # Additional tests for GAFF registration
    transitional_forcefield_4 = generate_transitional_forcefield(
        "amber14-all.xml", "tip3p.xml", True
    )
    assert isinstance(transitional_forcefield_4, app.ForceField)

import pytest
import simtk.openmm.app as app
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from openmmdl.openmmdl_simulation.scripts.forcefield_water import ff_selection, water_forcefield_selection, water_model_selection, generate_forcefield, generate_transitional_forcefield

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
    assert ff_selection('AMBER14') == 'amber14-all.xml'
    assert ff_selection('CHARMM36') == 'charmm36.xml'
    assert ff_selection('NonexistentFF') is None

def test_water_forcefield_selection():
    assert water_forcefield_selection('TIP3P', 'amber14-all.xml') == 'amber14/tip3p.xml'
    assert water_forcefield_selection('SPC/E', 'charmm36.xml') == 'charmm36/spce.xml'
    assert water_forcefield_selection('TIP3P', 'NonexistentFF') is None

def test_water_model_selection():
    assert water_model_selection('TIP3P', 'amber14-all.xml') == 'tip3p'
    assert water_model_selection('SPC/E', 'amber10.xml') == 'spce'
    assert water_model_selection('TIP3P-PME-B', 'charmm36.xml') == 'charmm'
    assert water_model_selection('SPC/E', 'amber14-all.xml') == 'spce'
    assert water_model_selection('TIP3P', 'NonexistentFF') is None

def test_generate_forcefield(sample_rdkit_molecule):
    forcefield = generate_forcefield('amber14-all.xml', 'amber14/tip3p.xml', False, sample_rdkit_molecule)
    assert isinstance(forcefield, app.ForceField)


def test_generate_forcefield_membrane_logic(sample_rdkit_molecule):
    forcefield_1 = generate_forcefield('amber10.xml', 'tip3p.xml', True, sample_rdkit_molecule)
    forcefield_2 = generate_forcefield('amber14-all.xml', 'amber14/tip3p.xml', True, sample_rdkit_molecule)
    forcefield_3 = generate_forcefield('amber14-all.xml', 'amber14/tip3p.xml', False, sample_rdkit_molecule)
    forcefield_4 = generate_forcefield('amber03.xml', 'tip3p.xml', False, sample_rdkit_molecule)
    assert isinstance(forcefield_1, app.ForceField)
    assert isinstance(forcefield_2, app.ForceField)
    assert isinstance(forcefield_3, app.ForceField)
    assert isinstance(forcefield_4, app.ForceField)

def test_generate_transitional_forcefield(sample_rdkit_molecule):
    transitional_forcefield = generate_transitional_forcefield('amber14-all.xml', 'tip3p.xml', True, sample_rdkit_molecule)
    assert isinstance(transitional_forcefield, app.ForceField)

import simtk.openmm.app as app
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator


def ff_selection(ff):
    """
    Selects the required xml forcefield file.
    Parameters
    ----------
    ff : str
        Input of forcefield
        
    Returns
    -------
    forcefield_selection: str
        Selected xml forcefield file.
    """
    forcefield_dict = {
        'AMBER14': 'amber14-all.xml',
        'AMBER99SB': 'amber99sb.xml',
        'AMBER99SB-ILDN': 'amber99sbildn.xml',
        'AMBER03': 'amber03.xml',
        'AMBER10': 'amber10.xml',
        'CHARMM36': 'charmm36.xml',
    }

    return forcefield_dict.get(ff, None)

def water_forcefield_selection(water, forcefield_selection):
    """
    Select the XML filename for water force field parameters based on the chosen force field and water model.

    Parameters
    ----------
    water : str 
        The chosen water model.
    forcefield_selection : str 
        The selected force field.

    Returns
    -------
    water_model : str
        The XML filename of the water forcefield.
    """
    old_amber = {'amber99sb.xml', 'amber99sbildn.xml', 'amber03.xml', 'amber10.xml'}

    # Define a dictionary to map water models
    water_model_mapping = {
        'TIP3P': 'tip3p.xml',
        'TIP3P-FB': 'tip3pfb.xml',
        'SPC/E': 'spce.xml',
        'TIP4P-Ew': 'tip4pew.xml',
        'TIP4P-FB': 'tip4pfb.xml',
        'TIP5P': 'tip5p.xml'
    }

    if forcefield_selection in old_amber:
        water_model = water_model_mapping.get(water, None)
    else:
        water_forcefields = {
            'amber14-all.xml': {
                'TIP3P': 'amber14/tip3p.xml',
                'TIP3P-FB': 'amber14/tip3pfb.xml',
                'SPC/E': 'amber14/spce.xml',
                'TIP4P-Ew': 'amber14/tip4pew.xml',
                'TIP4P-FB': 'amber14/tip4pfb.xml',
            },
            'charmm36.xml': {
                'CHARMM default': 'charmm36/water.xml',
                'TIP3P-PME-B': 'charmm36/tip3p-pme-b.xml',
                'TIP3P-PME-F': 'charmm36/tip3p-pme-f.xml',
                'SPC/E': 'charmm36/spce.xml',
                'TIP4P-Ew': 'charmm36/tip4pew.xml',
                'TIP4P-2005': 'charmm36/tip4p2005.xml',
                'TIP5P': 'charmm36/tip5p.xml',
                'TIP5P-Ew': 'charmm36/tip5pew.xml',
            }
        }
        water_model = water_forcefields.get(forcefield_selection, {}).get(water, None)
        
    return water_model

def water_model_selection(water, forcefield_selection):
    """
    Selects the required Water model forcefield .xml file according to water selection and previous force field selection.
    Parameters
    ----------
    water : str
    	Input of water model
    forcefield_selection : str
    	Input of selected forcefield .xml File
    Returns
    -------
    water_model: str
    	Water model forcefield .xml File.
    """
    old_amber = {'amber99sb.xml', 'amber99sbildn.xml', 'amber03.xml', 'amber10.xml'}
    
    if forcefield_selection in old_amber:
        if water == 'TIP3P':
            water_model = 'tip3p'
        elif water == 'TIP3P-FB':
            water_model = 'tip3pfb'
        elif water == "SPC/E":
            water_model = "spce"
        elif water == "TIP4P-Ew":
            water_model = "tip4pew"
        elif water == "TIP4P-FB":
            water_model = "tip4pfb"
        elif water == "TIP5P":
            water_model = "tip5p"
        else:
            return None
    
    elif forcefield_selection == 'amber14-all.xml':
        if water == 'TIP3P':
            water_model = 'tip3p'
        elif water == "TIP3P-FB":
            water_model = 'tip3pfb'
        elif water == "SPC/E":
            water_model = "spce"
        elif water == "TIP4P-Ew":
            water_model = "tip4pew"
        elif water == "TIP4P-FB":
            water_model = "tip4pfb"
        else:
            return None
    
    elif forcefield_selection == 'charmm36.xml':
        if water == 'CHARMM default':
            water_model = "charmm"
        elif water == "TIP3P-PME-B":
            water_model = "charmm"
        elif water == "TIP3P-PME-F":
            water_model = "charmm"
        elif water == "SPC/E":
            water_model = "charmm"
        elif water == "TIP4P-Ew":
            water_model = "tip4pew"
        elif water == "TIP4P-2005":
            water_model = "tip4pew"
        elif water == "TIP5P":
            water_model = "tip5p"
        elif water == "TIP5P-Ew":
            water_model = "tip5p"
        else:
            return None
    else:
        return None

    return water_model
    
def generate_forcefield(protein_ff, solvent_ff, add_membrane, rdkit_mol=None):
    """
    Generate an OpenMM Forcefield object and register a small molecule.

    Parameters
    ----------
    protein_ff: str
        Input of selected forcefield .xml File
    solvent_ff: str
        Input of selected water model forcefield .xml File
    add_membrane = bool
    	Selection if the system should be built with a membrane
    rdkit_mol: rdkit.Chem.rdchem.Mol
        Small molecule to register in the force field.
    Returns
    -------
    forcefield: simtk.openmm.app.Forcefield
        Forcefield with a registered small molecule.
    """
    old_amber = {'amber99sb.xml', 'amber99sbildn.xml', 'amber03.xml', 'amber10.xml'}
    
    # For older amber forcefields the additional lipid17.xml is required for templates
    if add_membrane == True:
        if protein_ff in old_amber:
            forcefield  = app.ForceField(protein_ff, solvent_ff, 'amber14/lipid17.xml')
        else:
            forcefield = app.ForceField(protein_ff,solvent_ff)
    else:
        forcefield = app.ForceField(protein_ff,solvent_ff)
    
    # If a ligand is present, a Forcefield with GAFF will be created for the ligand 
    if rdkit_mol is not None:
        gaff = GAFFTemplateGenerator(
            molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True), forcefield = 'gaff-2.11'
        )
        forcefield.registerTemplateGenerator(gaff.generator)

    return forcefield
    
def generate_transitional_forcefield(protein_ff, solvent_ff, add_membrane, rdkit_mol=None):
    """
    Generate an OpenMM transitional forcefield object with TIP3P water model for membrane building and register a small molecule.

    Parameters
    ----------
    protein_ff: str
        Name of the force field in xml format.
    solvent_ff: str
        Name of the water model force field in xml format.
    add_membrane = bool
    	Selection if the system should be build with a membrane
    rdkit_mol: rdkit.Chem.rdchem.Mol
        Small molecule to register in the force field.    
    Returns
    -------
    transitional_forcefield: simtk.openmm.app.Forcefield
        A transitional forcefield with tip3p water and a registered small molecule.
    """
    old_amber = {'amber99sb.xml', 'amber99sbildn.xml', 'amber03.xml', 'amber10.xml'}
    
    # For older amber forcefields the additional lipid17.xml is required for templates
    if add_membrane == True:
        if protein_ff in old_amber:
            transitional_forcefield  = app.ForceField(protein_ff, 'tip3p.xml', 'amber14/lipid17.xml')
        else:
            transitional_forcefield = app.ForceField(protein_ff,'amber14/tip3p.xml')
    else:
        transitional_forcefield = app.ForceField(protein_ff,solvent_ff)
    
    # If a ligand is present, a Forcefield with GAFF will be created for the ligand 
    if rdkit_mol is not None:
        gaff = GAFFTemplateGenerator(
            molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True), forcefield = 'gaff-2.11'
        )
        transitional_forcefield.registerTemplateGenerator(gaff.generator)

    return transitional_forcefield

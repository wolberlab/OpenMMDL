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
    if ff == 'AMBER14':
        forcefield_selection = 'amber14-all.xml'
    elif ff == 'AMBER99SB':
        forcefield_selection = 'amber99sb.xml'
    elif ff == 'AMBER99SB-ILDN':
        forcefield_selection = 'amber99sbildn.xml'
    elif ff == 'AMBER03':
        forcefield_selection = 'amber03.xml'
    elif ff == 'AMBER10':
        forcefield_selection = 'amber10.xml'
    elif ff == 'CHARMM36':
        forcefield_selection = 'charmm36.xml'
    else:
        return None

    return forcefield_selection

def water_forecfield_selection(water, forcefield_selection):
    """
    Selects the required Water model forcefield .xml file according to water selection and previous force field selection.
    Parameters
    ----------
    water : str
    	Input of water model
    force_selection : str
    	Input of Selected xml forcefield file
    Returns
    -------
    water_forcefield: str
    	Water model forcefield .xml File.
    """    
    old_amber = {'amber99sb.xml', 'amber99sbildn.xml', 'amber03.xml', 'amber10.xml'}
    
    if forcefield_selection in old_amber:
        if water == 'TIP3P':
            water_forcefield = "tip3p.xml"
        elif water == 'TIP3P-FB':
            water_forcefield = "tip3pfb.xml"
        elif water == "SPC/E":
            water_forcefield = "spce.xml"
        elif water == "TIP4P-Ew":
            water_forcefield = "tip4pew.xml"
        elif water == "TIP4P-FB":
            water_forcefield = "tip4pfb.xml"
        elif water == "TIP5P":
            water_forcefield = "tip5p.xml"
        else:
            return None
            
    elif forcefield_selection == 'amber14-all.xml':
        if water == 'TIP3P':
            water_forcefield = "amber14/tip3p.xml"
        elif water == "TIP3P-FB":
            water_forcefield = "amber14/tip3pfb.xml"
        elif water == "SPC/E":
            water_forcefield = "amber14/spce.xml"
        elif water == "TIP4P-Ew":
            water_forcefield = "amber14/tip4pew.xml"
        elif water == "TIP4P-FB":
            water_forcefield = "amber14/tip4pfb.xml"
        else:
            return None
    
    elif forcefield_selection == 'charmm36.xml':
        if water == 'CHARMM default':
            water_forcefield = "charmm36/water.xml"
        elif water == "TIP3P-PME-B":
            water_forcefield = "charmm36/tip3p-pme-b.xml"
        elif water == "TIP3P-PME-F":
            water_forcefield = "charmm36/tip3p-pme-f.xml"
        elif water == "SPC/E":
            water_forcefield = "charmm36/spce.xml"
        elif water == "TIP4P-Ew":
            water_forcefield = "charmm36/tip4pew.xml"
        elif water == "TIP4P-2005":
            water_forcefield = "charmm36/tip4p2005.xml"
        elif water == "TIP5P":
            water_forcefield = "charmm36/tip5p.xml"
        elif water == "TIP5P-Ew":
            water_forcefield = "charmm36/tip5pew.xml"
        else:
            return None
    else:
        return None
    
    return water_forcefield

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

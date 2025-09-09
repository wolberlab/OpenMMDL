import simtk.openmm.app as app
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator, SMIRNOFFTemplateGenerator
from openmm import app


def ff_selection(ff):
    """
    Selects the required XML forcefield file.

    Args:
        ff (str): Input forcefield.

    Returns:
        str: Selected XML forcefield file.
    """
    forcefield_dict = {
        "AMBER19": "amber19-all.xml",
        "AMBER14": "amber14-all.xml",
        "AMBER99SB": "amber99sb.xml",
        "AMBER99SB-ILDN": "amber99sbildn.xml",
        "AMBER03": "amber03.xml",
        "AMBER10": "amber10.xml",
        "CHARMM36": "charmm36.xml",
        "CHARMM2024" "charmm36_2024.xml",
    }

    return forcefield_dict.get(ff, None)


def water_forcefield_selection(water, forcefield_selection):
    """
    Selects the XML filename for water force field parameters based on the chosen force field and water model.

    Args:
        water (str): The chosen water model.
        forcefield_selection (str): The selected force field.

    Returns:
        str: The XML filename of the water forcefield.
    """
    old_amber = {"amber99sb.xml", "amber99sbildn.xml", "amber03.xml", "amber10.xml"}

    # Define a dictionary to map water models
    water_model_mapping = {
        "TIP3P": "tip3p.xml",
        "TIP3P-FB": "tip3pfb.xml",
        "SPC/E": "spce.xml",
        "TIP4P-Ew": "tip4pew.xml",
        "TIP4P-FB": "tip4pfb.xml",
        "TIP5P": "tip5p.xml",
        "OPC": "opc.xml",
        "OPC3": "opc3.xml",
    }

    if forcefield_selection in old_amber:
        water_model = water_model_mapping.get(water, None)
    else:
        water_forcefields = {
            "amber14-all.xml": {
                "TIP3P": "amber14/tip3p.xml",
                "TIP3P-FB": "amber14/tip3pfb.xml",
                "SPC/E": "amber14/spce.xml",
                "TIP4P-Ew": "amber14/tip4pew.xml",
                "TIP4P-FB": "amber14/tip4pfb.xml",
                "OPC": "amber14/opc.xml",
                "OPC3": "amber14/opc3.xml",
            },
            "amber19-all.xml": {
                "TIP3P": "amber19/tip3p.xml",
                "TIP3P-FB": "amber19/tip3pfb.xml",
                "SPC/E": "amber19/spce.xml",
                "TIP4P-Ew": "amber19/tip4pew.xml",
                "TIP4P-FB": "amber19/tip4pfb.xml",
                "OPC": "amber19/opc.xml",
                "OPC3": "amber19/opc3.xml",
            },
            "charmm36.xml": {
                "CHARMM default": "charmm36/water.xml",
                "TIP3P-PME-B": "charmm36/tip3p-pme-b.xml",
                "TIP3P-PME-F": "charmm36/tip3p-pme-f.xml",
                "SPC/E": "charmm36/spce.xml",
                "TIP4P-Ew": "charmm36/tip4pew.xml",
                "TIP4P-2005": "charmm36/tip4p2005.xml",
                "TIP5P": "charmm36/tip5p.xml",
                "TIP5P-Ew": "charmm36/tip5pew.xml",
            },
            "charmm36_2024.xml": {
                "CHARMM default": "charmm36_2024/water.xml",
                "TIP3P-PME-B": "charmm36_2024/tip3p-pme-b.xml",
                "TIP3P-PME-F": "charmm36_2024/tip3p-pme-f.xml",
                "SPC/E": "charmm36_2024/spce.xml",
                "TIP4P-Ew": "charmm36_2024/tip4pew.xml",
                "TIP4P-2005": "charmm36_2024/tip4p2005.xml",
                "TIP5P": "charmm36_2024/tip5p.xml",
                "TIP5P-Ew": "charmm36_2024/tip5pew.xml",
            },
        }
        water_model = water_forcefields.get(forcefield_selection, {}).get(water, None)

    return water_model


def water_model_selection(water, forcefield_selection):
    """
    Selects the required water model forcefield XML file according to water selection and previous force field selection.

    Args:
        water (str): Water model input.
        forcefield_selection (str): Input of selected forcefield XML file.

    Returns:
        str: Water model forcefield XML file.
    """
    old_amber = {"amber99sb.xml", "amber99sbildn.xml", "amber03.xml", "amber10.xml"}

    water_model_mapping = {
        "TIP3P": "tip3p",
        "TIP3P-FB": "tip3pfb",
        "SPC/E": "spce",
        "TIP4P-Ew": "tip4pew",
        "TIP4P-FB": "tip4pfb",
        "OPC": "opc",
        "OPC3": "opc3",
    }

    if forcefield_selection in old_amber:
        water_model = water_model_mapping.get(water)
    elif forcefield_selection == "amber14-all.xml":
        if water == "TIP5P":
            return None  # 'TIP5P' is not available in 'amber14-all.xml'
        water_model = water_model_mapping.get(water)
    elif forcefield_selection == "amber19-all.xml":
        if water == "TIP5P":
            return None  # 'TIP5P' is not available in 'amber19-all.xml'
        water_model = water_model_mapping.get(water)
    elif forcefield_selection == "charmm36.xml" or forcefield_selection == "charmm36_2024.xml":
        charmm_water_mapping = {
            "CHARMM default": "charmm",
            "TIP3P-PME-B": "charmm",
            "TIP3P-PME-F": "charmm",
            "SPC/E": "charmm",
            "TIP4P-Ew": "tip4pew",
            "TIP4P-2005": "tip4pew",
            "TIP5P": "tip5p",
            "TIP5P-Ew": "tip5p",
        }
        water_model = charmm_water_mapping.get(water)
    else:
        return None

    return water_model


def generate_forcefield(protein_ff, solvent_ff, add_membrane, smallMoleculeForceField=None, rdkit_mol=None):
    """
    Generate an OpenMM Forcefield object and register a small molecule.

    Args:
        protein_ff (str): Input of selected forcefield XML File.
        solvent_ff (str): Input of selected water model forcefield XML File.
        add_membrane (bool): Selection if the system should be built with a membrane.
        smallMoleculeForceField (str): Name of the force field applied to the small molecule.
        rdkit_mol (rdkit.Chem.rdchem.Mol): Small molecule to register in the force field.

    Returns:
        simtk.openmm.app.Forcefield: Forcefield with a registered small molecule.
    """
    old_amber = {"amber99sb.xml", "amber99sbildn.xml", "amber03.xml", "amber10.xml"}

    # For older amber forcefields, the additional lipid17.xml is required for templates
    if add_membrane == True:
        if protein_ff in old_amber:
            forcefield = app.ForceField(protein_ff, solvent_ff, "amber14/lipid17.xml")
        else:
            forcefield = app.ForceField(protein_ff, solvent_ff)
    else:
        forcefield = app.ForceField(protein_ff, solvent_ff)
    # If a ligand is present, a Forcefield with GAFF or SMIRNOFF will be created for the ligand
    if rdkit_mol is not None:
        if smallMoleculeForceField == "gaff":
            gaff = GAFFTemplateGenerator(
                molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True),
                forcefield="gaff-2.11",
            )
            forcefield.registerTemplateGenerator(gaff.generator)
        elif smallMoleculeForceField == "smirnoff":

            smirnoff = SMIRNOFFTemplateGenerator(
                molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True),
                forcefield="openff-2.2.1",
            )
            forcefield.registerTemplateGenerator(smirnoff.generator)
            
    return forcefield


def generate_transitional_forcefield(
    protein_ff, solvent_ff, add_membrane, smallMoleculeForceField=None, rdkit_mol=None
):
    """
    Generate an OpenMM transitional forcefield object with TIP3P water model for membrane building and register a small molecule.

    Args:
        protein_ff (str): Name of the force field in XML format.
        solvent_ff (str): Name of the water model force field in XML format.
        add_membrane (bool): Selection if the system should be built with a membrane.
        smallMoleculeForceField (str): Name of the force field applied to the small molecule.
        rdkit_mol (rdkit.Chem.rdchem.Mol): Small molecule to register in the force field.

    Returns:
        simtk.openmm.app.Forcefield: A transitional forcefield with TIP3P water and a registered small molecule.
    """
    old_amber = {"amber99sb.xml", "amber99sbildn.xml", "amber03.xml", "amber10.xml"}

    # For older amber forcefields, the additional lipid17.xml is required for templates
    if add_membrane == True:
        if protein_ff in old_amber:
            transitional_forcefield = app.ForceField(
                protein_ff, "tip3p.xml", "amber14/lipid17.xml"
            )
        else:
            transitional_forcefield = app.ForceField(protein_ff, "amber14/tip3p.xml")
    else:
        transitional_forcefield = app.ForceField(protein_ff, solvent_ff)

    # If a ligand is present, a Forcefield with GAFF or SMIRNOFF will be created for the ligand
    if rdkit_mol is not None:
        if smallMoleculeForceField == "gaff":
            gaff = GAFFTemplateGenerator(
                molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True),
                forcefield="gaff-2.11",
            )
            transitional_forcefield.registerTemplateGenerator(gaff.generator)
        elif smallMoleculeForceField == "smirnoff":
            smirnoff = SMIRNOFFTemplateGenerator(
                molecules=Molecule.from_rdkit(rdkit_mol, allow_undefined_stereo=True),
                forcefield="openff-2.2.1",
            )
            transitional_forcefield.registerTemplateGenerator(smirnoff.generator)

    return transitional_forcefield

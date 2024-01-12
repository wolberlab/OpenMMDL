from rdkit import Chem
from openff.toolkit.topology import Molecule, Topology
import simtk.openmm.app as app
from simtk.openmm.app import PDBFile, Modeller
from simtk.openmm import unit
from simtk.openmm import Vec3
import mdtraj as md
import numpy as np

    
def prepare_ligand(ligand_file, minimize_molecule=True):
    """
    Reads an SDF File into RDKit, adds hydrogens to the structure, minimizes it if selected and creates an openforcefield Molecule object.

    Parameters
    ----------
    ligand_file: str
        User input of SDF or MOL File
    minimize_molecule: bool
        Minimization of ligand.
        
    Returns
    -------
    rdkitmolh: rdkit.Chem.rdchem.Mol
    	The prepared and converted ligand.
    """
    # Reading of SDF File, converting to rdkit.
    file_name = ligand_file.lower()
    if file_name.endswith(".sdf"):
        rdkit_mol = Chem.SDMolSupplier(ligand_file, sanitize = False)
        for mol in rdkit_mol:
            rdkit_mol = mol
    elif file_name.endswith(".mol") and not file_name.endswith(".mol2"):
        print(ligand_file)
        rdkit_mol = Chem.rdmolfiles.MolFromMolFile(ligand_file, sanitize = False)
    elif file_name.endswith(".mol2"):
        rdkit_mol = Chem.rdmolfiles.MolFromMol2File(ligand_file, sanitize = False)
    # Adding of hydrogens and assigning chrial tags from the structure.
    print('Adding hydrogens')
    rdkitmolh = Chem.AddHs(rdkit_mol, addCoords=True)
    Chem.AssignAtomChiralTagsFromStructure(rdkitmolh)
    
    # Minimizes the molecule with the MMFF94s Forcefield if selected. 
    if minimize_molecule:
        Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(mol=rdkitmolh, mmffVariant='MMFF94s',maxIters=2000)
        
    # Converting of the ligand from rdkit to an opeenforcefield Molecule object.
    Molecule(rdkitmolh)

    return rdkitmolh
    
def rdkit_to_openmm(rdkit_mol, name):
    """
    Convert an RDKit molecule to an OpenMM molecule.
    Inspired by @hannahbrucemcdonald and @glass-w.

    Parameters
    ----------
    rdkit_mol: rdkit.Chem.rdchem.Mol
        RDKit molecule to convert.
    name: str
        Molecule name.
        
    Returns
    -------
    omm_molecule: simtk.openmm.app.Modeller
        OpenMM modeller object holding the molecule of interest.
    """
    # convert RDKit to OpenFF
    off_mol = Molecule.from_rdkit(rdkit_mol)

    # add name for molecule
    off_mol.name = name

    # add names for atoms
    element_counter_dict = {}
    for off_atom, rdkit_atom in zip(off_mol.atoms, rdkit_mol.GetAtoms()):
        element = rdkit_atom.GetSymbol()
        if element in element_counter_dict.keys():
            element_counter_dict[element] += 1
        else:
            element_counter_dict[element] = 1
        off_atom.name = element + str(element_counter_dict[element])
        
    # convert from OpenFF to OpenMM
    off_mol_topology = off_mol.to_topology()
    mol_topology = off_mol_topology.to_openmm()
    mol_positions = off_mol.conformers[0];
    new_mol_positions = []

    # convert units from Ångström to Nanometers
    for mol_position in off_mol.conformers[0]:
        new_mol_positions.append(mol_position.magnitude/10.0);

    # combine topology and positions in modeller object
    omm_mol = app.Modeller(mol_topology, new_mol_positions * unit.nanometers)

    return omm_mol
    
def merge_protein_and_ligand(protein, ligand):
    """
    Merge two OpenMM objects.

    Parameters
    ----------
    protein: pdbfixer.pdbfixer.PDBFixer
        Protein to merge.
    ligand: simtk.openmm.app.Modeller
        Ligand to merge.

    Returns
    -------
    complex_topology: simtk.openmm.app.topology.Topology
        The merged topology.
    complex_positions: simtk.unit.quantity.Quantity
        The merged positions.
    """
    # combine topologies
    md_protein_topology = md.Topology.from_openmm(protein.topology)  # using mdtraj for protein top
    md_ligand_topology = md.Topology.from_openmm(ligand.topology)  # using mdtraj for ligand top
    md_complex_topology = md_protein_topology.join(md_ligand_topology)  # add them together
    
    complex_topology = md_complex_topology.to_openmm()

    # combine positions
    total_atoms = len(protein.positions) + len(ligand.positions)

    # create an array for storing all atom positions as tupels containing a value and a unit
    # called OpenMM Quantities
    complex_positions = unit.Quantity(np.zeros([total_atoms, 3]), unit=unit.nanometers)
    complex_positions[: len(protein.positions)] = protein.positions  # add protein positions
    complex_positions[len(protein.positions) :] = ligand.positions  # add ligand positions

    return complex_topology, complex_positions

def water_padding_solvent_builder(model_water, forcefield, water_padding_distance, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein_name):
    """
    Build a Solvent Box with padding values.

    Parameters
    ----------
    model_water: str
        Selected water model.
    forcefield: openmm.app.forcefield.ForceField
        Selected Forcefield.
    water_padding_distance: float
        Solvent padding distance.
    protein_pdb: pdbfixer.pdbfixer.PDBFixer
        Protein as a pdbfixer object.
    modeller: openmm.app.modeller.Modeller
        Complex as a modeller object.
    water_positive_ion: str
        Positive ion.
    water_negative_ion: str
        Negative ion.
    water_ionicstrength: float
        Ionic strength.
    protein_name: str
        Protein file name.
    Returns
    -------
    modeller: openmm.app.modeller.Modeller
        The complex with solvent.
    """
    # Writing out the protein without solvent 
    with open(f'prepared_no_solvent_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(protein_pdb.topology, protein_pdb.positions, outfile)
    
    # Adds solvent to the selected protein
    if model_water == 'charmm' or model_water == 'tip3pfb' or model_water == 'tip3':
        modeller.addSolvent(forcefield, padding=water_padding_distance * unit.nanometers, positiveIon=water_positive_ion, negativeIon=water_negative_ion, ionicStrength=water_ionicstrength * unit.molar)
    elif model_water == 'charmm_tip4pew':
        protein_pdb.addSolvent(padding=water_padding_distance * unit.nanometers, positiveIon=water_positive_ion, negativeIon=water_negative_ion, ionicStrength=water_ionicstrength * unit.molar)
    else:
        if model_water == 'tip4pfb':
            model_water = 'tip4pew'
        modeller.addSolvent(forcefield,model=model_water, padding=water_padding_distance * unit.nanometers, positiveIon=water_positive_ion, negativeIon=water_negative_ion, ionicStrength=water_ionicstrength * unit.molar)

    # Writing out the protein with padding solvent
    with open(f'solvent_padding_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)     
    print(f"Protein with buffer solvent prepared")
    
    return modeller

def water_absolute_solvent_builder(model_water, forcefield, water_box_x, water_box_y, water_box_z, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein_name):
    """
    Build a Solvent Box with absolute values.

    Parameters
    ----------
    model_water: str
        Selected water model.
    forcefield: openmm.app.forcefield.ForceField
        Selected Forcefield.
    water_box_x: float
        Vector x of the solvent box.
    water_box_y: float
        Vector y of the solvent box.
    water_box_z: float
        Vector z of the solvent box.
    protein_pdb: pdbfixer.pdbfixer.PDBFixer
        protein as a pdbfixer object.
    modeller: openmm.app.modeller.Modeller
        complex as a modeller object.
    water_positive_ion: str
        Positive ion.
    water_negative_ion: str
        Negative ion.
    water_ionicstrength: float
        Ionic strength.
    protein_name: str
        Protein file name.

    Returns
    -------
    modeller: openmm.app.modeller.Modeller
        The complex with solvent.
    """
    # Writing out the protein without solvent 
    with open(f'prepared_no_solvent_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(protein_pdb.topology, protein_pdb.positions, outfile)
    
    # Adds solvent to the selected protein
    if model_water == 'charmm' or model_water == 'tip3pfb' or model_water == 'tip3':
        modeller.addSolvent(forcefield, boxSize=Vec3(water_box_x, water_box_y, water_box_z) * unit.nanometers, positiveIon=water_positive_ion, negativeIon=water_negative_ion, ionicStrength=water_ionicstrength * unit.molar)
    elif model_water == 'charmm_tip4pew':
        protein_pdb.addSolvent(boxSize=Vec3(water_box_x, water_box_y, water_box_z) * unit.nanometers, positiveIon=water_positive_ion, negativeIon=water_negative_ion, ionicStrength=water_ionicstrength * unit.molar)
    else:
        if model_water == 'tip4pfb':
            model_water = 'tip4pew'
        modeller.addSolvent(forcefield, model=model_water, boxSize=Vec3(water_box_x, water_box_y, water_box_z) * unit.nanometers, positiveIon=water_positive_ion, negativeIon=water_negative_ion, ionicStrength=water_ionicstrength * unit.molar)
        
    # Writing out the protein with absolute solvent    
    with open(f'solvent_absolute_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)     
    print(f"Protein with absolute solvent prepared")
    
    return modeller

def membrane_builder(ff, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, membrane_lipid_type, membrane_padding, membrane_positive_ion, membrane_negative_ion, membrane_ionicstrength, protein_name):
    """
    Build a membrane with minimum padding.

    Parameters
    ----------
    ff: str
        Selected Forcefield as a string.    
    model_water: str
        Selected water model.
    forcefield: openmm.app.forcefield.ForceField
        Selected Forcefield.
    transitional_forcefield: openmm.app.forcefield.ForceField
        Transitional forcefield for specific water models.    
    protein_pdb: pdbfixer.pdbfixer.PDBFixer
        protein as a pdbfixer object.    
    membrane_lipid_type: str
        Lipid type.
   membrane_padding: float
        Minimum membrane padding distance.
    membrane_positive_ion: str
        Positive ion.
    membrane_negative_ion: str
        Negative ion.
    membrane_ionicstrength: float
        Ionic strength.
    protein_name: str
        Protein file name.

    Returns
    -------
    modeller: openmm.app.modeller.Modeller
        The complex with solvent.
    """ 
    # Writing out the protein without solvent 
    with open(f'prepared_no_solvent_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(protein_pdb.topology, protein_pdb.positions, outfile)
    
    # Adds a membrane to the selected protein
    # The Water Models TIP4P and TIP5P require an transitional forcefield
    if ff == 'CHARMM36':
        protein_pdb.addMembrane(lipidType= membrane_lipid_type, minimumPadding= membrane_padding * unit.nanometer, positiveIon=membrane_positive_ion, negativeIon=membrane_negative_ion, ionicStrength=membrane_ionicstrength * unit.molar)
        modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
    else:
        if model_water == 'charmm':
            modeller.addMembrane(forcefield, lipidType= membrane_lipid_type, minimumPadding= membrane_padding * unit.nanometer, positiveIon=membrane_positive_ion, negativeIon=membrane_negative_ion, ionicStrength=membrane_ionicstrength * unit.molar)
        else:
            if model_water == 'tip4pew' or model_water == 'tip5p':
                modeller.addMembrane(transitional_forcefield, lipidType= membrane_lipid_type, minimumPadding= membrane_padding * unit.nanometer, positiveIon=membrane_positive_ion, negativeIon=membrane_negative_ion, ionicStrength=membrane_ionicstrength * unit.molar)
            else:
                modeller.addMembrane(forcefield, lipidType= membrane_lipid_type, minimumPadding= membrane_padding * unit.nanometer, positiveIon=membrane_positive_ion, negativeIon=membrane_negative_ion, ionicStrength=membrane_ionicstrength * unit.molar)

    with open(f'membrane_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    
    print(f"Protein with Membrane {membrane_lipid_type} prepared")
                    
    return modeller

def water_conversion(model_water, modeller_pre_conversion, protein_name):
    """
    Merge two OpenMM objects.

    Parameters
    ----------
    model_water: str
        The name of the prefered converted Water model.
    modeller_pre_conversion: pdbfixer.pdbfixer.PDBFixer
        The object that will be converted.
    protein_name: str
        Name of the Protein pdb file.

    Returns
    -------
    modeller: openmm.app.modeller.Modeller
        The converted object.
    """
    # Writes out the PDB of the preconverted pdb
    with open(f'pre_converted_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(modeller_pre_conversion.topology, modeller_pre_conversion.positions, outfile)
    
    # Converts the water from TIP3P to the required Water Model 
    modeller_pre_conversion.convertWater(model_water)
    modeller = modeller_pre_conversion
    
    # Writes out the pdb of the postconverted pdb
    with open(f'converted_{protein_name}', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    
    return modeller

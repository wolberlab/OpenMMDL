import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

def increase_ring_indices(ring, lig_index):
    """
    Increases the atom indices in a ring of the ligand obtained from the ligand to fit the atom indices present in the protein-ligand complex.

    Parameters
    ----------
    ring : str
        A list of atom indices belonging to a ring that need to be modified.
    lig_index : str
        An integer that is the first number of the ligand atom indices obtained from the protein-ligand, which is used to modify the ring indices

    Returns
    -------
    list of int :
        A new list with modified atom indicies.
    """
    return [atom_idx + lig_index for atom_idx in ring]

def convert_ligand_to_smiles(input_sdf, output_smi):
    """
    Converts ligand structures from an SDF file to SMILES :) format

    Parameters
    ----------
    input_sdf : str
        Path to the SDF file with the ligand that wll be converted.
    output_smi : str
        Path to the output SMILES files.

    Returns
    -------
    None
    """
    # Create a molecule supplier from an SDF file
    mol_supplier = Chem.SDMolSupplier(input_sdf)

    # Open the output SMILES file for writing
    with open(output_smi, "w") as output_file:
        # Iterate through molecules and convert each to SMILES
        for mol in mol_supplier:
            if mol is not None:  # Make sure the molecule was successfully read
                smiles = Chem.MolToSmiles(mol)
                print(smiles)
                output_file.write(smiles + "\n")
            else:
            	print("nono")


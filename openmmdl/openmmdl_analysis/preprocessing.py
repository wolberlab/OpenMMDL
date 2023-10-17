import MDAnalysis as mda
import subprocess
import os

def process_pdb_file(input_pdb_filename):
    """Process a PDB file to make it compatible with the openmmdl_analysis package.

    Args:
        input_pdb_filename (str): path to the input PDB file
    """
    # Load the PDB file
    u = mda.Universe(input_pdb_filename)

    # Iterate through the topology to modify residue names
    for atom in u.atoms:
        resname = atom.resname

        # Check if the residue name is one of the specified water names
        if resname in ["SPC", "TIP3", "TIP4", "WAT", "T3P", "T4P", "T5P"]:
            atom.residue.resname = "HOH"
        elif resname == "*":
            atom.residue.resname = "UNK"

    # Save the modified topology to a new PDB file
    u.atoms.write(input_pdb_filename)

# def extract_and_save_ligand_as_pdb(input_pdb_filename, output_pdb_filename, target_resname):
#     """Extract and save the ligand from the receptor ligand complex PDB file into a new PDB file by itself .

#     Args:
#         input_pdb_filename (str): name of the input PDB file
#         output_pdb_filename (str): name of the output PDB file
#         target_resname (str): resname of the ligand in the original PDB file
#     """
#     # Load the PDB file using MDAnalysis
#     u = mda.Universe(input_pdb_filename)

#     # Find the ligand by its residue name
#     ligand_atoms = u.select_atoms(f"resname {target_resname}")

#     if len(ligand_atoms) == 0:
#         print(f"No ligand with residue name '{target_resname}' found in the PDB file.")
#         return

#     # Create a new Universe with only the ligand
#     ligand_universe = mda.Merge(ligand_atoms)

#     # Save the ligand Universe as a PDB file
#     ligand_universe.atoms.write(output_pdb_filename)

def convert_pdb_to_sdf(input_pdb_filename, output_sdf_filename):
    """Convert ligand PDB file to SDF file for analysis using Open Babel.

    Args:
        input_pdb_filename (str): name of the input PDB file
        output_sdf_filename (str): name of the output SDF file
    """
    # Use subprocess to call Open Babel for the file format conversion
    try:
        subprocess.run(["obabel", input_pdb_filename, "-O", output_sdf_filename], check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error converting PDB to SDF: {e}")
        return

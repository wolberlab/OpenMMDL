import MDAnalysis as mda
import subprocess
import os
import re

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

def renumber_atoms_in_residues(input_pdb_file, output_pdb_file, lig_name):
    # Read the input PDB file
    with open(input_pdb_file, 'r') as f:
        pdb_lines = f.readlines()

    new_pdb_lines = []
    lig_residue_elements = {}

    for line in pdb_lines:
        if line.startswith("ATOM"):
            # Extract information from the ATOM line
            atom_serial = int(line[6:11])
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()

            # Check if the residue is residue_name
            if residue_name == lig_name:
                # Extract the element from the atom name (assuming it starts with a capital letter)
                element_match = re.match(r"([A-Z]+)", atom_name)
                if element_match:
                    element = element_match.group(1)
                else:
                    element = atom_name

                # Increment the count for the current element in the current residue_name residue
                lig_residue_elements[element] = lig_residue_elements.get(element, 0) + 1

                # Update the atom name based on the element count
                new_atom_name = f"{element}{lig_residue_elements[element]}"

                # Replace the old atom name with the new one
                line = line[:12] + f"{new_atom_name:4}" + line[16:]

        # Append the line to the new PDB lines
        new_pdb_lines.append(line)

    # Write the modified PDB lines to the output file
    with open(output_pdb_file, 'w') as f:
        f.writelines(new_pdb_lines)


def replace_atom_type(data):
    lines = data.split('\n')
    for i, line in enumerate(lines):
        if ' LIG  X' in line:
            # Extract the last column which contains the atom type (O/N/H)
            atom_type = line[12:13].strip()
            # Replace 'X' with the correct atom type
            lines[i] = line.replace(' LIG  X', f' LIG  {atom_type}')
    return '\n'.join(lines)

def process_pdb(input_file, output_file):
    with open(input_file, 'r') as f:
        pdb_data = f.read()

    modified_data = replace_atom_type(pdb_data)

    with open(output_file, 'w') as f:
        f.write(modified_data)

def move_hydrogens_to_end(structure, target_residue_name):
    # Counter for atom numbering within each residue
    atom_counter = 1

    # Iterate over all models in the structure
    for model in structure:
        # Iterate over all chains in the model
        for chain in model:
            # Iterate over all residues in the chain
            for residue in chain:
                # Check if the residue name matches the target residue name
                if residue.resname == target_residue_name:
                    # Collect hydrogen atoms in the residue
                    hydrogen_atoms = [atom for atom in residue if atom.element == 'H']

                    # Remove hydrogen atoms from the residue
                    for hydrogen_atom in hydrogen_atoms:
                        residue.detach_child(hydrogen_atom.id)

                    # Add hydrogen atoms to the end of the residue with renumbering
                    for hydrogen_atom in hydrogen_atoms:
                        hydrogen_atom.serial_number = atom_counter
                        atom_counter += 1
                        # Change the residue name to avoid conflicts
                        hydrogen_atom.name = (f'H{atom_counter}')
                        residue.add(hydrogen_atom)

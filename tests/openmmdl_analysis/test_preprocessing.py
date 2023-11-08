import os
import pytest
import numpy as np
import MDAnalysis as mda
from openmmdl.openmmdl_analysis.preprocessing import process_pdb_file

pdb_file_path = 'openmmdl/tests/data/in/0_unk_hoh.pdb' 

@pytest.fixture
def temp_pdb_file(tmp_path):
    input_pdb_filename = tmp_path / "test_input.pdb"
    # Copy the content of the provided PDB file to the temporary test file
    with open(pdb_file_path, "r") as src_pdb, open(input_pdb_filename, "w") as dest_pdb:
        dest_pdb.write(src_pdb.read())
    return input_pdb_filename

def test_process_pdb_file(temp_pdb_file):
    input_pdb_filename = temp_pdb_file
    # Check if the function modifies the PDB file as expected
    process_pdb_file(input_pdb_filename)

    # Load the modified PDB file
    u = mda.Universe(input_pdb_filename)

    # Check that the residue name "*" is changed to "UNK"
    for atom in u.atoms:
        if atom.residue.resname == "UNK":
            assert atom.residue.resname == "UNK"

    # You can add additional assertions as needed.

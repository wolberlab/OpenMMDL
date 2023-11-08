import pytest
import os
from pathlib import Path
import mdtraj as md

from openmmdl.openmmdl_simulation.scripts.post_md_conversions import mdtraj_conversion, MDanalysis_conversion



test_data_directory = Path("openmmdl/tests/data/in")
pdb_file = "0_unk_hoh.pdb"
dcd_file = "trajectory.dcd"
ligand_name = 'UNK'

def test_MDanalysis_conversion():
    original_cwd = os.getcwd()
    os.chdir(test_data_directory)
    # Create temporary directories to save the output files
    output_file_dcd = "centered_old_coordinates.dcd"
    output_file_xtc = 'centered_old_coordinates.xtc'
    output_file_pdb = 'centered_old_coordinates_top.pdb'
    output_file_gro = 'centered_old_coordinates_top.gro'
  
    mdtraj_conversion(pdb_file, "gro_xtc")
    mdtraj_conversion(pdb_file, "pdb_dcd")
    
    assert output_file_dcd is not None
    assert output_file_xtc is not None
    assert output_file_pdb is not None
    assert output_file_gro is not None
    os.chdir(original_cwd)

def test_mdtraj_conversion():
    original_cwd = os.getcwd()
    os.chdir(test_data_directory)
    
    # Create temporary directories to save the output files
    all_file_dcd = "centered_traj.dcd"
    all_file_pdb = 'centered_top.pdb'
    prot_lig_file_dcd = 'prot_lig_traj.dcd'
    prot_lig_file_pdb = 'prot_lig_top.pdb'
    all_file_xtc = "centered_traj.xtc"
    all_file_gro = "centered_top.gro"
    prot_lig_file_xtc = 'prot_lig_traj.xtc'
    prot_lig_file_gro = 'prot_lig_top.gro'

    MDanalysis_conversion(pdb_file, dcd_file, ligand_name, "pdb_dcd_gro_xtc", "mda_prot_lig_all")

    assert all_file_dcd is not None
    assert all_file_pdb is not None
    assert prot_lig_file_dcd is not None
    assert prot_lig_file_pdb is not None
    assert all_file_xtc is not None
    assert all_file_gro is not None
    assert prot_lig_file_xtc is not None
    assert prot_lig_file_gro is not None

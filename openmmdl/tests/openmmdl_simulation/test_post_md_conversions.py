import pytest
import os
import shutil
from pathlib import Path
import mdtraj as md

from openmmdl.openmmdl_simulation.scripts.post_md_conversions import mdtraj_conversion, MDanalysis_conversion

test_data_directory = Path("openmmdl/tests/data/in")
pdb_file = "0_unk_hoh.pdb"
dcd_file = "trajectory.dcd"
ligand_name = 'UNK'

def test_mdtraj_conversion():
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

def test_mdanalysis_conversion():
    original_cwd = Path(os.getcwd())
    test_data_directory = Path("openmmdl/tests/data/in")
    post_mdtraj_pdb_file = test_data_directory / "centered_old_coordinates_top.pdb"
    post_mdtraj_dcd_file = test_data_directory / "centered_old_coordinates.dcd"
    
    # Create temporary directories to save the output files
    all_file_dcd = "centered_traj.dcd"
    all_file_dcd_unaligned = "centered_traj_unaligned.dcd"
    all_file_pdb = 'centered_top.pdb'
    prot_lig_file_dcd = 'prot_lig_traj.dcd'
    prot_lig_file_dcd_unaligned = 'prot_lig_traj_unaligned.dcd'
    prot_lig_file_pdb = 'prot_lig_top.pdb'
    all_file_xtc = "centered_traj.xtc"
    all_file_xtc_unaligned = "centered_traj_unaligned.xtc"
    all_file_gro = "centered_top.gro"
    prot_lig_file_xtc = 'prot_lig_traj.xtc'
    prot_lig_file_xtc_unaligned = 'prot_lig_traj_unaligned.xtc'
    prot_lig_file_gro = 'prot_lig_top.gro'

    shutil.copy(str(post_mdtraj_pdb_file), '.')
    shutil.copy(str(post_mdtraj_dcd_file), '.')

    post_mdtraj_pdb_file = "centered_old_coordinates_top.pdb"
    post_mdtraj_dcd_file = "centered_old_coordinates.dcd"
    ligand_name = "UNK"
    mda_output = "pdb_dcd_gro_xtc"
    output_selection = "mda_prot_lig_all"
    

    #MDanalysis_conversion(pdb_file, dcd_file, ligand_name, "pdb_dcd_gro_xtc", "mda_prot_lig_all")
    MDanalysis_conversion(post_mdtraj_pdb_file, post_mdtraj_dcd_file, mda_output, output_selection, ligand_name)

    assert all_file_dcd is not None
    assert all_file_dcd_unaligned is not None
    assert all_file_pdb is not None
    assert prot_lig_file_dcd is not None
    assert prot_lig_file_dcd_unaligned is not None
    assert prot_lig_file_pdb is not None
    assert all_file_xtc is not None
    assert all_file_xtc_unaligned is not None
    assert all_file_gro is not None
    assert prot_lig_file_xtc is not None
    assert prot_lig_file_xtc_unaligned is not None
    assert prot_lig_file_gro is not None

    # Assertions or checks to verify the correctness of the results
    if "pdb" in mda_output:
        if output_selection != "mda_all":
            # Check if the expected PDB file exists
            pdb_file_path = original_cwd / "prot_lig_top.pdb"
            assert pdb_file_path.is_file()

            # Check if the expected DCD file exists
            dcd_file_path = original_cwd / "prot_lig_traj.dcd"
            assert dcd_file_path.is_file()

            # Check if the DCD file is not empty
            traj = md.load(dcd_file_path, top=pdb_file_path)
            assert traj.n_frames > 0

    if "gro" in mda_output:
        if output_selection != "mda_all":
            # Check if the expected GRO file exists
            gro_file_path = original_cwd / "prot_lig_top.gro"
            assert gro_file_path.is_file()

            # Check if the expected XTC file exists
            xtc_file_path = original_cwd / "prot_lig_traj.xtc"
            assert xtc_file_path.is_file()

            # Check if the XTC file is not empty
            traj = md.load(xtc_file_path, top=gro_file_path)
            assert traj.n_frames > 0

import pytest
from unittest.mock import MagicMock, patch
import MDAnalysis as mda
from pathlib import Path
# Import the class to be tested
from openmmdl.openmmdl_analysis.trajectory_saving import TrajectorySaver

test_data_directory = Path("openmmdl/tests/data/in")
pdb_file = test_data_directory / "interacting_waters.pdb"
dcd_file = test_data_directory / "interacting_waters.dcd"

def test_init():
    pdb_md = mda.Universe(pdb_file, dcd_file)
    saver = TrajectorySaver(pdb_md, "LIG", "HEM", nucleic=False)
    assert saver.pdb_md == pdb_md
    assert saver.ligname == "LIG"
    assert saver.special == "HEM"
    assert saver.nucleic is False

def test_save_interacting_waters_trajectory(tmp_path):
    pdb_md = mda.Universe(pdb_file, dcd_file)
    saver = TrajectorySaver(pdb_md, "UNK", None, nucleic=False)
    water_id = [7536]

    out_dir = tmp_path / "output"
    out_dir.mkdir()

    saver.save_interacting_waters_trajectory(interacting_waters=water_id, outputpath=str(out_dir) + "/")

    # Filepaths
    pdb_out = out_dir / "interacting_waters.pdb"
    dcd_out = out_dir / "interacting_waters.dcd"

    # Assertions
    assert pdb_out.exists(), "PDB file was not created."
    assert dcd_out.exists(), "DCD file was not created."

    assert pdb_out.stat().st_size > 0, "PDB file is empty."
    assert dcd_out.stat().st_size > 0, "DCD file is empty."

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
    assert saver.pdb_md == mock_universe
    assert saver.ligname == "LIG"
    assert saver.special == "HEM"
    assert saver.nucleic is False

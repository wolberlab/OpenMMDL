import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pytest
from unittest.mock import patch, MagicMock
import tempfile
import shutil

# Import the class to be tested
from openmmdl.openmmdl_analysis.rmsd_calculation import RMSDAnalyzer

@pytest.fixture
def analyzer():
    """Create and return an RMSDAnalyzer instance with mocked topology and trajectory files."""
    # Use simple file paths instead of MDAnalysis test data
    mock_top = "mock_topology.pdb"
    mock_traj = "mock_trajectory.dcd"
    
    # Patch the mda.Universe to avoid actually loading files
    with patch('MDAnalysis.Universe') as mock_universe:
        # Configure the mock universe
        universe_instance = mock_universe.return_value
        universe_instance.trajectory = MagicMock()
        
        # Create the analyzer with mocked files
        rmsd_analyzer = RMSDAnalyzer(mock_top, mock_traj)
        
        # Now replace the universe that was created during initialization
        # with our controlled mock
        rmsd_analyzer.universe = universe_instance
        
        yield rmsd_analyzer

@pytest.fixture
def test_dir():
    """Create a temporary directory for test outputs."""
    temp_dir = tempfile.mkdtemp()
    original_dir = os.getcwd()
    os.chdir(temp_dir)
    yield temp_dir
    os.chdir(original_dir)
    shutil.rmtree(temp_dir)

def test_initialization():
    """Test that the RMSDAnalyzer initializes correctly."""
    # Use basic mocking to avoid actual file loading
    with patch('MDAnalysis.Universe') as mock_universe:
        mock_top = "mock_topology.pdb"
        mock_traj = "mock_trajectory.dcd"
        analyzer = RMSDAnalyzer(mock_top, mock_traj)
        
        assert analyzer.prot_lig_top_file == mock_top
        assert analyzer.prot_lig_traj_file == mock_traj
        mock_universe.assert_called_once_with(mock_top, mock_traj)

@patch('os.makedirs')
@patch('pandas.DataFrame.to_csv')
@patch('matplotlib.pyplot.savefig')
def test_rmsd_for_atomgroups_protein_only(mock_savefig, mock_to_csv, mock_makedirs, analyzer, test_dir):
    """Test RMSD calculation for protein only."""
    with patch('MDAnalysis.analysis.rms.RMSD') as mock_rmsd:
        # Configure the mock
        mock_rmsd_instance = mock_rmsd.return_value
        mock_rmsd_instance.run.return_value = None
        mock_rmsd_instance.rmsd = np.array([[0, 0, 1.0], [1, 0, 1.5], [2, 0, 2.0]])
        
        # Run the method
        result = analyzer.rmsd_for_atomgroups("png", "protein")
        
        # Check that the method was called correctly
        mock_rmsd.assert_called_once()
        mock_rmsd_instance.run.assert_called_once()
        
        # Check the result
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3
        assert result.columns.tolist() == ["protein"]
        
        # Check that output directory was created
        mock_makedirs.assert_called_with("./RMSD/", exist_ok=True)
        
        # Check that files were saved
        mock_to_csv.assert_called_once()
        mock_savefig.assert_called_once()

@patch('os.makedirs')
@patch('pandas.DataFrame.to_csv')
@patch('matplotlib.pyplot.savefig')
def test_rmsd_for_atomgroups_with_multiple_selections(mock_savefig, mock_to_csv, mock_makedirs, analyzer, test_dir):
    """Test RMSD calculation with protein and ligand selections."""
    with patch('MDAnalysis.analysis.rms.RMSD') as mock_rmsd:
        # Configure the mock
        mock_rmsd_instance = mock_rmsd.return_value
        mock_rmsd_instance.run.return_value = None
        mock_rmsd_instance.rmsd = np.array([[0, 0, 1.0, 2.0], [1, 0, 1.5, 2.5], [2, 0, 2.0, 3.0]])
        
        # Run the method
        result = analyzer.rmsd_for_atomgroups("png", "protein", ["resname LIG"])
        
        # Check the result
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 3
        assert result.columns.tolist() == ["protein", "resname LIG"]


def test_calc_rmsd_2frames(analyzer):
    """Test RMSD calculation between two frames."""
    # Create two test frames
    ref = np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]])
    frame = np.array([[1.1, 0.9, 1.0], [2.1, 1.9, 2.0], [3.1, 2.9, 3.0]])
    
    # Calculate RMSD
    rmsd = analyzer.calc_rmsd_2frames(ref, frame)
    
    # Calculate expected RMSD by hand
    expected = np.sqrt(((0.1**2 + 0.1**2 + 0**2) + (0.1**2 + 0.1**2 + 0**2) + (0.1**2 + 0.1**2 + 0**2)) / 3)
    
    # Check the result
    assert rmsd == pytest.approx(expected, abs=1e-6)

@patch('tqdm.tqdm')
def test_calculate_distance_matrix(mock_tqdm, analyzer):
    """Test distance matrix calculation."""
    # Create mock trajectory with 3 frames
    analyzer.universe.trajectory = MagicMock()
    analyzer.universe.trajectory.__len__.return_value = 3
    
    # Mock the select_atoms method
    mock_atoms = MagicMock()
    analyzer.universe.select_atoms.return_value = mock_atoms
    
    # Set positions for each frame
    positions = [
        np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]),
        np.array([[1.1, 0.9, 1.0], [2.1, 1.9, 2.0]]),
        np.array([[1.5, 0.5, 1.0], [2.5, 1.5, 2.0]])
    ]
    
    def set_positions(frame_idx):
        mock_atoms.positions = positions[frame_idx]
    
    analyzer.universe.trajectory.__getitem__.side_effect = set_positions
    
    # Mock the RMSD calculation
    with patch.object(analyzer, 'calc_rmsd_2frames') as mock_calc_rmsd:
        # Define RMSD values
        rmsd_values = {
            (0, 1): 0.1,
            (0, 2): 0.5,
            (1, 2): 0.4
        }
        
        def mock_rmsd(ref, frame):
            # Find which positions are being compared
            for i, pos_i in enumerate(positions):
                if np.array_equal(ref, pos_i):
                    for j, pos_j in enumerate(positions):
                        if np.array_equal(frame, pos_j):
                            return rmsd_values.get((min(i, j), max(i, j)), 0.0)
            return 0.0
        
        mock_calc_rmsd.side_effect = mock_rmsd
        
        # Calculate the distance matrix
        result = analyzer.calculate_distance_matrix("protein")
        
        # Check the result
        expected = np.array([
            [0.0, 0.1, 0.5],
            [0.1, 0.0, 0.4],
            [0.5, 0.4, 0.0]
        ])
        np.testing.assert_array_almost_equal(result, expected, decimal=6)

def test_calculate_representative_frame(analyzer):
    """Test calculation of representative frame."""
    # Create a distance matrix
    dm = np.array([
        [0.0, 1.0, 2.0, 3.0],
        [1.0, 0.0, 0.5, 1.5],
        [2.0, 0.5, 0.0, 1.0],
        [3.0, 1.5, 1.0, 0.0]
    ])
    
    # Define frames in binding mode
    bmode_frames = [1, 2, 3, 4]  # 1-indexed as in the function
    
    # Calculate representative frame
    rep_frame = analyzer.calculate_representative_frame(bmode_frames, dm)
    
    # Frame 2 should have the lowest average RMSD to other frames
    # (1.0 + 0.5 + 1.5) / 3 = 1.0
    assert rep_frame == 2
    
    # Test with a subset of frames
    bmode_frames = [1, 3, 4]  # 1-indexed
    rep_frame = analyzer.calculate_representative_frame(bmode_frames, dm)
    
    # Frame 3 should have the lowest average RMSD to other frames
    # (2.0 + 1.0) / 2 = 1.5
    assert rep_frame == 3
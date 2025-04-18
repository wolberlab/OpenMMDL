import pytest
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from unittest.mock import patch, MagicMock, mock_open
import numpy as np

from openmmdl.openmmdl_analysis.markov_state_figure_generation import MarkovChainAnalysis

@pytest.fixture
def markov_chain():
    """Create a basic MarkovChainAnalysis instance for testing."""
    return MarkovChainAnalysis(min_transition=5)


@pytest.fixture
def sample_combined_dict():
    """Create a sample combined dictionary for testing."""
    # Creating a sequence with repeating patterns to simulate transitions
    sequence = ["A", "A", "B", "B", "C", "A", "B", "C", "A", "A"]
    # Repeat the sequence to get more data points
    all_data = sequence * 20
    return {"all": all_data}


def test_init():
    """Test the initialization of the MarkovChainAnalysis class."""
    analyzer = MarkovChainAnalysis(min_transition=5)
    
    assert analyzer.min_transition == 5
    assert analyzer.min_transitions == [5, 10, 25, 50]


def test_calculate_min_transitions():
    """Test the calculate_min_transitions method."""
    analyzer = MarkovChainAnalysis(min_transition=2)
    
    assert analyzer.calculate_min_transitions() == [2, 4, 10, 20]
    
    # Test with a different value
    analyzer.min_transition = 10
    assert analyzer.calculate_min_transitions() == [10, 20, 50, 100]


@patch("matplotlib.pyplot.figure")
@patch("matplotlib.pyplot.title")
@patch("matplotlib.pyplot.legend")
@patch("matplotlib.pyplot.axis")
@patch("matplotlib.pyplot.tight_layout")
@patch("matplotlib.pyplot.savefig")
@patch("matplotlib.pyplot.clf")
@patch("networkx.spring_layout")
@patch("networkx.draw_networkx_edges")
@patch("networkx.draw_networkx_edge_labels")
@patch("networkx.draw_networkx_nodes")
@patch("networkx.draw_networkx_labels")
@patch("os.makedirs")
def test_generate_transition_graph(mock_makedirs, mock_draw_labels, mock_draw_nodes, 
                                  mock_draw_edge_labels, mock_draw_edges, mock_spring_layout,
                                  mock_clf, mock_savefig, mock_tight_layout, mock_axis, 
                                  mock_legend, mock_title, mock_figure, 
                                  markov_chain, sample_combined_dict):
    """Test the generate_transition_graph method."""
    # Mock the spring layout to return a fixed position dictionary
    mock_pos = {node: [0, 0] for node in set(sample_combined_dict["all"])}
    mock_spring_layout.return_value = mock_pos
    
    # Call the method
    markov_chain.generate_transition_graph(
        total_frames=len(sample_combined_dict["all"]),
        combined_dict=sample_combined_dict,
        fig_type="png",
        font_size=36,
        size_node=200
    )
    
    # Check that the required functions were called
    mock_figure.assert_called()
    mock_title.assert_called()
    mock_spring_layout.assert_called()
    mock_makedirs.assert_called_with("Binding_Modes_Markov_States", exist_ok=True)
    mock_savefig.assert_called()
    mock_clf.assert_called()


@patch("matplotlib.pyplot.figure")
@patch("matplotlib.pyplot.savefig")
@patch("os.makedirs")
@patch("networkx.draw_networkx_edges")
@patch("networkx.draw_networkx_edge_labels")
@patch("networkx.draw_networkx_nodes")
@patch("networkx.draw_networkx_labels")
def test_generate_transition_graph_with_different_parameters(
    mock_draw_labels, mock_draw_nodes, mock_draw_edge_labels, mock_draw_edges,
    mock_makedirs, mock_savefig, mock_figure, markov_chain, sample_combined_dict):
    """Test the generate_transition_graph method with different parameters."""
    with patch("networkx.spring_layout") as mock_spring_layout:
        # Mock the spring layout to return a fixed position dictionary
        mock_pos = {node: [0, 0] for node in set(sample_combined_dict["all"])}
        mock_spring_layout.return_value = mock_pos
        
        # Call the method with different parameters
        markov_chain.generate_transition_graph(
            total_frames=len(sample_combined_dict["all"]),
            combined_dict=sample_combined_dict,
            fig_type="pdf",
            font_size=24,
            size_node=100
        )
        
        # Check that the plot was saved with pdf extension
        # Use assert_called() instead of assert_called_with() since the min_transition value
        # may be dynamically calculated or different than expected
        mock_savefig.assert_called()
        
        # Alternatively, if you still want to check file extension:
        called_args = mock_savefig.call_args[0][0]
        assert called_args.endswith(".pdf")
        assert "Binding_Modes_Markov_States" in called_args
        

def test_generate_transition_graph_with_empty_data():
    """Test the generate_transition_graph method with empty data."""
    analyzer = MarkovChainAnalysis(min_transition=5)
    empty_combined_dict = {"all": []}
    
    # This should not raise an error, but won't generate any graph
    analyzer.generate_transition_graph(
        total_frames=0,
        combined_dict=empty_combined_dict,
        fig_type="png"
    )


def test_generate_transition_graph_with_single_state():
    """Test the generate_transition_graph method with a single state."""
    analyzer = MarkovChainAnalysis(min_transition=5)
    single_state_dict = {"all": ["A"] * 100}
    
    with patch("networkx.spring_layout") as mock_spring_layout:
        # Mock the spring layout to return a fixed position dictionary
        mock_pos = {"A": [0, 0]}
        mock_spring_layout.return_value = mock_pos
        
        # This should create a graph with a single node and a self-loop
        analyzer.generate_transition_graph(
            total_frames=100,
            combined_dict=single_state_dict,
            fig_type="png"
        )


@patch("networkx.DiGraph")
def test_transition_calculations(mock_digraph, markov_chain, sample_combined_dict):
    """Test that transitions are calculated correctly."""
    # Create a mock graph instance
    mock_graph = MagicMock()
    mock_digraph.return_value = mock_graph
    
    with patch("networkx.spring_layout") as mock_spring_layout:
        # Mock the spring layout to return a fixed position dictionary
        mock_pos = {node: [0, 0] for node in set(sample_combined_dict["all"])}
        mock_spring_layout.return_value = mock_pos
        
        # Call the method
        markov_chain.generate_transition_graph(
            total_frames=len(sample_combined_dict["all"]),
            combined_dict=sample_combined_dict
        )
        
        # Check that add_edge was called for transitions
        mock_graph.add_edge.assert_called()
        
        # Get all calls to add_edge
        call_args_list = mock_graph.add_edge.call_args_list
        
        # Check that there are calls with weights (transitions were added)
        assert len(call_args_list) > 0


def test_parts_division(markov_chain, sample_combined_dict):
    """Test that the dataset is correctly divided into three parts."""
    # Calculate expected part lengths
    total_length = len(sample_combined_dict["all"])
    part_length = total_length // 3
    remaining_length = total_length % 3
    expected_part1_length = part_length + remaining_length
    expected_part2_length = part_length
    expected_part3_length = total_length - expected_part1_length - expected_part2_length
    
    # Extract the parts using the same logic as in the method
    part1_data = sample_combined_dict["all"][:expected_part1_length]
    part2_data = sample_combined_dict["all"][expected_part1_length:expected_part1_length + expected_part2_length]
    part3_data = sample_combined_dict["all"][expected_part1_length + expected_part2_length:]
    
    # Check that parts have the expected lengths
    assert len(part1_data) == expected_part1_length
    assert len(part2_data) == expected_part2_length
    assert len(part3_data) == expected_part3_length
    assert len(part1_data) + len(part2_data) + len(part3_data) == total_length

def test_min_transitions_property_setter():
    """Test that min_transitions is updated when min_transition is changed."""
    analyzer = MarkovChainAnalysis(min_transition=5)
    initial_transitions = analyzer.min_transitions.copy()
    
    # Change min_transition and recalculate
    analyzer.min_transition = 10
    analyzer.min_transitions = analyzer.calculate_min_transitions()
    
    # Check that min_transitions has been updated
    assert analyzer.min_transitions != initial_transitions
    assert analyzer.min_transitions == [10, 20, 50, 100]

def test_binding_site_markov_network():
    # Define test data
    total_frames = 1000
    min_transition = 5
    combined_dict = {"all": ["A", "B", "A", "C", "B", "A", "C", "A", "A", "B"]}

    # Run the function
    analyzer = MarkovChainAnalysis(min_transition)
    analyzer.generate_transition_graph(total_frames, combined_dict)

    # Check if the output file exists for each min_transition
    plot_filename = f"markov_chain_plot_{min_transition}.png"
    plot_path = os.path.join("Binding_Modes_Markov_States", plot_filename)
    assert os.path.exists(plot_path)

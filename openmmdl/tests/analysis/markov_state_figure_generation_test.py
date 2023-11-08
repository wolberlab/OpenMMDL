import networkx as nx
import os
from openmmdl.openmmdl_analysis.markov_state_figure_generation import min_transition_calculation, binding_site_markov_network

# Create a test for min_transition_calculation
def test_min_transition_calculation():
    min_transition = 10
    expected_output = [10, 20, 50, 100]
    result = min_transition_calculation(min_transition)
    assert result == expected_output

# Create a test for binding_site_markov_network
def test_binding_site_markov_network():
    # Define test data
    total_frames = 1000
    min_transitions = [5, 10]
    combined_dict = {
        'all': ['A', 'B', 'A', 'C', 'B', 'A', 'C', 'A', 'A', 'B']
    }

    # Run the function
    binding_site_markov_network(total_frames, min_transitions, combined_dict)

    # Check if the output file exists for each min_transition
    for min_transition_percent in min_transitions:
        plot_filename = f"markov_chain_plot_{min_transition_percent}.png"
        plot_path = os.path.join("Binding_Modes_Markov_States", plot_filename)
        assert os.path.exists(plot_path)

# Optionally, you can include more test cases to cover different scenarios and edge cases.

import networkx as nx
import matplotlib.pyplot as plt
import os

def min_transition_calculation(min_transition):
    """
    Calculates a list based on the minimum transition time provided values and returns it in factors 1, 2, 5, 10.

    Parameters
    ----------
    min_transition : float or int
        The minimum tranisiton time input for the generation of the factors.

    Returns
    -------
    list :
        List with the minimum transition time with factors 1, 2, 5, 10.
    """
    min_transitions = [min_transition, min_transition*2, min_transition*5, min_transition*10]
    
    return min_transitions

def binding_site_markov_network(total_frames, min_transitions, combined_dict, font_size=None, size_node=None):
    """
    Generate Markov Chain plots based on transition probabilities.

    Parameters
    ----------
    total_frames : int
        The number of frames in the protein-ligand MD simulation.
    min_transitions : list of int or float
        list of transition tresholds in %. A Markov Chain plot will be generated for each of the tresholds.
    combined_dict : dict 
        A dictionary with the information of the Binding Modes and their order of appearance during the simulation for all frames.
    font_size (optional) : int
        The font size for the node labels. The default value is set to 12.
    size_node (optional) : int
        The size of the nodes in the Markov Chain plot. the default value is set to 200.

    Returns
    -------
    None
    """
    font_size = 12 if font_size is None else font_size
    size_node = 200 if size_node is None else size_node

    # Calculate the number of elements in each part
    total_length = len(combined_dict['all'])
    part_length = total_length // 3
    remaining_length = total_length % 3

    # Divide the 'all_data' into three parts
    part1_length = part_length + remaining_length
    part2_length = part_length
    part1_data = combined_dict['all'][:part1_length]
    part2_data = combined_dict['all'][part1_length:part1_length + part2_length]
    part3_data = combined_dict['all'][part1_length + part2_length:]

    # Count the occurrences of each node in each part
    part1_node_occurrences = {node: part1_data.count(node) for node in set(part1_data)}
    part2_node_occurrences = {node: part2_data.count(node) for node in set(part2_data)}
    part3_node_occurrences = {node: part3_data.count(node) for node in set(part3_data)}

    # Get the top 10 nodes with the most occurrences
    node_occurrences = {node: combined_dict['all'].count(node) for node in set(combined_dict['all'])}
    top_10_nodes = sorted(node_occurrences, key=node_occurrences.get, reverse=True)[:10]

    for min_transition_percent in min_transitions:
        min_prob = min_transition_percent / 100  # Convert percentage to probability

        # Create a directed graph
        G = nx.DiGraph()

        # Count the occurrences of each transition and self-loop
        transitions = {}
        self_loops = {}
        for i in range(len(combined_dict['all']) - 1):
            current_state = combined_dict['all'][i]
            next_state = combined_dict['all'][i + 1]

            if current_state == next_state:  # Check for self-loop
                self_loop = (current_state, next_state)
                self_loops[self_loop] = self_loops.get(self_loop, 0) + 1
            else:
                transition = (current_state, next_state)
                transitions[transition] = transitions.get(transition, 0) + 1

        # Add edges to the graph with their probabilities
        for transition, count in transitions.items():
            current_state, next_state = transition
            probability = count / len(combined_dict['all']) * 100  # Convert probability to percentage
            if probability >= min_transition_percent:
                G.add_edge(current_state, next_state, weight=probability)

        # Add self-loops to the graph with their probabilities
        for self_loop, count in self_loops.items():
            state = self_loop[0]
            probability = count / len(combined_dict['all']) * 100  # Convert probability to percentage
            if probability >= min_transition_percent:
                G.add_edge(state, state, weight=probability)

        # Calculate transition probabilities for each direction (excluding self-loops)
        transition_probabilities_forward = {}
        transition_probabilities_backward = {}

        for transition, count in transitions.items():
            start_state, end_state = transition
            forward_transition = (start_state, end_state)
            backward_transition = (end_state, start_state)

            transition_probabilities_forward[forward_transition] = count / node_occurrences[start_state]
            transition_probabilities_backward[backward_transition] = count / node_occurrences[end_state]

        # Calculate self-loop probabilities
        self_loop_probabilities = {}
        for self_loop, count in self_loops.items():
            state = self_loop[0]
            self_loop_probabilities[state] = count / node_occurrences[state]

        # Generate the Markov Chain plot
        plt.figure(figsize=(30, 30))  # Increased figure size
        plt.title(f"Markov Chain Plot {min_transition_percent}% Frames Transition", fontsize=35)

        # Draw nodes and edges
        pos = nx.spring_layout(G, k=2, seed=42)  # Increased distance between nodes (k=2)
        edge_colors = []

        for u, v, data in G.edges(data=True):
            weight = data['weight']

            if u == v:  # Check if it is a self-loop
                edge_colors.append('green')  # Set green color for self-loop arrows
                width = 0.1  # Make self-loop arrows smaller
                connection_style = 'arc3,rad=-0.1'  # Make the self-loops more curved and compact
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=0.2, edge_color=edge_colors[-1],
                                       connectionstyle=connection_style)
            elif weight >= min_transition_percent:
                edge_colors.append('black')  # Highlight significant transitions in red
                width = 4.0
                edge_label = f"{transition_probabilities_forward.get((u, v), 0):.2f}% of Frames →\n{transition_probabilities_backward.get((v, u), 0):.2f}% of Frames ←"
                connection_style = 'arc3,rad=-0.1'
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=0.7, edge_color=edge_colors[-1], connectionstyle=connection_style)
                nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): edge_label}, font_size=15)
            else:
                edge_colors.append('grey')  # Use black for non-significant transitions
                width = 0.5
                edge_label = f"{transition_probabilities_forward.get((u, v), 0):.2f}% of Frames →\n{transition_probabilities_backward.get((v, u), 0):.2f}% of Frames ←"
                connection_style = 'arc3,rad=-0.1'
                nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=0.7, edge_color=edge_colors[-1], connectionstyle=connection_style)
                nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): edge_label}, font_size=12)

        # Update the node colors based on their appearance percentages in each part
        node_colors = []
        for node in G.nodes():
            if node in top_10_nodes:
                part1_percentage = part1_node_occurrences.get(node, 0) / node_occurrences[node]
                part2_percentage = part2_node_occurrences.get(node, 0) / node_occurrences[node]
                part3_percentage = part3_node_occurrences.get(node, 0) / node_occurrences[node]

                if part1_percentage > 0.5:
                    node_colors.append('green')
                elif part2_percentage > 0.5:
                    node_colors.append('orange')
                elif part3_percentage > 0.5:
                    node_colors.append('red')
                else:
                    node_colors.append('yellow')
            else:
                node_colors.append('skyblue')

        # Draw nodes with sizes correlated to occurrences and color top 10 nodes
        node_size = [size_node * node_occurrences[node] for node in G.nodes()]
        nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_colors, alpha=0.8)

        # Draw node labels with occurrence percentage and self-loop probability for top 10 nodes
        node_labels = {}
        for node in G.nodes():
            if node in top_10_nodes:
                node_occurrence_percentage = node_occurrences[node] / len(combined_dict['all']) * 100
                self_loop_probability = self_loop_probabilities.get(node, 0) * 100
                node_label = f"{node}\nOccurrences: {node_occurrence_percentage:.2f}%\nSelf-Loop Probability: {self_loop_probability:.2f}%"
            else:
                node_label = node
            node_labels[node] = node_label

        nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=font_size, font_color='black', verticalalignment="center")

        plt.axis('off')
        plt.tight_layout()

        # Save the plot as a PNG file
        plot_filename = f"markov_chain_plot_{min_transition_percent}.png"
        plot_path = os.path.join("Binding_Modes_Markov_States", plot_filename)
        os.makedirs("Binding_Modes_Markov_States", exist_ok=True)  # Create the folder if it doesn't exist

        plt.savefig(plot_path, dpi=300)
        plt.clf()  # Clear the current figure for the next iteration

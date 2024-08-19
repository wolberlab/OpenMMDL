import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
from typing import Dict, List, Tuple, Union

class MarkovChainAnalysis:
    def __init__(self, min_transition: float):
        self.min_transition = min_transition
        self.min_transitions: List[float] = self.calculate_min_transitions()

    def calculate_min_transitions(self) -> List[float]:
        """Calculates a list based on the minimum transition time provided values and returns it in factors 1, 2, 5, 10.

        Returns:
            List[float]: List with the minimum transition time with factors 1, 2, 5, 10.
        """
        min_transitions = [
            self.min_transition,
            self.min_transition * 2,
            self.min_transition * 5,
            self.min_transition * 10,
        ]
        return min_transitions

    def generate_transition_graph(
        self,
        total_frames: int,
        combined_dict: Dict[str, List[str]],
        fig_type: str = "png",
        font_size: int = 36,
        size_node: int = 200
    ) -> None:
        """Generate Markov Chain plots based on transition probabilities.

        Args:
            total_frames (int): The number of frames in the protein-ligand MD simulation.
            combined_dict (Dict[str, List[str]]): A dictionary with the information of the Binding Modes and their order of appearance during the simulation for all frames.
            fig_type (str, optional): File type for the output figures. Default is 'png'.
            font_size (int, optional): The font size for the node labels. The default value is set to 36.
            size_node (int, optional): The size of the nodes in the Markov Chain plot. The default value is set to 200.
        """
        # Calculate the number of elements in each part
        total_length: int = len(combined_dict["all"])
        part_length: int = total_length // 3
        remaining_length: int = total_length % 3

        # Divide the 'all_data' into three parts
        part1_length: int = part_length + remaining_length
        part2_length: int = part_length
        part1_data: List[str] = combined_dict["all"][:part1_length]
        part2_data: List[str] = combined_dict["all"][part1_length: part1_length + part2_length]
        part3_data: List[str] = combined_dict["all"][part1_length + part2_length:]

        # Count the occurrences of each node in each part
        part1_node_occurrences: Dict[str, int] = {
            node: part1_data.count(node) for node in set(part1_data)
        }
        part2_node_occurrences: Dict[str, int] = {
            node: part2_data.count(node) for node in set(part2_data)
        }
        part3_node_occurrences: Dict[str, int] = {
            node: part3_data.count(node) for node in set(part3_data)
        }

        # Create the legend
        legend_labels: Dict[str, str] = {
            "Blue": "Binding mode not in top 10 occurrence",
            "Green": "Binding Mode occurrence mostly in first third of frames",
            "Orange": "Binding Mode occurrence mostly in second third of frames",
            "Red": "Binding Mode occurrence mostly in third third of frames",
            "Yellow": "Binding Mode occurs throughout all trajectory equally",
        }

        legend_colors: List[str] = ["skyblue", "green", "orange", "red", "yellow"]

        legend_handles: List[Patch] = [
            Patch(color=color, label=label)
            for color, label in zip(legend_colors, legend_labels.values())
        ]

        # Get the top 10 nodes with the most occurrences
        node_occurrences: Dict[str, int] = {
            node: combined_dict["all"].count(node) for node in set(combined_dict["all"])
        }
        top_10_nodes: List[str] = sorted(node_occurrences, key=node_occurrences.get, reverse=True)[:10]

        for min_transition_percent in self.min_transitions:
            min_prob: float = min_transition_percent / 100  # Convert percentage to probability

            # Create a directed graph
            G: nx.DiGraph = nx.DiGraph()

            # Count the occurrences of each transition and self-loop
            transitions: Dict[Tuple[str, str], int] = {}
            self_loops: Dict[Tuple[str, str], int] = {}
            for i in range(len(combined_dict["all"]) - 1):
                current_state: str = combined_dict["all"][i]
                next_state: str = combined_dict["all"][i + 1]

                if current_state == next_state:  # Check for self-loop
                    self_loop: Tuple[str, str] = (current_state, next_state)
                    self_loops[self_loop] = self_loops.get(self_loop, 0) + 1
                else:
                    transition: Tuple[str, str] = (current_state, next_state)
                    transitions[transition] = transitions.get(transition, 0) + 1

            # Add edges to the graph with their probabilities
            for transition, count in transitions.items():
                current_state, next_state = transition
                probability: float = (count / len(combined_dict["all"])) * 100  # Convert probability to percentage
                if probability >= min_transition_percent:
                    G.add_edge(current_state, next_state, weight=probability)
                    # Include the reverse transition with a different color
                    reverse_transition: Tuple[str, str] = (next_state, current_state)
                    reverse_count: int = transitions.get(reverse_transition, 0)  # Use the correct count for the reverse transition
                    reverse_probability: float = (reverse_count / len(combined_dict["all"])) * 100
                    G.add_edge(next_state, current_state, weight=reverse_probability, reverse=True)

            # Add self-loops to the graph with their probabilities
            for self_loop, count in self_loops.items():
                state: str = self_loop[0]
                probability: float = (count / len(combined_dict["all"])) * 100  # Convert probability to percentage
                if probability >= min_transition_percent:
                    G.add_edge(state, state, weight=probability)

            # Calculate transition probabilities for each direction (excluding self-loops)
            transition_probabilities_forward: Dict[Tuple[str, str], float] = {}
            transition_probabilities_backward: Dict[Tuple[str, str], float] = {}
            transition_occurrences_forward: Dict[Tuple[str, str], float] = {}
            transition_occurrences_backward: Dict[Tuple[str, str], float] = {}

            for transition, count in transitions.items():
                start_state, end_state = transition
                forward_transition: Tuple[str, str] = (start_state, end_state)
                backward_transition: Tuple[str, str] = (end_state, start_state)

                # Separate counts for forward and backward transitions
                forward_count: int = transitions.get(forward_transition, 0)
                backward_count: int = transitions.get(backward_transition, 0)

                transition_probabilities_forward[forward_transition] = (forward_count / node_occurrences[start_state]) * 100
                transition_occurrences_forward[forward_transition] = (forward_count / len(combined_dict["all"])) * 100

                transition_probabilities_backward[backward_transition] = (backward_count / node_occurrences[end_state]) * 100
                transition_occurrences_backward[backward_transition] = (backward_count / len(combined_dict["all"])) * 100

            # Calculate self-loop probabilities
            self_loop_probabilities: Dict[str, float] = {}
            self_loop_occurences: Dict[str, float] = {}
            for self_loop, count in self_loops.items():
                state: str = self_loop[0]
                self_loop_probabilities[state] = count / node_occurrences[state]
                self_loop_occurences[state] = count / len(combined_dict["all"]) * 100

            # Generate the Markov Chain plot
            plt.figure(figsize=(60, 60))  # Increased figure size
            plt.title(f"Markov Chain Plot {min_transition_percent}% Frames Transition", fontsize=72)

            # Draw nodes and edges
            pos: Dict[str, Tuple[float, float]] = nx.spring_layout(G, k=2, seed=42)  # Increased distance between nodes (k=2)
            edge_colors: List[str] = []

            for u, v, data in G.edges(data=True):
                weight: float = data["weight"]

                if u == v:  # Check if it is a self-loop
                    edge_colors.append("green")  # Set green color for self-loop arrows
                    width: float = 0.1  # Make self-loop arrows smaller
                    connection_style: str = "arc3,rad=-0.1"  # Make the self-loops more curved and compact
                    nx.draw_networkx_edges(
                        G,
                        pos,
                        edgelist=[(u, v)],
                        width=width,
                        alpha=0.2,
                        edge_color=edge_colors[-1],
                        connectionstyle=connection_style,
                    )
                elif weight >= min_transition_percent:
                    edge_colors.append("black")  # Highlight significant transitions in black

                    # Check if both nodes are present before adding labels
                    if G.has_node(u) and G.has_node(v):
                        width: float = 4.0
                        forward_label: str = f"{transition_occurrences_forward.get((v, u), 0):.2f}% of Frames →, {transition_probabilities_forward.get((v, u), 0):.2f}% probability"
                        backward_label: str = f"{transition_occurrences_backward.get((u, v), 0):.2f}% of Frames ←, {transition_probabilities_backward.get((u, v), 0):.2f}% probability"
                        edge_label: str = f"{forward_label}\n{backward_label}"

                        connection_style: str = "arc3,rad=-0.1"
                        nx.draw_networkx_edges(
                            G,
                            pos,
                            edgelist=[(u, v)],
                            width=width,
                            alpha=0.7,
                            edge_color=edge_colors[-1],
                            connectionstyle=connection_style,
                        )
                        nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): edge_label}, font_size=26)
                else:
                    edge_colors.append("grey")  # Use grey for non-significant transitions
                    width: float = 0.5

                    # Check if both nodes are present before adding labels
                    if G.has_node(u) and G.has_node(v):
                        forward_label: str = f"{transition_occurrences_forward.get((v, u), 0):.2f}% of Frames →, {transition_probabilities_forward.get((v, u), 0):.2f}% probability"
                        backward_label: str = f"{transition_occurrences_backward.get((u, v), 0):.2f}% of Frames ←, {transition_probabilities_backward.get((u, v), 0):.2f}% probability"
                        edge_label: str = f"{forward_label}\n{backward_label}"

                        connection_style: str = "arc3,rad=-0.1"
                        nx.draw_networkx_edges(
                            G,
                            pos,
                            edgelist=[(u, v)],
                            width=width,
                            alpha=0.7,
                            edge_color=edge_colors[-1],
                            connectionstyle=connection_style,
                        )
                        nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): edge_label}, font_size=36)

            # Update the node colors based on their appearance percentages in each part
            node_colors: List[str] = []
            for node in G.nodes():
                if node in top_10_nodes:
                    part1_percentage: float = (part1_node_occurrences.get(node, 0) / node_occurrences[node])
                    part2_percentage: float = (part2_node_occurrences.get(node, 0) / node_occurrences[node])
                    part3_percentage: float = (part3_node_occurrences.get(node, 0) / node_occurrences[node])

                    if part1_percentage > 0.5:
                        node_colors.append("green")
                    elif part2_percentage > 0.5:
                        node_colors.append("orange")
                    elif part3_percentage > 0.5:
                        node_colors.append("red")
                    else:
                        node_colors.append("yellow")
                else:
                    node_colors.append("skyblue")

            # Draw nodes with sizes correlated to occurrences and color top 10 nodes
            node_size: List[int] = [size_node * node_occurrences[node] for node in G.nodes()]
            nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_colors, alpha=0.8)

            # Draw node labels with occurrence percentage and self-loop probability for nodes with edges
            node_labels: Dict[str, str] = {}

            for node in G.nodes():
                if G.degree(node) > 0:  # Check if the node has at least one edge
                    edges_with_node: List[Tuple[str, str, float]] = [
                        (u, v, data["weight"])
                        for u, v, data in G.edges(data=True)
                        if u == node or v == node
                    ]
                    relevant_edges: List[Tuple[str, str, float]] = [
                        edge
                        for edge in edges_with_node
                        if edge[2] >= min_transition_percent
                    ]

                    if relevant_edges:
                        if node in top_10_nodes:
                            node_occurrence_percentage: float = (node_occurrences[node] / len(combined_dict["all"])) * 100
                            self_loop_probability: float = self_loop_probabilities.get(node, 0) * 100
                            self_loop_occurence: float = self_loop_occurences.get(node, 0)
                            node_label: str = f"{node}\nOccurrences: {node_occurrence_percentage:.2f}%\nSelf-Loop Probability: {self_loop_probability:.2f}% \nSelf-Loop Occurrence: {self_loop_occurence:.2f}%"
                            node_labels[node] = node_label
                        else:
                            node_labels[node] = node

            nx.draw_networkx_labels(
                G,
                pos,
                labels=node_labels,
                font_size=font_size,
                font_color="black",
                verticalalignment="center",
            )

            # Add the legend to the plot
            plt.legend(handles=legend_handles, loc="upper right", fontsize=48)

            plt.axis("off")
            plt.tight_layout()

            # Save the plot
            plot_filename: str = f"markov_chain_plot_{min_transition_percent}.{fig_type}"
            plot_path: str = os.path.join("Binding_Modes_Markov_States", plot_filename)
            os.makedirs("Binding_Modes_Markov_States", exist_ok=True)  # Create the folder if it doesn't exist

            plt.savefig(plot_path, dpi=300)
            plt.clf()  # Clear the current figure for the next iteration

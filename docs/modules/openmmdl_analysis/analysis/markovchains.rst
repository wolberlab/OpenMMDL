API Documentation for markovchains
==================================

.. py:class:: MarkovChainAnalysis(min_transition)

    Class for analyzing binding mode transitions using Markov chains.

    This class provides functionality to analyze the transitions of binding modes over time in MD simulations.
    It allows the generation of Markov chain plots based on transition probabilities between the binding modes
    and highlights key binding modes based on their occurency in the trajectory.

    :ivar float min_transition: Minimum transition percentage threshold for transitions in the graph.
    
    :ivar list min_transitions: List of transition thresholds [x1, x2, x5, x10] derived from min_transition.


    .. py:method:: _calculate_min_transitions()

        Calculates a list based on the minimum transition time provided values and returns it in factors 1, 2, 5, 10.

        :returns: List with the minimum transition time with factors 1, 2, 5, 10.
        :rtype: list


    .. py:method:: generate_transition_graph(total_frames, combined_dict, fig_type='png', font_size=36, size_node=200)

        Generates Markov chain graphs showing binding mode transitions across a simulation trajectory.

        Nodes represent binding modes, with size proportional to their frequency. Node colors reflect their dominant
        occurrence in one-third of the trajectory. Arrows represent transitions. Dual direction arrows show reversible behavior.

        For each transition threshold in `min_transitions`, a separate plot is generated and saved.

        :param int total_frames: The number of frames in the protein-ligand MD simulation.
        :param dict combined_dict: A dictionary with the information of the Binding Modes and their order of appearance during the simulation for all frames.
        :param str fig_type: File type for the output figures. Default is 'png'.
        :param int font_size: The font size for the node labels. The default value is set to 36.
        :param int size_node: The size of the nodes in the Markov Chain plot. The default value is set to 200.

        :returns: None. Saves one plot per transition threshold in the directory `Binding_Modes_Markov_States`.
        :rtype: None

API Documentation for markovchains
==================================

.. py:class:: MarkovChainAnalysis(min_transition)

    Analyzes binding mode transitions in molecular dynamics simulations using Markov chain models.

    This class enables the construction of Markov chain graphs to visualize transition probabilities and trends
    in binding mode behavior over the course of a simulation.

    :param float min_transition: Minimum transition percentage threshold for transitions in the graph.
    
    :ivar list min_transitions: List of transition thresholds [x1, x2, x5, x10] derived from min_transition.


    .. py:method:: calculate_min_transitions()

        Calculates a list based on the minimum transition time provided values and returns it in factors 1, 2, 5, 10.

        :returns: List with the minimum transition time with factors 1, 2, 5, 10.
        :rtype: list


    .. py:method:: generate_transition_graph(total_frames, combined_dict, fig_type='png', font_size=36, size_node=200)

        Generates Markov chain graphs showing binding mode transitions across a simulation trajectory.

        Nodes represent binding modes, with size proportional to their frequency. Node colors reflect their dominant
        occurrence in one-third of the trajectory. Arrows represent transitions; dual direction arrows show reversible behavior.
        Self-loops are drawn for modes that persist over time and transition in themselves.

        For each transition threshold in `min_transitions`, a separate plot is generated and saved.

        :param int total_frames: Total number of frames in the simulation.
        :param dict combined_dict: Dictionary containing the frame-by-frame binding mode assignments. Expects a key `'all'`.
        :param str fig_type: File format for the output plots. Default is `'png'`.
        :param int font_size: Font size for node labels. Default is `36`.
        :param int size_node: Base size for graph nodes, scaled by mode occurrence. Default is `200`.

        :returns: Saves one plot per transition threshold in the directory `Binding_Modes_Markov_States`.
        :rtype: None

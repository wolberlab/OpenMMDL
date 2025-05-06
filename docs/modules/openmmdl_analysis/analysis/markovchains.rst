API Documentation for markovchains
============================


.. py:class:: MarkovChainAnalysis(min_transition)

    Performs Markov Chain analysis of binding mode transitions in molecular dynamics simulations.

    :param int min_transition: Minimum transition percentage used to threshold transitions for plotting.

    :ivar list min_transitions: A list of transition thresholds derived from the input using scaling factors (1×, 2×, 5×, 10×).


    .. py:method:: calculate_min_transitions()

        Calculates a list of minimum transition thresholds using multipliers (1, 2, 5, 10) on the `min_transition` value.

        :returns: List of threshold percentages to be used in Markov Chain analysis.
        :rtype: list


    .. py:method:: generate_transition_graph(total_frames, combined_dict, fig_type='png', font_size=36, size_node=200)

        Generates and saves Markov Chain visualizations showing state transitions and self-loops across a molecular dynamics trajectory.

        :param int total_frames: Total number of frames in the MD simulation.
        :param dict combined_dict: Dictionary with binding mode data, where `combined_dict["all"]` contains the frame-by-frame mode sequence.
        :param str fig_type: File format for the output plots. Defaults to `'png'`.
        :param int font_size: Font size for the graph node labels. Defaults to `36`.
        :param int size_node: Base node size in the plot. Scaled by number of occurrences. Defaults to `200`.

        :returns: Saves figures to disk in the directory `Binding_Modes_Markov_States` with one image per transition threshold.
        :rtype: None

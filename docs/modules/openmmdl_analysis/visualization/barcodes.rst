API Documentation for barcodes
============================

.. py:class:: BarcodeGenerator(df)

    Generates binary barcodes representing the presence of interactions across MD frames.

    :param pd.DataFrame df: DataFrame containing all interactions extracted from PLIP analysis.
    :ivar dict interactions: Dictionary mapping interaction types to their corresponding columns in the DataFrame.

    .. py:method:: generate_barcode(interaction)

        Generates barcodes for a given interaction.

        :param str interaction: Name of the interaction to generate a barcode for.
        :return: Binary array with 1 representing the interaction is present in the corresponding frame.
        :rtype: np.ndarray

    .. py:method:: interacting_water_ids(waterbridge_interactions)

        Generates a list of all water ids that form water bridge interactions.

        :param list waterbridge_interactions: list containing the names of all water bridge interactions.
        :return: Unique water IDs that form waterbridge interactions.
        :rtype: list

    .. py:method:: _generate_waterids_barcode(interaction)

        Generates a barcode containing corresponding water ids for a given interaction.

        :param str interaction: Name of the interaction to generate a barcode for.
        :return: List with water IDs for frames where interaction is present, 0 otherwise.
        :rtype: list

        .. note::
           Water IDs are obtained from the "WATER_IDX" column in the DataFrame.

    .. py:method:: _gather_interactions()

        Gathers interaction column names grouped by the corresponding interaction type.

        :return: Dictionary where the keys are interaction types and values are lists of corresponding column names.
        :rtype: dict

.. py:class:: BarcodePlotter(df_all)

    Visualizes interaction barcodes and waterbridge statistics using bar plots and pie charts.

    :param pd.DataFrame df_all: Full interaction dataframe passed for plotting.
    :ivar BarcodeGenerator barcode_gen: Instance of BarcodeGenerator used to compute barcodes from interaction data.

    .. py:method:: plot_waterbridge_piechart(waterbridge_barcodes, waterbridge_interactions, fig_type)

        Generates and saves pie charts showing the frequency of each of the water IDs participating in waterbridge interactions.

        :param dict waterbridge_barcodes: Dictionary of waterbridge interaction barcodes.
        :param list waterbridge_interactions: List of interaction column names related to waterbridge interactions.
        :param str fig_type: Image file format for saving (e.g., 'png', 'svg').
        :return: None. This function writes out a figure and does not return anything.
        :rtype: None

    .. py:method:: plot_barcodes_grouped(interactions, interaction_type, fig_type)

        Groups barcodes by ligand atom, plots individual and grouped barcodes, and saves them.

        :param list interactions: List of interaction names to be grouped and visualized.
        :param str interaction_type: The type of interaction (e.g., 'donor', 'acceptor', 'waterbridge').
        :param str fig_type: Image file format for saving (e.g., 'png', 'svg').
        :return: None. This function writes out a figure and does not return anything.
        :rtype: None

    .. py:method:: _plot_barcodes(barcodes, save_path)

        Plots barcodes of the interactions depending on the presence of the interaction.

        :param dict barcodes: Dictionary where keys are interaction names and values are 1D numpy arrays (barcodes).
        :param str save_path: Path to save the generated barcode plot image.
        :return: None. This function writes out a figure and does not return anything.
        :rtype: None

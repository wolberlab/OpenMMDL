API Documentation for figures
============================


.. py:class:: FigureMerger(binding_mode, occurrence_percent, split_data, merged_image_paths)

    Handles the creation and merging of binding mode figures with corresponding legends.

    :param str binding_mode: Name of the binding mode for which figures are created.
    :param float occurrence_percent: Occurrence percentage of the binding mode.
    :param list split_data: Interaction descriptors used for generating the figure legend.
    :param list merged_image_paths: List storing paths to the output merged images.

    .. py:method:: create_and_merge_images()

        Create and merge images to generate a legend for binding modes.

        :returns: Updated list of paths to the merged images.
        :rtype: list of str


.. py:class:: FigureArranger(merged_image_paths, output_path)

    Arranges multiple merged binding mode figures into a single image.

    :param list merged_image_paths: List of file paths to pre merged figures.
    :param str output_path: File path where the final arranged figure should be saved.

    .. py:method:: arranged_figure_generation()

        Generate an arranged figure by arranging merged images in rows and columns.

        :returns: None. This function writes out a figure and does not return anything.
        :rtype: None

API Documentation for FigureMerger
============================


.. py:class:: FigureMerger(binding_mode, occurrence_percent, split_data, merged_image_paths)

    A class to generate and merge visualization images representing binding modes and their associated legends.

    :param str binding_mode: Label or identifier for the binding mode.
    :param float occurrence_percent: Percentage occurrence of the binding mode.
    :param list split_data: Raw data describing interactions to be visualized.
    :param list merged_image_paths: List to store paths to the merged images.

    .. py:method:: create_and_merge_images()

        Generates line plots for binding mode interactions, creates a corresponding legend, merges them side by side, and saves the result.

        Returns:
            list: List of file paths to the saved merged images.

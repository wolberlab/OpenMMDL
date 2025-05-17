API Documentation for visualization
============================

.. py:class:: Visualizer(md, cloud_path, ligname, special)

    Visualizes the OpenMMDL Analysis output of the trajectories with annotated interaction clouds and water molecule interactions.

    The class provides tools to render 3D molecular scenes using NGLview, including trajectory
    visualization, ligand highlighting, interacting water molecules, and various interaction cloud types
    (e.g., hydrogen bonds, hydrophobic regions, etc.).


    :param str cloud_path: Path to the location where the JSON file containing interaction cloud coordinates is.
    :ivar mda.Universe md: The loaded trajectory and topology data used for visualization.
    :ivar dict cloud: Dictionary parsed from the JSON file representing interaction clouds (e.g., hydrophobic, donor, etc.).
    :ivar str ligname: Ligand selection string used in NGLview for highlighting.
    :ivar str special: Additional selection string for visual emphasis (can be None).

    .. py:method:: _load_cloud(cloud_path)

        Loads interaction clouds from a JSON file.

        :param str cloud_path: Path to the JSON file containing the interaction clouds.
        :return: Dictionary representing different types of interaction clouds and their properties.
        :rtype: dict

    .. py:method:: visualize(receptor_type='protein or nucleic', height='1000px', width='1000px')

        Generates visualization of the trajectory with the interacting waters and interaction clouds.

        :param str receptor_type: Type of receptor. Defaults to 'protein or nucleic'.
        :param str height: Height of the visualization. Defaults to '1000px'.
        :param str width: Width of the visualization. Defaults to '1000px'.
        :return: Returns an nglview.widget object containing the visualization
        :rtype: nglview.NGLWidget


.. py:function:: run_visualization()

    Runs the visualization notebook in the current directory. 
    The visualization notebook is copied from the package directory to the current directory and automaticaly started.

    :return: None
    :rtype: None

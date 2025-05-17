API Documentation for pharmacophore
============================


.. py:class:: PharmacophoreGenerator(df_all, ligand_name)

    A class for generating pharmacophore features from molecular dynamics (MD) interaction data.

    This class processes interaction data from MD simulations and generates pharmacophore features 
    including hydrogen bond donors/acceptors, hydrophobic features etc. exporting them in .pml format.

    :param pd.DataFrame df_all: DataFrame storing the input interaction data.
    :param str ligand_name: Name of the ligand.
    :ivar str complex_name: Name of the complex consisting of ligand name and "_complex".
    :param re.Pattern coord_pattern: Regular expression pattern for extracting 3D coordinates from strings.
    :param dict clouds: Dictionary containing interaction types and associated 3D coordinates with visualization metadata.

    .. py:method:: to_dict()

        Export the interaction cloud as a dictionary.

        :return: Dictionary representation of the interaction cloud.
        :rtype: dict

    .. py:method:: generate_md_pharmacophore_cloudcenters(output_filename, id_num=0)

        Generates pharmacophore from all interactions formed in the MD simulation.
        A feature is generated for each interaction at the center of all its ocurrences.

        :param str output_filename: Name the of the output .pml file.
        :param int id_num: ID number as an identifier in the PML file. Defaults to 0.
        :return: None. This function writes output directly to a .pml XML file and does not return anything.
        :rtype: None

    .. py:method:: generate_point_cloud_pml(outname)

        Generates pharmacophore point cloud and writes it to a .pml file.

        :param str outname: Name of the output .pml file.
        :return: None. This function writes output directly to a .pml file and does not return anything.
        :rtype: None

    .. py:method:: generate_bindingmode_pharmacophore(dict_bindingmode, outname, id_num=0)

        Generates pharmacophore from a binding mode and writes it to a .pml file.

        :param dict dict_bindingmode: Dictionary containing all interactions of the bindingmode and their coresponding ligand and protein coordinates.
        :param str outname: Name of the output .pml file.
        :param int id_num: if multiple id number can enumerate the diferent bindingmodes. Defaults to 0.
        :return: None. This function writes output directly to a .pml file and does not return anything.
        :rtype: None

    .. py:method:: _generate_clouds()
    
        Process the dataframe with the interactions to extract and categorize ligand/protein interactions.
    
        :return: A dict containing interaction types as keys and their 3D coordinates.
        :rtype: dict

    .. py:method:: _format_clouds(interaction_coords)
    
        Add visualization properties (color, radius) to the interactions.
    
        :param dict interaction_coords: Dictionary of raw 3D coordinates grouped by interaction type.
        :return: Dictionary formatted with interaction type as key, and a dictionary of coordinates, color, and radius as value.
        :rtype: dict

    .. py:method:: _generate_pharmacophore_centers(interactions)
    
        Generates pharmacophore points for interactions that are points such as hydrophobic and ionic interactions.
    
        :param list interactions: List of interactions to generate pharmacophore from.
        :return: Dict of interactions from which pharmacophore is generated as key and list of coordinates as value.
        :rtype: dict

    .. py:method:: _generate_pharmacophore_vectors(interactions)
    
        Generates pharmacophore points for interactions that are vectors such as hydrogen bond donors or acceptors.
    
        :param list interactions: List of interactions to generate pharmacophore vectors from.
        :return: Dict of interactions from which pharmacophore is generated as key and list of coordinates as value (first coords are ligand side, second are protein side).
        :rtype: dict

    .. py:method:: _generate_pharmacophore_centers_all_points(interactions)
    
        Generates pharmacophore points for all interactions to generate point cloud.
    
        :param list interactions: List of interactions to generate pharmacophore from.
        :return: Dict of interactions from which pharmacophore is generated as key and list of coordinates as value.
        :rtype: dict

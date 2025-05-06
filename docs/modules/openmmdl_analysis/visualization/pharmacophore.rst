API Documentation for MarkovChainAnalysis
============================


.. py:class:: PharmacophoreGenerator(df_all, ligand_name)

    A class to extract and organize pharmacophoric features from protein-ligand interaction data (e.g., PLIP analysis results).

    :param pandas.DataFrame df_all: DataFrame containing interaction data.
    :param str ligand_name: Name of the ligand for labeling outputs.

    .. py:method:: to_dict()

        Returns the pharmacophore clouds dictionary with interaction features and their spatial metadata.

        :return: Dictionary of pharmacophore interaction features.
        :rtype: dict

    .. py:method:: generate_pharmacophore_centers(interactions)

        Generates pharmacophore points (centroids) for interactions characterized by a single spatial location, such as hydrophobic or ionic contacts.

        :param list interactions: List of interaction types to process.
        :return: Dictionary of interaction types with corresponding centroid coordinates.
        :rtype: dict

    .. py:method:: generate_pharmacophore_vectors(interactions)

        Generates pharmacophore vector features for directional interactions, such as hydrogen bonds.

        :param list interactions: List of interaction types to process.
        :return: Dictionary mapping interaction types to ligand and protein coordinates.
        :rtype: dict

    .. py:method:: generate_md_pharmacophore_cloudcenters(output_filename, id_num=0)

        Generates and writes a pharmacophore model from MD-derived interaction data to a `.pml` XML file, with support for vector, point, and planar features.

        :param str output_filename: Path (without extension) to save the output `.pml` file.
        :param int id_num: Optional ID to tag the pharmacophore set. Default is 0.

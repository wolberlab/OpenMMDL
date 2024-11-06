import re
import pandas as pd
import xml.etree.ElementTree as ET
import numpy as np
from typing import Dict, List, Optional, Union


class PharmacophoreGenerator:
    def __init__(self, df_all: pd.DataFrame, ligand_name: str):
        self.df_all = df_all
        self.ligand_name = ligand_name
        self.complex_name = f"{ligand_name}_complex"
        self.coord_pattern = re.compile(r"\(([\d.-]+), ([\d.-]+), ([\d.-]+)\)")
        self.clouds = self._generate_clouds()

    def _generate_clouds(
        self,
    ) -> Dict[str, Dict[str, Union[List[List[float]], List[float], float]]]:
        interaction_coords = {
            "hydrophobic": [],
            "acceptor": [],
            "donor": [],
            "waterbridge": [],
            "negative_ionizable": [],
            "positive_ionizable": [],
            "pistacking": [],
            "pication": [],
            "halogen": [],
            "metal": [],
        }

        for index, row in self.df_all.iterrows():
            if row["LIGCOO"] != 0:
                coord_match = self.coord_pattern.match(row["LIGCOO"])
                if coord_match:
                    x, y, z = map(float, coord_match.groups())
                    x, y, z = round(x, 3), round(y, 3), round(z, 3)
                    interaction = row["INTERACTION"]
                    if interaction == "hbond":
                        interaction = (
                            "donor" if row["PROTISDON"] == "False" else "acceptor"
                        )
                    if interaction == "saltbridge":
                        interaction = (
                            "negative_ionizable"
                            if row["PROTISPOS"] == "True"
                            else "positive_ionizable"
                        )
                    if interaction in interaction_coords:
                        interaction_coords[interaction].append([x, y, z])

        for index, row in self.df_all.iterrows():
            if row["TARGETCOO"] != 0:
                coord_match = self.coord_pattern.match(row["TARGETCOO"])
                if coord_match:
                    x, y, z = map(float, coord_match.groups())
                    x, y, z = round(x, 3), round(y, 3), round(z, 3)
                    if row["INTERACTION"] == "metal":
                        interaction_coords["metal"].append([x, y, z])

        return self._format_clouds(interaction_coords)

    def _format_clouds(
        self, interaction_coords: Dict[str, List[List[float]]]
    ) -> Dict[str, Dict[str, Union[List[List[float]], List[float], float]]]:
        color_mapping = {
            "hydrophobic": [1.0, 1.0, 0.0],
            "acceptor": [1.0, 0.0, 0.0],
            "donor": [0.0, 1.0, 0.0],
            "waterbridge": [0.0, 1.0, 0.9],
            "negative_ionizable": [0.0, 0.0, 1.0],
            "positive_ionizable": [1.0, 0.0, 0.0],
            "pistacking": [0.0, 0.0, 1.0],
            "pication": [0.0, 0.0, 1.0],
            "halogen": [1.0, 0.0, 0.9],
            "metal": [1.0, 0.6, 0.0],
        }
        return {
            interaction: {
                "coordinates": coords,
                "color": color_mapping[interaction],
                "radius": 0.1,
            }
            for interaction, coords in interaction_coords.items()
        }

    def to_dict(
        self,
    ) -> Dict[str, Dict[str, Union[List[List[float]], List[float], float]]]:
        return self.clouds

    def generate_pharmacophore_centers(
        self, interactions: List[str]
    ) -> Dict[str, List[float]]:
        """Generates pharmacophore points for interactions that are points such as hydrophobic and ionic interactions.

        Args:
            interactions (List[str]): List of interactions to generate pharmacophore from.

        Returns:
            Dict[str, List[float]]: Dict of interactions from which pharmacophore is generated as key and list of coordinates as value.
        """
        coord_pattern = re.compile(r"\(([\d.-]+), ([\d.-]+), ([\d.-]+)\)")
        pharmacophore = {}
        for interaction in interactions:
            counter = 0
            sum_x, sum_y, sum_z = 0, 0, 0
            for index, row in self.df_all.iterrows():
                if row[interaction] == 1:
                    coord_match = coord_pattern.match(row["LIGCOO"])
                    if coord_match:
                        x, y, z = map(float, coord_match.groups())
                        sum_x += x
                        sum_y += y
                        sum_z += z
                        counter += 1

            if counter > 0:
                center_x = round((sum_x / counter), 3)
                center_y = round((sum_y / counter), 3)
                center_z = round((sum_z / counter), 3)
                pharmacophore[interaction] = [center_x, center_y, center_z]
        return pharmacophore

    def generate_pharmacophore_vectors(
        self, interactions: List[str]
    ) -> Dict[str, List[List[float]]]:
        """Generates pharmacophore points for interactions that are vectors such as hydrogen bond donors or acceptors.

        Args:
            interactions (List[str]): List of interactions to generate pharmacophore from.

        Returns:
            Dict[str, List[List[float]]]: Dict of interactions from which pharmacophore is generated as key and list of coordinates as value (first coords are ligand side, second are protein side).
        """
        coord_pattern = re.compile(r"\(([\d.-]+), ([\d.-]+), ([\d.-]+)\)")
        pharmacophore = {}
        for interaction in interactions:
            counter = 0
            sum_x, sum_y, sum_z = 0, 0, 0
            sum_a, sum_b, sum_c = 0, 0, 0
            for index, row in self.df_all.iterrows():
                if row[interaction] == 1:
                    coord_match = coord_pattern.match(row["LIGCOO"])
                    if coord_match:
                        x, y, z = map(float, coord_match.groups())
                        sum_x += x
                        sum_y += y
                        sum_z += z
                    coord_match = coord_pattern.match(row["PROTCOO"])
                    if coord_match:
                        a, b, c = map(float, coord_match.groups())
                        sum_a += a
                        sum_b += b
                        sum_c += c
                        counter += 1

            if counter > 0:
                center_x = round((sum_x / counter), 3)
                center_y = round((sum_y / counter), 3)
                center_z = round((sum_z / counter), 3)
                center_a = round((sum_a / counter), 3)
                center_b = round((sum_b / counter), 3)
                center_c = round((sum_c / counter), 3)
                pharmacophore[interaction] = [
                    [center_x, center_y, center_z],
                    [center_a, center_b, center_c],
                ]
        return pharmacophore

    def generate_md_pharmacophore_cloudcenters(
        self, output_filename: str, id_num: int = 0
    ) -> None:
        """Generates pharmacophore from all interactions formed in the MD simulation.
        A feature is generated for each interaction at the center of all its occurrences.

        Args:
            output_filename (str): Name of the output .pml file.
            id_num (int, optional): ID number as an identifier in the PML file. Defaults to 0.
        """
        feature_id_counter = 0

        root = ET.Element(
            "MolecularEnvironment",
            version="0.0",
            id=f"OpennMMDL_Analysis{id_num}",
            name=self.complex_name,
        )
        pharmacophore = ET.SubElement(
            root,
            "pharmacophore",
            name=self.complex_name,
            id=f"pharmacophore{id_num}",
            pharmacophoreType="LIGAND_SCOUT",
        )

        for interaction in [
            "Acceptor_hbond",
            "Donor_hbond",
            "pistacking",
            "hydrophobic",
            "PI_saltbridge",
            "NI_saltbridge",
        ]:
            feature_types = {
                "Acceptor_hbond": "HBA",
                "Donor_hbond": "HBD",
                "pistacking": "AR",
                "hydrophobic": "H",
                "PI_saltbridge": "PI",
                "NI_saltbridge": "NI",
            }
            if interaction not in feature_types:
                continue

            feature_type = feature_types[interaction]
            if interaction in ["Acceptor_hbond", "Donor_hbond"]:
                pharm = self.generate_pharmacophore_vectors(
                    self.df_all.filter(regex=interaction).columns
                )
                for feature_name, position in pharm.items():
                    feature_id_counter += 1
                    lig_loc = position[0]
                    prot_loc = position[1]
                    points_to_lig = "true" if feature_type == "HBA" else "false"
                    hasSyntheticProjectedPoint = "false"
                    vector = ET.SubElement(
                        pharmacophore,
                        "vector",
                        name=feature_type,
                        featureId=feature_name,
                        pointsToLigand=points_to_lig,
                        hasSyntheticProjectedPoint=hasSyntheticProjectedPoint,
                        optional="false",
                        disabled="false",
                        weight="1.0",
                        coreCompound=self.ligand_name,
                        id=f"feature{str(feature_id_counter)}",
                    )
                    if feature_type == "HBA":
                        origin = ET.SubElement(
                            vector,
                            "origin",
                            x3=str(prot_loc[0]),
                            y3=str(prot_loc[1]),
                            z3=str(prot_loc[2]),
                            tolerance="1.9499999",
                        )
                        target = ET.SubElement(
                            vector,
                            "target",
                            x3=str(lig_loc[0]),
                            y3=str(lig_loc[1]),
                            z3=str(lig_loc[2]),
                            tolerance="1.5",
                        )
                    if feature_type == "HBD":
                        origin = ET.SubElement(
                            vector,
                            "origin",
                            x3=str(lig_loc[0]),
                            y3=str(lig_loc[1]),
                            z3=str(lig_loc[2]),
                            tolerance="1.9499999",
                        )
                        target = ET.SubElement(
                            vector,
                            "target",
                            x3=str(prot_loc[0]),
                            y3=str(prot_loc[1]),
                            z3=str(prot_loc[2]),
                            tolerance="1.5",
                        )
            elif interaction in ["hydrophobic", "PI_saltbridge", "NI_saltbridge"]:
                pharm = self.generate_pharmacophore_centers(
                    self.df_all.filter(regex=interaction).columns
                )
                for feature_name, position in pharm.items():
                    feature_id_counter += 1
                    point = ET.SubElement(
                        pharmacophore,
                        "point",
                        name=feature_type,
                        featureId=feature_name,
                        optional="false",
                        disabled="false",
                        weight="1.0",
                        coreCompound=self.ligand_name,
                        id=f"feature{str(feature_id_counter)}",
                    )
                    location = ET.SubElement(
                        point,
                        "position",
                        x3=str(position[0]),
                        y3=str(position[1]),
                        z3=str(position[2]),
                        tolerance="1.5",
                    )
            elif interaction == "pistacking":
                pharm = self.generate_pharmacophore_vectors(
                    self.df_all.filter(regex=interaction).columns
                )
                feature_id_counter += 1
                lig_loc = position[0]
                prot_loc = position[1]

                vector = np.array(lig_loc) - np.array(prot_loc)
                normal_vector = vector / np.linalg.norm(vector)
                x, y, z = normal_vector

                plane = ET.SubElement(
                    pharmacophore,
                    "plane",
                    name=feature_type,
                    featureId=interaction,
                    optional="false",
                    disabled="false",
                    weight="1.0",
                    coreCompound=self.ligand_name,
                    id=f"feature{str(feature_id_counter)}",
                )
                position = ET.SubElement(
                    plane,
                    "position",
                    x3=str(lig_loc[0]),
                    y3=str(lig_loc[1]),
                    z3=str(lig_loc[2]),
                    tolerance="0.9",
                )
                normal = ET.SubElement(
                    plane,
                    "normal",
                    x3=str(x),
                    y3=str(y),
                    z3=str(z),
                    tolerance="0.43633232",
                )

        tree = ET.ElementTree(root)
        tree.write(f"{output_filename}.pml", encoding="UTF-8", xml_declaration=True)

    def generate_bindingmode_pharmacophore(
        self,
        dict_bindingmode: Dict[str, Dict[str, List[List[float]]]],
        outname: str,
        id_num: int = 0,
    ) -> None:
        """Generates pharmacophore from a binding mode and writes it to a .pml file.

        Args:
            dict_bindingmode (Dict[str, Dict[str, List[List[float]]]]): Dictionary containing all interactions of the binding mode and their corresponding ligand and protein coordinates.
            outname (str): Name of the output .pml file.
            id_num (int, optional): ID number for enumerating different binding modes. Defaults to 0.
        """
        feature_types = {
            "Acceptor_hbond": "HBA",
            "Donor_hbond": "HBD",
            "pistacking": "AR",
            "hydrophobic": "H",
            "PI_saltbridge": "PI",
            "NI_saltbridge": "NI",
        }
        feature_id_counter = 0
        root = ET.Element(
            "MolecularEnvironment",
            version="0.0",
            id=f"OpennMMDL_Analysis{id_num}",
            name=self.complex_name,
        )
        pharmacophore = ET.SubElement(
            root,
            "pharmacophore",
            name=self.complex_name,
            id=f"pharmacophore{id_num}",
            pharmacophoreType="LIGAND_SCOUT",
        )

        for interaction in dict_bindingmode.keys():
            # get feature type
            for interactiontype in feature_types.keys():
                if interactiontype in interaction:
                    feature_type = feature_types[interactiontype]
                    break

            # Generate vector features
            if feature_type in ["HBA", "HBD"]:
                if feature_type == "HBA":
                    orig_loc = dict_bindingmode[interaction]["PROTCOO"][0]
                    targ_loc = dict_bindingmode[interaction]["LIGCOO"][0]
                elif feature_type == "HBD":
                    orig_loc = dict_bindingmode[interaction]["LIGCOO"][0]
                    targ_loc = dict_bindingmode[interaction]["PROTCOO"][0]
                feature_id_counter += 1
                points_to_lig = "true" if feature_type == "HBA" else "false"
                hasSyntheticProjectedPoint = "false"
                vector = ET.SubElement(
                    pharmacophore,
                    "vector",
                    name=feature_type,
                    featureId=interaction,
                    pointsToLigand=points_to_lig,
                    hasSyntheticProjectedPoint=hasSyntheticProjectedPoint,
                    optional="false",
                    disabled="false",
                    weight="1.0",
                    coreCompound=self.complex_name,
                    id=f"feature{str(feature_id_counter)}",
                )
                origin = ET.SubElement(
                    vector,
                    "origin",
                    x3=str(orig_loc[0]),
                    y3=str(orig_loc[1]),
                    z3=str(orig_loc[2]),
                    tolerance="1.9499999",
                )
                target = ET.SubElement(
                    vector,
                    "target",
                    x3=str(targ_loc[0]),
                    y3=str(targ_loc[1]),
                    z3=str(targ_loc[2]),
                    tolerance="1.5",
                )
            # Generate point features
            elif feature_type in ["H", "PI", "NI"]:
                position = dict_bindingmode[interaction]["LIGCOO"][0]
                feature_id_counter += 1
                point = ET.SubElement(
                    pharmacophore,
                    "point",
                    name=feature_type,
                    featureId=interaction,
                    optional="false",
                    disabled="false",
                    weight="1.0",
                    coreCompound=self.complex_name,
                    id=f"feature{str(feature_id_counter)}",
                )
                position = ET.SubElement(
                    point,
                    "position",
                    x3=str(position[0]),
                    y3=str(position[1]),
                    z3=str(position[2]),
                    tolerance="1.5",
                )
            # Generate plane features
            elif feature_type == "AR":
                feature_id_counter += 1
                lig_loc = dict_bindingmode[interaction]["LIGCOO"][0]
                prot_loc = dict_bindingmode[interaction]["PROTCOO"][0]

                # Normalize vector of plane
                vector = np.array(lig_loc) - np.array(prot_loc)
                normal_vector = vector / np.linalg.norm(vector)
                x, y, z = normal_vector

                plane = ET.SubElement(
                    pharmacophore,
                    "plane",
                    name=feature_type,
                    featureId=interaction,
                    optional="false",
                    disabled="false",
                    weight="1.0",
                    coreCompound=self.complex_name,
                    id=f"feature{str(feature_id_counter)}",
                )
                position = ET.SubElement(
                    plane,
                    "position",
                    x3=str(lig_loc[0]),
                    y3=str(lig_loc[1]),
                    z3=str(lig_loc[2]),
                    tolerance="0.9",
                )
                normal = ET.SubElement(
                    plane,
                    "normal",
                    x3=str(x),
                    y3=str(y),
                    z3=str(z),
                    tolerance="0.43633232",
                )

        tree = ET.ElementTree(root)
        tree.write(
            f"./Binding_Modes_Markov_States/{outname}.pml",
            encoding="UTF-8",
            xml_declaration=True,
        )

    def generate_pharmacophore_centers_all_points(
        self, interactions: List[str]
    ) -> Dict[str, List[List[float]]]:
        """Generates pharmacophore points for all interactions to generate point cloud.

        Args:
            interactions (List[str]): List of interactions to generate pharmacophore from.

        Returns:
            Dict[str, List[List[float]]]: Dict of interactions with pharmacophore points as values.
        """
        coord_pattern = re.compile(r"\(([\d.-]+), ([\d.-]+), ([\d.-]+)\)")
        pharmacophore = {}
        for interaction in interactions:
            pharmacophore_points = []
            for index, row in self.df_all.iterrows():
                if row[interaction] == 1:
                    coord_match = coord_pattern.match(row["LIGCOO"])
                    if coord_match:
                        x, y, z = map(float, coord_match.groups())
                        pharmacophore_points.append([x, y, z])

            if pharmacophore_points:
                pharmacophore[interaction] = pharmacophore_points
        return pharmacophore

    def generate_point_cloud_pml(self, outname: str) -> None:
        """Generates pharmacophore point cloud and writes it to a .pml file.

        Args:
            outname (str): Name of the output .pml file.
        """
        cloud_dict = {}
        cloud_dict["H"] = self.generate_pharmacophore_centers_all_points(
            self.df_all.filter(regex="hydrophobic").columns
        )
        cloud_dict["HBA"] = self.generate_pharmacophore_centers_all_points(
            self.df_all.filter(regex="Acceptor_hbond").columns
        )
        cloud_dict["HBD"] = self.generate_pharmacophore_centers_all_points(
            self.df_all.filter(regex="Donor_hbond").columns
        )
        cloud_dict["AR"] = self.generate_pharmacophore_centers_all_points(
            self.df_all.filter(regex="pistacking").columns
        )
        cloud_dict["PI"] = self.generate_pharmacophore_centers_all_points(
            self.df_all.filter(regex="PI_saltbridge").columns
        )
        cloud_dict["NI"] = self.generate_pharmacophore_centers_all_points(
            self.df_all.filter(regex="NI_saltbridge").columns
        )
        cloud_dict["M"] = self.generate_pharmacophore_centers_all_points(
            self.df_all.filter(regex="metal").columns
        )

        pharmacophore = ET.Element(
            "pharmacophore",
            name=f"{self.complex_name}_pointcloud",
            id=f"pharmacophore0",
            pharmacophoreType="LIGAND_SCOUT",
        )
        feature_id_counter = 0
        for feature_type in cloud_dict.keys():
            for interaction in cloud_dict[feature_type].keys():
                if len(cloud_dict[feature_type][interaction]) > 1:
                    feature_cloud = ET.SubElement(
                        pharmacophore,
                        "featureCloud",
                        name=feature_type,
                        featureId=interaction,
                        optional="false",
                        disabled="false",
                        weight="1.0",
                        id=f"feature{str(feature_id_counter)}",
                    )
                    position = ET.SubElement(
                        feature_cloud,
                        "position",
                        x3=str(cloud_dict[feature_type][interaction][0][0]),
                        y3=str(cloud_dict[feature_type][interaction][0][1]),
                        z3=str(cloud_dict[feature_type][interaction][0][2]),
                    )
                    for additional_point in cloud_dict[feature_type][interaction][1:]:
                        additional_point = ET.SubElement(
                            feature_cloud,
                            "additionalPoint",
                            x3=str(round(additional_point[0], 2)),
                            y3=str(round(additional_point[1], 2)),
                            z3=str(round(additional_point[2], 2)),
                            weight="1.0",
                        )
                    feature_id_counter += 1

        tree = ET.ElementTree(pharmacophore)
        tree.write(f"{outname}.pml", encoding="UTF-8", xml_declaration=True)


# TODO if json cloudgen breaks here is original function
# def cloud_json_generation(df_all):
#     """generates dict for visualization of interaction clouds. Later saved as .json file.

#     Args:
#         df_all (pandas dataframe): dataframe containing all interactions from plip analysis (typicaly df_all)

#     Returns:
#         dict: dict containing all interaction clouds
#     """
#     coord_pattern = re.compile(r"\(([\d.-]+), ([\d.-]+), ([\d.-]+)\)")
#     hydrophobe_coords = []
#     acceptor_ccords = []
#     donor_coords = []
#     waterbridge_coords = []
#     negative_ionizable_coords = []
#     positive_ionizable_coords = []
#     pistacking_coords = []
#     pication_coords = []
#     halogen_coords = []
#     metal_coords = []

#     for index, row in df_all.iterrows():
#         if row["LIGCOO"] != 0:
#             coord_match = coord_pattern.match(row["LIGCOO"])
#             if coord_match:
#                 x, y, z = map(float, coord_match.groups())
#                 x, y, z = round(x, 3), round(y, 3), round(z, 3)
#                 interaction = row["INTERACTION"]
#                 if interaction == "hbond":
#                     if row["PROTISDON"] == "False":
#                         interaction = "donor"
#                     else:
#                         interaction = "acceptor"
#                 if interaction == "saltbridge":
#                     if row["PROTISPOS"] == "True":
#                         interaction = "negative_ionizable"
#                     else:
#                         interaction = "positive_ionizable"
#                 if interaction == "hydrophobic":
#                     hydrophobe_coords.append([x, y, z])
#                 if interaction == "acceptor":
#                     acceptor_ccords.append([x, y, z])
#                 if interaction == "donor":
#                     donor_coords.append([x, y, z])
#                 if interaction == "waterbridge":
#                     waterbridge_coords.append([x, y, z])
#                 if interaction == "negative_ionizable":
#                     negative_ionizable_coords.append([x, y, z])
#                 if interaction == "positive_ionizable":
#                     positive_ionizable_coords.append([x, y, z])
#                 if interaction == "pistacking":
#                     pistacking_coords.append([x, y, z])
#                 if interaction == "pication":
#                     pication_coords.append([x, y, z])
#                 if interaction == "halogen":
#                     halogen_coords.append([x, y, z])

#     for index, row in df_all.iterrows():
#         if row["TARGETCOO"] != 0:
#             coord_match = coord_pattern.match(row["TARGETCOO"])
#             if coord_match:
#                 x, y, z = map(float, coord_match.groups())
#                 x, y, z = round(x, 3), round(y, 3), round(z, 3)
#                 interaction = row["INTERACTION"]
#                 if interaction == "metal":
#                     metal_coords.append([x, y, z])

#     clouds = {}
#     clouds["hydrophobic"] = {
#         "coordinates": hydrophobe_coords,
#         "color": [1.0, 1.0, 0.0],
#         "radius": 0.1,
#     }
#     clouds["acceptor"] = {
#         "coordinates": acceptor_ccords,
#         "color": [1.0, 0.0, 0.0],
#         "radius": 0.1,
#     }
#     clouds["donor"] = {
#         "coordinates": donor_coords,
#         "color": [0.0, 1.0, 0.0],
#         "radius": 0.1,
#     }
#     clouds["waterbridge"] = {
#         "coordinates": waterbridge_coords,
#         "color": [0.0, 1.0, 0.9],
#         "radius": 0.1,
#     }
#     clouds["negative_ionizable"] = {
#         "coordinates": negative_ionizable_coords,
#         "color": [0.0, 0.0, 1.0],
#         "radius": 0.1,
#     }
#     clouds["positive_ionizable"] = {
#         "coordinates": positive_ionizable_coords,
#         "color": [1.0, 0.0, 0.0],
#         "radius": 0.1,
#     }
#     clouds["pistacking"] = {
#         "coordinates": pistacking_coords,
#         "color": [0.0, 0.0, 1.0],
#         "radius": 0.1,
#     }
#     clouds["pication"] = {
#         "coordinates": pication_coords,
#         "color": [0.0, 0.0, 1.0],
#         "radius": 0.1,
#     }
#     clouds["halogen"] = {
#         "coordinates": halogen_coords,
#         "color": [1.0, 0.0, 0.9],
#         "radius": 0.1,
#     }
#     clouds["metal"] = {
#         "coordinates": metal_coords,
#         "color": [1.0, 0.6, 0.0],
#         "radius": 0.1,
#     }

#     return clouds

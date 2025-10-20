from openmmdl.openmmdl_analysis.core.utils import (
    combine_subdict_values,
    remove_duplicate_values,
)


class BindingModeProcesser:
    """
    A class that processes protein-ligand interaction data for a given molecular dynamics (MD) simulation trajectory.
    The class performs multiple analyses on interaction data, such as filtering interactions based on frequency,
    generating fingerprints according to the combination of interactions present in frames creating specific binding modes.

    Attributes
    ----------
    pdb_md : str
        The path to the molecular dynamics (MD) PDB file for the protein-ligand complex.
    ligand : str
        The ligand information (e.g., name or structure).
    peptide : str or None
        The peptide chain ID, if applicable, or None if not considered.
    special : str
        Special handling or unique identifiers related to the ligand or protein.
    ligand_rings : list
        List of ligand ring information used for hydrophobic interaction detection.
    interaction_list : pd.DataFrame
        DataFrame containing the interactions for the trajectory frames.
    threshold : float
        The threshold for filtering interactions based on their occurrence frequency.
    total_frames : int
        The total number of frames in the molecular dynamics simulation trajectory.
    unique_columns_rings_grouped : dict
        A dictionary containing the grouped interactions based on frames.
    interactions_all : pd.DataFrame
        DataFrame containing all interaction data processed.
    unique_data_all : dict
        Dictionary of unique data generated for all interactions across all frames.
    unique_data : dict
        Dictionary of unique data generated based on filtered interactions and threshold.
    """

    def __init__(
        self,
        pdb_md,
        ligand,
        peptide,
        special,
        ligand_rings,
        interaction_list,
        threshold,
        total_frames,
    ):
        self.pdb_md = pdb_md
        self.ligand = ligand
        self.peptide = peptide
        self.special = special
        self.threshold = threshold
        self.ligand_rings = ligand_rings
        self.total_frames = total_frames
        self.unique_columns_rings_grouped = self._gather_interactions(interaction_list)
        self.interaction_list, self.unique_data = self._process_interaction_wraper(interaction_list, (threshold / 100))
        self.interactions_all, self.unique_data_all = self._process_interaction_wraper(interaction_list.copy(), 0.00001)

    def _gather_interactions(self, df):
        """
        Process a DataFrame with the protein-ligand interaction and generate column names for each unique interaction.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame that contains the interaction data for the whole trajectory.

        Returns
        -------
        dict
            A dictionary with the keys being 'FRAME' numbers and values being dictionaries containing
            row indices and their corresponding unique column names for interactions.
        """
        unique_columns_rings = {}
        unique_columns_rings_grouped = {}

        # Iterate through the rows of the DataFrame
        if self.peptide is None:
            for index, row in df.iterrows():
                # Check if the 'INTERACTION' is 'hydrophobic'
                if row["INTERACTION"] == "hydrophobic":
                    # Get the values from the current row
                    prot_partner = row["Prot_partner"]
                    ligcarbonidx = row["LIGCARBONIDX"]
                    interaction = row["INTERACTION"]
                    ring_found = False
                    # Concatenate the values to form a unique column name
                    for ligand_ring in self.ligand_rings:
                        if ligcarbonidx in ligand_ring:
                            numbers_as_strings = [str(ligcarbonidx) for ligcarbonidx in ligand_ring]
                            # Create the name with numbers separated by underscores
                            name_with_numbers = "_".join(numbers_as_strings)
                            col_name = f"{prot_partner}_{name_with_numbers}_{interaction}"
                            ring_found = True
                            break
                    if not ring_found:
                        ligcarbonidx = int(row["LIGCARBONIDX"])
                        col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
                elif row["INTERACTION"] == "hbond":
                    if row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        ligcarbonidx = int(row["ACCEPTORIDX"])
                        interaction = row["INTERACTION"]
                        type = "Acceptor"
                    elif not row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        ligcarbonidx = int(row["DONORIDX"])
                        interaction = row["INTERACTION"]
                        type = "Donor"
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
                elif row["INTERACTION"] == "halogen":
                    prot_partner = row["Prot_partner"]
                    ligcarbonidx = int(row["DON_IDX"])
                    halogen = row["DONORTYPE"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{ligcarbonidx}_{halogen}_{interaction}"
                elif row["INTERACTION"] == "waterbridge":
                    if row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        ligcarbonidx = int(row["ACCEPTOR_IDX"])
                        interaction = row["INTERACTION"]
                        type = "Acceptor"
                        # Concatenate the values to form a unique column name
                        col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
                    elif not row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        ligcarbonidx = int(row["DONOR_IDX"])
                        interaction = row["INTERACTION"]
                        type = "Donor"
                        # Concatenate the values to form a unique column name
                        col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
                elif row["INTERACTION"] == "pistacking":
                    prot_partner = row["Prot_partner"]
                    ligcarbonidx = row["LIG_IDX_LIST"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
                elif row["INTERACTION"] == "pication":
                    prot_partner = row["Prot_partner"]
                    ligidx = row["LIG_IDX_LIST"]
                    ligtype = row["LIG_GROUP"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
                    col_name = col_name.replace(",", "_")
                elif row["INTERACTION"] == "saltbridge":
                    prot_partner = row["Prot_partner"]
                    ligidx = row["LIG_IDX_LIST"]
                    lig_group = row["LIG_GROUP"]
                    interaction = row["INTERACTION"]
                    if row["PROTISPOS"]:
                        type = "NI"
                        # Concatenate the values to form a unique column name
                        col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                    elif not row["PROTISPOS"]:
                        type = "PI"
                        col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                elif row["INTERACTION"] == "metal":
                    special_ligand = row["RESTYPE_LIG"]
                    ligcarbonidx = int(row["TARGET_IDX"])
                    metal_type = row["METAL_TYPE"]
                    coordination = row["COORDINATION"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{special_ligand}_{ligcarbonidx}_{metal_type}_{coordination}_{interaction}"
                frame_value = row["FRAME"]
                if frame_value not in unique_columns_rings_grouped:
                    unique_columns_rings_grouped[frame_value] = {}
                if row["INTERACTION"] != "skip":
                    unique_columns_rings_grouped[frame_value][index] = col_name
                    # Add the column name and its value to the dictionary
                    unique_columns_rings[index] = col_name
        if self.peptide is not None:
            for index, row in df.iterrows():
                # Check if the 'INTERACTION' is 'hydrophobic'
                if row["INTERACTION"] == "hydrophobic":
                    # Get the values from the current row
                    prot_partner = row["Prot_partner"]
                    peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    ring_found = False
                    col_name = f"{prot_partner}_{peptide_partner}_{interaction}"
                elif row["INTERACTION"] == "hbond":
                    if row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                        interaction = row["INTERACTION"]
                        type = "Acceptor"
                    elif not row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                        interaction = row["INTERACTION"]
                        type = "Donor"
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{peptide_partner}_{type}_{interaction}"
                elif row["INTERACTION"] == "halogen":
                    prot_partner = row["Prot_partner"]
                    peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    halogen = row["DONORTYPE"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{peptide_partner}_{halogen}_{interaction}"
                elif row["INTERACTION"] == "waterbridge":
                    if row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                        interaction = row["INTERACTION"]
                        type = "Acceptor"
                        # Concatenate the values to form a unique column name
                        col_name = f"{prot_partner}_{peptide_partner}_{type}_{interaction}"
                    elif not row["PROTISDON"]:
                        prot_partner = row["Prot_partner"]
                        peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                        interaction = row["INTERACTION"]
                        type = "Donor"
                        # Concatenate the values to form a unique column name
                        col_name = f"{prot_partner}_{peptide_partner}_{type}_{interaction}"
                elif row["INTERACTION"] == "pistacking":
                    prot_partner = row["Prot_partner"]
                    peptide_partner = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{peptide_partner}_{interaction}"
                elif row["INTERACTION"] == "pication":
                    prot_partner = row["Prot_partner"]
                    ligidx = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    ligtype = row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
                    col_name = col_name.replace(",", "_")
                elif row["INTERACTION"] == "saltbridge":
                    prot_partner = row["Prot_partner"]
                    ligidx = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    lig_group = row["RESTYPE_LIG"]
                    interaction = row["INTERACTION"]
                    if row["PROTISPOS"]:
                        type = "NI"
                        # Concatenate the values to form a unique column name
                        col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                    elif not row["PROTISPOS"]:
                        type = "PI"
                        col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
                elif row["INTERACTION"] == "metal":
                    special_ligand = row["RESTYPE_LIG"]
                    ligcarbonidx = str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]
                    metal_type = row["METAL_TYPE"]
                    coordination = row["COORDINATION"]
                    interaction = row["INTERACTION"]
                    # Concatenate the values to form a unique column name
                    col_name = f"{special_ligand}_{ligcarbonidx}_{metal_type}_{coordination}_{interaction}"
                frame_value = row["FRAME"]
                if frame_value not in unique_columns_rings_grouped:
                    unique_columns_rings_grouped[frame_value] = {}
                if row["INTERACTION"] != "skip":
                    unique_columns_rings_grouped[frame_value][index] = col_name
                    # Add the column name and its value to the dictionary
                    unique_columns_rings[index] = col_name
        print("\033[1minteraction partners generated\033[0m")

        return unique_columns_rings_grouped

    def _process_interaction_wraper(self, interaction_list, threshold):
        """
        Apply filtering and interaction enumeration to an interaction DataFrame.

        Parameters
        ----------
        interaction_list : pd.DataFrame
            Interaction data obtained from the MD simulation.
        threshold : float
            Threshold for interaction occurrence (as a fraction of total frames).

        Returns
        -------
        interaction_list : pd.DataFrame
            Modified DataFrame including new interaction columns that contain the filtered values.
        unique_data : dict
            Dictionary containing unique filtered interaction names.
        """
        filtered_values = self._filtering_values(threshold, interaction_list)
        interaction_list.fillna(0, inplace=True)
        unique_data = self._unique_data_generation(filtered_values)
        self._df_iteration_numbering(interaction_list, unique_data)

        return interaction_list, unique_data

    def _filtering_values(self, threshold, df):
        """
        Filter and append values (interactions) to a DataFrame based on occurrence counts.

        Parameters
        ----------
        threshold : float
            A threshold value that is used for filtering of the values (interactions) based upon the occurence count.
        df : pd.DataFrame
            DataFrame to which the filtered values (interactions) will be added.

        Returns
        -------
        list of str
            A list of values, with unique values and their corresponding occurence counts.
        """
        # Call the function to remove duplicate keys
        unique_data = remove_duplicate_values(self.unique_columns_rings_grouped)

        # Call the function to combine sub-dictionary values
        unique_colums_rings_all = combine_subdict_values(unique_data)

        # Flatten the list of values
        all_values = unique_colums_rings_all["all"]

        # Count the occurrences of each value
        occurrences = {}
        for value in all_values:
            if value in occurrences:
                occurrences[value] += 1
            else:
                occurrences[value] = 1

        # Calculate the threshold (20% of 1000)
        threshold = threshold * self.total_frames

        # Filter out values that appear in less than 20% of 1000
        filtered_values = [value for value, count in occurrences.items() if count >= threshold]

        # Append the filtered values as columns to the DataFrame
        for value in filtered_values:
            df[value] = None

        return filtered_values

    def _unique_data_generation(self, filtered_values):
        """
        Generate a dictionary conataing the unique interactions from a list of filtered values obtained by filtering_values.


        Parameters
        ----------
        filtered_values : list of str
            A list of values, where the unique interactions are extracted from.

        Returns
        -------
        dict
            A dictionary containing the filtered unique interactions.
        """
        # Create a new dictionary to store unique values
        unique_data = {}

        # Loop through the elements of the data_list
        for element in filtered_values:
            # Check if the element is not already in the unique_data dictionary keys
            if element not in unique_data:
                # If not, add the element as both key and value
                unique_data[element] = element

        return unique_data

    def _df_iteration_numbering(self, df, unique_data):
        """
        Loop through the DataFrame and assign the values 1 and 0 to the rows, depending if the corresponding interaction from unique data is present.

        Parameters
        ----------
        df : pd.DataFrame
            DataFrame which has the interaction data for all of the frames.
        unique_data : dict
            Dictionary that contains the unique interactions obtained from unique_data_generation.

        Returns
        -------
        None
            Modifies the input DataFrame in-place by appending columns corresponding to recurring interactions.
        """
        if self.peptide is None:
            for index, row in df.iterrows():
                if row["INTERACTION"] == "hydrophobic":
                    for col in unique_data.values():
                        if "hydrophobic" in col:
                            ligcarbonidx_check = str(int(row["LIGCARBONIDX"]))
                            if ligcarbonidx_check in col:
                                parts = col.split("_")
                                prot_partner = parts[0]
                                interaction = parts[-1]
                                condition = (row["Prot_partner"] == prot_partner) & (row["INTERACTION"] == interaction)
                                df.at[index, col] = 1 if condition else 0
                            else:
                                continue
                elif row["INTERACTION"] == "hbond":
                    if row["PROTISDON"]:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                ligcarbonidx = int(ligcarbonidx)
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["ACCEPTORIDX"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                    elif not row["PROTISDON"]:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                ligcarbonidx = int(ligcarbonidx)
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["DONORIDX"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "halogen":
                    for col in unique_data.values():
                        if "halogen" in col:
                            prot_partner, ligcarbonidx, halogen, interaction = col.split("_")
                            ligcarbonidx = int(ligcarbonidx)
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (int(row["DON_IDX"]) == ligcarbonidx)
                                & (row["DONORTYPE"] == halogen)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pistacking":
                    for col in unique_data.values():
                        if "pistacking" in col:
                            prot_partner, ligcarbonidx, interaction = col.split("_")
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["LIG_IDX_LIST"] == ligcarbonidx)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "waterbridge":
                    for col in unique_data.values():
                        if "waterbridge" in col:
                            if row["PROTISDON"]:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["ACCEPTOR_IDX"]) == int(ligcarbonidx))
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                            elif not row["PROTISDON"]:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & (int(row["DONOR_IDX"]) == int(ligcarbonidx))
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pication":
                    for col in unique_data.values():
                        if "pication" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1:-2]
                            ligidx = ",".join(ligidx)
                            ligtype = parts[-2]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["LIG_IDX_LIST"] == ligidx)
                                & (row["LIG_GROUP"] == ligtype)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "saltbridge":
                    for col in unique_data.values():
                        if "saltbridge" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1:-3]
                            ligidx = ",".join(ligidx)
                            lig_group = parts[-3]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & (row["LIG_IDX_LIST"] == ligidx)
                                & (row["LIG_GROUP"] == lig_group)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "metal":
                    for col in unique_data.values():
                        if "metal" in col:
                            parts = col.split("_")
                            special_ligand = parts[0]
                            ligidx = int(parts[1])
                            metal_type = parts[2]
                            interaction = parts[-1]
                            condition = (
                                (row["RESTYPE_LIG"] == special_ligand)
                                & (int(row["TARGET_IDX"]) == ligidx)
                                & (row["METAL_TYPE"] == metal_type)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
        if self.peptide is not None:
            for index, row in df.iterrows():
                if row["INTERACTION"] == "hydrophobic":
                    for col in unique_data.values():
                        if "hydrophobic" in col:
                            ligcarbonidx_check = str(int(row["RESNR_LIG"]))
                            if ligcarbonidx_check in col:
                                parts = col.split("_")
                                prot_partner = parts[0]
                                interaction = parts[-1]
                                condition = (row["Prot_partner"] == prot_partner) & (row["INTERACTION"] == interaction)
                                df.at[index, col] = 1 if condition else 0
                            else:
                                continue
                elif row["INTERACTION"] == "hbond":
                    if row["PROTISDON"]:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                    elif not row["PROTISDON"]:
                        for col in unique_data.values():
                            if "hbond" in col:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "halogen":
                    for col in unique_data.values():
                        if "halogen" in col:
                            prot_partner, ligcarbonidx, halogen, interaction = col.split("_")
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligcarbonidx)
                                & (row["DONORTYPE"] == halogen)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pistacking":
                    for col in unique_data.values():
                        if "pistacking" in col:
                            prot_partner, ligcarbonidx, interaction = col.split("_")
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligcarbonidx)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "waterbridge":
                    for col in unique_data.values():
                        if "waterbridge" in col:
                            if row["PROTISDON"]:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                            elif not row["PROTISDON"]:
                                prot_partner, ligcarbonidx, type, interaction = col.split("_")
                                condition = (
                                    (row["Prot_partner"] == prot_partner)
                                    & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligcarbonidx)
                                    & (row["INTERACTION"] == interaction)
                                )
                                df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "pication":
                    for col in unique_data.values():
                        if "pication" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1]
                            ligtype = parts[-2]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligidx)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "saltbridge":
                    for col in unique_data.values():
                        if "saltbridge" in col:
                            parts = col.split("_")
                            prot_partner = parts[0]
                            ligidx = parts[1]
                            lig_group = parts[-3]
                            interaction = parts[-1]
                            condition = (
                                (row["Prot_partner"] == prot_partner)
                                & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligidx)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0
                elif row["INTERACTION"] == "metal":
                    for col in unique_data.values():
                        if "metal" in col:
                            parts = col.split("_")
                            special_ligand = parts[0]
                            ligidx = parts[1]
                            metal_type = parts[2]
                            interaction = parts[-1]
                            condition = (
                                (row["RESTYPE_LIG"] == special_ligand)
                                & ((str(row["RESNR_LIG"]) + row["RESTYPE_LIG"]) == ligidx)
                                & (row["METAL_TYPE"] == metal_type)
                                & (row["INTERACTION"] == interaction)
                            )
                            df.at[index, col] = 1 if condition else 0

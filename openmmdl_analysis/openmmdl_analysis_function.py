# Import libraries
from pathlib import Path
import time
import warnings

warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import rdkit
import MDAnalysis as mda
import matplotlib.pyplot as plt
from matplotlib import colors
import networkx as nx
from plip.structure.preparation import PDBComplex, LigandFinder, Mol, PLInteraction
from plip.exchange.report import BindingSiteReport

from collections import Counter
import itertools

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import cairosvg
import pylab
from PIL import Image
import os

def characterize_complex(pdb_file: str, binding_site_id: str) -> PLInteraction:
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    for ligand in pdb_complex.ligands:
        if ':'.join([ligand.hetid, ligand.chain, str(ligand.position)]) == binding_site_id:
            pdb_complex.characterize_complex(ligand)
    return pdb_complex.interaction_sets[binding_site_id]


def retrieve_plip_interactions(pdb_file):
    """
    Retrieves the interactions from PLIP.

    Parameters
    ----------
    pdb_file :
        The PDB file of the complex.

    Returns
    -------
    dict :
        A dictionary of the binding sites and the interactions.
    """
    protlig = PDBComplex()
    protlig.load_pdb(pdb_file)  # load the pdb file
    for ligand in protlig.ligands:
        protlig.characterize_complex(ligand)  # find ligands and analyze interactions
    sites = {}
    # loop over binding sites
    for key, site in sorted(protlig.interaction_sets.items()):
        binding_site = BindingSiteReport(site)  # collect data about interactions
        # tuples of *_features and *_info will be converted to pandas DataFrame
        keys = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )
        # interactions is a dictionary which contains relevant information for each
        # of the possible interactions: hydrophobic, hbond, etc. in the considered
        # binding site. Each interaction contains a list with
        # 1. the features of that interaction, e.g. for hydrophobic:
        # ('RESNR', 'RESTYPE', ..., 'LIGCOO', 'PROTCOO')
        # 2. information for each of these features, e.g. for hydrophobic
        # (residue nb, residue type,..., ligand atom 3D coord., protein atom 3D coord.)
        interactions = {
            k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info")
            for k in keys
        }
        sites[key] = interactions
    return sites

def create_df_from_binding_site(selected_site_interactions, interaction_type="hbond"):
    """
    Creates a data frame from a binding site and interaction type.

    Parameters
    ----------
    selected_site_interactions : dict
        Precaluclated interactions from PLIP for the selected site
    interaction_type : str
        The interaction type of interest (default set to hydrogen bond).

    Returns
    -------
    pd.DataFrame :
        DataFrame with information retrieved from PLIP.
    """

    # check if interaction type is valid:
    valid_types = [
        "hydrophobic",
        "hbond",
        "waterbridge",
        "saltbridge",
        "pistacking",
        "pication",
        "halogen",
        "metal",
    ]

    if interaction_type not in valid_types:
        print("!!! Wrong interaction type specified. Hbond is chosen by default!!!\n")
        interaction_type = "hbond"

    df = pd.DataFrame.from_records(
        # data is stored AFTER the column names
        selected_site_interactions[interaction_type][1:],
        # column names are always the first element
        columns=selected_site_interactions[interaction_type][0],
    )
    return df

def increase_ring_indices(ring, lig_index):
    return [atom_idx + lig_index for atom_idx in ring]

def combine_subdict_values(data):
    combined_data = {'all': []}
    for sub_dict in data.values():
        combined_data['all'].extend(sub_dict.values())
    return combined_data

def remove_duplicate_values(data):
    unique_data = {}

    for key, sub_dict in data.items():
        unique_sub_dict = {}
        seen_values = set()

        for sub_key, value in sub_dict.items():
            if value not in seen_values:
                unique_sub_dict[sub_key] = value
                seen_values.add(value)

        if unique_sub_dict:
            unique_data[key] = unique_sub_dict

    return unique_data


def gather_interactions(df, ligand_rings):
    # Create a dictionary to store the unique column names and their values
    unique_columns_rings = {}
    unique_columns_rings_grouped = {}

    # Iterate through the rows of the DataFrame
    for index, row in df.iterrows():
        # Check if the 'INTERACTION' is 'hydrophobic'
        if row['INTERACTION'] == 'hydrophobic':
            # Get the values from the current row
            prot_partner = row['Prot_partner']
            ligcarbonidx = row['LIGCARBONIDX']
            interaction = row['INTERACTION']
            ring_found = False
            # Concatenate the values to form a unique column name
            for ligand_ring in ligand_rings:
                if ligcarbonidx in ligand_ring:
                    numbers_as_strings = [str(ligcarbonidx) for ligcarbonidx in ligand_ring]
                    # Create the name with numbers separated by underscores
                    name_with_numbers = '_'.join(numbers_as_strings)
                    col_name = f"{prot_partner}_{name_with_numbers}_{interaction}"
                    ring_found = True
                    break
            if not ring_found:
                ligcarbonidx = int(row['LIGCARBONIDX'])
                col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
        elif row['INTERACTION'] == 'hbond':
            if row['PROTISDON'] == True:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['ACCEPTORIDX'])
                interaction = row['INTERACTION']
                type = "Acceptor"
            elif row['PROTISDON'] == False:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['DONORIDX'])
                interaction = row['INTERACTION']
                type = "Donor"
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
        elif row['INTERACTION'] == 'halogen':
            prot_partner = row['Prot_partner']
            ligcarbonidx = int(row['DON_IDX'])
            halogen = row['DONORTYPE']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{halogen}_{interaction}"
        elif row['INTERACTION'] == 'waterbridge':
            if row['PROTISDON'] == True:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['ACCEPTOR_IDX'])
                interaction = row['INTERACTION']
                type = "Acceptor"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
            elif row['PROTISDON'] == False:
                prot_partner = row['Prot_partner']
                ligcarbonidx = int(row['DONOR_IDX'])
                interaction = row['INTERACTION']
                type = "Donor"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligcarbonidx}_{type}_{interaction}"
        elif row['INTERACTION'] == 'pistacking':
            prot_partner = row['Prot_partner']
            ligcarbonidx = row['LIG_IDX_LIST']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{interaction}"
        elif row['INTERACTION'] == 'pication':
            prot_partner = row['Prot_partner']
            ligidx = row['LIG_IDX_LIST']
            ligtype = row['LIG_GROUP']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligidx}_{ligtype}_{interaction}"
            col_name = col_name.replace(',', '_')
        elif row['INTERACTION'] == 'saltbridge':
            prot_partner = row['Prot_partner']
            ligidx = row['LIG_IDX_LIST']
            lig_group = row['LIG_GROUP']
            interaction = row['INTERACTION']
            if row['PROTISPOS'] == True:
                type = "NI"
                # Concatenate the values to form a unique column name
                col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
            elif row['PROTISPOS'] == False:
                type = "PI"
                col_name = f"{prot_partner}_{ligidx}_{lig_group}_{type}_{interaction}"
        elif row['INTERACTION'] == 'metal':
            prot_partner = row['Prot_partner']
            ligcarbonidx = row['METAL_IDX']
            metal_type = row['METAL_TYPE']
            location = row['LOCATION']
            interaction = row['INTERACTION']
            # Concatenate the values to form a unique column name
            col_name = f"{prot_partner}_{ligcarbonidx}_{metal_type}_{location}_{interaction}"
        frame_value = row['FRAME']
        if frame_value not in unique_columns_rings_grouped:
            unique_columns_rings_grouped[frame_value] = {}
        unique_columns_rings_grouped[frame_value][index] = col_name
            
        # Add the column name and its value to the dictionary
        unique_columns_rings[index] = col_name
    print("interactions gathered")
    return unique_columns_rings_grouped


def filtering_values(threshold, frames, df, unique_columns_rings_grouped):
    # Call the function to remove duplicate keys
    unique_data = remove_duplicate_values(unique_columns_rings_grouped)

    # Call the function to combine sub-dictionary values
    unique_colums_rings_all = combine_subdict_values(unique_data)

    # Flatten the list of values
    all_values = unique_colums_rings_all['all']

    # Count the occurrences of each value
    occurrences = {}
    for value in all_values:
        if value in occurrences:
            occurrences[value] += 1
        else:
            occurrences[value] = 1

    # Calculate the threshold (20% of 1000)
    threshold = threshold * frames

    # Filter out values that appear in less than 20% of 1000
    filtered_values = [value for value, count in occurrences.items() if count >= threshold]

    # Append the filtered values as columns to the DataFrame
    for value in filtered_values:
        df[value] = None
    
    return filtered_values


def unique_data_generation(filtered_values):
    # Create a new dictionary to store unique values
    unique_data = {}

    # Loop through the elements of the data_list
    for element in filtered_values:
        # Check if the element is not already in the unique_data dictionary keys
        if element not in unique_data:
            # If not, add the element as both key and value
            unique_data[element] = element
    return unique_data


def df_iteration_numbering(df,unique_data):
    for index, row in df.iterrows():
        if row['INTERACTION'] == "hydrophobic":
            for col in unique_data.values():
                if "hydrophobic" in col:
                    ligcarbonidx_check = str(int(row['LIGCARBONIDX']))
                    if ligcarbonidx_check in col:
                        parts = col.split('_')
                        prot_partner = parts[0]
                        interaction = parts[-1]
                        condition = (row['Prot_partner'] == prot_partner) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
                    else:
                        continue
        elif row['INTERACTION'] == "hbond":
            if row['PROTISDON'] == True:
                for col in unique_data.values():
                    if "hbond" in col:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        ligcarbonidx = int(ligcarbonidx)
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['ACCEPTORIDX']) == ligcarbonidx) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
            elif row['PROTISDON'] == False:
                for col in unique_data.values():
                    if "hbond" in col:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        ligcarbonidx = int(ligcarbonidx)
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['DONORIDX']) == ligcarbonidx) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "halogen":
            for col in unique_data.values():
                if "halogen" in col:
                    prot_partner, ligcarbonidx, halogen, interaction = col.split('_')
                    ligcarbonidx = int(ligcarbonidx)
                    condition = (row['Prot_partner'] == prot_partner) & (int(row['DON_IDX']) == ligcarbonidx) & (row['DONORTYPE'] == halogen) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "pistacking":
            for col in unique_data.values():
                if "pistacking" in col:
                    prot_partner, ligcarbonidx, interaction = col.split('_')
                    condition = (row['Prot_partner'] == prot_partner) & (row['LIG_IDX_LIST'] == ligcarbonidx) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "waterbridge":
            for col in unique_data.values():
                if "waterbridge" in col:
                    if row['PROTISDON'] == True:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['ACCEPTOR_IDX']) == int(ligcarbonidx)) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
                    elif row['PROTISDON'] == False:
                        prot_partner, ligcarbonidx, type, interaction = col.split('_')
                        condition = (row['Prot_partner'] == prot_partner) & (int(row['DONOR_IDX']) == int(ligcarbonidx)) & (row['INTERACTION'] == interaction)
                        df.at[index, col] = 1 if condition else 0
        elif row['INTERACTION'] == "pication":
            for col in unique_data.values():
                if "pication" in col:
                    parts = col.split('_')
                    prot_partner = parts[0]
                    ligidx = parts[1:-2]
                    ligidx = ','.join(ligidx)
                    ligtype = parts[-2]
                    interaction = parts[-1]
                    condition = (row['Prot_partner'] == prot_partner) & (row['LIG_IDX_LIST'] == ligidx) & (row['LIG_GROUP'] == ligtype) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0

        elif row['INTERACTION'] == "saltbridge":
            for col in unique_data.values():
                if "saltbridge" in col:
                    parts = col.split('_')
                    prot_partner = parts[0]
                    ligidx = parts[1:-3]
                    ligidx = ','.join(ligidx)
                    lig_group = parts[-3]
                    type = parts[-2]
                    interaction = parts[-1]
                    condition = (row['Prot_partner'] == prot_partner) & (row['LIG_IDX_LIST'] == ligidx) & (row['LIG_GROUP'] == lig_group) & (row['INTERACTION'] == interaction)
                    df.at[index, col] = 1 if condition else 0

    print("done")

def update_values(df, new, unique_data):
    for idx, row in df.iterrows():
        frame_value = row['FRAME']
        values_to_update = new.loc[frame_value, list(unique_data.values())]
        df.loc[idx, list(unique_data.values())] = values_to_update


def binding_site_markov_network2(total_frames, min_transition, combined_dict, font_size=None, size_node=None):

    font_size=12 if font_size is None else font_size
    size_node=200 if size_node is None else size_node

    min_prob = min_transition / total_frames  # Set the minimum probability threshold (adjust this value as needed)

    # Calculate the number of elements in each part
    total_length = len(combined_dict['all'])
    part_length = total_length // 3
    remaining_length = total_length % 3

    # Divide the 'all_data' into three parts
    part1_length = part_length + remaining_length
    part2_length = part_length
    part1_data = combined_dict['all'][:part1_length]
    part2_data = combined_dict['all'][part1_length:part1_length + part2_length]
    part3_data = combined_dict['all'][part1_length + part2_length:]

    # Count the occurrences of each node in each part
    part1_node_occurrences = {node: part1_data.count(node) for node in set(part1_data)}
    part2_node_occurrences = {node: part2_data.count(node) for node in set(part2_data)}
    part3_node_occurrences = {node: part3_data.count(node) for node in set(part3_data)}

    # Get the top 10 nodes with the most occurrences
    node_occurrences = {node: combined_dict['all'].count(node) for node in set(combined_dict['all'])}
    top_10_nodes = sorted(node_occurrences, key=node_occurrences.get, reverse=True)[:10]

    # Create a directed graph
    G = nx.DiGraph()

    # Count the occurrences of each transition and self-loop
    transitions = {}
    self_loops = {}
    for i in range(len(combined_dict['all']) - 1):
        current_state = combined_dict['all'][i]
        next_state = combined_dict['all'][i + 1]

        if current_state == next_state:  # Check for self-loop
            self_loop = (current_state, next_state)
            self_loops[self_loop] = self_loops.get(self_loop, 0) + 1
        else:
            transition = (current_state, next_state)
            transitions[transition] = transitions.get(transition, 0) + 1

    # Add edges to the graph with their probabilities
    for transition, count in transitions.items():
        current_state, next_state = transition
        probability = count / len(combined_dict['all']) * 100  # Convert probability to percentage
        if probability >= 100 * min_prob:
            G.add_edge(current_state, next_state, weight=probability)

    # Add self-loops to the graph with their probabilities
    for self_loop, count in self_loops.items():
        state = self_loop[0]
        probability = count / len(combined_dict['all']) * 100  # Convert probability to percentage
        if probability >= 100 * min_prob:
            G.add_edge(state, state, weight=probability)

    # Calculate transition probabilities for each direction (excluding self-loops)
    transition_probabilities_forward = {}
    transition_probabilities_backward = {}

    for transition, count in transitions.items():
        start_state, end_state = transition
        forward_transition = (start_state, end_state)
        backward_transition = (end_state, start_state)

        transition_probabilities_forward[forward_transition] = count / node_occurrences[start_state]
        transition_probabilities_backward[backward_transition] = count / node_occurrences[end_state]

    # Calculate self-loop probabilities
    self_loop_probabilities = {}
    for self_loop, count in self_loops.items():
        state = self_loop[0]
        self_loop_probabilities[state] = count / node_occurrences[state]

    # Generate the Markov Chain plot
    plt.figure(figsize=(30, 30))  # Increased figure size
    plt.title(f"Markov Chain Plot {min_transition} Frames Transition", fontsize=35)

    # Draw nodes and edges
    pos = nx.spring_layout(G, k=2, seed=42)  # Increased distance between nodes (k=2)
    edge_colors = []

    for u, v, data in G.edges(data=True):
        weight = data['weight']

        if u == v:  # Check if it is a self-loop
            edge_colors.append('green')  # Set green color for self-loop arrows
            width = 0.1  # Make self-loop arrows smaller
            connection_style = 'arc3,rad=-0.1'  # Make the self-loops more curved and compact
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=0.2, edge_color=edge_colors[-1],
                               connectionstyle=connection_style)
        elif weight >= 250 * min_prob:
            edge_colors.append('black')  # Highlight significant transitions in red
            width = 4.0
            edge_label = f"{transition_probabilities_forward.get((u, v), 0):.2f}% of Frames →\n{transition_probabilities_backward.get((v, u), 0):.2f}% of Frames ←"
            connection_style = 'arc3,rad=-0.1'
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=0.7, edge_color=edge_colors[-1], connectionstyle=connection_style)
            nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): edge_label}, font_size=15)
        else:
            edge_colors.append('grey')  # Use black for non-significant transitions
            width = 0.5
            edge_label = f"{transition_probabilities_forward.get((u, v), 0):.2f}% of Frames →\n{transition_probabilities_backward.get((v, u), 0):.2f}% of Frames ←"
            connection_style = 'arc3,rad=-0.1' 
            nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=width, alpha=0.7, edge_color=edge_colors[-1], connectionstyle=connection_style)
            nx.draw_networkx_edge_labels(G, pos, edge_labels={(u, v): edge_label}, font_size=12)

    # Update the node colors based on their appearance percentages in each part
    node_colors = []
    for node in G.nodes():
        if node in top_10_nodes:
            part1_percentage = part1_node_occurrences.get(node, 0) / node_occurrences[node]
            part2_percentage = part2_node_occurrences.get(node, 0) / node_occurrences[node]
            part3_percentage = part3_node_occurrences.get(node, 0) / node_occurrences[node]

            if part1_percentage > 0.5:
                node_colors.append('green')
            elif part2_percentage > 0.5:
                node_colors.append('orange')
            elif part3_percentage > 0.5:
                node_colors.append('red')
            else:
                node_colors.append('yellow')
        else:
            node_colors.append('skyblue')

    # Draw nodes with sizes correlated to occurrences and color top 10 nodes
    node_size = [size_node * node_occurrences[node] for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color=node_colors, alpha=0.8)

    # Draw node labels with occurrence percentage and self-loop probability for top 10 nodes
    node_labels = {}
    for node in G.nodes():
        if node in top_10_nodes:
            node_occurrence_percentage = node_occurrences[node] / len(combined_dict['all']) * 100
            self_loop_probability = self_loop_probabilities.get(node, 0) * 100
            node_label = f"{node}\nOccurrences: {node_occurrence_percentage:.2f}%\nSelf-Loop Probability: {self_loop_probability:.2f}%"
        else:
            node_label = node
        node_labels[node] = node_label

    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=font_size, font_color='black', verticalalignment="center")

    plt.axis('off')
    plt.tight_layout()

    # Save the plot as a PNG file
    plot_filename = f"markov_chain_plot_{min_transition}_frames.png"
    plot_path = os.path.join("Binding_Modes_Markov_States", plot_filename)
    os.makedirs("Binding_Modes_Markov_States", exist_ok=True)  # Create the folder if it doesn't exist

    plt.savefig(plot_path, dpi=300)
    plt.show()


def split_interaction_data(data):
    split_data = []
    for item in data:
        parts = item.split('_')
        entity_name = parts[0]
        numeric_codes = " ".join(parts[1:-1])  # Join numeric codes with spaces
        interaction_type = parts[-1]
        split_value = f"{entity_name} {numeric_codes} {interaction_type}"
        split_data.append(split_value)
    return split_data

def highlight_numbers(split_data, starting_idx):
    highlighted_hbond_acceptor = []
    highlighted_hbond_donor = []    
    highlighted_hydrophobic = []
    highlighted_hbond_both = []
    highlighted_waterbridge = []
    highlighted_pistacking = []
    highlighted_halogen = []
    highlighted_ni = []
    highlighted_pi = []
    highlighted_pication = []
    highlighted_metal = []
    
    
    for item in split_data:
        parts = item.split()
        entity_name = parts[0]
        numeric_codes = parts[1:-1]
        interaction_type = parts[-1]

        if interaction_type == 'hbond':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-2]
            type = parts[-2]
            interaction_type = parts[-1]
            for code in numeric_codes:
                atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
                if type == "Donor":
                    highlighted_hbond_donor.append(atom_index)
                elif type == "Acceptor":
                    highlighted_hbond_acceptor.append(atom_index)                    

        elif interaction_type == 'hydrophobic':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-1]
            interaction_type = parts[-1]
            for code in numeric_codes:
                atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
                highlighted_hydrophobic.append(atom_index)

        elif interaction_type == 'waterbridge':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-1]
            interaction_type = parts[-1]
            for code in numeric_codes:
                atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
                highlighted_waterbridge.append(atom_index)

        elif interaction_type == 'pistacking':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-1]
            interaction_type = parts[-1]
            for code in numeric_codes:
                atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
                highlighted_pistacking.append(atom_index)

        elif interaction_type == 'halogen':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-2]
            interaction_type = parts[-1]
            halogen_type = parts[-2]
            for code in numeric_codes:
                atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
                highlighted_halogen.append(atom_index)

        elif interaction_type == 'saltbridge':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-2]
            interaction_type = parts[-1]
            saltbridge_type = parts[-2] 
            if saltbridge_type == "NI":
                highlighted_ni.append(atom_index)
                for code in numeric_codes:
                    atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
            elif saltbridge_type == "PI":
                highlighted_pi.append(atom_index)
                for code in numeric_codes:
                    atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code

        elif interaction_type == 'pication':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-2]
            interaction_type = parts[-1]
            halogen_type = parts[-2]
            for code in numeric_codes:
                atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
                highlighted_pication.append(atom_index)

        elif interaction_type == 'metal':
            parts = item.split()
            entity_name = parts[0]
            numeric_codes = parts[1:-2]
            interaction_type = parts[-1]
            halogen_type = parts[-2]
            for code in numeric_codes:
                atom_index = int(code) - starting_idx  # Subtract starting_idx from the numeric code
                highlighted_metal.append(atom_index)
    
    for value in highlighted_hbond_donor[:]:  # Using a copy of the list to avoid modifying while iterating
        if value in highlighted_hbond_acceptor:
            highlighted_hbond_donor.remove(value)
            highlighted_hbond_acceptor.remove(value)
            highlighted_hbond_both.append(value)
    
    return highlighted_hbond_donor, highlighted_hbond_acceptor, highlighted_hbond_both, highlighted_hydrophobic, highlighted_waterbridge, highlighted_pistacking, highlighted_halogen, highlighted_ni, highlighted_pi, highlighted_pication, highlighted_metal

def barcodegeneration(df, interaction):
    barcode = []
    
    unique_frames = df['FRAME'].unique()
    
    for frame in unique_frames:
        frame_data = df[df['FRAME'] == frame]
        
        if 1 in frame_data[interaction].values:
            barcode.append(1)
            
            
        else:
            barcode.append(0)
    
    return np.array(barcode)

def water_id_finder(df, interaction):
    water_id_list = []
    for index, row in df.iterrows():
        if row[interaction] == 1:
            water_id_list.append(int(row['WATER_IDX']))
    return water_id_list

def waterids_barcode_generator(df, interaction):
    water_id_list = []
    waterid_barcode = []
    for index, row in df.iterrows():
        if row[interaction] == 1:
            water_id_list.append(int(row['WATER_IDX']))
    
    barcode = barcodegeneration(df, interaction)
    
    for value in barcode:
        if value == 1:
            waterid_barcode.append(water_id_list.pop(0))
        else:
            waterid_barcode.append(0)
    return waterid_barcode


def plot_hydrophobic_barcodes(hydrophobicinteraction_barcodes, save_path='Barcodes/hydrophobic_barcodes.png'):
    if not hydrophobicinteraction_barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(hydrophobicinteraction_barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(hydrophobicinteraction_barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

def plot_acceptor_barcodes(acceptor_barcodes, save_path='Barcodes/hbond_acceptor_barcodes.png'):
    if not acceptor_barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(acceptor_barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(acceptor_barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')


def plot_donor_barcodes(donor_barcodes, save_path='Barcodes/hbond_donor_barcodes.png'):
    if not donor_barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(donor_barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(donor_barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

def plot_pistacking_barcodes(pistacking_barcodes, save_path='Barcodes/pistacking_barcodes.png'):
    if not pistacking_barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(pistacking_barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(pistacking_barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

def plot_halogen_barcodes(halogen_barcodes, save_path='Barcodes/halogen_barcodes.png'):
    if not halogen_barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(halogen_barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(halogen_barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

def plot_waterbridge_barcodes(waterbridge_barcodes, save_path='Barcodes/waterbridge_barcodes.png'):
    if not waterbridge_barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(waterbridge_barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(waterbridge_barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')



def plot_waterbridge_piechart(df_all, waterbridge_barcodes, waterbridge_interactions):
	if not waterbridge_barcodes:
	        print("No Piecharts to plot.")
	        return

	os.makedirs('Barcodes/Waterbridge_Piecharts', exist_ok=True)
	plt.figure(figsize=(6, 6))
	for waterbridge_interaction in waterbridge_interactions:
	    plt.clf()
	    waterid_barcode = waterids_barcode_generator(df_all, waterbridge_interaction)
	    waters_count = {}
	    
	    for waterid in waterid_barcode:
	        if waterid != 0:
	            if waterid in waters_count:
	                waters_count[waterid] += 1
	            else:
	                waters_count[waterid] = 1
	    
	    labels = [f'ID {id}' for id in waters_count.keys()]
	    values = waters_count.values()
	    
		    # Combine small categories into "Other" category
	    threshold = 7  # You can adjust this threshold. It is the percentage of the pie chart, not the total number
	    total_second_values = sum(value for _, value in waters_count.items())
	    small_ids = [id for id, value in waters_count.items() if (value / total_second_values) * 100 < threshold]
	
	    if small_ids:
	        small_count = sum(count for id, count in waters_count.items() if id in small_ids)
	        values = [count if id not in small_ids else small_count for id, count in waters_count.items()]
	        labels = [f'ID {id}' if id not in small_ids else '' for id in waters_count.keys()]
	    def autopct_func(pct):
	        total = sum(waters_count.values())
	        value = int(round(pct/100.*total))
	        return f'{pct:.1f}%\n({value})'
	    plt.pie(values, labels=labels, autopct=autopct_func, shadow=False, startangle=140)
	    plt.axis('equal')
	    plt.title(str(waterbridge_interaction), fontweight='bold')
	    # Manually create the legend with the correct labels
	    legend_labels = [f'ID {id}' for id in waters_count.keys()]
	    legend = plt.legend(legend_labels, loc="upper right", bbox_to_anchor=(1.2, 1))
	    plt.setp(legend.get_texts(), fontsize='small')  # Adjust font size for legend
	    plt.text(0.5, 0, f"Total frames with waterbridge: {round(((sum(1 for val in waterid_barcode if val != 0) / len(waterid_barcode)) * 100), 2)}%", size=12, ha="center", 	transform=plt.gcf().transFigure)
	    # Adjust the position of the subplots within the figure
	    plt.subplots_adjust(top=0.99, bottom=0.01)  # You can change the value as needed
	    plt.savefig(f'Barcodes/Waterbridge_Piecharts/{waterbridge_interaction}.png', bbox_inches='tight', dpi=300)


def plot_pication_barcodes(pication_barcodes, save_path='Barcodes/pication_barcodes.png'):
    if not pication_barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(pication_barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(pication_barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

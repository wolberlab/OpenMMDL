"""
mmdl_simulation.py
Perform Simulations of Protein-ligand complexes with OpenMM
"""
import argparse
import sys
import os
import shutil
import argparse
import matplotlib
import pickle
import json
from openmmdl_analysis_function import *
from tqdm import tqdm
from visualization_functions import interacting_water_ids, save_interacting_waters_trajectory, cloud_json_generation
parser = argparse.ArgumentParser()


logo = '\n'.join(["     ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      ",
                  "   .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      ",
                  "  / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      ",
                  " ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    ",
                  " |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    ",
                  " : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    ",
                  "  \ `_/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  ",
                  "   '. \_/``'.'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ ",
                  "     '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` ",
                  "              Prepare and Perform OpenMM Protein-Ligand MD Simulations                                 ",
                  "                                     Alpha Version                                                     "])




if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='openmmdl', description=logo, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-t', dest='topology', help='Topology File after MD Simulation', required=True)
    parser.add_argument('-d', dest='trajectory', help='Trajectory File in DCD Format', required=True)
    parser.add_argument('-l', dest='ligand_sdf', help='Ligand in SDF Format', required=True)       
    parser.add_argument('-n', dest='ligand_name', help='Ligand Name (3 Letter Code in PDB)', required=True)
    parser.add_argument('-b', dest='binding', help='Binding Mode Treshold for Binding Mode in %', default=40)   
    parser.add_argument('-df', dest='dataframe', help='Dataframe (use if the interactions were already calculated, default name would be "df_all.csv")', default=None)
    parser.add_argument('-m', dest='min_transition', help='Minimal Transition % for Markov State Model', default=1)

    input_formats = ['.pdb', '.dcd', '.sdf', '.csv'] 
    args = parser.parse_args()
    if input_formats[0] not in args.topology:
        print("PDB is missing, try the absolute path")
    if input_formats[1] not in args.trajectory:
        print("DCD is missing, try the absolute path")
    if input_formats[2] not in args.ligand_sdf:
        print("SDF is missing, try the absolute path")
    if args.ligand_name == None:
        print("Ligand Name is Missing, Add Ligand Name")

    topology = args.topology
    trajectory = args.trajectory
    ligand_sdf = args.ligand_sdf
    ligand = args.ligand_name
    treshold = int(args.binding)
    dataframe = args.dataframe
    min_transition = args.min_transition

    pdb_md = mda.Universe(topology, trajectory)
    print(len(pdb_md.trajectory)-1)

    complex = pdb_md.select_atoms("protein or resname UNK or (resname HOH and around 10 resname UNK)")
    complex.write("complex.pdb")
    ligand_complex = pdb_md.select_atoms("resname UNK")
    ligand_complex.write("lig.pdb")

    lig_rd = rdkit.Chem.rdmolfiles.MolFromPDBFile("lig.pdb")

    lig_rd_ring = lig_rd.GetRingInfo()


    novel = mda.Universe("complex.pdb")
    novel_lig = novel.select_atoms("resname UNK")

    for atom in novel_lig:
        lig_index = atom.id
        break
    print(lig_index)

    ligand_rings = []

    # Iterate through each ring, increase indices by 1, and print the updated rings
    for atom_ring in lig_rd_ring.AtomRings():
        updated_ring = increase_ring_indices(atom_ring, lig_index)
        print(updated_ring)
        ligand_rings.append(updated_ring)

    # Create a molecule supplier from an SDF file
    mol_supplier = Chem.SDMolSupplier(f"{ligand_sdf}")

    # Open the output SMILES file for writing
    with open("lig.smi", "w") as output_file:
        # Iterate through molecules and convert each to SMILES
        for mol in mol_supplier:
            if mol is not None:  # Make sure the molecule was successfully read
                smiles = Chem.MolToSmiles(mol)
                output_file.write(smiles + "\n")


    # Assuming lig_rd is your RDKit molecule
    m = Chem.MolToSmiles(lig_rd, kekuleSmiles=True)

    # Create an RDKit molecule from the SMILES string
    mol = Chem.MolFromSmiles(m)

    # Check if the conversion was successful
    if mol is not None:
        w1 = Chem.SmilesWriter("./lig2.smi")
        w1.write(mol)  # Write the molecule
        w1.close()
    else:
        print("Failed to convert SMILES string to molecule.")
    
    # Read the content of the lig2.smi file
    with open("./lig2.smi", "r") as file:
        lines = file.readlines()

    # Remove the first line (if there are any lines)
    if lines:
        lines = lines[1:]

    # Write the modified content back to the file
    with open("./lig2.smi", "w") as file:
        file.writelines(lines)
    
    interaction_list = pd.DataFrame(columns=["RESNR", "RESTYPE", "RESCHAIN", "RESNR_LIG", "RESTYPE_LIG", "RESCHAIN_LIG", "DIST", "LIGCARBONIDX", "PROTCARBONIDX", "LIGCOO", "PROTCOO"])

    total_frames = len(pdb_md.trajectory[1:])
    print(dataframe)

    if dataframe == None:
        for num, frame in tqdm(enumerate(pdb_md.trajectory[1:], 1), total=total_frames, desc="Processing frames"):
            sussy = pdb_md.select_atoms("protein or resname UNK or (resname HOH and around 10 resname UNK)")
            sussy.write("md.pdb")
            interactions_by_site = retrieve_plip_interactions("md.pdb")
            index_of_selected_site = -1
            selected_site = list(interactions_by_site.keys())[index_of_selected_site]

    
            interaction_type = "hydrophobic"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)

    
            interaction_type = "hbond"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)
    
            interaction_type = "waterbridge"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            for index, row in tmp_interaction.iterrows():
                nope = int(row['WATER_IDX'])
                if nope != '':
                    water_pdb = mda.Universe("md.pdb")
                    water_index = water_pdb.select_atoms(f"id {nope}")
                    water_index_residues = water_index.resids
                    for water_residue in water_index_residues:
                        water_residue = water_residue
                    water_residues_list = water_residue.tolist()
                    tmp_interaction.loc[index, ['WATER_IDX']] = water_residues_list
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)
    
            interaction_type = "saltbridge"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)

            interaction_type = "pistacking"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)

            interaction_type = "pication"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)

            interaction_type = "halogen"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)

            interaction_type = "metal"
            tmp_interaction= create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
            tmp_interaction['FRAME'] = int(num)
            tmp_interaction['INTERACTION'] = interaction_type
            frames = [interaction_list, tmp_interaction]
            interaction_list = pd.concat(frames)
            interaction_list.to_csv("interactions_gathered.csv")
    elif dataframe != None:
        interaction_tmp = pd.read_csv(f"{dataframe}")
        interaction_list = interaction_tmp.drop(interaction_tmp.columns[0],axis=1)

    interaction_list["Prot_partner"] = interaction_list["RESNR"].astype(str) + interaction_list["RESTYPE"] + interaction_list["RESCHAIN"]

    interaction_list = interaction_list.reset_index(drop=True)

    interaction_list.to_csv("preparation.csv")

    unique_columns_rings_grouped = gather_interactions(interaction_list, ligand_rings)

    interactions_all = interaction_list.copy()

    # Add Frames + Treshold by user
    filtered_values = filtering_values(threshold=treshold/100, frames=len(pdb_md.trajectory)-1, df=interaction_list, unique_columns_rings_grouped=unique_columns_rings_grouped)

    filtering_all = filtering_values(threshold=0.00001, frames=len(pdb_md.trajectory)-1, df=interactions_all, unique_columns_rings_grouped=unique_columns_rings_grouped)

    # Replace NaN values with 0 in the entire DataFrame
    interaction_list.fillna(0, inplace=True)
    interactions_all.fillna(0, inplace=True)

    unique_data = unique_data_generation(filtered_values)
    unique_data_all = unique_data_generation(filtering_all)

    df_iteration_numbering(interaction_list,unique_data)

    df_iteration_numbering(interactions_all,unique_data_all)

    interactions_all.to_csv("df_all.csv")

    # Group by 'FRAME' and transform each group to set all values to 1 if there is at least one 1 in each column
    grouped_frames_treshold = interaction_list.groupby('FRAME', as_index=False)[list(unique_data.values())].max()

    grouped_frames_treshold = grouped_frames_treshold.set_index('FRAME', drop=False)

    update_values(interaction_list, grouped_frames_treshold, unique_data)

    grouped_frames_treshold['FRAME'] = grouped_frames_treshold['FRAME'].astype(int)

    # Extract all columns except 'FRAME' and the index column
    selected_columns = grouped_frames_treshold.columns[1:-1]

    # Create a list of lists with the values from selected columns for each row
    treshold_result_list = [row[selected_columns].values.tolist() for _, row in grouped_frames_treshold.iterrows()]

    # Calculate the occurrences of each list in the result_list
    treshold_occurrences = Counter(tuple(lst) for lst in treshold_result_list)

    # Find the top 50 duplicate occurrences
    top_50_treshold_duplicates = treshold_occurrences.most_common(50)

    # Create a new column 'fingerprint' in the DataFrame
    grouped_frames_treshold['fingerprint'] = None

    # Set the 'fingerprint' column values based on the corresponding index in result_list
    for index, fingerprint_value in enumerate(treshold_result_list,1):
        grouped_frames_treshold.at[index, 'fingerprint'] = fingerprint_value

    # Assuming your original DataFrame is named 'df'
    # First, we'll create a new column 'Binding_fingerprint_hbond'
    grouped_frames_treshold['Binding_fingerprint_treshold'] = ''

    # Dictionary to keep track of encountered fingerprints and their corresponding labels
    treshold_fingerprint_dict = {}

    # Counter to generate the labels (Hbond_Binding_1, Hbond_Binding_2, etc.)
    label_counter = 1

    # Iterate through the rows and process the 'fingerprint' column
    for index, row in grouped_frames_treshold.iterrows():
        fingerprint = tuple(row['fingerprint'])
    
        # Check if the fingerprint has been encountered before
        if fingerprint in treshold_fingerprint_dict:
            grouped_frames_treshold.at[index, 'Binding_fingerprint_treshold'] = treshold_fingerprint_dict[fingerprint]
        else:
            # Assign a new label if the fingerprint is new
            label = f'Binding_Mode_{label_counter}'
            treshold_fingerprint_dict[fingerprint] = label
            grouped_frames_treshold.at[index, 'Binding_fingerprint_treshold'] = label
            label_counter += 1

    # Display the updated DataFrame with the new column
    print(grouped_frames_treshold)

    # Group the DataFrame by the 'Binding_fingerprint_hbond' column and create the dictionary
    fingerprint_dict = grouped_frames_treshold['Binding_fingerprint_treshold'].to_dict()

    combined_dict = {'all': []}
    for key, value in fingerprint_dict.items():
        combined_dict['all'].append(value)

    print(combined_dict)

    transition_1_percent = 0.01 * total_frames
    transition_2_percent = 0.02 * total_frames
    transition_5_percent = 0.05 * total_frames
    transition_10_percent = 0.1 * total_frames

    binding_site_markov_network2(total_frames,min_transition=transition_1_percent,combined_dict=combined_dict)

    binding_site_markov_network2(total_frames,min_transition=transition_2_percent,combined_dict=combined_dict)

    binding_site_markov_network2(total_frames,min_transition=transition_5_percent,combined_dict=combined_dict)

    binding_site_markov_network2(total_frames,min_transition=transition_10_percent,combined_dict=combined_dict)

    # Get the top 10 nodes with the most occurrences
    node_occurrences = {node: combined_dict['all'].count(node) for node in set(combined_dict['all'])}
    top_10_nodes = sorted(node_occurrences, key=node_occurrences.get, reverse=True)[:10]

    top_10_nodes_with_occurrences = {node: node_occurrences[node] for node in top_10_nodes}

    # The provided dictionary with Binding_fingerprint_treshold
    binding_fingerprint_treshold = top_10_nodes_with_occurrences

    columns_with_value_1 = {}  # Initialize an empty dictionary to store the result

    for treshold, row_count in binding_fingerprint_treshold.items():
        for i in range(1, row_count + 1):
            # Get the row corresponding to the treshold value
            row = grouped_frames_treshold.loc[grouped_frames_treshold['Binding_fingerprint_treshold'] == treshold].iloc[i - 1]

            # Extract column names with value 1 in that row
            columns_with_1 = row[row == 1].index.tolist()

            # Convert the list to a set to remove duplicates
            columns_set = set(columns_with_1)

            # Add the columns to the dictionary under the corresponding threshold
            if treshold not in columns_with_value_1:
                columns_with_value_1[treshold] = set()
            columns_with_value_1[treshold].update(columns_set)



    matplotlib.use("Agg") 
    binding_site = {}
    merged_image_paths = []
    warnings.resetwarnings()
    for key, values in columns_with_value_1.items():
        binding_site[key] = values
        occurrence_count = top_10_nodes_with_occurrences[key]
        occurrence_percent = 100* occurrence_count / total_frames
        binding_mode = key
        with open("lig.smi", "r") as file:
            reference_smiles = file.read().strip()  # Read the SMILES from the file and remove any leading/trailing whitespace
        reference_mol = Chem.MolFromSmiles(reference_smiles)
        prepared_ligand = AllChem.AssignBondOrdersFromTemplate(reference_mol, lig_rd)



        # Generate 2D coordinates for the molecule
        AllChem.Compute2DCoords(prepared_ligand)

        split_data = split_interaction_data(values)

        # Get the highlighted atom indices based on interaction type
        highlighted_hbond_donor, highlighted_hbond_acceptor, highlighted_hbond_both, highlighted_hydrophobic, highlighted_waterbridge, highlighted_pistacking, highlighted_halogen, highlighted_ni, highlighted_pi, highlighted_pication, highlighted_metal = highlight_numbers(split_data, starting_idx=lig_index)
        highlighted_hbond_donor = [int(x) for x in highlighted_hbond_donor]
        highlighted_hbond_acceptor = [int(x) for x in highlighted_hbond_acceptor]
        highlighted_hbond_both = [int(x) for x in highlighted_hbond_both]
        highlighted_hydrophobic = [int(x) for x in highlighted_hydrophobic]
        highlighted_waterbridge = [int(x) for x in highlighted_waterbridge]
        highlighted_pistacking = [int(x) for x in highlighted_pistacking]
        highlighted_halogen = [int(x) for x in highlighted_halogen]
        highlighted_ni = [int(x) for x in highlighted_ni]
        highlighted_pi = [int(x) for x in highlighted_pi]
        highlighted_pication = [int(x) for x in highlighted_pication]
        highlighted_metal = [int(x) for x in highlighted_metal]
    
    
        # Generate a dictionary for hydrogen bond acceptors
        hbond_acceptor_dict = {key: (1.0,0.6,0.6) for key in highlighted_hbond_acceptor}

        # Generate a dictionary for hydrogen bond acceptors and donors
        hbond_both_dict = {key: (0.6,0.0,0.5) for key in highlighted_hbond_both}

        # Generate a dictionary for hydrogen bond donors
        hbond_donor_dict = {key: (0.3,0.5,1) for key in highlighted_hbond_donor}

        # Generate a dictionary for hydrophobic features
        hydrophobic_dict = {key: (1.0, 1.0, 0.0) for key in highlighted_hydrophobic}

        # Generate a dictionary for water bridge interactions
        waterbridge_dict = {key: (0.0, 1.0, 0.9) for key in highlighted_waterbridge}

        # Generate a dictionary for pistacking
        pistacking_dict = {key: (0.0, 0.0, 1.0) for key in highlighted_pistacking}

        # Generate a dictionary for halogen interactions
        halogen_dict = {key: (1.0, 0.0, 0.9) for key in highlighted_halogen}

        # Generate a dictionary for negative ionizables
        ni_dict = {key: (0.0, 0.0, 1.0) for key in highlighted_ni}

        # Generate a dictionary for negative ionizables
        pi_dict = {key: (1.0, 0.0, 0.0) for key in highlighted_pi}

        # Generate a dictionary for pication
        pication_dict = {key: (0.0, 0.0, 1.0) for key in highlighted_pication}

        # Generate a dictionary for metal interactions
        metal_dict = {key: (1.0, 0.6, 0.0) for key in highlighted_metal}

        hbond_donor_dict.update((k, v) for k, v in hbond_acceptor_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in hydrophobic_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in hbond_both_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in waterbridge_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in pistacking_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in halogen_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in ni_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in pi_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in pication_dict.items() if k not in hbond_donor_dict)
        hbond_donor_dict.update((k, v) for k, v in metal_dict.items() if k not in hbond_donor_dict)
    
        higlight_atoms = highlighted_hbond_donor + highlighted_hbond_acceptor + highlighted_hbond_both + highlighted_hydrophobic + highlighted_waterbridge + highlighted_pistacking + highlighted_halogen + highlighted_ni + highlighted_pi + highlighted_pication + highlighted_metal 

        higlight_atoms = list(set(higlight_atoms))
    
        # Convert the RDKit molecule to SVG format with atom highlights
        drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)

        #drawer.DrawMolecule(prepared_ligand, highlightAtoms=highlighted_hydrophobic,  highlightAtomColors={k: (1.0, 1.0, 0.0)for k in highlighted_hydrophobic})
        drawer.DrawMolecule(prepared_ligand, highlightAtoms=higlight_atoms, highlightAtomColors=hbond_donor_dict)

        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:', '')

        # Save the SVG to a file
        with open(f'{binding_mode}.svg', 'w') as f:
            f.write(svg)

        cairosvg.svg2png(url=f'{binding_mode}.svg', write_to=f'{binding_mode}.png')

    
        # Create the main figure and axis (unchanged)
        fig = pylab.figure()
        ax = fig.add_subplot(111)

        # Data for the x-axis and random data for demonstration (unchanged)
        x = range(10)
        data_points = [pylab.randn(10) for _ in range(len(split_data))]


        # Plot lines on the same axis and collect them into a list (unchanged)
        lines = []
        for i, data in enumerate(split_data):
            y = data_points[i]
            label = data.split()[-1]
            print(data)
            if label == 'hydrophobic':
                line, = ax.plot(x, y, label=data, color=(0.8, 1, 0), linewidth=5.0)
            elif label == 'hbond':
                line, = ax.plot(x, y, label=data, color=(1.0, 0.6, 0.6), linewidth=5.0)
            else:
                line, = ax.plot(x, y, label=data)
            lines.append(line)

        # Create a separate figure for the legend (unchanged)
        figlegend = pylab.figure(figsize=(8, 6))

        # Add a legend to the subplot (ax) using the lines and full entries as labels (unchanged)
        legend = figlegend.legend(lines, split_data, loc='center')

        # Set the linewidth of the legend lines to be thicker (unchanged)
        for line in legend.get_lines():
            line.set_linewidth(5.0)

        # Add the text "Binding_Mode_9" above the legend and closer to it
        figlegend.text(0.5, 0.9, f"{binding_mode}", ha='center', fontsize=12, weight='bold')
        figlegend.text(0.5, 0.85, f"Occurence {occurrence_percent}%", ha='center', fontsize=12, weight='bold')


        # Save the legend figure to a file (unchanged)
        figlegend.savefig(f'{binding_mode}_legend.png')

        merged_image_paths.append(f"{binding_mode}_merged.png")


        #Read the two images
        image1 = Image.open(f'{binding_mode}.png')
        image2 = Image.open(f'{binding_mode}_legend.png')
        #resize, first image
        image1_size = image1.size
        image2_size = image2.size
        total_width = image1_size[0] + image2_size[0]
        new_image = Image.new('RGB', (total_width, image1_size[1]))
        new_image.paste(image1,(0,0))
        new_image.paste(image2,(image1_size[0],0))
        new_image.save(f"{binding_mode}_merged.png","PNG")
        os.remove(f'{binding_mode}.png')
        os.remove(f'{binding_mode}_legend.png')
        os.remove(f'{binding_mode}.svg')

    merged_images = [Image.open(path) for path in merged_image_paths]

    # Calculate the maximum width and height for the images
    max_width = max(image.size[0] for image in merged_images)
    max_height = max(image.size[1] for image in merged_images)

    # Determine the number of images per row (in your case, 2 images per row)
    images_per_row = 2

    # Calculate the number of rows and columns required
    num_rows = (len(merged_images) + images_per_row - 1) // images_per_row
    total_width = max_width * images_per_row
    total_height = max_height * num_rows

    # Create a new image with the calculated width and height
    big_figure = Image.new('RGB', (total_width, total_height), (255, 255, 255))  # Set background to white

    x_offset = 0
    y_offset = 0

    for image in merged_images:
        # Paste the image onto the big figure
        big_figure.paste(image, (x_offset, y_offset))
    
        # Update offsets
        x_offset += max_width
    
        # Move to the next row if necessary
        if x_offset >= total_width:
            x_offset = 0
            y_offset += max_height

    # Save the big figure
    big_figure.save("all_binding_modes_arranged.png", "PNG")
    os.rename("all_binding_modes_arranged.png", "Binding_Modes_Markov_States/all_binding_modes_arranged.png")
    for path in merged_image_paths:
        os.remove(path)


    df_all = pd.read_csv('df_all.csv')
    hydrophobic_interactions = df_all.filter(regex='hydrophobic').columns
    acceptor_interactions = df_all.filter(regex='Acceptor_hbond').columns
    donor_interactions = df_all.filter(regex='Donor_hbond').columns
    pistacking_interactions = df_all.filter(regex='pistacking').columns
    halogen_interactions = df_all.filter(regex='halogen').columns
    waterbridge_interactions = df_all.filter(regex='waterbridge').columns
    pication_interactions = df_all.filter(regex='pication').columns
    saltbridge_ni_interactions = df_all.filter(regex='NI').columns
    saltbridge_pi_interactions = df_all.filter(regex='PI').columns

    hydrophobicinteraction_barcodes = {}
    for hydrophobic_interaction in hydrophobic_interactions:
        barcode = barcodegeneration(df_all, hydrophobic_interaction)
        hydrophobicinteraction_barcodes[hydrophobic_interaction] = barcode

    acceptor_barcodes = {}
    for acceptor_interaction in acceptor_interactions:
        barcode = barcodegeneration(df_all, acceptor_interaction)
        acceptor_barcodes[acceptor_interaction] = barcode

    donor_barcodes = {}
    for donor_interaction in donor_interactions:
        barcode = barcodegeneration(df_all, donor_interaction)
        donor_barcodes[donor_interaction] = barcode

    pistacking_barcodes = {}
    for pistacking_interaction in pistacking_interactions:
        barcode = barcodegeneration(df_all, pistacking_interaction)
        pistacking_barcodes[pistacking_interaction] = barcode

    halogen_barcodes = {}
    for halogen_interaction in halogen_interactions:
        barcode = barcodegeneration(df_all, halogen_interaction)
        halogen_barcodes[halogen_interaction] = barcode

    waterbridge_barcodes = {}
    for waterbridge_interaction in waterbridge_interactions:
        barcode = barcodegeneration(df_all, waterbridge_interaction)
        waterbridge_barcodes[waterbridge_interaction] = barcode

    pication_barcodes = {}
    for pication_interaction in pication_interactions:
        barcode = barcodegeneration(df_all, pication_interaction)
        pication_barcodes[pication_interaction] = barcode

    plot_hydrophobic_barcodes(hydrophobicinteraction_barcodes)
    plot_acceptor_barcodes(acceptor_barcodes)
    plot_donor_barcodes(donor_barcodes)
    plot_pistacking_barcodes(pistacking_barcodes)
    plot_halogen_barcodes(halogen_barcodes)
    plot_waterbridge_barcodes(waterbridge_barcodes)
    plot_pication_barcodes(pication_barcodes)

    plot_waterbridge_piechart(df_all, waterbridge_barcodes, waterbridge_interactions)

  
    interacting_water_id_list = interacting_water_ids(df_all, waterbridge_interactions)

    # dump interacting waters for visualization
    with open('interacting_waters', 'wb') as f:
        pickle.dump(interacting_water_id_list, f)
    save_interacting_waters_trajectory(topology, trajectory, interacting_water_id_list)

    # save clouds for visualization
    with open('clouds.json', 'w') as f:
        json.dump(cloud_json_generation(df_all), f)
    

    

import os
import numpy as np
import matplotlib.pyplot as plt


def barcodegeneration(df, interaction):
    """Generates barcodes for a given interaction  .

    Args:
        df (pandas dataframe): Dataframe containing all interactions from plip analysis (typicaly df_all)
        interaction (str): name of the interaction to generate a barcode for

    Returns:
        numpy array: returns an binary array of wit 1 representing the interaction is present in the corresponding frame
    """
    barcode = []
    
    unique_frames = df['FRAME'].unique()
    
    for frame in unique_frames:
        frame_data = df[df['FRAME'] == frame]
        
        if 1 in frame_data[interaction].values:
            barcode.append(1)
            
            
        else:
            barcode.append(0)
    
    return np.array(barcode)


def waterids_barcode_generator(df, interaction):
    """Generates a barcode containing coresponding water ids for a given interaction.

    Args:
        df (pandas dataframe): dataframe containing all interactions from plip analysis (typicaly df_all)
        interaction (str): name of the interaction to generate a barcode for

    Returns:
        list: returns a list of waterids for the frames where the interaction is present 0 if no interaction present
    """
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


def plot_barcodes(barcodes, save_path):
    """Generates picture of barcodes for interactions of a specific type.

    Args:
        barcodes (list): list of np arrays containing the barcodes for each interaction
        save_path (str): name of the file to save the picture to
    """     
    if not barcodes:
        print("No barcodes to plot.")
        return
    
    num_plots = len(barcodes)
    num_cols = 1
    num_rows = (num_plots + num_cols - 1) // num_cols

    fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))
    
    # If only one row, axs is a single Axes object, not an array
    if num_rows == 1:
        axs = [axs]

    for i, (title, barcode) in enumerate(barcodes.items()):
        ax = axs[i]
        ax.set_axis_off()
        im = ax.imshow(barcode.reshape(1, -1), cmap='binary', aspect='auto', interpolation='nearest')

        percent_occurrence = (barcode.sum() / len(barcode)) * 100
        ax.text(1.05, 0.5, f"{percent_occurrence:.2f}%", transform=ax.transAxes, va='center', fontsize=8)

        ax.set_title(title, fontweight='bold', fontsize=8)

    os.makedirs(os.path.dirname("./Barcodes/"), exist_ok=True)
    plt.tight_layout()
    plt.savefig(f"./Barcodes/{save_path}", dpi=300, bbox_inches='tight')


def plot_waterbridge_piechart(df_all, waterbridge_barcodes, waterbridge_interactions):
    """Generates piecharts for each waterbridge interaction with the water ids of the interacting waters.

    Args:
        df_all (pandas dataframe): dataframe contaning all interactions (typicaly df_all)
        waterbridge_barcodes (list): list of np arrays containing the barcodes for each interaction
        waterbridge_interactions (list): list of strings containing the names of the waterbridge interactions
    """
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
        
        plt.pie(values, labels=labels,
                autopct=lambda pct: f'{pct:.1f}%\n({int(round(pct/100.0 * sum(values)))})',
                shadow=False,
                startangle=140)
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

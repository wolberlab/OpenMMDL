import os
import numpy as np
import matplotlib.pyplot as plt


class BarcodeGenerator:
    """
    Generates binary barcodes representing the presence of interactions across MD frames.

    Attributes
    ----------
    df : pd.DataFrame
        DataFrame containing all interactions extracted from PLIP analysis.
    interactions : dict
        Dictionary mapping interaction types to their corresponding columns in the DataFrame.
    """
    def __init__(self, df):
        self.df = df
        self.interactions = self._gather_interactions()

    def _gather_interactions(self):
        """
        Gathers interaction column names grouped by the corresponding interaction type.
    
        Returns
        -------
        dict
            Dictionary where the keys are interaction types and values are lists of corresponding column names.
        """
        hydrophobic_interactions = self.df.filter(regex="hydrophobic").columns
        acceptor_interactions = self.df.filter(regex="Acceptor_hbond").columns
        donor_interactions = self.df.filter(regex="Donor_hbond").columns
        pistacking_interactions = self.df.filter(regex="pistacking").columns
        halogen_interactions = self.df.filter(regex="halogen").columns
        waterbridge_interactions = self.df.filter(regex="waterbridge").columns
        pication_interactions = self.df.filter(regex="pication").columns
        saltbridge_ni_interactions = self.df.filter(regex="NI_saltbridge").columns
        saltbridge_pi_interactions = self.df.filter(regex="PI_saltbridge").columns
        metal_interactions = self.df.filter(regex="metal").columns

        return {
            "hydrophobic": hydrophobic_interactions,
            "acceptor": acceptor_interactions,
            "donor": donor_interactions,
            "pistacking": pistacking_interactions,
            "halogen": halogen_interactions,
            "waterbridge": waterbridge_interactions,
            "pication": pication_interactions,
            "saltbridge_ni": saltbridge_ni_interactions,
            "saltbridge_pi": saltbridge_pi_interactions,
            "metal": metal_interactions,
        }

    def generate_barcode(self, interaction):
        """
        Generates barcodes for a given interaction.

        Parameters
        ----------
        interaction : str
            Name of the interaction to generate a barcode for

        Returns
        -------
        np.ndarray
            Binary array with 1 representing the interaction is present in the corresponding frame
        """
        barcode = []
        unique_frames = self.df["FRAME"].unique()

        for frame in unique_frames:
            frame_data = self.df[self.df["FRAME"] == frame]
            if 1 in frame_data[interaction].values:
                barcode.append(1)
            else:
                barcode.append(0)

        return np.array(barcode)

    def _generate_waterids_barcode(self, interaction):
        """
        Generates a barcode containing corresponding water ids for a given interaction.

        Parameters
        ----------
        interaction : str
            Name of the interaction to generate a barcode for.

        Returns
        -------
        list of int
            List with water IDs for frames where interaction is present, 0 otherwise.
    
        Notes
        -----
        Water IDs are obtained from the "WATER_IDX" column in the DataFrame.
        """
        water_id_list = []
        waterid_barcode = []

        for index, row in self.df.iterrows():
            if row[interaction] == 1:
                water_id_list.append(int(float(row["WATER_IDX"])))

        barcode = self.generate_barcode(interaction)

        for value in barcode:
            if value == 1:
                waterid_barcode.append(water_id_list.pop(0))
            else:
                waterid_barcode.append(0)

        return waterid_barcode

    def interacting_water_ids(self, waterbridge_interactions):
        """
        Generates a list of all water ids that form water bridge interactions.

        Parameters
        ----------
        waterbridge_interactions : list of str
            list containing the names of all water bridge interactions.

        Returns
        -------
        list of int
            Unique water IDs that form waterbridge interactions.
        """
        interacting_waters = []
        for waterbridge_interaction in waterbridge_interactions:
            waterid_barcode = self._generate_waterids_barcode(waterbridge_interaction)
            for waterid in waterid_barcode:
                if waterid != 0:
                    interacting_waters.append(waterid)

        return list(set(interacting_waters))


class BarcodePlotter:
    """
    Visualizes interaction barcodes and waterbridge statistics using bar plots and pie charts.

    Attributes
    ----------
    df_all : pd.DataFrame
        Full PLIP interaction dataframe passed for plotting.
    barcode_gen : BarcodeGenerator
        Instance of BarcodeGenerator used to compute barcodes from interaction data.
    """
    def __init__(self, df_all):
        self.df_all = df_all
        self.barcode_gen = BarcodeGenerator(df_all)

    def _plot_barcodes(self, barcodes, save_path):
        """
        Plots barcodes of the interactions depending on the presence of the interaction.
    
        Parameters
        ----------
        barcodes : dict
            Dictionary where keys are interaction names and values are 1D numpy arrays (barcodes).
        save_path : str
            Path to save the generated barcode plot image.
    
        Returns
        -------
        None
        """
        if not barcodes:
            print("No barcodes to plot.")
            return

        num_plots = len(barcodes)
        num_cols = 1
        num_rows = (num_plots + num_cols - 1) // num_cols

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(8.50, num_rows * 1))

        if num_rows == 1:
            axs = [axs]

        for i, (title, barcode) in enumerate(barcodes.items()):
            ax = axs[i]
            ax.set_axis_off()
            ax.imshow(
                barcode.reshape(1, -1),
                cmap="binary",
                aspect="auto",
                interpolation="nearest",
                vmin=0,
                vmax=1,
            )

            percent_occurrence = (barcode.sum() / len(barcode)) * 100
            ax.text(
                1.05,
                0.5,
                f"{percent_occurrence:.2f}%",
                transform=ax.transAxes,
                va="center",
                fontsize=8,
            )
            ax.set_title(title, fontweight="bold", fontsize=8)

        save_dir = os.path.dirname(save_path)
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    def plot_waterbridge_piechart(
        self, waterbridge_barcodes, waterbridge_interactions, fig_type
    ):
        """
        Generates and saves pie charts showing the frequency of each of the water IDs participating in waterbridge interactions.
    
        Parameters
        ----------
        waterbridge_barcodes : dict
            Dictionary of waterbridge interaction barcodes.
        waterbridge_interactions : list of str
            List of interaction column names related to waterbridge interactions.
        fig_type : str
            Image file format for saving (e.g., 'png', 'svg').
    
        Returns
        -------
        None
        """
        if not waterbridge_barcodes:
            print("No Piecharts to plot.")
            return

        os.makedirs("Barcodes/Waterbridge_Piecharts", exist_ok=True)
        plt.figure(figsize=(6, 6))

        for waterbridge_interaction in waterbridge_interactions:
            plt.clf()
            waterid_barcode = self.barcode_gen._generate_waterids_barcode(
                waterbridge_interaction
            )
            waters_count = {}

            for waterid in waterid_barcode:
                if waterid != 0:
                    waters_count[waterid] = waters_count.get(waterid, 0) + 1

            labels = [f"ID {id}" for id in waters_count.keys()]
            values = list(waters_count.values())

            threshold = 7
            total_values = sum(values)
            small_ids = [
                id
                for id, value in waters_count.items()
                if (value / total_values) * 100 < threshold
            ]

            if small_ids:
                small_count = sum(
                    count for id, count in waters_count.items() if id in small_ids
                )
                values = [
                    count if id not in small_ids else small_count
                    for id, count in waters_count.items()
                ]
                labels = [
                    f"ID {id}" if id not in small_ids else ""
                    for id in waters_count.keys()
                ]

            plt.pie(
                values,
                labels=labels,
                autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct/100.0 * sum(values)))})",
                shadow=False,
                startangle=140,
            )
            plt.axis("equal")
            plt.title(str(waterbridge_interaction), fontweight="bold")
            legend_labels = [f"ID {id}" for id in waters_count.keys()]
            legend = plt.legend(
                legend_labels, loc="upper right", bbox_to_anchor=(1.2, 1)
            )
            plt.setp(legend.get_texts(), fontsize="small")
            plt.text(
                0.5,
                0,
                f"Total frames with waterbridge: {round(((sum(1 for val in waterid_barcode if val != 0) / len(waterid_barcode)) * 100), 2)}%",
                size=12,
                ha="center",
                transform=plt.gcf().transFigure,
            )
            plt.subplots_adjust(top=0.99, bottom=0.01)
            plt.savefig(
                f"Barcodes/Waterbridge_Piecharts/{waterbridge_interaction}.{fig_type}",
                bbox_inches="tight",
                dpi=300,
            )

    def plot_barcodes_grouped(self, interactions, interaction_type, fig_type):
        """
        Groups barcodes by ligand atom, plots individual and grouped barcodes, and saves them.
    
        Parameters
        ----------
        interactions : list of str
            List of interaction names to be grouped and visualized.
        interaction_type : str 
            The type of interaction (e.g., 'donor', 'acceptor', 'waterbridge').
        fig_type : str 
            Image file format for saving (e.g., 'png', 'svg').
    
        Returns
        -------
        None
        """
        ligatoms_dict = {}
        for interaction in interactions:
            ligatom = interaction.split("_")
            ligatom.pop(0)
            ligatom.pop(-1)
            if interaction_type in [
                "acceptor",
                "donor",
                "waterbridge",
                "saltbridge_ni",
                "saltbridge_pi",
            ]:
                ligatom.pop(-1)
                if interaction_type in ["saltbridge_ni", "saltbridge_pi"]:
                    ligatom.pop(-1)
            ligatom = "_".join(ligatom)
            ligatoms_dict.setdefault(ligatom, []).append(interaction)

        total_interactions = {}
        for ligatom in ligatoms_dict:
            ligatom_interaction_barcodes = {}
            for interaction in ligatoms_dict[ligatom]:
                barcode = self.barcode_gen.generate_barcode(interaction)
                ligatom_interaction_barcodes[interaction] = barcode
            os.makedirs(f"./Barcodes/{ligatom}", exist_ok=True)
            self._plot_barcodes(
                ligatom_interaction_barcodes,
                f"./Barcodes/{ligatom}/{ligatom}_{interaction_type}_barcodes.{fig_type}",
            )

            barcodes_list = list(ligatom_interaction_barcodes.values())
            grouped_array = np.logical_or.reduce(barcodes_list)
            grouped_array[np.all(np.vstack(barcodes_list) == 0, axis=0)] = 0
            grouped_array = grouped_array.astype(int)
            total_interactions[ligatom] = grouped_array

        self._plot_barcodes(
            total_interactions, f"./Barcodes/{interaction_type}_interactions.{fig_type}"
        )

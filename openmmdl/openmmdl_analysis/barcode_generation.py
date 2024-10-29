import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import List, Dict, Union


class BarcodeGenerator:
    def __init__(self, df: pd.DataFrame):
        """
        Initializes the BarcodeGenerator with a dataframe.

        Args:
            df (pd.DataFrame): Dataframe containing all interactions from plip analysis (typically df_all).
        """
        self.df = df
        self.interactions = self.gather_interactions()

    def gather_interactions(self) -> Dict[str, pd.Index]:
        """
        Gathers columns related to different types of interactions.

        Returns:
            Dict[str, pd.Index]: Dictionary where keys are interaction types and values are columns corresponding to those interactions.
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

    def generate_barcode(self, interaction: str) -> np.ndarray:
        """
        Generates barcodes for a given interaction.

        Args:
            interaction (str): Name of the interaction to generate a barcode for.

        Returns:
            np.ndarray: Binary array with 1 representing the interaction is present in the corresponding frame.
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

    def generate_waterids_barcode(self, interaction: str) -> List[Union[int, int]]:
        """
        Generates a barcode containing corresponding water ids for a given interaction.

        Args:
            interaction (str): Name of the interaction to generate a barcode for.

        Returns:
            List[Union[int, int]]: List of water ids for the frames where the interaction is present, 0 if no interaction present.
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

    def interacting_water_ids(self, waterbridge_interactions: List[str]) -> List[int]:
        """
        Generates a list of all water ids that form water bridge interactions.

        Args:
            waterbridge_interactions (List[str]): List of strings containing the names of all water bridge interactions.

        Returns:
            List[int]: List of all unique water ids that form water bridge interactions.
        """
        interacting_waters = []
        for waterbridge_interaction in waterbridge_interactions:
            waterid_barcode = self.generate_waterids_barcode(waterbridge_interaction)
            for waterid in waterid_barcode:
                if waterid != 0:
                    interacting_waters.append(waterid)
        return list(set(interacting_waters))


class BarcodePlotter:
    def __init__(self, df_all: pd.DataFrame):
        """
        Initializes the BarcodePlotter with a dataframe and BarcodeGenerator instance.

        Args:
            df_all (pd.DataFrame): Dataframe containing all interactions from PLIP analysis.
        """
        self.df_all = df_all
        self.barcode_gen = BarcodeGenerator(df_all)

    def plot_barcodes(self, barcodes: Dict[str, np.ndarray], save_path: str) -> None:
        """
        Plots barcodes and saves the figure to the specified path.

        Args:
            barcodes (Dict[str, np.ndarray]): Dictionary where keys are interaction names and values are binary barcode arrays.
            save_path (str): Path to save the plotted figure.
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

        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches="tight")

    def plot_waterbridge_piechart(
        self,
        waterbridge_barcodes: Dict[str, np.ndarray],
        waterbridge_interactions: List[str],
        fig_type: str,
    ) -> None:
        """
        Plots pie charts for waterbridge interactions and saves them to files.

        Args:
            waterbridge_barcodes (Dict[str, np.ndarray]): Dictionary where keys are interaction names and values are binary barcode arrays.
            waterbridge_interactions (List[str]): List of water bridge interaction names.
            fig_type (str): File extension for the saved figures (e.g., 'png', 'svg').
        """
        if not waterbridge_barcodes:
            print("No Piecharts to plot.")
            return

        os.makedirs("Barcodes/Waterbridge_Piecharts", exist_ok=True)
        plt.figure(figsize=(6, 6))

        for waterbridge_interaction in waterbridge_interactions:
            plt.clf()
            waterid_barcode = self.barcode_gen.generate_waterids_barcode(
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
                autopct=lambda pct: f"{pct:.1f}%\n({int(round(pct / 100.0 * sum(values)))})",
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

    def plot_barcodes_grouped(
        self, interactions: List[str], interaction_type: str, fig_type: str
    ) -> None:
        """
        Plots grouped barcodes for interactions and saves the figure to a file.

        Args:
            interactions (List[str]): List of interaction names.
            interaction_type (str): Type of interaction for grouping.
            fig_type (str): File extension for the saved figure (e.g., 'png', 'pdf').
        """
        ligatoms_dict: Dict[str, List[str]] = {}
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
            self.plot_barcodes(
                ligatom_interaction_barcodes,
                f"./Barcodes/{ligatom}/{ligatom}_{interaction_type}_barcodes.{fig_type}",
            )

            barcodes_list = list(ligatom_interaction_barcodes.values())
            grouped_array = np.logical_or.reduce(barcodes_list)
            grouped_array[np.all(np.vstack(barcodes_list) == 0, axis=0)] = 0
            grouped_array = grouped_array.astype(int)
            total_interactions[ligatom] = grouped_array

        self.plot_barcodes(
            total_interactions, f"./Barcodes/{interaction_type}_interactions.{fig_type}"
        )

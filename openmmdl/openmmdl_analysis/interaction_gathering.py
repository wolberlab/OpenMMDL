import os
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from plip.basic import config
from plip.structure.preparation import PDBComplex, PLInteraction
from plip.exchange.report import BindingSiteReport
from multiprocessing import Pool
from functools import partial
from typing import List, Dict, Tuple, Optional, Union

config.KEEPMOD = True


class InteractionAnalyzer:
    def __init__(
        self, 
        pdb_md: mda.Universe, 
        dataframe: Optional[str], 
        num_processes: int, 
        lig_name: str, 
        special_ligand: Optional[str], 
        peptide: Optional[str]
    ) -> None:
        self.pdb_md = pdb_md
        self.dataframe = dataframe
        self.num_processes = num_processes
        self.lig_name = lig_name
        self.special = special_ligand
        self.peptide = peptide
        self.interaction_list = self.process_trajectory()

    def characterize_complex(
        self, pdb_file: str, binding_site_id: str
    ) -> Optional[PLInteraction]:
        """Characterize the protein-ligand complex and return their interaction set.

        Args:
            pdb_file (str): Path to the PDB file.
            binding_site_id (str): Identifier of the binding site.

        Returns:
            Optional[PLInteraction]: Interactions if found; otherwise, None.
        """
        pdb_complex = PDBComplex()
        pdb_complex.load_pdb(pdb_file)
        for ligand in pdb_complex.ligands:
            if (
                ":".join([ligand.hetid, ligand.chain, str(ligand.position)])
                == binding_site_id
            ):
                pdb_complex.characterize_complex(ligand)
                return pdb_complex.interaction_sets[binding_site_id]

        return None

    def retrieve_plip_interactions(
        self, pdb_file: str, lig_name: str
    ) -> Dict[str, Dict[str, List]]:
        """Retrieve interactions from PLIP.

        Args:
            pdb_file (str): Path to the PDB file of the complex.
            lig_name (str): Name of the ligand in the complex topology.

        Returns:
            Dict[str, Dict[str, List]]: Dictionary of binding sites and interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)
        for ligand in protlig.ligands:
            if str(ligand.longname) == lig_name:
                protlig.characterize_complex(ligand)
        sites: Dict[str, Dict[str, List]] = {}
        for key, site in sorted(protlig.interaction_sets.items()):
            binding_site = BindingSiteReport(site)
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
            interactions = {
                k: [getattr(binding_site, k + "_features")]
                + getattr(binding_site, k + "_info")
                for k in keys
            }
            sites[key] = interactions

        return sites

    def retrieve_plip_interactions_peptide(
        self, pdb_file: str
    ) -> Dict[str, Dict[str, List]]:
        """Retrieve interactions from PLIP for a peptide.

        Args:
            pdb_file (str): Path to the PDB file of the complex.

        Returns:
            Dict[str, Dict[str, List]]: Dictionary of binding sites and interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)
        protlig.characterize_complex(protlig.ligands[-1])
        sites: Dict[str, Dict[str, List]] = {}
        for key, site in sorted(protlig.interaction_sets.items()):
            binding_site = BindingSiteReport(site)
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
            interactions = {
                k: [getattr(binding_site, k + "_features")]
                + getattr(binding_site, k + "_info")
                for k in keys
            }
            sites[key] = interactions

        return sites

    def create_df_from_binding_site(
        self, selected_site_interactions: Dict[str, List], interaction_type: str = "hbond"
    ) -> pd.DataFrame:
        """Create a DataFrame from a binding site and interaction type.

        Args:
            selected_site_interactions (Dict[str, List]): Precalculated interactions from PLIP for the selected site.
            interaction_type (str, optional): Interaction type of interest (default is "hbond").

        Returns:
            pd.DataFrame: DataFrame with information from PLIP.
        """
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
            print(
                "\033[1m!!! Wrong interaction type specified. Hbond is chosen by default!!!\033[0m\n"
            )
            interaction_type = "hbond"

        df = pd.DataFrame.from_records(
            selected_site_interactions[interaction_type][1:],
            columns=selected_site_interactions[interaction_type][0],
        )
        return df

    def change_lig_to_residue(
        self, file_path: str, new_residue_name: str
    ) -> None:
        """Reformat the topology file to change the ligand to a residue.

        Args:
            file_path (str): Path to the topology file.
            new_residue_name (str): New residue name.
        """
        with open(file_path, "r") as file:
            lines = file.readlines()

        with open(file_path, "w") as file:
            for line in lines:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()

                    if residue_name == self.lig_name:
                        modified_line = line[:17] + new_residue_name + line[20:]
                        file.write(modified_line)
                    else:
                        file.write(line)
                else:
                    file.write(line)

    def process_frame_special(self, frame: int) -> pd.DataFrame:
        """Process a single frame of MD simulation for special ligands.

        Args:
            frame (int): Frame number to be processed.

        Returns:
            pd.DataFrame: DataFrame containing interaction data for the processed frame.
        """
        atoms_selected = self.pdb_md.select_atoms(
            f"protein or nucleic or resname {self.lig_name} or (resname HOH and around 10 resname {self.lig_name}) or resname {self.special}"
        )
        for num in self.pdb_md.trajectory[(frame):(frame + 1)]:
            atoms_selected.write(f"processing_frame_{frame}.pdb")
        interaction_list = pd.DataFrame()
        interactions_by_site = self.retrieve_plip_interactions(
            f"processing_frame_{frame}.pdb", self.lig_name
        )
        index_of_selected_site = -1
        selected_site = list(interactions_by_site.keys())[index_of_selected_site]

        interaction_types = [
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        ]

        for interaction_type in interaction_types:
            tmp_interaction = self.create_df_from_binding_site(
                interactions_by_site[selected_site],
                interaction_type=interaction_type,
            )
            tmp_interaction["frame"] = frame
            tmp_interaction["interaction_type"] = interaction_type
            interaction_list = pd.concat([interaction_list, tmp_interaction])
        
        os.remove(f"processing_frame_{frame}.pdb")
        return interaction_list

    def process_frame(self, frame: int) -> pd.DataFrame:
        """Process a single frame of MD simulation.

        Args:
            frame (int): Frame number to be processed.

        Returns:
            pd.DataFrame: DataFrame containing interaction data for the processed frame.
        """
        atoms_selected = self.pdb_md.select_atoms(
            f"protein or nucleic or resname {self.lig_name} or (resname HOH and around 10 resname {self.lig_name}) or resname {self.special}"
        )
        for num in self.pdb_md.trajectory[(frame):(frame + 1)]:
            atoms_selected.write(f"processing_frame_{frame}.pdb")
        interaction_list = pd.DataFrame()
        if self.peptide is None:
            interactions_by_site = self.retrieve_plip_interactions(
                f"processing_frame_{frame}.pdb", self.lig_name
            )
            index_of_selected_site = -1
            selected_site = list(interactions_by_site.keys())[index_of_selected_site]

            interaction_types = [
                "hydrophobic",
                "hbond",
                "waterbridge",
                "saltbridge",
                "pistacking",
                "pication",
                "halogen",
                "metal",
            ]

            for interaction_type in interaction_types:
                tmp_interaction = self.create_df_from_binding_site(
                    interactions_by_site[selected_site],
                    interaction_type=interaction_type,
                )
                tmp_interaction["frame"] = frame
                tmp_interaction["interaction_type"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])
        else:
            interactions_by_site = self.retrieve_plip_interactions_peptide(
                f"processing_frame_{frame}.pdb"
            )
            selected_site = list(interactions_by_site.keys())[0]
            interaction_types = [
                "hydrophobic",
                "hbond",
                "waterbridge",
                "saltbridge",
                "pistacking",
                "pication",
                "halogen",
                "metal",
            ]
            for interaction_type in interaction_types:
                tmp_interaction = self.create_df_from_binding_site(
                    interactions_by_site[selected_site],
                    interaction_type=interaction_type,
                )
                tmp_interaction["frame"] = frame
                tmp_interaction["interaction_type"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])
        
        os.remove(f"processing_frame_{frame}.pdb")
        return interaction_list

    def process_frame_wrapper(self, frame: int) -> pd.DataFrame:
        """Wrapper function to handle processing of a frame, including special cases.

        Args:
            frame (int): Frame number to be processed.

        Returns:
            pd.DataFrame: DataFrame containing interaction data for the processed frame.
        """
        if self.special:
            return self.process_frame_special(frame)
        return self.process_frame(frame)

    def fill_missing_frames(self, df: pd.DataFrame) -> pd.DataFrame:
        """Fill in missing frames in the DataFrame with NaNs.

        Args:
            df (pd.DataFrame): Original DataFrame with potentially missing frames.

        Returns:
            pd.DataFrame: DataFrame with missing frames filled.
        """
        all_frames = set(range(len(self.pdb_md.trajectory)))
        existing_frames = set(df["frame"].unique())
        missing_frames = all_frames - existing_frames
        missing_data = []

        for frame in missing_frames:
            for interaction_type in df["interaction_type"].unique():
                missing_data.append({
                    "frame": frame,
                    "interaction_type": interaction_type,
                    # Add other necessary columns with NaNs or empty values as appropriate
                })

        missing_df = pd.DataFrame(missing_data)
        filled_df = pd.concat([df, missing_df], ignore_index=True)
        return filled_df

    def process_trajectory(self) -> pd.DataFrame:
        """Main method for processing the trajectory.

        Returns:
            pd.DataFrame: A DataFrame containing all interaction data for the trajectory.
        """
        pool = Pool(processes=self.num_processes)
        with tqdm(total=len(self.pdb_md.trajectory)) as pbar:
            interaction_list = pd.concat(
                list(
                    tqdm(
                        pool.imap(self.process_frame_wrapper, range(len(self.pdb_md.trajectory))),
                        total=len(self.pdb_md.trajectory),
                    )
                )
            )
        pool.close()
        pool.join()

        if self.dataframe:
            interaction_list.to_csv(self.dataframe)
        
        # Ensure all frames are represented in the final DataFrame
        interaction_list = self.fill_missing_frames(interaction_list)
        
        return interaction_list

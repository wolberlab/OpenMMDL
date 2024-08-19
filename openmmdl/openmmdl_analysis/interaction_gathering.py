import os
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from plip.basic import config
from plip.structure.preparation import PDBComplex, PLInteraction
from plip.exchange.report import BindingSiteReport
from multiprocessing import Pool
from typing import Optional, List, Dict, Tuple

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
    ):
        self.pdb_md = pdb_md
        self.dataframe = dataframe
        self.num_processes = num_processes
        self.lig_name = lig_name
        self.special = special_ligand
        self.peptide = peptide
        self.ineraction_list = self.process_trajectory()

    def characterize_complex(
        self, pdb_file: str, binding_site_id: str
    ) -> Optional[PLInteraction]:
        """Characterize the protein-ligand complex and return their interaction set.

        Args:
            pdb_file (str): Path to the PDB file.
            binding_site_id (str): Identifier of the binding site.

        Returns:
            Optional[PLInteraction]: Interaction object or None if the binding site is not found.
        """
        pdb_complex = PDBComplex()
        pdb_complex.load_pdb(pdb_file)
        for ligand in pdb_complex.ligands:
            if (
                ":".join([ligand.hetid, ligand.chain, str(ligand.position)])
                == binding_site_id
            ):
                pdb_complex.characterize_complex(ligand)
                return pdb_complex.interaction_sets.get(binding_site_id)
        return None

    def retrieve_plip_interactions(
        self, pdb_file: str, lig_name: str
    ) -> Dict[str, Dict[str, List]]:
        """Retrieves the interactions from PLIP.

        Args:
            pdb_file (str): Path of the PDB file of the complex.
            lig_name (str): Name of the ligand in the complex.

        Returns:
            Dict[str, Dict[str, List]]: Dictionary of binding sites and their interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)
        for ligand in protlig.ligands:
            if str(ligand.longname) == lig_name:
                protlig.characterize_complex(ligand)
        
        sites = {}
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
        """Retrieves the interactions from PLIP for a peptide.

        Args:
            pdb_file (str): Path of the PDB file of the complex.

        Returns:
            Dict[str, Dict[str, List]]: Dictionary of binding sites and their interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)
        protlig.characterize_complex(protlig.ligands[-1])
        
        sites = {}
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
        """Creates a DataFrame from a binding site and interaction type.

        Args:
            selected_site_interactions (Dict[str, List]): Pre-calculated interactions from PLIP.
            interaction_type (str, optional): The interaction type of interest. Defaults to "hbond".

        Returns:
            pd.DataFrame: DataFrame with information retrieved from PLIP.
        """
        valid_types = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )

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
        """Reformats the topology file to change the ligand to a residue.

        Args:
            file_path (str): Filepath of the topology file.
            new_residue_name (str): New residue name of the ligand.
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

    def process_frame(self, frame: int) -> pd.DataFrame:
        """Process a single frame of MD simulation.

        Args:
            frame (int): The number of the frame that will be processed.

        Returns:
            pd.DataFrame: DataFrame containing the interaction data for the processed frame.
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
                tmp_interaction["FRAME"] = int(frame)
                tmp_interaction["INTERACTION"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])

            if os.path.exists(f"processing_frame_{frame}.pdb"):
                os.remove(f"processing_frame_{frame}.pdb")

            if self.special is not None:
                combi_lig_special = mda.Universe("ligand_special.pdb")
                complex = mda.Universe("complex.pdb")
                complex_all = complex.select_atoms("all")
                result = self.process_frame_special(frame)
                results_df = pd.concat(result, ignore_index=True)
                results_df = results_df[results_df["LOCATION"] == "protein.sidechain"]
                results_df["RESTYPE"] = results_df["RESTYPE"].replace(
                    ["HIS", "SER", "CYS"], self.lig_name
                )
                results_df["LOCATION"] = results_df["LOCATION"].replace(
                    "protein.sidechain", "ligand"
                )
                updated_target_idx = []

                for index, row in results_df.iterrows():
                    ligand_special_int_nr = int(row["TARGET_IDX"])
                    ligand_special_int_nr_atom = combi_lig_special.select_atoms(
                        f"id {ligand_special_int_nr}"
                    )
                    for atom in ligand_special_int_nr_atom:
                        atom_name = atom.name
                        if atom_name in ["N", "C", "O", "S"]:
                            atom_name = f"{atom_name}1"
                        else:
                            base_name, atom_number = atom_name[:-1], int(atom_name[-1])
                            new_atom_number = atom_number + 1
                            atom_name = f"{base_name}{new_atom_number}"
                        for complex_atom in complex_all:
                            complex_atom_name = complex_atom.name
                            if complex_atom_name == atom_name:
                                updated_target_idx.append(row["TARGET_IDX"])
                                break

                results_df["TARGET_IDX"] = updated_target_idx
                results_df = results_df[
                    (results_df["INTERACTION"] == "hbond") & (results_df["TARGET_IDX"].notna())
                ]
                interaction_list = pd.concat([interaction_list, results_df])

        return interaction_list

    def process_trajectory(self) -> pd.DataFrame:
        """Process the entire trajectory of the MD simulation.

        Returns:
            pd.DataFrame: DataFrame containing interaction data across all frames.
        """
        total_frames = len(self.pdb_md.trajectory)
        frame_range = range(total_frames)
        with Pool(processes=self.num_processes) as pool:
            results = list(tqdm(pool.imap(self.process_frame, frame_range), total_frames))
        
        all_interactions = pd.concat(results, ignore_index=True)
        return all_interactions

    def process_frame_special(self, frame: int) -> List[pd.DataFrame]:
        """Placeholder method for processing special frame interactions.

        Args:
            frame (int): The number of the frame to process.

        Returns:
            List[pd.DataFrame]: List of DataFrames with interaction data.
        """
        # Implement specific frame processing for special cases if needed.
        return []

# Example of usage:
# pdb_md = mda.Universe('your_pdb_file.pdb')
# analyzer = InteractionAnalyzer(
#     pdb_md=pdb_md,
#     dataframe='your_dataframe.csv',
#     num_processes=4,
#     lig_name='LIG',
#     special_ligand='SPL',
#     peptide='PEP'
# )
# df = analyzer.ineraction_list
# print(df)

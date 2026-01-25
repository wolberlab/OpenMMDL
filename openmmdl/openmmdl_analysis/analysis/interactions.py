import os
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from plip.basic import config
from plip.structure.preparation import PDBComplex
from plip.exchange.report import BindingSiteReport
from multiprocessing import Pool
from typing import Dict, List, Tuple, Optional, Any

config.KEEPMOD = True


from openmmdl.openmmdl_analysis.core.utils import coord_str


class InteractionAnalyzer:
    """
    Analyzes molecular interactions between a protein and a ligand/peptide
    throughout an MD trajectory using PLIP (Protein-Ligand Interaction Profiler).

    Attributes
    ----------
    pdb_md : mda.Universe
        MDAnalysis Universe object representing the topology and trajectory.
    dataframe : str or None
        Path to an existing interaction CSV file. If None, the trajectory will be processed anew.
    num_processes : int
        Number of CPU cores to use for parallel frame analysis.
    lig_name : str
        Residue name of the ligand in the complex.
    special_ligand : str
        Residue name for special ligands like metal ions (optional).
    peptide : str
        Chain ID of the peptide ligand (optional).
    md_len : int
        Number of frames in the trajectory.
    interaction_list : pd.DataFrame
        DataFrame storing the extracted interactions across the trajectory.
    interaction_package : str
        The interaction package used for calulating the interactions
    """

    def __init__(
        self,
        pdb_md,
        dataframe,
        num_processes,
        lig_name,
        special_ligand,
        peptide,
        md_len,
        interaction_package,
    ):
        """
        Initialize an interaction analyzer and (optionally) process the trajectory.

        Depending on `interaction_package`, this initializer triggers the corresponding
        trajectory processing routine and populates `interaction_list`:

        - ``interaction_package == "plip"``: interactions are computed via PLIP on a
          per-frame PDB export.
        - ``interaction_package == "prolif"``: interactions are computed via ProLIF
          directly on the MDAnalysis trajectory selections.

        If `dataframe` is provided and `interaction_package` is ``"plip"``, interactions
        are loaded from the CSV instead of recomputing them.

        Parameters
        ----------
        pdb_md : MDAnalysis.Universe
            Universe containing topology and trajectory for the protein–ligand system.
        dataframe : str or None
            Path to a previously generated interaction CSV file. If None, interactions
            will be computed from the trajectory (subject to `interaction_package`).
        num_processes : int
            Number of processes/CPU cores used for parallel analysis.
        lig_name : str
            Residue name of the ligand (e.g., ``"LIG"``) used for selections and/or
            PLIP ligand identification.
        special_ligand : str or None
            Residue name for special ligands (e.g., metal ions) that require additional
            handling in the PLIP workflow. Not supported for ProLIF in this implementation.
        peptide : str or None
            Chain identifier for peptide ligands. If provided, the peptide chain is treated
            as the ligand selection and the remaining protein is treated as the pocket.
        md_len : int
            Number of frames in the trajectory to consider. This implementation typically
            skips frame 0 and processes frames ``1..md_len-1``.
        interaction_package : {"plip", "prolif"}
            Backend used to compute interactions.
        """
        self.pdb_md = pdb_md
        self.dataframe = dataframe
        self.num_processes = num_processes
        self.lig_name = lig_name
        self.special = special_ligand
        self.peptide = peptide
        self.md_len = md_len
        self.interaction_package = interaction_package
        if self.interaction_package == "plip":
            self.interaction_list = self._process_trajectory_plip()
        elif self.interaction_package == "prolif":
            self.interaction_list = self._process_trajectory_prolif()
        else:
            raise ValueError(f"Unknown interaction_package={self.interaction_package!r}. Expected 'plip' or 'prolif'.")

    def _retrieve_plip_interactions(self, pdb_file, lig_name):
        """
        Retrieves the interactions from PLIP.

        Parameters
        ----------
        pdb_file : str
            The path of the PDB file of the complex.
        lig_name : str
            Name of the Ligand in the complex topology that will be analyzed.

        Returns
        -------
        dict
            A dictionary of the binding sites and the interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)  # load the pdb file
        for ligand in protlig.ligands:
            if str(ligand.longname) == lig_name:
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
            # binding site.
            interactions = {
                k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info") for k in keys
            }
            sites[key] = interactions

        return sites

    def _retrieve_plip_interactions_peptide(self, pdb_file):
        """
        Retrives the interactions from PLIP for a peptide.

        Parameters
        ----------
        pdb_file : str
            The path of the PDB file of the complex.

        Returns
        -------
        dict
            A dictionary of the binding sites and the interactions.
        """
        protlig = PDBComplex()
        protlig.load_pdb(pdb_file)  # load the pdb file
        protlig.characterize_complex(protlig.ligands[-1])  # find ligands and analyze interactions
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
            # binding site.
            interactions = {
                k: [getattr(binding_site, k + "_features")] + getattr(binding_site, k + "_info") for k in keys
            }
            sites[key] = interactions

        return sites

    def _create_df_from_binding_site(self, selected_site_interactions, interaction_type="hbond"):
        """
        Creates a data frame from a binding site and interaction type.

        Parameters
        ----------
        selected_site_interactions : dict
            Precaluclated interactions from PLIP for the selected site.
        interaction_type : str, optional
            The interaction type of interest (default set to hydrogen bond). Defaults to "hbond".

        Returns
        -------
        pd.DataFrame
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
            print("\033[1m!!! Wrong interaction type specified. Hbond is chosen by default!!!\033[0m\n")
            interaction_type = "hbond"

        df = pd.DataFrame.from_records(
            # data is stored AFTER the column names
            selected_site_interactions[interaction_type][1:],
            # column names are always the first element
            columns=selected_site_interactions[interaction_type][0],
        )
        return df

    def _change_lig_to_residue(self, file_path, new_residue_name):
        """
        Reformats the topology file to change the ligand to a residue. This is needed for interactions with special ligands such as metal ions.

        Parameters
        ----------
        file_path : str
            Filepath of the topology file.
        new_residue_name : str
            New residue name of the ligand now changed to mimic an amino acid residue.

        Returns
        -------
        None
            Modifies and writes out new topology file.
        """
        with open(file_path, "r") as file:
            lines = file.readlines()

        with open(file_path, "w") as file:
            for line in lines:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    # Assuming the standard PDB format for simplicity
                    # You may need to adapt this part based on your specific PDB file
                    residue_name = line[17:20].strip()

                    # Check if the residue name matches the one to be changed
                    if residue_name == self.lig_name:
                        # Change the residue name to the new one
                        modified_line = line[:17] + new_residue_name + line[20:]
                        file.write(modified_line)
                    else:
                        file.write(line)
                else:
                    file.write(line)

    def _process_frame(self, frame):
        """
        Process a single frame of MD simulation.

        Parameters
        ----------
        frame : int
            The number of the frame that will be processed.

        Returns
        -------
        pd.DataFrame
            A dataframe conatining the interaction data for the processed frame.
        """
        atoms_selected = self.pdb_md.select_atoms(
            f"protein or nucleic or resname {self.lig_name} or (resname HOH and around 10 resname {self.lig_name}) or resname {self.special}"
        )
        for num in self.pdb_md.trajectory[(frame) : (frame + 1)]:
            atoms_selected.write(f"processing_frame_{frame}.pdb")
        if self.peptide is None:
            interactions_by_site = self._retrieve_plip_interactions(f"processing_frame_{frame}.pdb", self.lig_name)
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

            interaction_list = pd.DataFrame()
            for interaction_type in interaction_types:
                tmp_interaction = self._create_df_from_binding_site(
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
                result = self._process_frame_special(frame)
                results_df = pd.concat(result, ignore_index=True)
                results_df = results_df[results_df["LOCATION"] == "protein.sidechain"]
                results_df["RESTYPE"] = results_df["RESTYPE"].replace(["HIS", "SER", "CYS"], self.lig_name)
                results_df["LOCATION"] = results_df["LOCATION"].replace("protein.sidechain", "ligand")
                updated_target_idx = []

                for index, row in results_df.iterrows():
                    ligand_special_int_nr = int(row["TARGET_IDX"])
                    ligand_special_int_nr_atom = combi_lig_special.select_atoms(f"id {ligand_special_int_nr}")
                    for atom in ligand_special_int_nr_atom:
                        atom_name = atom.name
                        # Adjust atom_name based on the specified conditions
                        if atom_name in ["N", "C", "O", "S"]:
                            atom_name = f"{atom_name}1"
                        else:
                            # Assuming the format is a single letter followed by a number
                            base_name, atom_number = atom_name[:-1], int(atom_name[-1])
                            new_atom_number = atom_number + 1
                            atom_name = f"{base_name}{new_atom_number}"
                        for complex_atom in complex_all:
                            complex_atom_name = complex_atom.name
                            if atom_name == complex_atom_name:
                                true_number = complex_atom.id
                                break  # Exit the loop once a match is found
                        updated_target_idx.append(true_number)

                # Update 'TARGET_IDX' in interaction_list
                results_df["TARGET_IDX"] = updated_target_idx
                interaction_list["TARGET_IDX"] = interaction_list["TARGET_IDX"]

                # Concatenate the updated results_df to interaction_list
                interaction_list = pd.concat([interaction_list, results_df])
        if self.peptide is not None:
            interactions_by_site = self._retrieve_plip_interactions_peptide(f"processing_frame_{frame}.pdb")
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

            interaction_list = pd.DataFrame()
            for interaction_type in interaction_types:
                tmp_interaction = self._create_df_from_binding_site(
                    interactions_by_site[selected_site],
                    interaction_type=interaction_type,
                )
                tmp_interaction["FRAME"] = int(frame)
                tmp_interaction["INTERACTION"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])
            if os.path.exists(f"processing_frame_{frame}.pdb"):
                os.remove(f"processing_frame_{frame}.pdb")

        return interaction_list

    def _process_frame_special(self, frame):
        """
        Function extension of process_frame to process special ligands.

        Parameters
        ----------
        frame : int
            Number of the frame that will be processed.

        Returns
        -------
        list of pd.DataFrame
            List of dataframes containing the interaction data for the processed frame with the special ligand.
        """
        res_renaming = ["HIS", "SER", "CYS"]
        interaction_dfs = []
        for res in res_renaming:
            self.pdb_md.trajectory[frame]
            atoms_selected = self.pdb_md.select_atoms(f"resname {self.lig_name} or resname {self.special}")
            atoms_selected.write(f"processing_frame_{frame}.pdb")
            self._change_lig_to_residue(f"processing_frame_{frame}.pdb", res)
            interactions_by_site = self._retrieve_plip_interactions(f"processing_frame_{frame}.pdb", self.special)
            index_of_selected_site = -1
            selected_site = list(interactions_by_site.keys())[index_of_selected_site]
            interaction_types = ["metal"]
            interaction_list = pd.DataFrame()
            for interaction_type in interaction_types:
                tmp_interaction = self._create_df_from_binding_site(
                    interactions_by_site[selected_site],
                    interaction_type=interaction_type,
                )
                tmp_interaction["FRAME"] = int(frame)
                tmp_interaction["INTERACTION"] = interaction_type
                interaction_list = pd.concat([interaction_list, tmp_interaction])
            interaction_dfs.append(interaction_list)
            os.remove(f"processing_frame_{frame}.pdb")

        return interaction_dfs

    def _process_frame_wrapper(self, args):
        """
        Wrapper for the MD Trajectory procession.

        Parameters
        ----------
        args : tuple
            Tuple containing (frame_idx: int - number of the frame to be processed).

        Returns
        -------
        tuple
            Tuple containing the frame index and the result of from the process_frame function.
        """
        frame_idx, pdb_md, lig_name, special_ligand, peptide = args

        return frame_idx, self._process_frame(frame_idx)

    def _fill_missing_frames(self, df):
        """
        Fills the frames with no interactions in the DataFrame with placeholder values.

        Parameters
        ----------
        df : pd.DataFrame
            The input DataFrame with frames that have no interactions.

        Returns
        -------
        pd.DataFrame
            DataFrame with placeholder values in the frames with no interactions.
        """

        # Create a set containing all unique values in the 'FRAME' column
        existing_frames = set(df["FRAME"])

        # Create a list to store new rows for missing numbers
        missing_rows = []

        # Iterate through numbers from 0 to md_len
        for frame_number in range(1, self.md_len):
            if frame_number not in existing_frames:
                # Create a new row with 'FRAME' set to the missing number and other columns set to "skip"
                missing_row = {"FRAME": frame_number}
                for col in df.columns:
                    if col != "FRAME":
                        missing_row[col] = "skip"
                missing_rows.append(missing_row)

        # Concatenate the missing rows with the original DataFrame
        df = pd.concat([df, pd.DataFrame(missing_rows)], ignore_index=True)

        # Sort the DataFrame by the 'FRAME' column
        df.sort_values(by="FRAME", inplace=True)

        return df

    def _process_trajectory_plip(self):
        """
        Process protein-ligand trajectory with multiple CPUs in parallel.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
        """
        if self.dataframe is None:
            print("\033[1mProcessing protein-ligand trajectory\033[0m")
            print(f"\033[1mUsing {self.num_processes} CPUs\033[0m")

            with Pool(processes=self.num_processes) as pool:
                frame_args = [
                    (i, self.pdb_md, self.lig_name, self.special, self.peptide) for i in range(1, self.md_len)
                ]

                # Initialize the progress bar with the total number of frames
                pbar = tqdm(
                    total=self.md_len - 1,
                    ascii=True,
                    desc="\033[1mAnalyzing frames\033[0m",
                )

                results = []
                for result in pool.imap(self._process_frame_wrapper, frame_args):
                    results.append(result)
                    pbar.update(1)  # Update the progress manually

            # Close the progress bar
            pbar.close()

            # Extract the results and sort them by frame index
            results.sort(key=lambda x: x[0])
            interaction_lists = [result[1] for result in results]

            interaction_list = pd.concat(interaction_lists)

            interaction_list.to_csv("interactions_gathered.csv")

        elif self.dataframe is not None:
            print(f"\033[1mGathering data from {self.dataframe}\033[0m")
            interaction_tmp = pd.read_csv(self.dataframe)
            interaction_list = interaction_tmp.drop(interaction_tmp.columns[0], axis=1)

        interaction_list["Prot_partner"] = (
            interaction_list["RESNR"].astype(str) + interaction_list["RESTYPE"] + interaction_list["RESCHAIN"]
        )

        interaction_list = self._fill_missing_frames(
            interaction_list,
        )

        print("\033[1mProtein-ligand trajectory processed\033[0m")

        return interaction_list

    # -------------------------
    # ProLIF helpers (NEW)
    # -------------------------
    _PROLIF_BASE_COLUMNS = [
        "FRAME",
        "INTERACTION",
        "RESTYPE",
        "RESNR",
        "RESCHAIN",
        "RESTYPE_LIG",
        "RESNR_LIG",
        "PROTCOO",
        "LIGCOO",
        "LOCATION",
        "PROTISDON",
        "DONORIDX",
        "ACCEPTORIDX",
        "DONOR_IDX",
        "ACCEPTOR_IDX",
        "WATER_IDX",
        "PROTISPOS",
        "LIGCARBONIDX",
        "LIG_GROUP",
        "LIG_IDX_LIST",
        "PROT_IDX_LIST",
        "DON_IDX",
        "DONORTYPE",
        "TARGET_IDX",
        "TARGETCOO",
        "METAL_TYPE",
        "COORDINATION",
    ]

    @staticmethod
    def _element_upper(atom) -> str:
        """
        Return an atom element symbol in uppercase, with a name-based fallback.

        Parameters
        ----------
        atom : object
            Atom-like object with optional attributes `element` and `name`.

        Returns
        -------
        str
            Uppercase element symbol. If `atom.element` is missing/empty, the first
            character of `atom.name` is used; if neither is available, returns ``"X"``.
        """
        el = getattr(atom, "element", None)
        if el and str(el).strip():
            return str(el).strip().upper()
        name = getattr(atom, "name", "") or ""
        return name[0].upper() if name else "X"

    @classmethod
    def _pick_point_atom(cls, atoms, prefer_elements=None):
        """
        Select a representative atom for coordinate and index fields.

        The selection rules are applied in order:
        1. If `prefer_elements` is provided, return the first non-hydrogen atom whose
           element is in that set.
        2. Otherwise return the first non-hydrogen (heavy) atom.
        3. If only hydrogens exist, return the first atom.

        Parameters
        ----------
        atoms : MDAnalysis.core.groups.AtomGroup or sequence
            Collection of atoms from which to select a representative.
        prefer_elements : sequence of str or None, optional
            Element symbols to prefer (case-insensitive), e.g. ``{"F", "Cl", "Br", "I"}``.

        Returns
        -------
        atom : object or None
            A selected atom-like object, or None if `atoms` is empty.

        Notes
        -----
        This is primarily used to populate single-atom fields (e.g., `LIGCOO`, `PROTCOO`,
        donor/target indices) when an interaction involves a group of atoms.
        """
        if atoms is None or len(atoms) == 0:
            return None
        prefer_elements = set(e.upper() for e in (prefer_elements or []))

        # prefer specific elements first (non-H)
        if prefer_elements:
            for a in atoms:
                el = cls._element_upper(a)
                if el in prefer_elements and el != "H":
                    return a

        # otherwise first heavy atom
        for a in atoms:
            if cls._element_upper(a) != "H":
                return a

        # fallback
        return atoms[0]

    def _atoms_from_indices(self, idxs, selection: mda.AtomGroup) -> mda.AtomGroup:
        """
        Resolve ProLIF-provided indices into an MDAnalysis AtomGroup.

        ProLIF may provide indices that refer either to:
        - Universe-level atom indices (0-based), or
        - selection-local indices relative to the AtomGroup passed to ProLIF.

        This method attempts Universe-level indexing first and validates that the result
        is a subset of `selection`. If that fails, it falls back to selection-local indexing.

        Parameters
        ----------
        idxs : sequence of int or None
            Indices to resolve. If empty/None, an empty AtomGroup is returned.
        selection : MDAnalysis.core.groups.AtomGroup
            AtomGroup representing the intended endpoint domain (ligand or protein selection).

        Returns
        -------
        MDAnalysis.core.groups.AtomGroup
            Resolved AtomGroup. Returns an empty AtomGroup when indices cannot be interpreted.
        """
        if not idxs:
            return self.pdb_md.atoms[[]]

        try:
            idx_list = [int(i) for i in idxs]
        except Exception:
            return self.pdb_md.atoms[[]]

        # Try Universe indices first
        try:
            atoms_u = self.pdb_md.atoms[idx_list]
            # sanity: should be subset of selection for ligand/protein endpoints
            if set(atoms_u.indices).issubset(set(selection.indices)):
                return atoms_u
        except Exception:
            pass

        # Fallback: selection-local indexing
        try:
            return selection[idx_list]
        except Exception:
            return self.pdb_md.atoms[[]]

    def _ensure_topology_for_prolif(
        self, ligand_ag: mda.core.groups.AtomGroup, protein_ag: mda.core.groups.AtomGroup
    ) -> None:
        """
        Best-effort preparation of topology information required by ProLIF.

        ProLIF typically relies on MDAnalysis-to-RDKit conversion, which may require
        element information and bond topology. PDB-based or stripped topologies may
        lack these fields, leading to missing or empty interaction results.

        This method attempts to:
        - ensure `elements` are present on the Universe (best-effort guessing), and
        - guess bonds on the ligand/protein selections (and/or Universe), when absent.

        Parameters
        ----------
        ligand_ag : MDAnalysis.core.groups.AtomGroup
            Ligand AtomGroup passed to ProLIF.
        protein_ag : MDAnalysis.core.groups.AtomGroup
            Protein/pocket AtomGroup passed to ProLIF.

        Returns
        -------
        None
        """
        # Ensure elements are present (required by bond guessing and RDKit conversion in many cases)
        try:
            elems = getattr(self.pdb_md.atoms, "elements", None)
            if elems is None or any((e is None or str(e).strip() == "") for e in elems):
                from MDAnalysis.topology.guessers import guess_types

                self.pdb_md.add_TopologyAttr("elements", guess_types(self.pdb_md.atoms.names))
        except Exception:
            # Don't hard-fail: some topologies provide elements as a non-iterable proxy
            pass

        # Guess bonds on the specific selections if missing
        for ag, label in ((ligand_ag, "ligand"), (protein_ag, "protein")):
            try:
                # AtomGroup.bonds is present when Universe has bond topology
                if not hasattr(ag, "bonds") or len(ag.bonds) == 0:
                    ag.guess_bonds()
            except Exception as e:
                # Fallback: try Universe-level guessing (MDAnalysis versions differ)
                try:
                    self.pdb_md.guess_TopologyAttrs(to_guess=["bonds"], force_guess=["bonds"])
                except Exception:
                    print(
                        f"Warning: could not guess bonds for {label} selection ({e}). ProLIF may produce no interactions."
                    )

        try:
            self.pdb_md.atoms.guess_bonds()
        except Exception:
            pass

    def _infer_water_resnames(self) -> List[str]:
        """
        Infer candidate water residue names from the current Universe.

        This helper scans residues and identifies small residues consistent with water-like
        composition (presence of oxygen and optionally hydrogens). It is used as a fallback
        when MDAnalysis' built-in `water` selection returns no atoms.

        Returns
        -------
        list of str
            Sorted list of residue names that likely correspond to water in this system.

        Notes
        -----
        This is heuristic and may include false positives in unusual topologies. It is
        intended only as a fallback to make WaterBridge selection more robust.
        """
        candidates = set()
        for res in self.pdb_md.residues:
            try:
                if res.atoms.n_atoms > 6:
                    continue
                elems = [self._element_upper(a) for a in res.atoms]
                if "O" in elems and ("H" in elems or elems.count("O") == 1):
                    candidates.add(str(res.resname))
            except Exception:
                continue
        return sorted(candidates)

    def _select_water_ag(self, ligand_ag, pocket_ag, order: int) -> Optional[mda.AtomGroup]:
        """
        Select water atoms near the ligand/pocket for ProLIF WaterBridge interactions.

        Parameters
        ----------
        ligand_ag : MDAnalysis.core.groups.AtomGroup
            Ligand AtomGroup used as one endpoint of the proximity selection.
        pocket_ag : MDAnalysis.core.groups.AtomGroup
            Protein/pocket AtomGroup used as the other endpoint of the proximity selection.
        order : int
            WaterBridge order (number of water molecules in the bridge). Higher order
            typically requires a larger spatial cutoff.

        Returns
        -------
        MDAnalysis.core.groups.AtomGroup or None
            Water AtomGroup suitable for passing to ProLIF WaterBridge parameters, or
            None if no water atoms could be selected.
        """
        # ProLIF tutorial uses ~8 Å for order 1 and recommends increasing for higher order. :contentReference[oaicite:4]{index=4}
        water_cutoff = 8.0 + 2.0 * max(0, int(order) - 1)

        # 1) try MDAnalysis 'water' selection keyword
        sel = f"water and byres around {water_cutoff} (group ligand or group pocket)"
        try:
            wat = self.pdb_md.select_atoms(sel, ligand=ligand_ag, pocket=pocket_ag, updating=True)
            if len(wat) > 0:
                return wat
        except Exception:
            pass

        # 2) fallback: infer resnames
        resnames = self._infer_water_resnames()
        if not resnames:
            return None
        sel = f"resname {' '.join(resnames)} and byres around {water_cutoff} (group ligand or group pocket)"
        try:
            wat = self.pdb_md.select_atoms(sel, ligand=ligand_ag, pocket=pocket_ag, updating=True)
            if len(wat) > 0:
                return wat
        except Exception:
            pass
        try:
            print("Debug: unique resnames (first 30):", sorted(set(self.pdb_md.atoms.resnames))[:30])
        except Exception:
            pass
        return None

    def _residue_fields(self, residue_obj) -> Tuple[str, str, str]:
        """
        Extract residue identity fields from a ProLIF residue identifier.

        Parameters
        ----------
        residue_obj : object or None
            Residue identifier object provided by ProLIF (commonly a ResidueId-like object).
            If None, placeholder values are returned.

        Returns
        -------
        RESTYPE : str
            Residue name/type (e.g., ``"ASP"``) or ``"skip"`` if not available.
        RESNR : str
            Residue number/index (e.g., ``"145"``) or ``"skip"`` if not available.
        RESCHAIN : str
            Chain identifier (e.g., ``"A"``) or ``"skip"`` if not available.

        Notes
        -----
        This function is defensive against ProLIF version differences by trying multiple
        attribute names for each field (e.g., `name`/`resname`, `number`/`resid`, `chain`/`segid`).
        """
        if residue_obj is None:
            return "skip", "skip", "skip"

        # ProLIF ResidueId typically has name/number/chain
        for name_attr in ("name", "resname", "resn"):
            if hasattr(residue_obj, name_attr):
                restype = getattr(residue_obj, name_attr)
                break
        else:
            restype = "skip"

        for num_attr in ("number", "resid", "resnum", "resnr"):
            if hasattr(residue_obj, num_attr):
                resnr = getattr(residue_obj, num_attr)
                break
        else:
            resnr = "skip"

        for chain_attr in ("chain", "chainID", "segid"):
            if hasattr(residue_obj, chain_attr):
                reschain = getattr(residue_obj, chain_attr)
                break
        else:
            reschain = "skip"

        return str(restype), str(resnr), str(reschain)

    def _process_trajectory_prolif(self) -> pd.DataFrame:
        """
        Compute protein–ligand interactions over a trajectory using ProLIF.

        This method:
        1. Builds ligand and protein/pocket selections from `self.pdb_md`.
        2. Performs best-effort topology preparation required for ProLIF (elements/bonds).
        3. Runs a ProLIF fingerprint over frames ``1..self.md_len-1`` in parallel.
        4. Converts ProLIF interaction outputs into the OpenMMDL interaction row schema.

        Returns
        -------
        pandas.DataFrame
            Interaction table with columns defined by `_PROLIF_BASE_COLUMNS` plus a
            ``Prot_partner`` column constructed as ``RESNR + RESTYPE + RESCHAIN``.
            Frames without interactions are filled using `_fill_missing_frames`.

        Notes
        -----
        - WaterBridge interactions require an explicit water AtomGroup. If no water atoms are
          selected, WaterBridge is disabled and a warning is emitted.
        - Interaction names emitted by ProLIF are mapped to OpenMMDL interaction labels
          (e.g., ``"HBDonor"``/``"HBAcceptor"`` → ``"hbond"``).
        - Atom index fields are generally stored as `Atom.id` (often PDB serial IDs) to match
          OpenMMDL downstream expectations.
        """
        try:
            import prolif as plf
        except Exception as e:
            raise ImportError(
                "interaction_engine='prolif' requested but ProLIF is not installed. "
                "Install it with: pip install prolif==2.1.0"
            ) from e

        if self.special is not None:
            raise NotImplementedError(
                "interaction_engine='prolif' is currently not supported together with special_ligand. "
                "Use interaction_engine='plip' for special_ligand."
            )

        # ---- selections ----
        if self.peptide is None:
            lig_sel = f"resname {self.lig_name}"
            ligand_ag = self.pdb_md.select_atoms(lig_sel)

            # ProLIF tutorial recommends selecting a pocket around the ligand (no updating=True). :contentReference[oaicite:8]{index=8}
            base = "protein or nucleic" if getattr(config, "DNARECEPTOR", False) else "protein"
            protein_ag = self.pdb_md.select_atoms(
                f"({base}) and byres around 12 group ligand",
                ligand=ligand_ag,
            )
        else:
            # peptide treated as ligand chain; pocket = rest of protein
            ligand_ag = self.pdb_md.select_atoms(f"chainID {self.peptide}")
            base = "protein or nucleic" if getattr(config, "DNARECEPTOR", False) else "protein"
            protein_ag = self.pdb_md.select_atoms(f"({base}) and not chainID {self.peptide}")

        if len(ligand_ag) == 0:
            raise ValueError(f"ProLIF: ligand selection returned 0 atoms: {lig_sel}")
        if len(protein_ag) == 0:
            # fallback to full protein if pocket selection is empty
            base = "protein or nucleic" if getattr(config, "DNARECEPTOR", False) else "protein"
            protein_ag = self.pdb_md.select_atoms(base)

        self._ensure_topology_for_prolif(ligand_ag, protein_ag)

        # ---- interactions ----
        prolif_interactions = [
            "Hydrophobic",
            "HBDonor",
            "HBAcceptor",
            "PiStacking",
            "CationPi",
            "PiCation",
            "Cationic",
            "Anionic",
            "XBDonor",
            "XBAcceptor",
            "WaterBridge",
        ]

        # WaterBridge requires explicit water parameter, and order=3 is supported as per docs. :contentReference[oaicite:9]{index=9}
        WATER_BRIDGE_ORDER = 3
        water_ag = self._select_water_ag(ligand_ag, protein_ag, WATER_BRIDGE_ORDER)
        parameters = {}
        if water_ag is not None and len(water_ag) > 0:
            parameters["WaterBridge"] = {"water": water_ag, "order": WATER_BRIDGE_ORDER}
        else:
            print("\033[33mWarning: WaterBridge disabled (no water atoms selected).\033[0m")
            prolif_interactions = [x for x in prolif_interactions if x != "WaterBridge"]

        fp = plf.Fingerprint(prolif_interactions, count=True, parameters=(parameters or None))

        # OpenMMDL historically skips frame 0
        traj = self.pdb_md.trajectory[1 : self.md_len]
        fp.run(traj, ligand_ag, protein_ag, n_jobs=self.num_processes, progress=True)

        # ---- convert to OpenMMDL schema ----
        rows: List[Dict[str, Any]] = []

        # robust extraction: iterate InteractionData objects as in ProLIF docs. :contentReference[oaicite:10]{index=10}
        for frame_key, ifp in (getattr(fp, "ifp", {}) or {}).items():
            try:
                frame = int(frame_key)
            except Exception:
                continue
            if frame < 1 or frame >= self.md_len:
                continue

            # ensure coordinates correspond to this frame
            try:
                self.pdb_md.trajectory[frame]
            except Exception:
                pass

            for interaction_data in ifp.interactions():
                meta = interaction_data.metadata or {}
                idx_block = meta.get("parent_indices") or meta.get("indices") or {}

                lig_atoms = self._atoms_from_indices(idx_block.get("ligand"), ligand_ag)
                prot_atoms = self._atoms_from_indices(idx_block.get("protein"), protein_ag)

                lig_point = self._pick_point_atom(lig_atoms)
                prot_point = self._pick_point_atom(prot_atoms)

                restype, resnr, reschain = self._residue_fields(getattr(interaction_data, "protein", None))
                lig_restype, lig_resnr, _lig_chain = self._residue_fields(getattr(interaction_data, "ligand", None))

                row = {c: "skip" for c in self._PROLIF_BASE_COLUMNS}
                row["TARGETCOO"] = 0
                row.update(
                    {
                        "FRAME": frame,
                        "RESTYPE": restype,
                        "RESNR": int(resnr),
                        "RESCHAIN": str(reschain),
                        "RESTYPE_LIG": str(lig_restype),
                        "RESNR_LIG": str(lig_resnr),
                        "LIGCOO": coord_str(lig_point.position if lig_point is not None else None),
                        "PROTCOO": coord_str(prot_point.position if prot_point is not None else None),
                    }
                )

                name = str(interaction_data.interaction)

                # indices as PDB serial IDs (Atom.id), consistent with OpenMMDL downstream int() casts
                lig_ids = [int(a.id) for a in lig_atoms] if len(lig_atoms) else []
                prot_ids = [int(a.id) for a in prot_atoms] if len(prot_atoms) else []

                if name == "Hydrophobic":
                    row["INTERACTION"] = "hydrophobic"
                    if lig_point is not None:
                        row["LIGCARBONIDX"] = int(lig_point.id)

                elif name in ("HBDonor", "HBAcceptor"):
                    row["INTERACTION"] = "hbond"
                    if name == "HBAcceptor":
                        # ligand acceptor, protein donor
                        row["PROTISDON"] = True
                        if prot_point is not None:
                            row["DONORIDX"] = int(prot_point.id)
                        if lig_point is not None:
                            row["ACCEPTORIDX"] = int(lig_point.id)
                    else:
                        # ligand donor, protein acceptor
                        row["PROTISDON"] = False
                        if lig_point is not None:
                            row["DONORIDX"] = int(lig_point.id)
                        if prot_point is not None:
                            row["ACCEPTORIDX"] = int(prot_point.id)

                elif name == "PiStacking":
                    row["INTERACTION"] = "pistacking"
                    if lig_ids:
                        row["LIG_IDX_LIST"] = ",".join(map(str, lig_ids))
                    if prot_ids:
                        row["PROT_IDX_LIST"] = ",".join(map(str, prot_ids))

                elif name in ("CationPi", "PiCation"):
                    row["INTERACTION"] = "pication"
                    row["LIG_GROUP"] = "cation" if name == "CationPi" else "pi"
                    if lig_ids:
                        row["LIG_IDX_LIST"] = ",".join(map(str, lig_ids))
                    if prot_ids:
                        row["PROT_IDX_LIST"] = ",".join(map(str, prot_ids))

                elif name in ("Cationic", "Anionic"):
                    row["INTERACTION"] = "saltbridge"
                    row["PROTISPOS"] = True if name == "Anionic" else False
                    row["LIG_GROUP"] = "anion" if name == "Anionic" else "cation"
                    if lig_ids:
                        row["LIG_IDX_LIST"] = ",".join(map(str, lig_ids))
                    if prot_ids:
                        row["PROT_IDX_LIST"] = ",".join(map(str, prot_ids))

                elif name in ("XBDonor", "XBAcceptor"):
                    row["INTERACTION"] = "halogen"
                    # prefer actual halogen element for donor
                    hal = self._pick_point_atom(lig_atoms, prefer_elements={"F", "CL", "BR", "I"})
                    if hal is not None:
                        row["DON_IDX"] = int(hal.id)
                        row["DONORTYPE"] = self._element_upper(hal)
                    if prot_point is not None:
                        row["TARGET_IDX"] = int(prot_point.id)

                elif name == "WaterBridge":
                    row["INTERACTION"] = "waterbridge"

                    # metadata contains roles + water residue identifiers :contentReference[oaicite:11]{index=11}
                    protein_role = str(meta.get("protein_role", ""))
                    ligand_role = str(meta.get("ligand_role", ""))

                    # Map to the boolean expected by BindingModeProcesser (it assumes exactly one endpoint "donor"):
                    if protein_role == "HBDonor":
                        row["PROTISDON"] = True
                    elif ligand_role == "HBDonor":
                        row["PROTISDON"] = False
                    else:
                        # e.g. both acceptors: water is donor; choose acceptor labelling (PROTISDON=True)
                        row["PROTISDON"] = True

                    water_res = meta.get("water_residues") or ()
                    if water_res:
                        try:
                            row["WATER_IDX"] = int(getattr(water_res[0], "number", "skip"))
                        except Exception:
                            row["WATER_IDX"] = "skip"

                    # Endpoint indices for OpenMMDL naming
                    if row["PROTISDON"]:
                        if prot_point is not None:
                            row["DONOR_IDX"] = int(prot_point.id)
                        if lig_point is not None:
                            row["ACCEPTOR_IDX"] = int(lig_point.id)
                    else:
                        if lig_point is not None:
                            row["DONOR_IDX"] = int(lig_point.id)
                        if prot_point is not None:
                            row["ACCEPTOR_IDX"] = int(prot_point.id)

                else:
                    # not mapped → ignore
                    continue

                rows.append(row)

        # Always return stable schema
        interaction_list = pd.DataFrame.from_records(rows, columns=self._PROLIF_BASE_COLUMNS)

        # Ensure Prot_partner exists (BindingModeProcesser expects it in residue schema)
        # Mirror PLIP path
        if not interaction_list.empty:
            interaction_list["Prot_partner"] = (
                interaction_list["RESNR"].astype(str)
                + interaction_list["RESTYPE"].astype(str)
                + interaction_list["RESCHAIN"].astype(str)
            )
        else:
            # create the column even if empty
            interaction_list["Prot_partner"] = pd.Series(dtype="object")

        # Fill missing frames (will also fill Prot_partner with 'skip' for missing frames)
        interaction_list = self._fill_missing_frames(interaction_list)

        # keep OpenMMDL behavior: write gathered interactions
        interaction_list.to_csv("interactions_gathered.csv", index=False)
        return interaction_list

import os
import pandas as pd
import MDAnalysis as mda
from tqdm import tqdm
from plip.structure.preparation import PDBComplex, LigandFinder, Mol, PLInteraction
from plip.exchange.report import BindingSiteReport
from multiprocessing import Pool
from functools import partial


def characterize_complex(pdb_file: str, binding_site_id: str) -> PLInteraction:
    """
    Characterize the protein-ligand complex and return their interaction set

    Parameters
    ----------
    pdb_file : str
        A string, which represents the path to the PDB File
    binding_site_id : str
        A string that specifies the identifier of the binding site

    Returns
    -------
    PLInteraction :
        A object representing the interactions if. If Binding site is not found returns None
    """
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
        # binding site.
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


def process_frame(frame, pdb_md):
    """
    Process a single frame of MD simulation.

    Parameters
    ----------
    frame : int
        The number of the frame that is going to be processed.
    pdb_md : mda.Universe
        The MDAnalysis Universe class representation of the topology and the trajectory of the file that is being processed.

    Returns
    -------
    pd.DataFrame :
        A dataframe conatining the interaction data for the processed frame.
    """
    atoms_selected = pdb_md.select_atoms("protein or resname UNK or (resname HOH and around 10 resname UNK)")
    for num in pdb_md.trajectory[(frame):(frame+1)]:
        atoms_selected.write(f'{frame}.pdb')
    interactions_by_site = retrieve_plip_interactions(f"{frame}.pdb")
    index_of_selected_site = -1
    selected_site = list(interactions_by_site.keys())[index_of_selected_site]

    interaction_types = ["hydrophobic", "hbond", "waterbridge", "saltbridge", "pistacking", "pication", "halogen", "metal"]

    interaction_list = pd.DataFrame()
    for interaction_type in interaction_types:
        tmp_interaction = create_df_from_binding_site(interactions_by_site[selected_site], interaction_type=interaction_type)
        tmp_interaction['FRAME'] = int(frame)
        tmp_interaction['INTERACTION'] = interaction_type
        interaction_list = pd.concat([interaction_list, tmp_interaction])

    os.remove(f"{frame}.pdb")

    return interaction_list


def process_frame_wrapper(args):
    """
    Wrapper for the MD Trajectory procession.

    Parameters
    ----------
    args : tuple
        - frame_idx : int
            Integer representing the index of the processing frame.
        - pdb_md : mda.Universe
            The MDAnalysis Universe class representation of the topology and the trajectory of the file that is being processed.
        
    Returns
    -------
    tuple :
        tuple containing the frame index and the result of from `process_frame(frame_idx, pdb_md)`.
    """
    frame_idx, pdb_md = args

    return frame_idx, process_frame(frame_idx, pdb_md)


def process_trajectory(pdb_md, dataframe, num_processes=4):
    """
    Process protein-ligand trajectory with multiple CPUs in parallel.

    Parameters
    ----------
    pdb_md : mda.Universe
        MDAnalysis Universe object representing the protein-ligand topology and trajectory.
    dataframe : str
        Name of a CSV file as str, where the interaction data will be read from if not None.
    num_processes : int (optional)
        The number of CPUs that will be used for the processing of the protein-ligand trajectory
        
    Returns
    -------
    pd.DataFrame :
        A DataFrame containing all the protein-ligand interaction data from the whole trajectory.
    """
    if dataframe is None:
        print("\033[1mProcessing protein-ligand trajectory\033[0m")
        print(f"\033[1mUsing {num_processes} CPUs\033[0m")
        total_frames = len(pdb_md.trajectory) - 1

        with Pool(processes=num_processes) as pool:
            frame_args = [(i, pdb_md) for i in range(1, total_frames + 1)]
            
            # Initialize the progress bar with the total number of frames
            pbar = tqdm(total=total_frames, ascii=True, desc="Analyzing frames")
            
            results = []
            for result in pool.imap(process_frame_wrapper, frame_args):
                results.append(result)
                pbar.update(1)  # Update the progress manually

        # Close the progress bar
        pbar.close()

        # Extract the results and sort them by frame index
        results.sort(key=lambda x: x[0])
        interaction_lists = [result[1] for result in results]

        interaction_list = pd.concat(interaction_lists)

        interaction_list.to_csv("interactions_gathered.csv")
    elif dataframe is not None:
        print(f"\033[1mGathering data from {dataframe}\033[0m")
        interaction_tmp = pd.read_csv(dataframe)
        interaction_list = interaction_tmp.drop(interaction_tmp.columns[0], axis=1)

    print("\033[1mProtein-ligand trajectory processed\033[0m")
    
    return interaction_list

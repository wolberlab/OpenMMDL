import pytest
import subprocess
import os
import openmmdl

from pathlib import Path

from openmmdl.openmmdl_analysis.preprocessing import process_pdb_file, convert_pdb_to_sdf
from openmmdl.openmmdl_analysis.rmsd_calculation import rmsd_for_atomgroups, RMSD_dist_frames
from openmmdl.openmmdl_analysis.ligand_processing import increase_ring_indices, convert_ligand_to_smiles
from openmmdl.openmmdl_analysis.interaction_gathering import characterize_complex, retrieve_plip_interactions, create_df_from_binding_site, process_frame, process_trajectory, fill_missing_frames
from openmmdl.openmmdl_analysis.binding_mode_processing import gather_interactions, remove_duplicate_values, combine_subdict_values, filtering_values, unique_data_generation, df_iteration_numbering, update_values
from openmmdl.openmmdl_analysis.markov_state_figure_generation import min_transition_calculation, binding_site_markov_network
from openmmdl.openmmdl_analysis.rdkit_figure_generation import split_interaction_data, highlight_numbers, generate_interaction_dict, update_dict, create_and_merge_images, arranged_figure_generation
from openmmdl.openmmdl_analysis.barcode_generation import barcodegeneration,plot_barcodes,plot_waterbridge_piechart
from openmmdl.openmmdl_analysis.visualization_functions import interacting_water_ids, save_interacting_waters_trajectory, cloud_json_generation
from openmmdl.openmmdl_analysis.pml_writer import generate_md_pharmacophore_cloudcenters, generate_bindingmode_pharmacophore, generate_pharmacophore_centers_all_points, generate_point_cloud_pml

# Print current working directory
print("Current working directory:", os.getcwd())

# Print the full path to the input file
input_pdb_filename = "openmmdl/tests/data/in/0_unk_hoh.pdb"
print("Full path to input file:", os.path.abspath(input_pdb_filename))

test_data_directory = Path("openmmdl/tests/data/in")

@pytest.fixture(scope="session")
def test_data_dir(tmp_path_factory):
    data_dir = tmp_path_factory.mktemp("test_data")
    return data_dir


def test_script_execution(test_data_dir):
    # Define the root directory
    root_dir = "/home/runner/work/OpenMMDL/OpenMMDL"

    # Specify the relative path to the input data
    relative_input_path = "openmmdl/tests/data/in"

    # Combine the root directory and relative input path
    input_data_path = os.path.join(root_dir, relative_input_path)

    # Ensure that the script runs successfully without errors
    script_path = "openmmdlanalysis.py"
    topology_file = os.path.join(input_data_path, "0_unk_hoh.pdb")
    trajectory_file = os.path.join(root_dir, "openmmdl/tests/data/in/all_50.dcd")
    ligand_sdf_file = os.path.join(input_data_path, "lig.sdf")
    ligand_name = "UNK"
    
    cmd = f" openmmdl_analysis -t {topology_file} -d {trajectory_file} -l {ligand_sdf_file} -n {ligand_name} -b 10 -c 2"
    
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=test_data_dir)
    
    assert result.returncode == 0, f"Script execution failed with error:\n{result.stderr.decode()}"

    # Check that expected output files are generated
    assert os.path.exists(os.path.join(test_data_dir, "complex.pdb"))
    assert os.path.exists(os.path.join(test_data_dir, "lig.pdb"))
    assert os.path.exists(os.path.join(test_data_dir, "df_all.csv"))

    # Check for the presence of "Barcodes" folder
    barcodes_folder = os.path.join(test_data_dir, "Barcodes")
    assert os.path.exists(barcodes_folder), "The 'Barcodes' folder is missing."

    # Check for the existence of "hydrophobic_barcodes.png" inside "Barcodes"
    hydro_file_path = os.path.join(barcodes_folder, "hydrophobic_barcodes.png")
    assert os.path.exists(hydro_file_path), "The 'hydrophobic_barcodes.png' file is missing inside 'Barcodes'."


    # Check for the presence of "Binding_Modes_Markov_States" folder
    markov_states_folder = os.path.join(test_data_dir, "Binding_Modes_Markov_States")
    assert os.path.exists(markov_states_folder), "The 'Binding_Modes_Markov_States' folder is missing."

    # Check for the existence of "all_binding_modes_arranged.png" inside "Binding_Modes_Markov_States"
    png_file_path = os.path.join(markov_states_folder, "all_binding_modes_arranged.png")
    assert os.path.exists(png_file_path), "The 'all_binding_modes_arranged.png' file is missing inside 'Binding_Modes_Markov_States'."

    # Add more checks for other output files as needed

import os
import shutil
from typing import List


def cleanup(protein_name):
    """Cleans up the PDB Reporter Output File and MDTraj Files of the performed simulation.

    Args:
        protein_name (str): Name of the protein PDB.

    Returns:
        None
    """
    print("Cleaning Up :)")
    try:
        os.remove(f"output_{protein_name}")
        os.remove("centered_old_coordinates.pdb")
        os.remove("centered_old_coordinates.dcd")
    except FileNotFoundError:
        print("One or more files not found. Cleanup skipped.")
    print("Cleanup is done.")


def create_directory_if_not_exists(directory_path):
    """Create a directory if it doesn't exist, or overwrite it if already does.

    Args:
        directory_path (str): Path of the directory that you want to create.

    Returns:
        None
    """
    if not os.path.exists(directory_path):
        os.mkdir(directory_path)
    else:
        shutil.rmtree(directory_path)
        os.mkdir(directory_path)


def copy_file(src, dest):
    """Copy a file to the destination path.

    Args:
        src (str): Path of the file that needs to be copied.
        dest (str): Path of destination where the file needs to be copied to.

    Returns:
        None
    """
    if os.path.exists(src):
        shutil.copy(src, dest)


def organize_files(source, destination):
    """Organizes the files and moves them from the source to the destination directory.

    Args:
        source (str): Path of the file that needs to be moved.
        destination (str): Path of destination where the file needs to be moved to.

    Returns:
        None
    """
    for file in source:
        if os.path.exists(file):
            os.rename(file, os.path.join(destination, os.path.basename(file)))


def post_md_file_movement(protein_name: str, prmtop: str = None, inpcrd: str = None, ligands: List[str] = None):
    """Organizes and moves the files after the MD simulation to their respective directories.

    Args:
        protein_name (str): Name of the protein PDB.
        prmtop (str, optional): Path to the AMBER topology file.
        inpcrd (str, optional): Path to the AMBER coordinate file.
        ligands (List[str], optional): List of paths to the ligand files.

    Returns:
        None
    """
    # Create necessary directories
    create_directory_if_not_exists("Input_Files")
    create_directory_if_not_exists("MD_Files")
    create_directory_if_not_exists("MD_Files/Pre_MD")
    create_directory_if_not_exists("MD_Files/Minimization_Equilibration")
    create_directory_if_not_exists("MD_Files/MD_Output")
    create_directory_if_not_exists("MD_Postprocessing")
    create_directory_if_not_exists("Final_Output")
    create_directory_if_not_exists("Final_Output/All_Atoms")
    create_directory_if_not_exists("Final_Output/Prot_Lig")
    create_directory_if_not_exists("Checkpoints")
    # Move input files
    if ligands:
        for lig in ligands:
            copy_file(lig, "Input_Files")
            copy_file(lig, "Final_Output/All_Atoms")
            copy_file(lig, "Final_Output/Prot_Lig")

    copy_file(protein_name, "Input_Files")
    copy_file(prmtop, "Input_Files") if prmtop else None
    copy_file(inpcrd, "Input_Files") if inpcrd else None

    # Organize pre-MD files
    source_pre_md_files = [
        "prepared_no_solvent_",
        "solvent_padding_",
        "solvent_absolute_",
        "membrane_",
    ]
    destination_pre_md = "MD_Files/Pre_MD"
    organize_files(
        [f"{prefix}{protein_name}" for prefix in source_pre_md_files],
        destination_pre_md,
    )

    # Organize topology files after minimization and equilibration
    source_topology_files = ["Energyminimization_", "Equilibration_"]
    destination_topology = "MD_Files/Minimization_Equilibration"
    organize_files(
        [f"{prefix}{protein_name}" for prefix in source_topology_files],
        destination_topology,
    )

    # Organize simulation output files
    organize_files([f"output_{protein_name}", "trajectory.dcd"], "MD_Files/MD_Output")

    # Organize MDtraj and MDAnalysis files
    organize_files(
        [
            "centered_old_coordinates_top.pdb",
            "centered_old_coordinates.dcd",
            "centered_old_coordinates_top.gro",
            "centered_old_coordinates.xtc",
            "centered_traj_unaligned.dcd",
            "centered_traj_unaligned.xtc",
            "prot_lig_traj_unaligned.dcd",
            "prot_lig_traj_unaligned.xtc",
        ],
        "MD_Postprocessing",
    )
    organize_files(
        [
            "centered_top.pdb",
            "centered_traj.dcd",
            "centered_top.gro",
            "centered_traj.xtc",
        ],
        "Final_Output/All_Atoms",
    )
    organize_files(
        [
            "prot_lig_top.pdb",
            "prot_lig_traj.dcd",
            "prot_lig_top.gro",
            "prot_lig_traj.xtc",
        ],
        "Final_Output/Prot_Lig",
    )

    # Organize checkpoint files
    organize_files(
        ["checkpoint.chk", "10x_checkpoint.chk", "100x_checkpoint.chk"],
        "Checkpoints",
    )

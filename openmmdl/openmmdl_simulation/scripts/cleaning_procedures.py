import glob
import os
import shutil
import sys
import time
from typing import List


def _safe_rmtree(path: str, retries: int = 5, delay: float = 0.2):
    """Retry rmtree to tolerate transient EBUSY / directory-not-empty races."""
    for attempt in range(retries):
        if not os.path.exists(path):
            return
        try:
            shutil.rmtree(path)
            return
        except OSError as exc:
            if exc.errno not in (16, 39) or attempt == retries - 1:
                raise
            time.sleep(delay)


def close_reporters(simulation):
    """Close file-backed reporters without closing sys.stdout/sys.stderr."""
    for reporter in list(getattr(simulation, "reporters", [])):
        out = getattr(reporter, "_out", None)

        try:
            if out not in (sys.stdout, sys.stderr):
                reporter.close()
                continue
        except Exception:
            pass

        for attr in ("_out", "_file", "file"):
            handle = getattr(reporter, attr, None)
            if handle is not None and handle not in (sys.stdout, sys.stderr):
                try:
                    handle.close()
                except Exception:
                    pass

    simulation.reporters.clear()


def cleanup_post_md():
    """Remove intermediate folders after outputs have been organized.

    Keeps only the copied input files and final outputs.
    """
    print("Removing intermediate postprocessing folders :)")
    for directory in ("MD_Files", "MD_Postprocessing", "Checkpoints"):
        if os.path.exists(directory):
            _safe_rmtree(directory)
    print("Cleanup is done.")


def create_directory_if_not_exists(directory_path):
    """Create a directory if it doesn't exist, or overwrite it if it does.

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
    if src and os.path.exists(src):
        shutil.copy(src, dest)


def organize_files(source, destination):
    """Move files from source list into destination if they exist.

    Args:
        source (list[str]): Paths of files that need to be moved.
        destination (str): Path of destination where the files need to be moved.

    Returns:
        None
    """
    for file in source:
        if os.path.exists(file):
            os.rename(file, os.path.join(destination, os.path.basename(file)))


def post_md_file_movement(
    protein_name: str,
    prmtop: str = None,
    inpcrd: str = None,
    ligands: List[str] = None,
    mda_selection: str = "mda_prot_lig_all",
):
    """Organize and move the files after the MD simulation.

    Args:
        protein_name (str): Name of the protein PDB.
        prmtop (str, optional): Path to the AMBER topology file.
        inpcrd (str, optional): Path to the AMBER coordinate file.
        ligands (List[str], optional): List of paths to the ligand files.
        mda_selection (str, optional): Which MDAnalysis outputs were requested.
            Supported values:
                - "mda_all"
                - "mda_prot_lig"
                - "mda_prot_lig_all"

    Returns:
        None
    """
    write_all_atoms = mda_selection in ("mda_prot_lig_all", "mda_all")
    write_prot_lig = mda_selection in ("mda_prot_lig_all", "mda_prot_lig")

    # Create necessary directories
    create_directory_if_not_exists("Input_Files")
    create_directory_if_not_exists("MD_Files")
    create_directory_if_not_exists("MD_Files/Pre_MD")
    create_directory_if_not_exists("MD_Files/Minimization_Equilibration")
    create_directory_if_not_exists("MD_Files/MD_Output")
    create_directory_if_not_exists("MD_Postprocessing")
    create_directory_if_not_exists("Final_Output")

    if write_all_atoms:
        create_directory_if_not_exists("Final_Output/All_Atoms")
    if write_prot_lig:
        create_directory_if_not_exists("Final_Output/Prot_Lig")

    create_directory_if_not_exists("Checkpoints")

    # Move input ligand files
    if ligands:
        for lig in ligands:
            copy_file(lig, "Input_Files")
            if write_all_atoms:
                copy_file(lig, "Final_Output/All_Atoms")
            if write_prot_lig:
                copy_file(lig, "Final_Output/Prot_Lig")

    # Copy ligand partial-charge MOL2 exports (if generated)
    for mol2 in glob.glob("ligand_*_pc.mol2"):
        copy_file(mol2, "Input_Files")
        if write_all_atoms:
            copy_file(mol2, "Final_Output/All_Atoms")
        if write_prot_lig:
            copy_file(mol2, "Final_Output/Prot_Lig")

    # Copy core input files
    copy_file(protein_name, "Input_Files")
    if prmtop:
        copy_file(prmtop, "Input_Files")
    if inpcrd:
        copy_file(inpcrd, "Input_Files")

    # Organize pre-MD files
    source_pre_md_files = [
        "prepared_no_solvent_",
        "solvent_padding_",
        "solvent_absolute_",
        "membrane_",
    ]
    organize_files(
        [f"{prefix}{protein_name}" for prefix in source_pre_md_files],
        "MD_Files/Pre_MD",
    )

    # Organize topology files after minimization and equilibration
    source_topology_files = ["Energyminimization_", "Equilibration_"]
    organize_files(
        [f"{prefix}{protein_name}" for prefix in source_topology_files],
        "MD_Files/Minimization_Equilibration",
    )

    # Organize simulation output files
    organize_files(
        [f"output_{protein_name}", "trajectory.dcd"],
        "MD_Files/MD_Output",
    )

    # Organize intermediate MDtraj / MDAnalysis files
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

    # Organize final all-atom outputs
    if write_all_atoms:
        organize_files(
            [
                "centered_top.pdb",
                "centered_traj.dcd",
                "centered_top.gro",
                "centered_traj.xtc",
            ],
            "Final_Output/All_Atoms",
        )

    # Organize final protein-ligand outputs
    if write_prot_lig:
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

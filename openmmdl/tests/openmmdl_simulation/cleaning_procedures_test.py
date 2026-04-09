import os
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

from openmmdl.openmmdl_simulation.scripts.cleaning_procedures import (
    cleanup_post_md,
    copy_file,
    create_directory_if_not_exists,
    organize_files,
    post_md_file_movement,
)


@pytest.fixture
def test_directory_path():
    return "test_directory"


def test_cleanup_post_md(tmp_path):
    original_cwd = Path.cwd()
    try:
        os.chdir(tmp_path)

        Path("MD_Files").mkdir()
        Path("MD_Postprocessing").mkdir()
        Path("Checkpoints").mkdir()
        Path("Input_Files").mkdir()
        Path("Final_Output").mkdir()

        cleanup_post_md()

        assert not Path("MD_Files").exists()
        assert not Path("MD_Postprocessing").exists()
        assert not Path("Checkpoints").exists()
        assert Path("Input_Files").exists()
        assert Path("Final_Output").exists()
    finally:
        os.chdir(original_cwd)


def test_create_directory_if_not_exists(test_directory_path):
    create_directory_if_not_exists(test_directory_path)
    assert os.path.exists(test_directory_path)

    create_directory_if_not_exists(test_directory_path)
    assert os.path.exists(test_directory_path)

    shutil.rmtree(test_directory_path)
    assert not os.path.exists(test_directory_path)


@patch("os.path.exists")
@patch("shutil.copy")
def test_copy_file(mock_copy, mock_exists):
    src = "source_file.txt"
    dest = "destination_directory"

    mock_exists.return_value = True

    copy_file(src, dest)

    mock_exists.assert_called_with(src)
    mock_copy.assert_called_with(src, dest)


@patch("os.path.exists")
@patch("os.rename")
def test_organize_files(mock_rename, mock_exists):
    source = ["file1.txt", "file2.txt", "file3.txt"]
    destination = "destination_directory"

    mock_exists.side_effect = [True] * len(source)

    organize_files(source, destination)

    assert mock_rename.call_count == 3
    mock_rename.assert_any_call("file1.txt", os.path.join(destination, "file1.txt"))
    mock_rename.assert_any_call("file2.txt", os.path.join(destination, "file2.txt"))
    mock_rename.assert_any_call("file3.txt", os.path.join(destination, "file3.txt"))


def test_post_md_file_movement(tmp_path):
    original_cwd = Path.cwd()
    repo_root = Path(__file__).resolve().parents[2]
    test_data_directory = repo_root / "tests" / "data" / "in"

    ligand = test_data_directory / "CVV.sdf"
    protein_name = test_data_directory / "6b73.pdb"
    prmtop = test_data_directory / "6b73.prmtop"
    inpcrd = test_data_directory / "6b73.inpcrd"
    protein_no_solvent = test_data_directory / "prepared_no_solvent_6b73.pdb"
    protein_solvent = test_data_directory / "solvent_padding_6b73.pdb"
    protein_equilibration = test_data_directory / "Equilibration_6b73.pdb"
    protein_minimization = test_data_directory / "Energyminimization_6b73.pdb"
    output_pdb = test_data_directory / "output_6b73.pdb"
    mdtraj_top = test_data_directory / "centered_old_coordinates_top.pdb"
    prot_lig_top = test_data_directory / "prot_lig_top.pdb"
    checkpoint = test_data_directory / "checkpoint.chk"
    checkpoint_10x = test_data_directory / "10x_checkpoint.chk"

    assert ligand.exists()
    assert protein_name.exists()
    assert prmtop.exists()
    assert inpcrd.exists()
    assert protein_no_solvent.exists()

    try:
        os.chdir(tmp_path)

        shutil.copy(str(protein_no_solvent), ".")
        shutil.copy(str(protein_solvent), ".")
        shutil.copy(str(protein_equilibration), ".")
        shutil.copy(str(protein_minimization), ".")
        shutil.copy(str(output_pdb), ".")
        shutil.copy(str(mdtraj_top), ".")
        shutil.copy(str(prot_lig_top), ".")
        shutil.copy(str(checkpoint), ".")
        shutil.copy(str(checkpoint_10x), ".")
        shutil.copy(str(protein_name), ".")

        local_protein_name = "6b73.pdb"

        post_md_file_movement(
            str(local_protein_name),
            str(prmtop),
            str(inpcrd),
            [str(ligand)],
        )

        input_files_dir = Path("Input_Files")
        md_files_dir = Path("MD_Files")
        md_postprocessing_dir = Path("MD_Postprocessing")
        final_output_dir = Path("Final_Output")
        checkpoints_dir = Path("Checkpoints")

        assert input_files_dir.exists()
        assert (md_files_dir / "Pre_MD").exists()
        assert (md_files_dir / "Pre_MD" / "prepared_no_solvent_6b73.pdb").exists()
        assert (md_files_dir / "Pre_MD" / "solvent_padding_6b73.pdb").exists()
        assert (
            md_files_dir / "Minimization_Equilibration" / "Equilibration_6b73.pdb"
        ).exists()
        assert (
            md_files_dir
            / "Minimization_Equilibration"
            / "Energyminimization_6b73.pdb"
        ).exists()
        assert (md_files_dir / "MD_Output" / "output_6b73.pdb").exists()
        assert (md_postprocessing_dir / "centered_old_coordinates_top.pdb").exists()
        assert (final_output_dir / "Prot_Lig" / "prot_lig_top.pdb").exists()
        assert (checkpoints_dir / "checkpoint.chk").exists()
        assert (checkpoints_dir / "10x_checkpoint.chk").exists()
    finally:
        os.chdir(original_cwd)


if __name__ == "__main__":
    pytest.main()

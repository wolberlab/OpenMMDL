from pathlib import Path

from openmmdl.openmmdl_simulation.scripts.export import (
    PreparedSystemBuilder,
    PreparedSystemExporter,
)


TEST_DIR = Path(__file__).resolve().parent
DATA_DIR = TEST_DIR / "data/in"

PROTEIN_FILE = DATA_DIR / "6b73.pdb"
LIGAND_FILE = DATA_DIR / "CVV.sdf"


def protein_only_config():
    return {
        "protein": str(PROTEIN_FILE),
        "ligand": None,
        "ligand_name": None,
        "ligand_minimization": False,
        "ligand_sanitization": False,
        "forcefield": "AMBER19",
        "water_model": "TIP3P",
        "use_solvent": True,
        "add_membrane": False,
        "small_molecule_forcefield": None,
        "small_molecule_forcefield_version": None,
        "nonbonded_method": "PME",
        "nonbonded_cutoff": 1.0,
        "ewald_error_tolerance": 0.0005,
        "constraints": None,
        "rigid_water": False,
        "hydrogen_mass": None,
        "water_box_mode": "Buffer",
        "water_box_shape": "cube",
        "water_padding_distance": 1.0,
        "water_positive_ion": "Na+",
        "water_negative_ion": "Cl-",
        "water_ionicstrength": 0.15,
        "output_prefix": "prepared_system",
    }


def protein_ligand_config():
    return {
        "protein": str(PROTEIN_FILE),
        "ligand": str(LIGAND_FILE),
        "ligand_name": "UNK",
        "ligand_minimization": False,
        "ligand_sanitization": False,
        "forcefield": "AMBER19",
        "water_model": "TIP3P",
        "use_solvent": True,
        "add_membrane": False,
        "small_molecule_forcefield": "gaff",
        "small_molecule_forcefield_version": "gaff-2.2.20",
        "nonbonded_method": "PME",
        "nonbonded_cutoff": 1.0,
        "ewald_error_tolerance": 0.0005,
        "constraints": None,
        "rigid_water": False,
        "hydrogen_mass": None,
        "water_box_mode": "Buffer",
        "water_box_shape": "cube",
        "water_padding_distance": 1.0,
        "water_positive_ion": "Na+",
        "water_negative_ion": "Cl-",
        "water_ionicstrength": 0.15,
        "output_prefix": "prepared_system",
    }


def test_input_files_exist():
    assert PROTEIN_FILE.exists(), f"Missing test protein file: {PROTEIN_FILE}"
    assert LIGAND_FILE.exists(), f"Missing test ligand file: {LIGAND_FILE}"


def test_builder_builds_protein_only_system():
    builder = PreparedSystemBuilder(protein_only_config())
    topology, system, positions = builder.build()

    assert topology is not None
    assert system is not None
    assert positions is not None
    assert topology.getNumAtoms() > 0


def test_builder_builds_protein_ligand_system():
    builder = PreparedSystemBuilder(protein_ligand_config())
    topology, system, positions = builder.build()

    assert topology is not None
    assert system is not None
    assert positions is not None
    assert topology.getNumAtoms() > 0


def test_export_processed_pdb_and_xml(tmp_path):
    builder = PreparedSystemBuilder(protein_ligand_config())
    topology, system, positions = builder.build()

    exporter = PreparedSystemExporter(
        topology=topology,
        system=system,
        positions=positions,
        output_prefix="prepared_system",
        output_dir=tmp_path,
    )

    written = exporter.export_selected(["processed_pdb", "openmm_xml"])
    written_paths = {Path(p) for p in written}

    expected_files = {
        tmp_path / "prepared_system.pdb",
        tmp_path / "prepared_system.xml",
    }

    assert expected_files.issubset(written_paths)
    for path in expected_files:
        assert path.exists()
        assert path.stat().st_size > 0


def test_export_amber(tmp_path):
    builder = PreparedSystemBuilder(protein_ligand_config())
    topology, system, positions = builder.build()

    exporter = PreparedSystemExporter(
        topology=topology,
        system=system,
        positions=positions,
        output_prefix="prepared_system",
        output_dir=tmp_path,
    )

    written = exporter.export_selected(["amber"])
    written_paths = {Path(p) for p in written}

    expected_files = {
        tmp_path / "prepared_system.prmtop",
        tmp_path / "prepared_system.inpcrd",
    }

    assert expected_files.issubset(written_paths)
    for path in expected_files:
        assert path.exists()
        assert path.stat().st_size > 0


def test_export_gromacs(tmp_path):
    builder = PreparedSystemBuilder(protein_ligand_config())
    topology, system, positions = builder.build()

    exporter = PreparedSystemExporter(
        topology=topology,
        system=system,
        positions=positions,
        output_prefix="prepared_system",
        output_dir=tmp_path,
    )

    written = exporter.export_selected(["gromacs"])
    written_paths = {Path(p) for p in written}

    expected_files = {
        tmp_path / "prepared_system.top",
        tmp_path / "prepared_system.gro",
    }

    assert expected_files.issubset(written_paths)
    for path in expected_files:
        assert path.exists()
        assert path.stat().st_size > 0


def test_export_psf(tmp_path):
    builder = PreparedSystemBuilder(protein_ligand_config())
    topology, system, positions = builder.build()

    exporter = PreparedSystemExporter(
        topology=topology,
        system=system,
        positions=positions,
        output_prefix="prepared_system",
        output_dir=tmp_path,
    )

    written = exporter.export_selected(["psf"])
    written_paths = {Path(p) for p in written}

    expected_files = {
        tmp_path / "prepared_system.psf",
        tmp_path / "prepared_system.pdb",
    }

    assert expected_files.issubset(written_paths)
    for path in expected_files:
        assert path.exists()
        assert path.stat().st_size > 0
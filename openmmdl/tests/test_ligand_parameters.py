"""Bespoke ligand parameters: BespokeFit orchestration, OpenMM use and Amber export."""

import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import openmmdl.cli.cli as cli
from openmmdl.ligand_parameters import core
from openmmdl.ligand_parameters.cli import main as ligand_main
from openmmdl.openmmdl_setup.ligand_parameters import amber_prepare_command, read_parameter_settings


# --------------------------------------------------------------------------- BespokeFit command


def test_fit_command_passes_qc_spec_and_memory_per_core(tmp_path):
    options = core.FitOptions(qc_program="xtb", qc_method="gfn2xtb", qc_basis="none", cores=4, memory_gb=6)
    command = core.fit_command(tmp_path / "ligand.sdf", tmp_path, options)
    assert command[:3] == ["openff-bespoke", "executor", "run"]
    assert command[command.index("--default-qc-spec") + 1 : command.index("--default-qc-spec") + 4] == ["xtb", "gfn2xtb", "none"]
    assert command[command.index("--qc-compute-max-mem") + 1] == "1.5"
    assert command[command.index("--output-force-field") + 1] == str(tmp_path / "final.offxml")


@pytest.mark.parametrize("kwargs", [{"workers": 0}, {"cores": -1}, {"memory_gb": 0}])
def test_fit_options_reject_non_positive_resources(kwargs):
    with pytest.raises(ValueError):
        core.FitOptions(**kwargs)


def _fake_bespoke(tmp_path, body):
    """A stand-in for openff-bespoke: a python script that receives the executor arguments."""
    script = tmp_path / "openff-bespoke"
    script.write_text(f"#!{sys.executable}\nimport sys\nargs = sys.argv[1:]\n{body}\n")
    script.chmod(0o755)
    return str(script)


def test_run_fit_logs_output_and_returns_force_field(tmp_path):
    executable = _fake_bespoke(
        tmp_path, "print('fitting'); open(args[args.index('--output-force-field') + 1], 'w').write('<SMIRNOFF/>')"
    )
    (tmp_path / "fit").mkdir()
    final = core.run_fit(tmp_path / "ligand.sdf", tmp_path / "fit", core.FitOptions(executable=executable))
    assert final.read_text() == "<SMIRNOFF/>"
    assert "fitting" in (tmp_path / "fit" / "fit.log").read_text()


def test_run_fit_fails_without_force_field(tmp_path):
    executable = _fake_bespoke(tmp_path, "print('QC failed'); sys.exit(1)")
    (tmp_path / "fit").mkdir()
    with pytest.raises(RuntimeError, match="did not produce a force field"):
        core.run_fit(tmp_path / "ligand.sdf", tmp_path / "fit", core.FitOptions(executable=executable))
    with pytest.raises(RuntimeError, match="not found"):
        core.run_fit(tmp_path / "ligand.sdf", tmp_path / "fit", core.FitOptions(executable="no-such-bespoke"))


@pytest.mark.skipif(os.name != "posix", reason="process groups")
def test_run_fit_stops_executor_that_hangs_after_writing_the_result(tmp_path, monkeypatch):
    monkeypatch.setattr(core, "POLL_SECONDS", 0.2)
    monkeypatch.setattr(core, "SHUTDOWN_GRACE_SECONDS", 0.5)
    executable = _fake_bespoke(
        tmp_path,
        "import time; open(args[args.index('--output-force-field') + 1], 'w').write('<SMIRNOFF/>'); time.sleep(60)",
    )
    (tmp_path / "fit").mkdir()
    assert core.run_fit(tmp_path / "ligand.sdf", tmp_path / "fit", core.FitOptions(executable=executable)).is_file()


# --------------------------------------------------------------------------- setup form and Amber command


def test_setup_form_settings(tmp_path):
    assert read_parameter_settings({}, {}, "sdfFile") == {"mode": "standard"}
    with pytest.raises(ValueError, match="single complex"):
        read_parameter_settings({"ligandParameterMode": "fit"}, {}, "sdfFile", library=True)
    with pytest.raises(ValueError, match="Upload the ligand"):
        read_parameter_settings({"ligandParameterMode": "import"}, {}, "sdfFile")


def test_amber_prepare_command_is_shell_safe():
    command = shlex.split(amber_prepare_command({"mode": "import", "offxml": "my fit.offxml"}, "lig $(x).sdf"))
    assert command[:3] == ["openmmdl", "ligand", "import"]
    assert command[command.index("--ligand") + 1] == "lig $(x).sdf"
    assert command[command.index("--offxml") + 1] == "my fit.offxml"
    assert command[command.index("--target") + 1] == "amber"
    options = core.FitOptions(executable="/env/bin/openff-bespoke", cores=8)
    command = shlex.split(amber_prepare_command({"mode": "fit", "options": options.__dict__}, "lig.sdf"))
    assert command[command.index("--executable") + 1] == "/env/bin/openff-bespoke"
    assert command[command.index("--cores") + 1] == "8"


# --------------------------------------------------------------------------- CLI and launcher


def test_ligand_cli_is_registered_and_shows_help(capsys):
    assert cli.COMMANDS["ligand"][0] == "openmmdl.ligand_parameters.cli:main"
    assert cli.main(["ligand", "--help"]) == 0
    assert "fit" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        ligand_main(["fit", "--help"])
    assert "--qc-spec" in capsys.readouterr().out


def test_simulation_launcher_copies_parameters_and_keeps_fit_on_retry(tmp_path, monkeypatch):
    from openmmdl.openmmdl_simulation import openmmdlsimulation as runner
    from openmmdl.openmmdl_simulation.cli import build_parser

    monkeypatch.chdir(tmp_path)
    for name, content in [("run.py", "pass"), ("protein.pdb", "END"), ("custom.offxml", "<SMIRNOFF/>")]:
        (tmp_path / name).write_text(content)
    (tmp_path / "old_fit").mkdir()
    (tmp_path / "old_fit" / "final.offxml").write_text("<SMIRNOFF/>")
    calls = []

    class Process:
        def __init__(self, *args, **kwargs):
            assert Path("custom.offxml").is_file() and Path("old_fit/final.offxml").is_file()
            if not calls:
                Path("ligand_parameters").mkdir()
                Path("ligand_parameters/manifest.json").write_text("{}")
                self.returncode, self.stdout = 1, iter(["Particle coordinate is NaN\n"])
            else:
                assert Path("ligand_parameters/manifest.json").is_file()
                self.returncode, self.stdout = 0, iter([])
            calls.append(1)

        def wait(self):
            return self.returncode

    monkeypatch.setattr(runner.subprocess, "Popen", Process)
    args = build_parser().parse_args(
        ["-f", str(tmp_path / "run"), "-t", "protein.pdb", "-s", "run.py", "-p", "custom.offxml", "-p", "old_fit", "--failure-retries", "1"]
    )
    assert runner.run_simulation(args) == 0
    assert len(calls) == 2


# --------------------------------------------------------------------------- scientific checks


@pytest.fixture
def ligand_and_forcefield(tmp_path):
    toolkit = pytest.importorskip("openff.toolkit")
    pytest.importorskip("openff.interchange")
    pytest.importorskip("openmmforcefields")
    from openff.units import unit

    def make(smiles="CCO"):
        mol = toolkit.Molecule.from_smiles(smiles)
        mol.generate_conformers(n_conformers=1)
        # Molecule-specific charges and a changed torsion make the force field distinguishable
        # from the stock Sage parameters used by the standard SMIRNOFF option.
        charges = np.linspace(-0.2, 0.2, mol.n_atoms)
        mol.partial_charges = (charges - charges.mean()) * unit.elementary_charge
        ff = toolkit.ForceField("openff-2.2.0.offxml")
        ff["LibraryCharges"].add_parameter(
            {"smirks": mol.to_smiles(mapped=True), **{f"charge{i + 1}": q for i, q in enumerate(mol.partial_charges)}}
        )
        for parameter in ff["ProperTorsions"].parameters:
            parameter.k1 *= 1.25
        sdf, offxml = tmp_path / "ligand.sdf", tmp_path / "ligand.offxml"
        mol.partial_charges = None
        mol.to_file(str(sdf), "SDF")
        ff.to_file(str(offxml))
        return mol, ff, sdf, offxml

    return make


def _energy(system, positions):
    import openmm
    from openmm import unit

    context = openmm.Context(system, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    context.setPositions(positions)
    return context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)


def _charges(system):
    from openmm import NonbondedForce, unit

    force = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    return [force.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge) for i in range(force.getNumParticles())]


def test_load_ligand_and_validate_offxml(ligand_and_forcefield, tmp_path):
    mol, ff, sdf, offxml = ligand_and_forcefield()
    assert core.load_ligand(sdf).n_atoms == mol.n_atoms
    two = tmp_path / "two.sdf"
    two.write_text(sdf.read_text() * 2)
    with pytest.raises(ValueError, match="exactly one molecule"):
        core.load_ligand(two)
    with pytest.raises(ValueError, match="SDF or MOL"):
        core.load_ligand(tmp_path / "ligand.pdb")
    assert core.validate_offxml(offxml) is not None
    bad = tmp_path / "system.offxml"
    bad.write_text("<ForceField><AtomTypes/></ForceField>")
    with pytest.raises(ValueError, match="not a readable SMIRNOFF"):
        core.validate_offxml(bad)


def test_import_is_reused_but_not_mixed_with_other_runs(ligand_and_forcefield, tmp_path):
    mol, ff, sdf, offxml = ligand_and_forcefield()
    output = tmp_path / "parameters"
    final = core.prepare_parameters(sdf, output, offxml=offxml)
    assert final == output / "final.offxml" and (output / "manifest.json").is_file()
    assert core.prepare_parameters(sdf, output, offxml=offxml) == final
    with pytest.raises(ValueError, match="another or unfinished run"):
        core.prepare_parameters(sdf, output, offxml=offxml, target="amber")
    with pytest.raises(ValueError, match="either"):
        core.prepare_parameters(sdf, tmp_path / "x", offxml=offxml, options=core.FitOptions())


def test_fit_result_is_validated_and_cosmetic_attributes_dropped(ligand_and_forcefield, tmp_path, monkeypatch):
    from openff.toolkit import ForceField

    mol, ff, sdf, offxml = ligand_and_forcefield()
    text = offxml.read_text().replace("<Bonds ", '<Bonds parameterize="k" ', 1)
    assert "parameterize" in text

    def fake_fit(sdf, directory, options):
        assert Path(sdf).is_file()
        (directory / "final.offxml").write_text(text)

    monkeypatch.setattr(core, "run_fit", fake_fit)
    final = core.prepare_parameters(sdf, tmp_path / "fit", options=core.FitOptions())
    assert "parameterize" not in final.read_text()
    assert ForceField(str(final)) is not None
    assert json.loads((tmp_path / "fit" / "manifest.json").read_text())["options"]["qc_program"] == "psi4"


def test_import_cli_reports_uncovered_ligand(ligand_and_forcefield, tmp_path):
    from openff.toolkit import ForceField

    mol, ff, sdf, offxml = ligand_and_forcefield()
    incomplete = tmp_path / "incomplete.offxml"
    ff = ForceField(str(offxml))
    ff.deregister_parameter_handler("Bonds")
    ff.to_file(str(incomplete))
    with pytest.raises(Exception):
        ligand_main(["import", "--ligand", str(sdf), "--offxml", str(incomplete), "--output", str(tmp_path / "out")])
    assert ligand_main(["import", "--ligand", str(sdf), "--offxml", str(offxml), "--output", str(tmp_path / "ok")]) == 0
    assert (tmp_path / "ok" / "final.offxml").is_file()


def test_openmm_forcefield_uses_bespoke_parameters_for_the_first_ligand(ligand_and_forcefield):
    from openff.toolkit import Molecule
    from openmm import app, unit
    from openmmdl.openmmdl_simulation.scripts.forcefield_water import generate_forcefield

    mol, ff, sdf, offxml = ligand_and_forcefield("CCO")
    other = Molecule.from_smiles("CC(=O)O")
    other.generate_conformers(n_conformers=1)
    forcefield = generate_forcefield(
        "amber14-all.xml", "amber14/tip3p.xml", False, "gaff", "gaff-2.11",
        [mol.to_rdkit(), other.to_rdkit()], ligand_offxml=str(offxml),
    )
    for molecule in (mol, other):
        topology = molecule.to_topology().to_openmm()
        system = forcefield.createSystem(topology, nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
        positions = molecule.conformers[0].to_openmm()
        if molecule is mol:
            direct = ff.create_openmm_system(mol.to_topology())
            for force in direct.getForces():
                if hasattr(force, "setNonbondedMethod"):
                    force.setNonbondedMethod(0)
            assert _charges(system) == pytest.approx(_charges(direct), abs=1e-10)
            assert _energy(system, positions) == pytest.approx(_energy(direct, positions), abs=1e-6)
        else:
            assert np.isfinite(_energy(system, positions))  # parameterized by GAFF


@pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1"])
def test_amber_export_and_restore_keep_the_fitted_parameters(ligand_and_forcefield, tmp_path, smiles):
    from openmm import app, unit
    from openmmdl.ligand_parameters.amber import restore_ligand_parameters

    tleap = shutil.which("tleap")
    if not tleap:
        pytest.skip("AmberTools tleap is required")
    mol, ff, sdf, offxml = ligand_and_forcefield(smiles)
    output = tmp_path / "bespoke_ligand"
    core.prepare_parameters(sdf, output, offxml=offxml, target="amber")
    (output / "tleap.in").write_text(
        "loadamberparams ligand.frcmod\nloadoff ligand.lib\nx = loadpdb ligand.pdb\nsaveamberparm x test.prmtop test.inpcrd\nquit\n"
    )
    run = subprocess.run([tleap, "-f", "tleap.in"], cwd=output, capture_output=True, text=True)
    assert run.returncode == 0 and "Errors = 0" in run.stdout, run.stdout
    restore_ligand_parameters(output / "test.prmtop", output / "test.inpcrd", output)
    settings = dict(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    reference = app.AmberPrmtopFile(str(output / "ligand.prmtop")).createSystem(**settings)
    restored = app.AmberPrmtopFile(str(output / "test.prmtop")).createSystem(**settings)
    coordinates = app.AmberInpcrdFile(str(output / "test.inpcrd")).positions
    assert _charges(restored) == pytest.approx(_charges(reference), abs=1e-6)
    # distorted geometry, so that the improper terms of the planar ring contribute
    distorted = np.array(coordinates.value_in_unit(unit.nanometer)) + np.random.default_rng(4).normal(scale=0.015, size=(mol.n_atoms, 3))
    assert _energy(restored, distorted * unit.nanometer) == pytest.approx(_energy(reference, distorted * unit.nanometer), abs=1e-3)


@pytest.mark.parametrize("water,box", [("tip3p", "TIP3PBOX"), ("opc", "OPCBOX")])
def test_amber_restore_keeps_the_solvated_receptor(ligand_and_forcefield, tmp_path, water, box):
    import parmed as pmd
    from openmm import app, unit
    from openmmdl.ligand_parameters.amber import restore_ligand_parameters

    tleap = shutil.which("tleap")
    if not tleap:
        pytest.skip("AmberTools tleap is required")
    mol, ff, sdf, offxml = ligand_and_forcefield("c1ccccc1")
    output = tmp_path / "bespoke_ligand"
    core.prepare_parameters(sdf, output, offxml=offxml, target="amber")
    (output / "tleap.in").write_text(
        f"source leaprc.protein.ff14SB\nsource leaprc.water.{water}\n"
        "loadamberparams ligand.frcmod\nloadoff ligand.lib\n"
        "lig = loadpdb ligand.pdb\nrcp = sequence { ACE ALA NME }\n"
        "translate rcp { 0 0 12 }\nx = combine { rcp lig }\n"
        f"solvatebox x {box} 5\nsaveamberparm x system.prmtop system.inpcrd\nquit\n"
    )
    run = subprocess.run([tleap, "-f", "tleap.in"], cwd=output, capture_output=True, text=True)
    assert run.returncode == 0 and "Errors = 0" in run.stdout, run.stdout
    original = pmd.load_file(str(output / "system.prmtop"), xyz=str(output / "system.inpcrd"))
    restore_ligand_parameters(output / "system.prmtop", output / "system.inpcrd", output)
    restored = pmd.load_file(str(output / "system.prmtop"), xyz=str(output / "system.inpcrd"))
    assert len(original.atoms) == len(restored.atoms)
    assert np.allclose(original.box, restored.box) and np.allclose(original.coordinates, restored.coordinates)
    settings = dict(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    before, after = original["!:BSP"], restored["!:BSP"]
    assert [(a.name, a.residue.name) for a in before.atoms] == [(a.name, a.residue.name) for a in after.atoms]
    assert _energy(before.createSystem(**settings), before.coordinates * unit.angstrom) == pytest.approx(
        _energy(after.createSystem(**settings), before.coordinates * unit.angstrom), abs=1e-5
    )
    reference = pmd.load_file(str(output / "ligand.prmtop"), xyz=str(output / "ligand.inpcrd"))
    ligand = restored[":BSP"]
    distorted = (ligand.coordinates + np.random.default_rng(14).normal(scale=0.15, size=ligand.coordinates.shape)) * unit.angstrom
    assert _energy(ligand.createSystem(**settings), distorted) == pytest.approx(
        _energy(reference.createSystem(**settings), distorted), abs=1e-5
    )

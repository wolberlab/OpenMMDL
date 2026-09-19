"""Bespoke ligand parameters with OpenFF BespokeFit.

Fits (or imports) a SMIRNOFF force field for the primary ligand and prepares it
for the OpenMM (PDB) path or, as a tleap library, for the Amber path. 
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import signal
import subprocess
import time
from copy import deepcopy
from dataclasses import asdict, dataclass
from pathlib import Path

# Residue code of the bespoke ligand in the Amber path.
RESNAME = "BSP"
# How often the running fit is checked, and how long a finished fit may keep its executor running.
POLL_SECONDS = 10
SHUTDOWN_GRACE_SECONDS = 120


@dataclass(frozen=True)
class FitOptions:
    """Settings passed on to `openff-bespoke executor run`."""

    executable: str = "openff-bespoke"  # e.g. a path into a separate conda environment
    qc_program: str = "psi4"
    qc_method: str = "b3lyp-d3bj"
    qc_basis: str = "dzvp"
    workers: int = 1  # QC compute workers
    cores: int = 2  # per worker
    memory_gb: int = 4  # per worker

    def __post_init__(self):
        if min(self.workers, self.cores, self.memory_gb) < 1:
            raise ValueError("QC workers, cores and memory must be at least 1")


def load_ligand(path):
    """Read a single ligand (SDF or MOL, explicit hydrogens) as an OpenFF Molecule."""
    from openff.toolkit import Molecule

    path = Path(path)
    if path.suffix.lower() not in {".sdf", ".mol"}:
        raise ValueError("Bespoke ligand parameters need the ligand as an SDF or MOL file")
    try:
        molecule = Molecule.from_file(str(path), allow_undefined_stereo=True)
    except Exception as exc:
        raise ValueError(f"Could not read the ligand {path.name}: {exc}") from exc
    if isinstance(molecule, list):
        raise ValueError(f"{path.name} must contain exactly one molecule")
    return molecule


def validate_offxml(path):
    """Return the SMIRNOFF force field stored in an .offxml file, or raise ValueError."""
    from openff.toolkit import ForceField

    try:
        return ForceField(str(path), allow_cosmetic_attributes=True)
    except Exception as exc:
        raise ValueError(f"{Path(path).name} is not a readable SMIRNOFF force field: {exc}") from exc


def fit_command(sdf, directory, options):
    """Command line for `openff-bespoke executor run` on `sdf`, writing into `directory`."""
    directory = Path(directory)
    return [
        options.executable, "executor", "run",
        "--file", str(sdf),
        "--workflow", "default",
        "--default-qc-spec", options.qc_program, options.qc_method, options.qc_basis,
        "--output", str(directory / "fit.json"),
        "--output-force-field", str(directory / "final.offxml"),
        "--directory", str(directory / "executor"),
        "--n-qc-compute-workers", str(options.workers),
        "--qc-compute-n-cores", str(options.cores),
        "--qc-compute-max-mem", str(options.memory_gb / options.cores),  # BespokeFit takes GB per core
    ]


def run_fit(sdf, directory, options):
    """Run BespokeFit in a subprocess (output in directory/fit.log) and return final.offxml."""
    directory = Path(directory)
    executable = shutil.which(options.executable)
    if executable is None:
        raise RuntimeError(
            f"BespokeFit executable not found: {options.executable}. Install openff-bespokefit "
            "(with psi4 or xtb) or give the full path to openff-bespoke."
        )
    command = [executable, *fit_command(sdf, directory, options)[1:]]
    # The executor starts redis and celery workers, which must come from the same environment.
    env = dict(os.environ, PATH=str(Path(executable).parent) + os.pathsep + os.environ.get("PATH", ""))
    env["PSI_SCRATCH"] = str(directory / "scratch")
    (directory / "scratch").mkdir(exist_ok=True)
    final = directory / "final.offxml"
    with (directory / "fit.log").open("w") as log:
        process = subprocess.Popen(
            command, cwd=directory, env=env, stdout=log, stderr=subprocess.STDOUT, start_new_session=True
        )
        completed_at = None
        while process.poll() is None:
            time.sleep(POLL_SECONDS)
            # The force field is only written after a successful fit. The executor can hang
            # while shutting down its workers afterwards, so stop it once the result is in.
            if final.is_file():
                completed_at = completed_at or time.monotonic()
                if time.monotonic() - completed_at > SHUTDOWN_GRACE_SECONDS:
                    os.killpg(process.pid, signal.SIGTERM)
                    process.wait()
    if not final.is_file():
        raise RuntimeError(
            f"BespokeFit did not produce a force field (exit code {process.returncode}), see {directory / 'fit.log'}"
        )
    return final


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_parameters(ligand, output_dir, *, offxml=None, options=None, target="openmm"):
    """Fit (`options`) or import (`offxml`) parameters for `ligand` and return the path of final.offxml.

    `output_dir` receives final.offxml, the fitting files and, for target "amber", the tleap
    library files written by `export_amber`. A directory that already holds the result for the
    same inputs is reused, so re-running a simulation script does not fit again.
    """
    if (offxml is None) == (options is None):
        raise ValueError("Give either an .offxml file to import or fitting options")
    if target not in {"openmm", "amber"}:
        raise ValueError("Target must be 'openmm' or 'amber'")
    output = Path(output_dir).resolve()  # the fit runs inside this directory
    final = output / "final.offxml"
    manifest = output / "manifest.json"
    inputs = {
        "ligand": _sha256(ligand),
        "offxml": _sha256(offxml) if offxml else None,
        "options": asdict(options) if options else None,
        "target": target,
    }
    if manifest.is_file() and json.loads(manifest.read_text()) == inputs:
        return final
    if output.is_dir() and any(output.iterdir()):
        raise ValueError(f"{output} contains files of another or unfinished run; remove it or use another directory")
    output.mkdir(parents=True, exist_ok=True)

    molecule = load_ligand(ligand)
    if offxml:
        validate_offxml(offxml)
        shutil.copyfile(offxml, final)
    else:
        sdf = output / "ligand.sdf"
        molecule.to_file(str(sdf), "SDF")
        run_fit(sdf, output, options)
    # Rewrite without ForceBalance's bookkeeping attributes so every SMIRNOFF reader accepts the file.
    forcefield = validate_offxml(final)
    forcefield.to_file(str(final), discard_cosmetic_attributes=True)

    # Check that the force field covers the ligand (raises on unassigned parameters).
    for i, atom in enumerate(molecule.atoms):
        atom.name = f"{atom.symbol}{i + 1}"
    interchange = forcefield.create_interchange(molecule.to_topology())
    if target == "amber":
        export_amber(interchange, output)
    manifest.write_text(json.dumps(inputs, indent=2) + "\n")
    return final


def export_amber(interchange, directory):
    """Write the parameterized ligand as a tleap library (ligand.lib, ligand.frcmod, ligand.pdb).

    tleap only keeps one of the three SMIRNOFF improper terms per center, so the exact
    ligand.prmtop/ligand.inpcrd are written too and put back after tleap by
    `openmmdl.ligand_parameters.amber.restore_ligand_parameters`.
    """
    import parmed as pmd

    directory = Path(directory)
    interchange.to_prmtop(directory / "ligand.prmtop")
    interchange.to_inpcrd(directory / "ligand.inpcrd")
    structure = pmd.load_file(str(directory / "ligand.prmtop"), xyz=str(directory / "ligand.inpcrd"))
    # Give every atom its own two-character type (digit + alphanumeric, unlike any Amber or
    # GAFF type) so the ligand parameters never overwrite receptor parameters in tleap.
    alphabet = "0123456789abcdefghijklmnopqrstuvwxyz"
    if len(structure.atoms) > 10 * len(alphabet):
        raise ValueError("The Amber export supports ligands with at most 360 atoms")
    for i, atom in enumerate(structure.atoms):
        atom.type = str(i // len(alphabet)) + alphabet[i % len(alphabet)]
        atom.atom_type = deepcopy(atom.atom_type)
        atom.atom_type.name = atom.type
        atom.residue.name = RESNAME
    pmd.amber.AmberParameterSet.from_structure(structure).write(str(directory / "ligand.frcmod"))
    residue = pmd.modeller.ResidueTemplate.from_residue(structure.residues[0])
    pmd.amber.AmberOFFLibrary.write({RESNAME: residue}, str(directory / "ligand.lib"))
    structure.save(str(directory / "ligand.pdb"), overwrite=True)

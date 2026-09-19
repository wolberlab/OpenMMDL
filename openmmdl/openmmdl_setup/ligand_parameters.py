"""Bespoke ligand parameter settings of the setup form, shared by the PDB and Amber pages."""

from __future__ import annotations

import shlex
import tempfile
from dataclasses import asdict
from pathlib import Path

from openmmdl.ligand_parameters.core import FitOptions, load_ligand, validate_offxml


def read_parameter_settings(form, uploads, ligand_key, library=False):
    """Return the settings dict stored in the session, or raise ValueError for invalid input.

    The dict holds "mode" ("standard", "fit" or "import") plus "options" (FitOptions fields)
    for a fit and "offxml" (uploaded filename) for an import.
    """
    mode = form.get("ligandParameterMode", "standard")
    if mode == "standard":
        return {"mode": mode}
    if mode not in {"fit", "import"}:
        raise ValueError("Unknown ligand parameter mode")
    if library:
        raise ValueError("Bespoke ligand parameters need the single complex mode; prepare each ligand separately")
    if ligand_key not in uploads:
        raise ValueError("Upload the ligand to fit or import bespoke parameters")
    _check_upload(uploads[ligand_key][0], load_ligand)
    if mode == "import":
        if "ligandOffxml" not in uploads:
            raise ValueError("Upload the SMIRNOFF force field (.offxml) to import")
        _check_upload(uploads["ligandOffxml"][0], validate_offxml)
        return {"mode": mode, "offxml": uploads["ligandOffxml"][0][1]}
    program, method, basis = form.get("bespokeQcSpec", "psi4 b3lyp-d3bj dzvp").split()
    options = FitOptions(
        executable=form.get("bespokeExecutable", "").strip() or "openff-bespoke",
        qc_program=program,
        qc_method=method,
        qc_basis=basis,
        workers=int(form.get("bespokeWorkers", 1)),
        cores=int(form.get("bespokeCores", 2)),
        memory_gb=int(form.get("bespokeMemory", 4)),
    )
    return {"mode": mode, "options": asdict(options)}


def _check_upload(upload, check):
    """Run `check` on a temporary copy of an uploaded (file object, filename) pair."""
    file, name = upload
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / name
        file.seek(0)
        path.write_bytes(file.read())
        file.seek(0)
        check(path)


def amber_prepare_command(settings, ligand):
    """Shell command of the Amber script that fits or imports the parameters into bespoke_ligand/."""
    command = ["openmmdl", "ligand", settings["mode"], "--ligand", ligand, "--output", "bespoke_ligand", "--target", "amber"]
    if settings["mode"] == "import":
        command += ["--offxml", settings["offxml"]]
    else:
        options = settings["options"]
        command += [
            "--executable", options["executable"],
            "--qc-spec", options["qc_program"], options["qc_method"], options["qc_basis"],
            "--workers", str(options["workers"]),
            "--cores", str(options["cores"]),
            "--memory-gb", str(options["memory_gb"]),
        ]
    return shlex.join(command)

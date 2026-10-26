"""Run options of ``openmmdl analysis``, read from the analysis argparse parser.

Flags, defaults, choices and help texts come from the real CLI parser, so the GUI form
cannot drift from the command line. ``FIELDS`` only adds what a form needs on top.
"""

from __future__ import annotations

import argparse
import os
import re
import shlex
import subprocess
from pathlib import Path

from openmmdl.openmmdl_analysis.cli import build_parser, parse_bool_flag

FIELDS = {
    "topology": {"label": "Topology", "group": "inputs", "kind": "file", "ext": [".pdb", ".tpr"]},
    "trajectory": {
        "label": "Trajectory",
        "group": "inputs",
        "kind": "file",
        "ext": [".dcd", ".xtc", ".trr"],
    },
    "ligand_sdf": {"label": "Ligand SDF", "group": "inputs", "kind": "file", "ext": [".sdf"]},
    "ligand_name": {"label": "Ligand residue name", "group": "ligand", "kind": "text"},
    "special_ligand": {"label": "Special ligand residue name", "group": "ligand", "kind": "text"},
    "peptide": {"label": "Peptide chain", "group": "ligand", "kind": "text"},
    "receptor_nucleic": {"label": "Nucleic acid receptor", "group": "ligand", "kind": "bool"},
    "interaction_package": {"label": "Interaction engine", "group": "analysis", "kind": "choice"},
    "binding": {
        "label": "Binding mode threshold (%)",
        "group": "analysis",
        "kind": "int",
        "min": 0,
        "max": 100,
    },
    "min_transition": {
        "label": "Minimal Markov transition (%)",
        "group": "analysis",
        "kind": "int",
        "min": 0,
        "max": 100,
    },
    "frames": {"label": "Last frame to analyze", "group": "analysis", "kind": "int", "min": 1},
    "cpu_count": {"label": "CPU cores", "group": "analysis", "kind": "int", "min": 1},
    "representative_frame": {
        "label": "Representative frame per binding mode",
        "group": "outputs",
        "kind": "bool",
    },
    "frame_rmsd": {"label": "RMSD between frames", "group": "outputs", "kind": "bool"},
    "generate_pml": {"label": "Pharmacophore .pml files", "group": "outputs", "kind": "bool"},
    "stable_water_analysis": {
        "label": "Stable water analysis",
        "group": "outputs",
        "kind": "bool",
    },
    "water_eps": {
        "label": "Water cluster size (Å)",
        "group": "outputs",
        "kind": "float",
        "min": 0.1,
    },
    "figure_type": {
        "label": "Figure format",
        "group": "outputs",
        "kind": "choice",
        "choices": ["png", "svg", "pdf", "jpg"],
    },
    "reference": {
        "label": "Reference PDB for residue numbering",
        "group": "advanced",
        "kind": "file",
        "ext": [".pdb"],
    },
    "dataframe": {
        "label": "Precomputed interactions CSV",
        "group": "advanced",
        "kind": "file",
        "ext": [".csv"],
    },
    "verbose": {"label": "Verbose log", "group": "advanced", "kind": "bool"},
}

# A checkbox needs no "(accepts true/false, yes/no, y/n)" hint.
_ACCEPTS_RE = re.compile(r"\s*\(accepts [^)]*\)")

GROUPS = [
    ("inputs", "Input files"),
    ("ligand", "Ligand"),
    ("analysis", "Analysis"),
    ("outputs", "Optional outputs"),
    ("advanced", "Advanced"),
]


def _actions():
    return {a.dest: a for a in build_parser()._actions if a.dest in FIELDS}


def _is_bool(action):
    return action.type is parse_bool_flag or isinstance(action, argparse._StoreTrueAction)


def option_schema():
    """Describe every form field: flag, default, help text and widget hints."""
    actions = _actions()
    fields = []
    for dest, meta in FIELDS.items():
        action = actions[dest]
        fields.append(
            {
                "dest": dest,
                "flag": action.option_strings[0],
                "default": action.default,
                "help": _ACCEPTS_RE.sub("", (action.help or "").replace("%%", "%")),
                "required": bool(action.required),
                "choices": meta.get("choices") or action.choices,
                **{k: v for k, v in meta.items() if k != "choices"},
            }
        )
    return {"fields": fields, "groups": [{"id": g, "label": label} for g, label in GROUPS]}


def _present(value):
    return value is not None and value != "" and value is not False


def _is_default(value, default):
    try:
        return float(value) == float(default)
    except (TypeError, ValueError):
        return str(value) == str(default)


def build_argv(values):
    """Turn form values into ``openmmdl analysis`` arguments, leaving defaults out."""
    actions = _actions()
    argv = []
    for dest in FIELDS:
        action, value = actions[dest], values.get(dest)
        if not _present(value):
            continue
        flag = action.option_strings[0]
        if isinstance(action, argparse._StoreTrueAction):
            argv.append(flag)
        elif _is_bool(action):
            if parse_bool_flag(value) and not action.default:
                argv += [flag, "True"]
        elif not _is_default(value, action.default):
            argv += [flag, str(value)]
    return argv


def format_command(argv):
    if os.name == "nt":
        return "openmmdl analysis " + subprocess.list2cmdline(argv)
    return "openmmdl analysis " + shlex.join(argv)


def _check_file(values, dest, errors, required=False):
    path = values.get(dest)
    if not path:
        if required:
            errors[dest] = "Required."
        return
    p = Path(path)
    if not p.is_file():
        errors[dest] = "File not found."
        return
    ext = FIELDS[dest]["ext"]
    if p.suffix.lower() not in ext:
        errors[dest] = f"Expected {' or '.join(ext)}."


def _check_number(values, dest, errors, cast):
    value = values.get(dest)
    if value is None or value == "":
        return
    try:
        number = cast(value)
    except (TypeError, ValueError):
        errors[dest] = "Not a number."
        return
    meta = FIELDS[dest]
    if "min" in meta and number < meta["min"]:
        errors[dest] = f"Must be at least {meta['min']}."
    elif "max" in meta and number > meta["max"]:
        errors[dest] = f"Must be at most {meta['max']}."


def _long_paths_enabled():
    if os.name != "nt":
        return True
    try:
        import winreg

        with winreg.OpenKey(
            winreg.HKEY_LOCAL_MACHINE, r"SYSTEM\CurrentControlSet\Control\FileSystem"
        ) as key:
            return winreg.QueryValueEx(key, "LongPathsEnabled")[0] == 1
    except OSError:
        return False


# The barcode figures land about 150 characters below the output folder
# (BindingModes_*/Barcodes/<ligand atoms>/<ligand atoms>_<type>_barcodes.png).
MAX_OUTPUT_DIR_LENGTH = 90


def validate(values, output_dir):
    """Return (errors per field, warnings) for a run request."""
    errors, warnings = {}, []
    _check_file(values, "topology", errors, required=True)
    _check_file(values, "trajectory", errors, required=True)
    for dest in ("ligand_sdf", "reference", "dataframe"):
        _check_file(values, dest, errors)
    for dest in ("binding", "min_transition", "frames", "cpu_count"):
        _check_number(values, dest, errors, int)
    _check_number(values, "water_eps", errors, float)

    if not values.get("ligand_name") and not values.get("peptide"):
        errors["ligand_name"] = "Give the ligand residue name, or a peptide chain."
    if values.get("interaction_package") == "prolif" and values.get("special_ligand"):
        errors["special_ligand"] = "Special ligands are not supported with ProLIF."
    if str(values.get("topology", "")).lower().endswith(".tpr"):
        warnings.append(
            "GROMACS inputs are converted to PDB/DCD by overwriting the input files. "
            "Run on copies of your .tpr and trajectory."
        )
    if values.get("stable_water_analysis"):
        warnings.append("Stable water analysis runs after the interaction analysis and adds time.")

    if not output_dir:
        errors["output_dir"] = "Choose an output folder."
    else:
        out = Path(output_dir)
        if len(str(out.resolve())) > MAX_OUTPUT_DIR_LENGTH and not _long_paths_enabled():
            warnings.append(
                "This output folder path is long. Windows limits paths to 260 characters and the "
                "analysis writes files up to about 150 characters deep inside it, so it can fail "
                "late in the run. Pick a shorter folder (for example C:\\openmmdl\\run1) or enable "
                "Windows long paths."
            )
        if out.exists() and not out.is_dir():
            errors["output_dir"] = "This is a file, not a folder."
        elif out.is_dir() and any(
            (out / name).exists()
            for name in ("missing_frames_filled.csv", "BindingModes_residue", "RMSD")
        ):
            warnings.append(
                "The output folder already holds analysis results; they will be overwritten."
            )
    return errors, warnings

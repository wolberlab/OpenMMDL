"""Local file browsing and input detection for the analysis GUI."""

from __future__ import annotations

import os
import shutil
from collections import Counter, OrderedDict
from pathlib import Path

# Where `openmmdl simulation` leaves analysis-ready topology/trajectory pairs, best first.
SIMULATION_LAYOUTS = [
    ("Final_Output/All_Atoms", "centered_top.pdb", "centered_traj.dcd", "All atoms, centered"),
    (
        "Final_Output/Prot_Lig",
        "prot_lig_top.pdb",
        "prot_lig_traj.dcd",
        "Protein and ligand, centered",
    ),
    (
        "MD_Postprocessing",
        "centered_old_coordinates_top.pdb",
        "centered_old_coordinates.dcd",
        "Postprocessing, centered",
    ),
]

AMINO_ACIDS = set(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
    "HID HIE HIP HSD HSE HSP CYX CYM ASH GLH LYN ACE NME NMA SEC PYL MSE".split()
)
NUCLEOTIDES = set("DA DT DG DC DU A U G C T RA RU RG RC DA5 DA3 DT5 DT3 DG5 DG3 DC5 DC3".split())
WATERS = set("HOH WAT SOL TIP3 TIP4 TIP5 T3P T4P T5P SPC OPC H2O".split())
IONS = set(
    "NA CL K MG CA ZN MN FE CU CO NI CD LI RB CS BR IOD F SOD CLA POT CAL MG2 ZN2 NA+ CL- K+".split()
)


def places():
    """Starting points for the file browser."""
    result = [{"label": "Home", "path": str(Path.home())}]
    cwd = os.getcwd()
    if Path(cwd) != Path.home():
        result.append({"label": "Working folder", "path": cwd})
    listdrives = getattr(os, "listdrives", None)
    if listdrives:
        result += [{"label": drive, "path": drive} for drive in listdrives()]
    elif os.name != "nt":
        result.append({"label": "/", "path": "/"})
    return result


def default_output_dir(name):
    """A fresh folder under ~/openmmdl_analysis, named after an input file."""
    base = Path.home() / "openmmdl_analysis"
    stem = Path(name.replace("\\", "/")).stem or "run"
    candidate, n = base / stem, 2
    while candidate.exists() and any(candidate.iterdir()):
        candidate, n = base / f"{stem}_{n}", n + 1
    return candidate


def save_upload(stream, name, output_dir):
    """Write a dropped file to <output_dir>/inputs; browsers never reveal its original path."""
    filename = Path(name.replace("\\", "/")).name
    if filename in ("", ".", ".."):
        raise ValueError("Missing file name.")
    target_dir = Path(output_dir).expanduser().resolve() / "inputs"
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / filename
    with open(target, "wb") as handle:
        shutil.copyfileobj(stream, handle, 1 << 20)
    return target


def list_dir(path):
    folder = Path(path).expanduser() if path else Path.home()
    folder = folder.resolve()
    if folder.is_file():
        folder = folder.parent
    entries = []
    for child in folder.iterdir():
        if child.name.startswith("."):
            continue
        try:
            is_dir = child.is_dir()
            size = None if is_dir else child.stat().st_size
        except OSError:
            continue
        entries.append({"name": child.name, "path": str(child), "dir": is_dir, "size": size})
    entries.sort(key=lambda e: (not e["dir"], e["name"].lower()))
    parent = folder.parent if folder.parent != folder else None
    return {
        "path": str(folder),
        "parent": str(parent) if parent else None,
        "entries": entries,
        "isResults": is_results_folder(folder),
    }


def is_results_folder(folder):
    folder = Path(folder)
    return any(
        (folder / name).exists()
        for name in (
            "missing_frames_filled.csv",
            "BindingModes_residue",
            "BindingModes_ligand",
            "df_all.csv",
        )
    )


def detect_simulation(folder):
    """Find topology/trajectory/ligand files in an `openmmdl simulation` output folder."""
    folder = Path(folder)
    candidates = []
    for sub, top, traj, label in SIMULATION_LAYOUTS:
        base = folder / sub
        if (base / top).is_file() and (base / traj).is_file():
            ligands = sorted(str(p) for p in base.glob("*.sdf"))
            if not ligands:
                ligands = sorted(str(p) for p in (folder / "Input_Files").glob("*.sdf"))
            candidates.append(
                {
                    "label": label,
                    "topology": str(base / top),
                    "trajectory": str(base / traj),
                    "ligands": ligands,
                }
            )
    return candidates


def _element(line):
    element = line[76:78].strip()
    if element:
        return element.upper()
    name = line[12:16].strip()
    return name.lstrip("0123456789")[:1].upper()


def inspect_topology(path):
    """Summarize a PDB topology: ligand candidates, chains and atom counts."""
    path = Path(path)
    if path.suffix.lower() != ".pdb":
        return {"supported": False, "ligands": [], "chains": []}

    residues = OrderedDict()
    with path.open("r", errors="replace") as handle:
        for line in handle:
            record = line[:6]
            if record.startswith("ENDMDL"):
                break
            if record not in ("ATOM  ", "HETATM"):
                continue
            key = (line[21], line[22:27].strip(), line[17:21].strip())
            entry = residues.setdefault(key, {"heavy": 0, "atoms": 0})
            entry["atoms"] += 1
            if _element(line) != "H":
                entry["heavy"] += 1

    chains = OrderedDict()
    ligands = OrderedDict()
    waters = ions = atoms = 0
    for (chain, _resid, resname), entry in residues.items():
        atoms += entry["atoms"]
        if resname in WATERS:
            waters += 1
        elif resname in IONS:
            ions += 1
        elif resname in AMINO_ACIDS or resname in NUCLEOTIDES:
            chains.setdefault(chain, Counter())["residues"] += 1
            if resname in NUCLEOTIDES:
                chains[chain]["nucleic"] += 1
        else:
            lig = ligands.setdefault(
                resname, {"name": resname, "copies": 0, "heavyAtoms": 0, "chains": set()}
            )
            lig["copies"] += 1
            lig["heavyAtoms"] = max(lig["heavyAtoms"], entry["heavy"])
            lig["chains"].add(chain.strip() or "-")

    ligand_list = sorted(ligands.values(), key=lambda l: -l["heavyAtoms"])
    for lig in ligand_list:
        lig["chains"] = sorted(lig["chains"])
    suggested = next((l["name"] for l in ligand_list if l["heavyAtoms"] >= 5), None)
    return {
        "supported": True,
        "atoms": atoms,
        "waters": waters,
        "ions": ions,
        "ligands": ligand_list,
        "suggestedLigand": suggested,
        "chains": [
            {"id": chain, "residues": c["residues"], "nucleic": c["nucleic"] > 0}
            for chain, c in chains.items()
        ],
        "sdfNearby": sorted(str(p) for p in path.parent.glob("*.sdf")),
    }

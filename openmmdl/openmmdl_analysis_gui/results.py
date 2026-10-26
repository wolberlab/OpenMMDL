"""Read an ``openmmdl analysis`` results folder into JSON-ready structures."""

from __future__ import annotations

import ast
import json
import pickle
import re
import struct
import tempfile
import threading
import zipfile
from collections import OrderedDict
from pathlib import Path

from openmmdl.openmmdl_analysis.openmmdlanalysis import parse_xyz
from openmmdl.openmmdl_analysis_gui.runner import RUN_FILE

SCHEMA_DIRS = [("residue", "BindingModes_residue"), ("ligand", "BindingModes_ligand")]
IMAGE_EXT = {".png", ".svg", ".jpg", ".jpeg", ".gif"}
# Bigger images (the Markov plots are ~18000 px wide) are linked, not shown inline.
INLINE_IMAGE_MAX_PX = 4000

_ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_lock = threading.Lock()
_frames_cache = OrderedDict()
_universe_cache = OrderedDict()


def _cached(cache, key, build, size=4):
    with _lock:
        if key in cache:
            cache.move_to_end(key)
            return cache[key]
    value = build()
    with _lock:
        cache[key] = value
        while len(cache) > size:
            cache.popitem(last=False)
    return value


def _rel(root, path):
    return Path(path).relative_to(root).as_posix()


def _png_size(path):
    with open(path, "rb") as handle:
        header = handle.read(24)
    if header[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">II", header[16:24])


def _read_run_file(root):
    path = root / RUN_FILE
    if path.is_file():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return None
    return None


def _first_hetatm_resname(pdb):
    if not pdb.is_file():
        return None
    with pdb.open("r", errors="replace") as handle:
        for line in handle:
            if line.startswith(("HETATM", "ATOM  ")):
                return line[17:21].strip()
    return None


def schema_dirs(root):
    found = [{"id": sid, "dir": name} for sid, name in SCHEMA_DIRS if (root / name).is_dir()]
    if not found and (root / "df_all.csv").is_file():
        found.append({"id": "default", "dir": ""})
    return found


def _read_rmsd(root):
    path = root / "RMSD" / "RMSD_over_time.csv"
    if not path.is_file():
        return None
    import pandas as pd

    df = pd.read_csv(path, sep=" ", index_col=0)
    return {
        "frames": [int(f) for f in df.index],
        "series": [
            {"name": str(name), "values": [None if v != v else float(v) for v in df[name]]}
            for name in df.columns
        ],
    }


class _PlainUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError(f"refusing to load {module}.{name}")


def _read_interacting_waters(path):
    try:
        with open(path, "rb") as handle:
            data = _PlainUnpickler(handle).load()
    except (OSError, pickle.UnpicklingError, EOFError):
        return []
    return sorted(int(x) for x in data) if isinstance(data, (list, tuple, set)) else []


def _stable_waters(root):
    result = []
    for folder in sorted(root.glob("stableWaters_clusterEps_*")):
        levels = []
        for level in sorted(folder.glob("clusterSize*"), key=lambda p: int(p.name[11:] or 0)):
            rep = level / "representative_waters.pdb"
            levels.append(
                {
                    "name": level.name,
                    "minFrames": int(level.name[11:] or 0),
                    "clusters": len(list(level.glob("cluster_*.pdb"))),
                    "representative": _rel(root, rep) if rep.is_file() else None,
                    "residues": (
                        _rel(root, level / "interacting_residues.csv")
                        if (level / "interacting_residues.csv").is_file()
                        else None
                    ),
                }
            )
        stable_csv = folder / "stable_waters.csv"
        result.append(
            {
                "dir": folder.name,
                "eps": folder.name.split("_")[-1],
                "levels": levels,
                "csv": _rel(root, stable_csv) if stable_csv.is_file() else None,
            }
        )
    return result


def _log_messages(root):
    path = root / "openmmdl.log"
    if not path.is_file():
        return None
    messages = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = _ANSI_RE.sub("", line).strip()
        level, _, text = line.partition(": ")
        if level in ("WARNING", "ERROR", "CRITICAL"):
            messages.append({"level": level.lower(), "text": text})
    return {"path": "openmmdl.log", "messages": messages}


def list_files(root):
    files = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name.startswith("."):
            continue
        rel = _rel(root, path)
        entry = {"path": rel, "size": path.stat().st_size, "ext": path.suffix.lower()}
        if entry["ext"] in IMAGE_EXT:
            size = _png_size(path) if entry["ext"] == ".png" else None
            entry["image"] = True
            entry["inline"] = size is None or max(size) <= INLINE_IMAGE_MAX_PX
            if size:
                entry["width"], entry["height"] = size
        files.append(entry)
    return files


FIGURE_EXT = IMAGE_EXT | {".pdf", ".eps", ".ps", ".tif", ".tiff"}


def zip_folder(root, sub="", figures_only=False):
    """ZIP of a results folder (or one of its subfolders) in a temporary file, rewound."""
    root = Path(root).resolve()
    base = (root / sub).resolve()
    if base != root and root not in base.parents:
        raise ValueError("Folder is outside the results folder.")
    if not base.is_dir():
        raise FileNotFoundError(str(base))
    handle = tempfile.TemporaryFile()
    with zipfile.ZipFile(handle, "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(base.rglob("*")):
            if path.is_file() and (not figures_only or path.suffix.lower() in FIGURE_EXT):
                archive.write(path, f"{root.name}/{_rel(root, path)}")
    handle.seek(0)
    return handle


def load_summary(root):
    root = Path(root).resolve()
    run = _read_run_file(root)
    values = (run or {}).get("values", {})
    schemas = schema_dirs(root)
    ligand = values.get("ligand_name") or _first_hetatm_resname(root / "lig.pdb")

    structure = None
    clouds = None
    waters = []
    for schema in schemas:
        vis = root / schema["dir"] / "Visualization"
        if (vis / "interacting_waters.pdb").is_file() and structure is None:
            traj = vis / "interacting_waters.dcd"
            structure = {
                "pdb": _rel(root, vis / "interacting_waters.pdb"),
                "traj": _rel(root, traj) if traj.is_file() else None,
            }
            if (vis / "interacting_waters.pkl").is_file():
                waters = _read_interacting_waters(vis / "interacting_waters.pkl")
        if (vis / "clouds.json").is_file() and clouds is None:
            clouds = _rel(root, vis / "clouds.json")
    if structure is None and (root / "complex.pdb").is_file():
        structure = {"pdb": "complex.pdb", "traj": None}
    if structure and structure["traj"]:
        structure["frames"] = trajectory_length(root, structure["pdb"], structure["traj"])

    rmsd_matrix = next(iter(sorted((root / "RMSD").glob("RMSD_between_the_frames.*"))), None)
    return {
        "root": str(root),
        "name": root.name,
        "run": run,
        "ligand": ligand,
        "peptide": values.get("peptide") or None,
        "threshold": int(values.get("binding") or 40),
        "minTransition": int(values.get("min_transition") or 1),
        "schemas": schemas,
        "rmsd": _read_rmsd(root),
        "rmsdMatrix": _rel(root, rmsd_matrix) if rmsd_matrix else None,
        "structure": structure,
        "clouds": clouds,
        "interactingWaters": waters,
        "stableWaters": _stable_waters(root),
        "log": _log_messages(root),
        "files": list_files(root),
    }


# --- binding-mode features (df_all.csv) ---------------------------------------------------

TYPE_ORDER = [
    "hydrophobic",
    "acceptor",
    "donor",
    "waterbridge",
    "pistacking",
    "pication",
    "saltbridge_ni",
    "saltbridge_pi",
    "halogen",
    "metal",
]


def _serials(tokens):
    serials = []
    for token in tokens:
        for part in token.split(","):
            if not part.isdigit():
                return []
            serials.append(int(part))
    return serials


def parse_feature_name(name, kind):
    """Split a df_all feature column into partner, ligand atom serials and group.

    Mirrors the naming in BindingModeProcesser and the ligand-atom grouping in
    BarcodePlotter.
    """
    tokens = name.split("_")
    group = None
    if kind == "hydrophobic":
        atoms = tokens[1:-1]
    elif kind in ("acceptor", "donor", "waterbridge"):
        atoms = tokens[1:-2]
    elif kind == "halogen":
        atoms, group = tokens[1:2], tokens[2] if len(tokens) > 3 else None
    elif kind == "pistacking":
        atoms = tokens[1:-1]
    elif kind == "pication":
        atoms, group = tokens[1:-2], tokens[-2]
    elif kind in ("saltbridge_ni", "saltbridge_pi"):
        atoms, group = tokens[1:-3], tokens[-3]
    elif kind == "metal":
        atoms, group = tokens[1:2], tokens[2] if len(tokens) > 3 else None
    else:
        atoms = []
    return {"partner": tokens[0], "ligandAtoms": _serials(atoms), "group": group}


def _runs(flags):
    """Run-length encode a 0/1 sequence as [start, length, start, length, ...]."""
    runs, start = [], None
    for i, flag in enumerate(flags):
        if flag and start is None:
            start = i
        elif not flag and start is not None:
            runs += [start, i - start]
            start = None
    if start is not None:
        runs += [start, len(flags) - start]
    return runs


def _int_or_none(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _load_frames(df_path):
    import pandas as pd

    from openmmdl.openmmdl_analysis.visualization.barcodes import BarcodeGenerator

    df = pd.read_csv(df_path, low_memory=False)
    df = df.drop(columns=[c for c in df.columns if c.startswith("Unnamed")])
    df["FRAME"] = pd.to_numeric(df["FRAME"], errors="coerce")
    df = df.dropna(subset=["FRAME"])
    df["FRAME"] = df["FRAME"].astype(int)

    kinds = {}
    for kind, columns in BarcodeGenerator(df).interactions.items():
        for column in columns:
            kinds.setdefault(column, kind)
    feature_cols = [c for c in df.columns if c in kinds]
    for column in feature_cols:
        df[column] = pd.to_numeric(df[column], errors="coerce").fillna(0).astype("int8")
    return df, feature_cols, kinds


def _frames_df(root, schema_dir):
    path = Path(root) / schema_dir / "df_all.csv"
    return _cached(_frames_cache, (str(path), path.stat().st_mtime), lambda: _load_frames(path))


def _top_modes(folder):
    path = folder / "top_10_binding_modes.csv"
    if not path.is_file():
        return None
    import pandas as pd

    df = pd.read_csv(path)
    return [
        {
            "mode": row["Binding Mode"],
            "firstFrame": _int_or_none(row["First Frame"]),
            "frames": [int(f) for f in ast.literal_eval(str(row["All Frames"]))],
            "representativeFrame": _int_or_none(row["Representative Frame"]),
            "percent": float(row["Percentage Occurrence"]),
        }
        for _, row in df.iterrows()
    ]


def _schema_figures(root, folder):
    states = folder / "Binding_Modes_Markov_States"
    individual = {}
    for path in sorted((states / "individual_figures").glob("*_merged.png")):
        individual[path.name[: -len("_merged.png")]] = _rel(root, path)
    markov = []
    for path in states.glob("markov_chain_plot_*.*"):
        threshold = path.stem.rsplit("_", 1)[-1]
        if threshold.isdigit():
            markov.append({"threshold": int(threshold), "path": _rel(root, path)})
    numbering = next(iter(sorted(folder.glob("ligand_numbering.*"))), None)
    arranged = states / "all_binding_modes_arranged.png"
    pdbs = {p.stem: _rel(root, p) for p in states.glob("Binding_Mode_*.pdb")}
    return {
        "arranged": _rel(root, arranged) if arranged.is_file() else None,
        "individual": individual,
        "markov": sorted(markov, key=lambda m: m["threshold"]),
        "ligandNumbering": _rel(root, numbering) if numbering else None,
        "modePdbs": pdbs,
    }


def load_schema(root, schema_dir):
    root = Path(root).resolve()
    folder = root / schema_dir
    df, feature_cols, kinds = _frames_df(root, schema_dir)
    frames = sorted(df["FRAME"].unique().tolist())
    presence = df.groupby("FRAME")[feature_cols].max().reindex(frames, fill_value=0)

    features = []
    for column in feature_cols:
        flags = presence[column].tolist()
        rows = df.loc[df[column] == 1]
        residues = []
        if {"RESNR", "RESTYPE", "RESCHAIN"} <= set(df.columns):
            seen = set()
            for resnr, resname, chain in rows[["RESNR", "RESTYPE", "RESCHAIN"]].itertuples(
                index=False
            ):
                key = (_int_or_none(resnr), str(resname), str(chain))
                if key[0] is None or key in seen:
                    continue
                seen.add(key)
                residues.append({"resnr": key[0], "resname": key[1], "chain": key[2]})
        features.append(
            {
                "id": column,
                "type": kinds[column],
                **parse_feature_name(column, kinds[column]),
                "residues": residues,
                "count": int(sum(flags)),
                "runs": _runs(flags),
            }
        )
    return {
        "frames": frames,
        "features": features,
        "topModes": _top_modes(folder),
        "figures": _schema_figures(root, folder),
    }


def frame_interactions(root, schema_dir, frame):
    """Interaction geometry of one frame: ligand/protein (and water) points per row."""
    df, feature_cols, _ = _frames_df(Path(root).resolve(), schema_dir)
    rows = df.loc[df["FRAME"] == int(frame)]
    result = []
    for _, row in rows.iterrows():
        active = [c for c in feature_cols if row[c] == 1]
        if not active:
            continue
        lig = parse_xyz(row.get("LIGCOO")) or parse_xyz(row.get("TARGETCOO"))
        prot = parse_xyz(row.get("PROTCOO")) or parse_xyz(row.get("METALCOO"))
        if lig is None or prot is None:
            continue
        result.append(
            {
                "features": active,
                "ligand": lig,
                "partner": prot,
                "water": parse_xyz(row.get("WATERCOO")),
            }
        )
    return result


# --- 2D ligand depiction ---------------------------------------------------------------------


def _ligand_mol(root, ligand):
    import MDAnalysis as mda
    from rdkit import Chem

    universe = mda.Universe(str(root / "complex.pdb"))
    atoms = universe.select_atoms(f"resname {ligand}")
    if not len(atoms):
        raise ValueError(f"No residue named {ligand} in complex.pdb")
    first = atoms.residues[0].atoms
    serials = [int(i) for i in first.ids]
    try:
        mol = first.convert_to("RDKIT")
    except Exception:
        block = "".join(_pdb_lines(root / "complex.pdb", set(serials)))
        mol = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
        if mol is None:
            raise
        for atom, serial in zip(mol.GetAtoms(), serials):
            atom.SetIntProp("_serial", serial)
        return mol
    for atom in mol.GetAtoms():
        atom.SetIntProp("_serial", serials[atom.GetIntProp("_MDAnalysis_index")])
    return mol


def _pdb_lines(pdb, serials):
    with pdb.open("r", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and _int_or_none(line[6:11]) in serials:
                yield line


def ligand_depiction(root, ligand, width=460, height=340):
    """SVG of the ligand plus the drawing position of every heavy atom, keyed by serial."""
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    root = Path(root).resolve()
    mol = _ligand_mol(root, ligand)
    heavy = Chem.RWMol(mol)
    for atom in sorted(heavy.GetAtoms(), key=lambda a: -a.GetIdx()):
        if atom.GetAtomicNum() == 1:
            heavy.RemoveAtom(atom.GetIdx())
    mol = heavy.GetMol()
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        mol.UpdatePropertyCache(strict=False)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.08
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    atoms = []
    for atom in mol.GetAtoms():
        point = drawer.GetDrawCoords(atom.GetIdx())
        atoms.append(
            {
                "serial": atom.GetIntProp("_serial"),
                "symbol": atom.GetSymbol(),
                "x": round(point.x, 2),
                "y": round(point.y, 2),
            }
        )
    return {
        "svg": drawer.GetDrawingText().replace("svg:", ""),
        "width": width,
        "height": height,
        "atoms": atoms,
    }


# --- ProLIF ligand network --------------------------------------------------------------------

PROLIF_NAMES = {
    "hydrophobic": "Hydrophobic",
    "acceptor": "HBAcceptor",
    "donor": "HBDonor",
    "waterbridge": "WaterBridge",
    "pistacking": "PiStacking",
    "saltbridge_ni": "Anionic",
    "saltbridge_pi": "Cationic",
    "halogen": "XBDonor",
    "metal": "MetalAcceptor",
}

# Must be in <head> so vis sizes its canvas to the frame from the start.
_NETWORK_STYLE = """
<style>
  html, body { margin: 0; height: 100%; overflow: hidden; background: #fff; }
  #mynetwork { width: 100% !important; height: 100vh !important; }
  #networklegend { display: none; }
</style>
"""

# Lets the GUI toggle interactions in place and hear clicks on residues and edges.
_NETWORK_BRIDGE = """
<script>
(function () {
  var net = ifp, E = net.body.data.edges, N = net.body.data.nodes;
  var OFF = { color: "#e2e8f0", opacity: 0.6 };
  function apply(active, emphasis) {
    var eu = [], nu = [], live = {};
    E.forEach(function (e) {
      if (!e.oaFeature) return;
      var on = active.has(e.oaFeature);
      if (on) live[e.to] = true;
      var w = emphasis.has(e.oaFeature) ? e.oaWidth * 1.8 + 2 : e.oaWidth;
      eu.push({ id: e.id, color: on ? { color: e.oaColor } : OFF, width: on ? w : 1,
                font: { color: on ? "#0f172a" : "#cbd5e1" } });
    });
    E.update(eu);
    N.forEach(function (n) {
      if (n.group !== "protein") return;
      var on = !!live[n.id];
      nu.push({ id: n.id, color: on ? n.oaColor : "#f1f5f9",
                font: { color: on ? n.oaFont : "#94a3b8" } });
    });
    N.update(nu);
  }
  E.forEach(function (e) {
    if (e.oaFeature) E.update({ id: e.id, oaColor: e.color, oaWidth: Math.max(e.width, 2) });
  });
  N.forEach(function (n) {
    if (n.group === "protein")
      N.update({ id: n.id, oaColor: n.color, oaFont: (n.font && n.font.color) || "black" });
  });
  // Stabilization can finish before this script runs; ProLIF then has already turned physics off.
  var settled = false;
  function settle() {
    if (settled) return;
    settled = true;
    net.setOptions({ physics: false });
    net.fit({ animation: false });
    parent.postMessage({ type: "oa-network-ready" }, "*");
  }
  net.once("stabilizationIterationsDone", settle);
  setTimeout(function () { if (net.physics && net.physics.physicsEnabled === false) settle(); }, 0);
  setTimeout(settle, 4000);
  window.addEventListener("resize", function () { net.fit({ animation: false }); });
  net.on("click", function (p) {
    var features = [];
    if (p.nodes.length) {
      E.forEach(function (e) { if (e.oaFeature && e.to === p.nodes[0]) features.push(e.oaFeature); });
    } else if (p.edges.length) {
      var e = E.get(p.edges[0]);
      if (e && e.oaFeature) features.push(e.oaFeature);
    }
    if (features.length) parent.postMessage({ type: "oa-network-toggle", features: features }, "*");
  });
  window.addEventListener("message", function (ev) {
    var d = ev.data || {};
    if (d.type === "oa-network-state") apply(new Set(d.active), new Set(d.emphasis || []));
    if (d.type === "oa-network-fit") net.fit({ animation: { duration: 300 } });
    if (d.type === "oa-network-export") {
      var url = null, error = null;
      try {
        var canvas = document.querySelector("#mynetwork canvas");
        url = canvas ? canvas.toDataURL("image/png") : null;
      } catch (e) {
        error = String(e);
      }
      parent.postMessage({ type: "oa-network-image", url: url, error: error }, "*");
    }
  });
})();
</script>
"""


def _prolif_interaction(kind, group):
    if kind == "pication":
        return "PiCation" if (group or "").lower() in ("aromatic", "pi") else "CationPi"
    return PROLIF_NAMES.get(kind)


def _mean_distance(rows):
    distances = []
    for lig, prot in rows[["LIGCOO", "PROTCOO"]].head(200).itertuples(index=False):
        a, b = parse_xyz(lig), parse_xyz(prot)
        if a and b:
            distances.append(sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5)
    return sum(distances) / len(distances) if distances else 0.0


_CDN_SCRIPT_RE = re.compile(
    r'<script[^>]*src="https://unpkg\.com/vis-network[^"]*"[^>]*></script>'
)
_CDN_STYLE_RE = re.compile(r'<link[^>]*href="https://unpkg\.com/vis-network[^"]*"[^>]*/?>')


def lignetwork_html(root, schema_dir, ligand, feature_ids, vis_url):
    """ProLIF LigNetwork of the given df_all features; edges carry their feature id.

    ``vis_url`` replaces ProLIF's CDN copy of vis-network, so the view also works offline.
    """
    import pandas as pd
    from prolif.plotting.network import LigNetwork
    from prolif.residue import ResidueId

    root = Path(root).resolve()
    mol = _ligand_mol(root, ligand)
    index_of = {atom.GetIntProp("_serial"): atom.GetIdx() for atom in mol.GetAtoms()}
    lig_resid = str(ResidueId.from_atom(mol.GetAtomWithIdx(0)))
    df, _, kinds = _frames_df(root, schema_dir)
    n_frames = df["FRAME"].nunique()

    records, edge_ids = [], []
    for fid in feature_ids:
        kind = kinds.get(fid)
        parsed = parse_feature_name(fid, kind) if kind else None
        interaction = _prolif_interaction(kind, parsed["group"]) if parsed else None
        atoms = (
            tuple(index_of[s] for s in parsed["ligandAtoms"] if s in index_of) if parsed else ()
        )
        if not interaction or not atoms:
            continue
        if kind in ("donor", "halogen"):
            # LigNetwork draws donors from the second atom of the tuple (the hydrogen).
            atoms = (atoms[0], atoms[0])
        rows = df.loc[df[fid] == 1]
        for (resnr, resname, chain), group in rows.groupby(
            ["RESNR", "RESTYPE", "RESCHAIN"], dropna=False
        ):
            number = _int_or_none(resnr)
            if number is None:
                continue
            chain = None if chain != chain or chain in ("", "skip") else str(chain)
            protein = str(ResidueId(str(resname), number, chain))
            records.append(
                {
                    "ligand": lig_resid,
                    "protein": protein,
                    "interaction": interaction,
                    "atoms": atoms,
                    "weight": group["FRAME"].nunique() / n_frames,
                    "distance": _mean_distance(group),
                    "components": "ligand_protein",
                }
            )
            edge_ids.append((f"{fid}|{protein}", fid))
    if not records:
        raise ValueError("None of these interactions can be drawn on the ligand.")

    class FeatureNetwork(LigNetwork):
        def _make_interactions(self, mass=2):
            super()._make_interactions(mass)
            for edge, (edge_id, fid) in zip(self.edges, edge_ids):
                edge["id"] = edge_id
                edge["oaFeature"] = fid

    table = pd.DataFrame(records).set_index(["ligand", "protein", "interaction", "atoms"])
    network = FeatureNetwork(table[["weight", "distance", "components"]], mol)
    html = network._get_html(width="100%", height="100%", fontsize=18, show_interaction_data=True)
    html = _CDN_SCRIPT_RE.sub(f'<script type="text/javascript" src="{vis_url}"></script>', html)
    html = _CDN_STYLE_RE.sub("", html)
    html = html.replace("</head>", _NETWORK_STYLE + "</head>")
    return html.replace("</body>", _NETWORK_BRIDGE + "</body>")


# --- trajectory frames ------------------------------------------------------------------------


def _universe(root, pdb_rel, traj_rel):
    import MDAnalysis as mda

    pdb, traj = root / pdb_rel, root / traj_rel
    key = (str(pdb), str(traj), pdb.stat().st_mtime, traj.stat().st_mtime)
    return _cached(_universe_cache, key, lambda: mda.Universe(str(pdb), str(traj)), size=2)


def trajectory_length(root, pdb_rel, traj_rel):
    return len(_universe(Path(root), pdb_rel, traj_rel).trajectory)


def frame_coordinates(root, pdb_rel, traj_rel, index):
    """Coordinates of one frame as little-endian float32 xyz triples."""
    universe = _universe(Path(root).resolve(), pdb_rel, traj_rel)
    with _lock:
        universe.trajectory[int(index)]
        positions = universe.atoms.positions.astype("<f4")
    return positions.tobytes()

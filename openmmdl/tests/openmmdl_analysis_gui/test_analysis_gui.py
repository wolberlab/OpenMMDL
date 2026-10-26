import json
import pickle
import threading
import time
import zipfile
from pathlib import Path

import pytest

from openmmdl.openmmdl_analysis.cli import build_parser
from openmmdl.openmmdl_analysis_gui import filesystem, options, results
from openmmdl.openmmdl_analysis_gui.runner import AnalysisRun
from openmmdl.openmmdl_analysis_gui.server import app

DATA_IN = Path(__file__).resolve().parents[1] / "data" / "in"

FEATURES = [
    "166ARGA_4220,4221_Carboxylate_NI_saltbridge",
    "165ASPA_4203_Acceptor_hbond",
    "60GLUA_4206_4207_4216_4217_4218_4205_hydrophobic",
]


def _pdb_line(record, serial, name, resname, chain, resid, element, xyz=(0.0, 0.0, 0.0)):
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:<3} {chain}{resid:>4}    "
        f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          {element:>2}\n"
    )


@pytest.fixture
def results_folder(tmp_path):
    """A minimal results folder: df_all.csv for one schema, RMSD, a run file and the log."""
    schema = tmp_path / "BindingModes_residue"
    schema.mkdir()
    rows = []
    presence = {
        1: [FEATURES[0], FEATURES[1]],
        2: [FEATURES[0]],
        3: [FEATURES[0], FEATURES[2]],
        4: [],
    }
    for frame, present in presence.items():
        if not present:
            rows.append({"FRAME": frame, "INTERACTION": "skip"})
        for feature in present:
            resnr, restype = (166, "ARG") if "ARG" in feature else (165, "ASP")
            if "GLU" in feature:
                resnr, restype = 60, "GLU"
            row = {
                "FRAME": frame,
                "INTERACTION": feature.rsplit("_", 1)[1],
                "RESNR": resnr,
                "RESTYPE": restype,
                "RESCHAIN": "A",
                "LIGCOO": "(1.0, 2.0, 3.0)",
                "PROTCOO": "(np.float64(4.0), np.float64(6.0), np.float64(3.0))",
            }
            row.update({f: int(f == feature) for f in FEATURES})
            rows.append(row)
    columns = ["FRAME", "INTERACTION", "RESNR", "RESTYPE", "RESCHAIN", "LIGCOO", "PROTCOO"]

    def csv_line(cells):
        return ",".join(f'"{c}"' if "," in c else c for c in cells) + "\n"

    with open(schema / "df_all.csv", "w") as handle:
        handle.write(csv_line([""] + columns + FEATURES))
        for i, row in enumerate(rows):
            handle.write(csv_line([str(i)] + [str(row.get(c, 0)) for c in columns + FEATURES]))

    (tmp_path / "RMSD").mkdir()
    (tmp_path / "RMSD" / "RMSD_over_time.csv").write_text(
        'frame backbone protein "resname UNK"\n0 0.0 0.0 0.0\n1 1.25 1.5 2.0\n'
    )
    (tmp_path / "missing_frames_filled.csv").write_text("FRAME\n1\n")
    (tmp_path / "openmmdl.log").write_text("INFO: fine\nWARNING: \x1b[1mcareful\x1b[0m\n")
    (tmp_path / "analysis_gui_run.json").write_text(
        json.dumps(
            {"command": "openmmdl analysis", "values": {"ligand_name": "UNK", "binding": 30}}
        )
    )
    return tmp_path


def test_option_schema_follows_the_cli_parser():
    schema = options.option_schema()
    fields = {f["dest"]: f for f in schema["fields"]}
    assert fields["topology"]["flag"] == "-t"
    assert fields["topology"]["required"] is True
    assert fields["binding"]["default"] == 40
    assert "accepts" not in fields["frame_rmsd"]["help"]
    assert {g["id"] for g in schema["groups"]} >= {"inputs", "ligand", "analysis"}


def test_build_argv_leaves_out_defaults_and_round_trips():
    values = {
        "topology": "top.pdb",
        "trajectory": "traj.dcd",
        "ligand_name": "UNK",
        "binding": 40,
        "water_eps": 1,
        "representative_frame": True,
        "frame_rmsd": False,
        "verbose": True,
        "min_transition": 5,
    }
    argv = options.build_argv(values)
    assert argv[:6] == ["-t", "top.pdb", "-d", "traj.dcd", "-n", "UNK"]
    assert "-b" not in argv and "--watereps" not in argv and "-r" not in argv
    parsed = build_parser().parse_args(argv)
    assert parsed.representative_frame is True
    assert parsed.verbose is True
    assert int(parsed.min_transition) == 5


def test_validate_reports_missing_inputs_and_conflicts(tmp_path):
    errors, _ = options.validate({}, "")
    assert {"topology", "trajectory", "ligand_name", "output_dir"} <= set(errors)

    top = tmp_path / "top.pdb"
    top.write_text("END\n")
    errors, _ = options.validate(
        {
            "topology": str(top),
            "trajectory": str(tmp_path / "missing.dcd"),
            "peptide": "B",
            "interaction_package": "prolif",
            "special_ligand": "HEM",
            "binding": 140,
        },
        str(tmp_path / "out"),
    )
    assert "topology" not in errors and "ligand_name" not in errors
    assert errors["trajectory"] == "File not found."
    assert "special_ligand" in errors and "binding" in errors


def test_validate_warns_about_long_output_paths(tmp_path, monkeypatch):
    monkeypatch.setattr(options, "_long_paths_enabled", lambda: False)
    _, warnings = options.validate({}, str(tmp_path / ("x" * 120)))
    assert any("260 characters" in w for w in warnings)


def test_inspect_topology_suggests_the_ligand(tmp_path):
    lines = [
        _pdb_line("ATOM", 1, "CA", "ALA", "A", 1, "C"),
        _pdb_line("ATOM", 2, "CA", "GLY", "A", 2, "C"),
    ]
    lines += [_pdb_line("HETATM", 3 + i, f"C{i}", "LIG", "A", 10, "C") for i in range(6)]
    lines.append(_pdb_line("HETATM", 9, "H1", "LIG", "A", 10, "H"))
    lines.append(_pdb_line("HETATM", 10, "O", "HOH", "W", 20, "O"))
    lines.append(_pdb_line("HETATM", 11, "NA", "NA", "I", 30, "NA"))
    pdb = tmp_path / "system.pdb"
    pdb.write_text("".join(lines) + "END\n")
    (tmp_path / "lig.sdf").write_text("")

    info = filesystem.inspect_topology(pdb)
    assert info["suggestedLigand"] == "LIG"
    assert info["ligands"][0]["heavyAtoms"] == 6
    assert info["waters"] == 1 and info["ions"] == 1
    assert info["chains"] == [{"id": "A", "residues": 2, "nucleic": False}]
    assert info["sdfNearby"] == [str(tmp_path / "lig.sdf")]


def test_detect_simulation_finds_final_output(tmp_path):
    base = tmp_path / "Final_Output" / "All_Atoms"
    base.mkdir(parents=True)
    for name in ("centered_top.pdb", "centered_traj.dcd", "UNK.sdf"):
        (base / name).write_text("")
    found = filesystem.detect_simulation(tmp_path)
    assert len(found) == 1
    assert found[0]["topology"].endswith("centered_top.pdb")
    assert found[0]["ligands"] == [str(base / "UNK.sdf")]


def test_default_output_dir_skips_used_folders(tmp_path, monkeypatch):
    monkeypatch.setattr(filesystem.Path, "home", lambda: tmp_path)
    first = filesystem.default_output_dir("C:\\data\\system.pdb")
    assert first == tmp_path / "openmmdl_analysis" / "system"
    (first / "inputs").mkdir(parents=True)
    assert filesystem.default_output_dir("system.pdb").name == "system_2"


def test_api_upload_writes_into_the_inputs_folder(client, tmp_path):
    reply = client.post(
        "/api/upload",
        query_string={"name": "..\\evil/../top.pdb", "outputDir": str(tmp_path)},
        data=b"ATOM\n",
    )
    body = reply.get_json()
    assert reply.status_code == 200
    assert Path(body["path"]) == tmp_path / "inputs" / "top.pdb"
    assert (tmp_path / "inputs" / "top.pdb").read_bytes() == b"ATOM\n"


def test_list_dir_hides_dotfiles_and_lists_folders_first(results_folder):
    (results_folder / ".hidden").write_text("")
    listing = filesystem.list_dir(str(results_folder))
    names = [e["name"] for e in listing["entries"]]
    assert ".hidden" not in names
    assert names.index("RMSD") < names.index("openmmdl.log")
    assert listing["isResults"] is True


@pytest.mark.parametrize(
    "name, kind, expected",
    [
        (FEATURES[0], "saltbridge_ni", ("166ARGA", [4220, 4221], "Carboxylate")),
        (FEATURES[1], "acceptor", ("165ASPA", [4203], None)),
        (FEATURES[2], "hydrophobic", ("60GLUA", [4206, 4207, 4216, 4217, 4218, 4205], None)),
        ("ligand_5_Acceptor_hbond", "acceptor", ("ligand", [5], None)),
        ("ligand_3GLY_Acceptor_hbond", "acceptor", ("ligand", [], None)),
        ("UNK_4210_ZN_4_metal", "metal", ("UNK", [4210], "ZN")),
    ],
)
def test_parse_feature_name(name, kind, expected):
    parsed = results.parse_feature_name(name, kind)
    assert (parsed["partner"], parsed["ligandAtoms"], parsed["group"]) == expected


def test_runs_encode_presence():
    assert results._runs([1, 1, 0, 1, 0, 0, 1]) == [0, 2, 3, 1, 6, 1]
    assert results._runs([0, 0]) == []


def test_load_schema_decodes_df_all(results_folder):
    data = results.load_schema(results_folder, "BindingModes_residue")
    assert data["frames"] == [1, 2, 3, 4]
    features = {f["id"]: f for f in data["features"]}
    salt = features[FEATURES[0]]
    assert salt["type"] == "saltbridge_ni"
    assert salt["count"] == 3 and salt["runs"] == [0, 3]
    assert salt["residues"] == [{"resnr": 166, "resname": "ARG", "chain": "A"}]
    assert features[FEATURES[1]]["runs"] == [0, 1]
    assert data["topModes"] is None


def test_frame_interactions_return_geometry(results_folder):
    rows = results.frame_interactions(results_folder, "BindingModes_residue", 1)
    assert len(rows) == 2
    assert rows[0]["ligand"] == [1.0, 2.0, 3.0] and rows[0]["partner"] == [4.0, 6.0, 3.0]
    assert results.frame_interactions(results_folder, "BindingModes_residue", 4) == []


def test_load_summary(results_folder):
    summary = results.load_summary(results_folder)
    assert summary["ligand"] == "UNK" and summary["threshold"] == 30
    assert summary["schemas"] == [{"id": "residue", "dir": "BindingModes_residue"}]
    assert [s["name"] for s in summary["rmsd"]["series"]] == ["backbone", "protein", "resname UNK"]
    assert summary["rmsd"]["series"][2]["values"] == [0.0, 2.0]
    assert summary["log"]["messages"] == [{"level": "warning", "text": "careful"}]


def test_interacting_waters_pickle_must_be_plain(tmp_path):
    plain = tmp_path / "plain.pkl"
    plain.write_bytes(pickle.dumps([5, 3]))
    assert results._read_interacting_waters(plain) == [3, 5]
    unsafe = tmp_path / "unsafe.pkl"
    unsafe.write_bytes(pickle.dumps([Path("x")]))
    assert results._read_interacting_waters(unsafe) == []


def test_ligand_depiction_maps_atoms_to_serials():
    depiction = results.ligand_depiction(DATA_IN, "UNK")
    serials = [a["serial"] for a in depiction["atoms"]]
    assert depiction["svg"].lstrip().startswith("<?xml") or "<svg" in depiction["svg"]
    assert serials and min(serials) >= 4187
    assert all(a["symbol"] != "H" for a in depiction["atoms"])


def test_frame_coordinates_are_float32_xyz():
    data = results.frame_coordinates(
        DATA_IN, "interacting_waters.pdb", "interacting_waters.dcd", 0
    )
    n_atoms = sum(
        1
        for line in open(DATA_IN / "interacting_waters.pdb")
        if line.startswith(("ATOM", "HETATM"))
    )
    assert len(data) == n_atoms * 3 * 4


def test_lignetwork_tags_edges_with_features(results_folder):
    pytest.importorskip("prolif")
    (results_folder / "complex.pdb").write_bytes((DATA_IN / "complex.pdb").read_bytes())
    serials = sorted(a["serial"] for a in results.ligand_depiction(DATA_IN, "UNK")["atoms"])
    df = results_folder / "BindingModes_residue" / "df_all.csv"
    df.write_text(df.read_text().replace("4220,4221", f"{serials[0]},{serials[1]}"))
    feature = FEATURES[0].replace("4220,4221", f"{serials[0]},{serials[1]}")
    html = results.lignetwork_html(
        results_folder, "BindingModes_residue", "UNK", [feature], "/static/vis-network.min.js"
    )
    assert '"oaFeature": "' + feature + '"' in html
    assert "unpkg.com" not in html and "/static/vis-network.min.js" in html


def test_zip_folder_packs_figures_or_everything(results_folder):
    (results_folder / "RMSD" / "RMSD_over_time.png").write_bytes(b"\x89PNG")
    with zipfile.ZipFile(results.zip_folder(results_folder, figures_only=True)) as archive:
        assert archive.namelist() == [f"{results_folder.name}/RMSD/RMSD_over_time.png"]
    with zipfile.ZipFile(results.zip_folder(results_folder, "RMSD")) as archive:
        assert sorted(archive.namelist()) == [
            f"{results_folder.name}/RMSD/RMSD_over_time.csv",
            f"{results_folder.name}/RMSD/RMSD_over_time.png",
        ]
    with pytest.raises(ValueError):
        results.zip_folder(results_folder, "..")


def test_api_zip_downloads_and_refuses_outside_paths(client, results_folder):
    ok = client.get("/api/results/zip", query_string={"root": str(results_folder), "kind": "all"})
    assert ok.status_code == 200 and ok.mimetype == "application/zip"
    escape = client.get(
        "/api/results/zip", query_string={"root": str(results_folder), "path": ".."}
    )
    assert escape.status_code == 403


def test_runner_splits_progress_bars_from_log_lines():
    run = AnalysisRun.__new__(AnalysisRun)
    run.lines, run.times, run.dropped, run.progress = [], [], 0, ""
    run.started = time.time()
    run._lock = threading.Lock()
    rest = run._consume(
        "INFO: \x1b[1mFiles are preprocessed\x1b[0m\n"
        "\rAnalyzing frames:   3%|3   | 3/99 [00:26<10:19,  6.45s/it]C:\\x.py:1: UserWarning: noisy\n"
        "\rAnalyzing frames:   4%|4   | 4/99 [00:30<10:00,  6.00s/it]\r"
        "\rAnalyzing frames:   5%"
    )
    assert run.lines == ["INFO: Files are preprocessed", "C:\\x.py:1: UserWarning: noisy"]
    assert run.progress.startswith("Analyzing frames:   4%")
    assert rest == "Analyzing frames:   5%"
    assert len(run.times) == len(run.lines)


@pytest.fixture
def client():
    app.config["TESTING"] = True
    return app.test_client()


def test_api_options_and_command(client):
    assert client.get("/api/options").get_json()["fields"]
    reply = client.post("/api/command", json={"values": {"ligand_name": "UNK"}, "outputDir": ""})
    body = reply.get_json()
    assert body["command"].startswith("openmmdl analysis") and "-n UNK" in body["command"]
    assert "topology" in body["errors"]


def test_api_rejects_folders_without_results(client, tmp_path):
    reply = client.get("/api/results/summary", query_string={"root": str(tmp_path)})
    assert reply.status_code == 400


def test_api_file_stays_inside_the_results_folder(client, results_folder):
    ok = client.get(
        "/api/results/file", query_string={"root": str(results_folder), "path": "openmmdl.log"}
    )
    assert ok.status_code == 200
    escape = client.get(
        "/api/results/file",
        query_string={"root": str(results_folder / "RMSD"), "path": "../openmmdl.log"},
    )
    assert escape.status_code == 403


def test_api_run_refuses_invalid_requests(client):
    reply = client.post("/api/runs", json={"values": {}, "outputDir": ""})
    assert reply.status_code == 400
    assert "topology" in reply.get_json()["errors"]

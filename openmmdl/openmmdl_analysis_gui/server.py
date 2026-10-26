"""Flask server behind ``openmmdl analysis-gui``."""

from __future__ import annotations

import os
import threading
import time
import webbrowser
from pathlib import Path
from urllib.parse import quote

from flask import Flask, Response, abort, jsonify, render_template, request, send_file, url_for

from openmmdl import __version__
from openmmdl.openmmdl_analysis_gui import filesystem, options, results
from openmmdl.openmmdl_analysis_gui.runner import AnalysisRun

SETUP_STATIC = Path(__file__).resolve().parent.parent / "openmmdl_setup" / "static"
# Reading results imports pyplot (via BarcodeGenerator); no GUI backend in a server thread.
os.environ.setdefault("MPLBACKEND", "Agg")

app = Flask(__name__)
runs = {}


class BadRequest(Exception):
    pass


@app.errorhandler(BadRequest)
def _bad_request(error):
    return jsonify({"error": str(error)}), 400


@app.errorhandler(FileNotFoundError)
def _not_found(error):
    return jsonify({"error": f"Not found: {error.filename or error}"}), 404


def _arg(name, required=True):
    value = request.args.get(name)
    if required and not value:
        raise BadRequest(f"Missing parameter '{name}'.")
    return value


def _results_root():
    root = Path(_arg("root")).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(str(root))
    return root


def _inside(root, rel):
    path = (root / rel).resolve()
    if root != path and root not in path.parents:
        abort(403)
    if not path.is_file():
        raise FileNotFoundError(rel)
    return path


@app.route("/")
def index():
    return render_template("index.html", version=__version__, initial=app.config.get("INITIAL"))


@app.route("/brand/<name>")
def brand(name):
    return send_file(_inside(SETUP_STATIC, name))


@app.route("/api/options")
def api_options():
    return jsonify(options.option_schema())


@app.route("/api/command", methods=["POST"])
def api_command():
    body = request.get_json(force=True)
    values, output_dir = body.get("values", {}), body.get("outputDir", "")
    errors, warnings = options.validate(values, output_dir)
    argv = options.build_argv(values)
    return jsonify(
        {
            "argv": argv,
            "command": options.format_command(argv),
            "errors": errors,
            "warnings": warnings,
        }
    )


@app.route("/api/fs/places")
def api_places():
    return jsonify(filesystem.places())


@app.route("/api/fs/list")
def api_list():
    path = request.args.get("path") or str(Path.home())
    try:
        return jsonify(filesystem.list_dir(path))
    except PermissionError:
        raise BadRequest(f"No permission to open {path}")


@app.route("/api/upload", methods=["POST"])
def api_upload():
    name = _arg("name")
    output_dir = request.args.get("outputDir") or str(filesystem.default_output_dir(name))
    try:
        target = filesystem.save_upload(request.stream, name, output_dir)
    except ValueError as error:
        raise BadRequest(str(error))
    return jsonify(
        {
            "path": str(target),
            "size": target.stat().st_size,
            "outputDir": str(Path(output_dir).expanduser().resolve()),
        }
    )


@app.route("/api/inspect/topology")
def api_inspect_topology():
    path = Path(_arg("path"))
    if not path.is_file():
        raise FileNotFoundError(str(path))
    return jsonify(filesystem.inspect_topology(path))


@app.route("/api/inspect/simulation")
def api_inspect_simulation():
    return jsonify(filesystem.detect_simulation(Path(_arg("path"))))


@app.route("/api/runs", methods=["GET", "POST"])
def api_runs():
    if request.method == "GET":
        return jsonify([dict(run.snapshot(), lines=[]) for run in runs.values()])
    body = request.get_json(force=True)
    values, output_dir = body.get("values", {}), body.get("outputDir", "")
    errors, warnings = options.validate(values, output_dir)
    if errors:
        return jsonify({"errors": errors, "warnings": warnings}), 400
    run = AnalysisRun(values, output_dir, copy_topology=body.get("copyTopology", True))
    runs[run.id] = run
    return jsonify(run.snapshot())


def _run(run_id):
    if run_id not in runs:
        raise FileNotFoundError(run_id)
    return runs[run_id]


@app.route("/api/runs/<run_id>")
def api_run(run_id):
    return jsonify(_run(run_id).snapshot(int(request.args.get("since", 0))))


@app.route("/api/runs/<run_id>/cancel", methods=["POST"])
def api_run_cancel(run_id):
    _run(run_id).cancel()
    return jsonify(_run(run_id).snapshot())


@app.route("/api/results/summary")
def api_results_summary():
    root = _results_root()
    if not filesystem.is_results_folder(root):
        raise BadRequest(f"{root} does not look like an OpenMMDL Analysis results folder.")
    return jsonify(results.load_summary(root))


@app.route("/api/results/schema")
def api_results_schema():
    return jsonify(results.load_schema(_results_root(), request.args.get("dir", "")))


@app.route("/api/results/frame-interactions")
def api_frame_interactions():
    root, frame = _results_root(), int(_arg("frame"))
    return jsonify(results.frame_interactions(root, request.args.get("dir", ""), frame))


@app.route("/api/results/ligand2d")
def api_ligand2d():
    try:
        return jsonify(results.ligand_depiction(_results_root(), _arg("ligand")))
    except ValueError as error:
        raise BadRequest(str(error))


@app.route("/api/results/lignetwork", methods=["POST"])
def api_lignetwork():
    body = request.get_json(force=True)
    root = Path(body.get("root", "")).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(str(root))
    try:
        html = results.lignetwork_html(
            root,
            body.get("dir", ""),
            body.get("ligand", ""),
            body.get("features", []),
            url_for("static", filename="vis-network.min.js", _external=True),
        )
    except ValueError as error:
        raise BadRequest(str(error))
    return Response(html, mimetype="text/html")


@app.route("/api/results/coords")
def api_coords():
    root = _results_root()
    data = results.frame_coordinates(root, _arg("pdb"), _arg("traj"), int(_arg("frame")))
    return Response(data, mimetype="application/octet-stream")


@app.route("/api/results/zip")
def api_zip():
    root = _results_root()
    sub = request.args.get("path", "").strip("/")
    figures = request.args.get("kind") == "figures"
    try:
        archive = results.zip_folder(root, sub, figures_only=figures)
    except ValueError:
        abort(403)
    name = "_".join(filter(None, [root.name, sub.replace("/", "_"), "figures" if figures else ""]))
    return send_file(
        archive, mimetype="application/zip", as_attachment=True, download_name=f"{name}.zip"
    )


@app.route("/api/results/file")
def api_file():
    root = _results_root()
    return send_file(_inside(root, _arg("path")), conditional=True)


def _warm_imports():
    # The first results request otherwise pays several seconds of cold imports.
    import MDAnalysis  # noqa: F401
    import matplotlib.pyplot  # noqa: F401
    import pandas  # noqa: F401
    from rdkit.Chem.Draw import rdMolDraw2D  # noqa: F401


def main(host="127.0.0.1", port=5050, open_browser=True, results_folder=None):
    threading.Thread(target=_warm_imports, daemon=True).start()
    url = f"http://{host}:{port}/"
    if results_folder:
        root = str(Path(results_folder).expanduser().resolve())
        app.config["INITIAL"] = {"results": root}
        url += "#/results?root=" + quote(root)
    if open_browser:

        def _open():
            time.sleep(1)
            webbrowser.open(url)

        threading.Thread(target=_open, daemon=True).start()
    print(f"OpenMMDL Analysis GUI running at {url}  (Ctrl+C to stop)")
    app.run(host=host, port=port, debug=False, threaded=True)

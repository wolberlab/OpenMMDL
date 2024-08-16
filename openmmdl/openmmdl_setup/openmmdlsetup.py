import openmm as mm
import openmm.unit as unit
from openmm.app import PDBFile, PDBxFile
from pdbfixer.pdbfixer import (
    PDBFixer,
    proteinResidues,
    dnaResidues,
    rnaResidues,
    _guessFileFormat,
)
from flask import (
    Flask,
    request,
    session,
    g,
    render_template,
    make_response,
    send_file,
    url_for,
)
from werkzeug.utils import secure_filename
from multiprocessing import Process, Pipe
import datetime
import os
import shutil
import signal
import sys
import tempfile
import threading
import time
import traceback
import webbrowser
import zipfile
from typing import Optional, Dict, Any, Tuple, List, Union

from openmmdl.openmmdl_setup.setup_options import (
    SetupOptionsConfigurator,
    RequestSessionManager,
)
from openmmdl.openmmdl_setup.configfile_creator import ConfigCreator, ConfigWriter
from openmmdl.openmmdl_setup.amberscript_creator import (
    AmberScriptGenerator,
    AmberBashScriptWriter,
)
from openmmdl.openmmdl_setup.file_operator import FileUploader

if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO


app = Flask(__name__)
app.config.from_object(__name__)
app.config.update({"SECRET_KEY": "development key"})
app.jinja_env.globals["mm"] = mm

# Typing for global variables
uploadedFiles: Dict[str, List[Tuple[tempfile.NamedTemporaryFile, str]]] = {}
fixer: Optional[PDBFixer] = None
scriptOutput: Optional[str] = None
simulationProcess: Optional[Process] = None


@app.route("/headerControls")
def headerControls() -> str:
    """
    Handle requests related to header controls such as restarting or stopping the server.
    """
    if "startOver" in request.args:
        return showSelectFileType()
    if "quit" in request.args:
        func = request.environ.get("werkzeug.server.shutdown")
        if func is None:
            raise RuntimeError("Not running with the Werkzeug Server")
        func()
        return "OpenMM Setup has stopped running. You can close this window."


@app.route("/")
def showSelectFileType() -> str:
    """
    Render the page to select the file type PDB or Amber for the simulation setup.
    """
    return render_template("selectFileType.html")


@app.route("/selectFiles")
def selectFiles() -> str:
    """
    Handle file type selection and configure files based on the selected type.
    """
    session["fileType"] = request.args.get("type", "")
    return showConfigureFiles()


def showConfigureFiles() -> str:
    """
    Render the appropriate configuration page based on the file type selected.
    """
    try:
        fileType = session["fileType"]
        if fileType == "pdb":
            return render_template("configurePdbFile.html")
        elif fileType == "amber":
            return render_template("configureAmberFiles.html")
    except:
        app.logger.error("Error displaying configure files page", exc_info=True)
    # The file type is invalid, so send them back to the select file type page.
    return showSelectFileType()


@app.route("/configureFiles", methods=["POST"])
def configureFiles() -> str:
    """
    Handles the configuration of files based on the file type and form data provided.
    """
    fileType = session["fileType"]
    if fileType == "pdb":
        if "file" not in request.files or request.files["file"].filename == "":
            # They didn't select a file. Send them back.
            return showConfigureFiles()

        FileUploader.save_uploaded_files(uploadedFiles, request)
        configurefile_request_session = RequestSessionManager(request.form)
        configurefile_request_session.configureFiles_add_forcefield_ligand_settings()
        session["sdfFile"] = uploadedFiles["sdfFile"][0][1]
        simulationoptions_session_manager = SetupOptionsConfigurator(session)
        simulationoptions_session_manager.configure_default_options()
        file, name = uploadedFiles["file"][0]
        file.seek(0, 0)
        session["pdbType"] = _guessFileFormat(file, name)
        if session["pdbType"] == "pdb":
            global fixer
            fixer = PDBFixer(pdbfile=file)
        return showSelectChains()
    elif fileType == "amber":
        session["has_files"] = request.form.get("has_files", "")
        has_files = session["has_files"]
        if has_files == "yes":
            if (
                "prmtopFile" not in request.files
                or request.files["prmtopFile"].filename == ""
                or "inpcrdFile" not in request.files
                or request.files["inpcrdFile"].filename == ""
            ):
                # if the user doesn't select prmtop or incprd file. Send them back.
                return showConfigureFiles()
            session["nmLig"] = "nmLig" in request.form
            session["spLig"] = "spLig" in request.form
            if session["nmLig"]:
                # If the user doesn't type the resname of the ligand. Send them back.
                if "nmLigName" not in request.form or request.form["nmLigName"] == "":
                    return showConfigureFiles()
                else:
                    session["nmLigName"] = request.form.get("nmLigName", "")
            if session["spLig"]:
                if "spLigName" not in request.form or request.form["spLigName"] == "":
                    return showConfigureFiles()
                else:
                    session["spLigName"] = request.form.get("spLigName", "")
            FileUploader.save_uploaded_files(uploadedFiles, request)
        elif has_files == "no":
            simulation_configurator = SetupOptionsConfigurator(session)
            simulation_configurator.configureDefaultAmberOptions()
            return showAmberOptions()

    simulation_configurator = SetupOptionsConfigurator(session)
    simulation_configurator.configure_default_options()
    return showSimulationOptions()


@app.route("/showAmberOptions")
def showAmberOptions() -> str:
    """
    Render the page for setting the Amber simulation options.
    """
    return render_template("AmberOptions.html")


@app.route("/setAmberOptions", methods=["POST"])
def setAmberOptions() -> str:
    """
    Handles the setting of Amber options and file uploads.
    """
    for key in request.form:
        session[key] = request.form[key]
    ######## Receptor ########
    form = request.form
    ######## Amber Receptor Sessions ########
    amber_session_manager = RequestSessionManager(form)
    amber_session_manager.setAmberOptions_rcp_session()
    # save uploaded pdb file for receptor
    rcpType = session["rcpType"]
    if rcpType == "protRcp":
        if "protFile" not in request.files or request.files["protFile"].filename == "":
            showAmberOptions()
        FileUploader.save_uploaded_files(uploadedFiles, request)
    elif rcpType == "dnaRcp":
        if "dnaFile" not in request.files or request.files["dnaFile"].filename == "":
            showAmberOptions()
        FileUploader.save_uploaded_files(uploadedFiles, request)
    elif rcpType == "rnaRcp":
        if "rnaFile" not in request.files or request.files["rnaFile"].filename == "":
            showAmberOptions()
        FileUploader.save_uploaded_files(uploadedFiles, request)
    elif rcpType == "carboRcp":
        if (
            "carboFile" not in request.files
            or request.files["carboFile"].filename == ""
        ):
            showAmberOptions()
        FileUploader.save_uploaded_files(uploadedFiles, request)

    ######## Ligand ########
    session["nmLig"] = (
        "nmLig" in request.form
    )  # store whether the nmLig checkbox is checked, e.g. True or False
    session["spLig"] = "spLig" in request.form
    # save uploaded pdb or sdf file for ligand
    ## for normal ligand
    if session["nmLig"]:
        if (
            "nmLigFile" not in request.files
            or request.files["nmLigFile"].filename == ""
        ):
            showAmberOptions()
        FileUploader.save_uploaded_files(uploadedFiles, request)

    ## for special ligand
    if session["spLig"]:
        if (
            "spLigFile" not in request.files
            or request.files["spLigFile"].filename == ""
            or "prepcFile" not in request.files
            or request.files["prepcFile"].filename == ""
            or "frcmodFile" not in request.files
            or request.files["frcmodFile"].filename == ""
        ):
            showAmberOptions()
        FileUploader.save_uploaded_files(uploadedFiles, request)

    ######## Add Water/Membrane ########

    amber_session_manager.setAmberOptions_water_membrane_session()
    amber_script_creator = AmberBashScriptWriter(session, uploadedFiles)
    return amber_script_creator.write_amber_script()


@app.route("/downloadAmberBashScript")
def downloadAmberBashScript() -> Any:
    """
    Downloads the amber bash script if it exists.
    """
    filename = "amber.bash"
    if os.path.exists(filename):
        return send_file(filename, as_attachment=True, download_name=filename)
    return "The amber.bash script does not exist."


@app.route("/downloadZipFile")
def downloadZipFile() -> Any:
    """
    Downloads the zip file containing the necessary simulation files.
    """
    if not os.path.exists("output_files.zip"):
        return "The zip file does not exist."
    return send_file("output_files.zip", as_attachment=True, download_name="output_files.zip")


@app.route("/submitSimulation")
def submitSimulation() -> str:
    """
    Starts the simulation process based on the configured options.
    """
    if simulationProcess:
        simulationProcess.terminate()

    # Wait for the termination
    time.sleep(1)

    # Create the simulation process
    pipe_parent, pipe_child = Pipe()
    global simulationProcess
    simulationProcess = Process(target=run_simulation, args=(pipe_child,))
    simulationProcess.start()

    # Wait for the process to finish and get the results
    results = pipe_parent.recv()
    if "success" in results:
        return "Simulation complete."
    return "Simulation failed."


def run_simulation(pipe_child) -> None:
    """
    Function to run the simulation process.
    """
    try:
        # Run the actual simulation logic here
        # ...
        pipe_child.send({"success": True})
    except Exception as e:
        pipe_child.send({"error": str(e)})


@app.route("/stopSimulation")
def stopSimulation() -> str:
    """
    Stops the ongoing simulation if it's running.
    """
    if simulationProcess:
        simulationProcess.terminate()
        return "Simulation stopped."
    return "No simulation is currently running."


if __name__ == "__main__":
    app.run(debug=True)

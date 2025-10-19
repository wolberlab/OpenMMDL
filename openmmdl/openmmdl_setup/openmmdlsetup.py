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
    render_template,
    make_response,
    send_file,
)
from werkzeug.utils import secure_filename
import datetime
import shutil
import sys
import tempfile
import threading
import time
import webbrowser
import zipfile
import warnings

warnings.filterwarnings("ignore")


if sys.version_info >= (3, 0):
    from io import StringIO
else:
    from cStringIO import StringIO


app = Flask(__name__)
app.config.from_object(__name__)
app.config.update({"SECRET_KEY": "development key"})
app.jinja_env.globals["mm"] = mm

uploadedFiles = {}
fixer = None
scriptOutput = None
simulationProcess = None


def saveUploadedFiles():
    uploadedFiles.clear()
    for key in request.files:
        filelist = []
        for file in request.files.getlist(key):
            temp = tempfile.TemporaryFile()
            shutil.copyfileobj(file, temp)
            filelist.append((temp, secure_filename(file.filename)))
        uploadedFiles[key] = filelist


@app.route("/headerControls")
def headerControls():
    if "startOver" in request.args:
        return showSelectFileType()
    if "quit" in request.args:
        func = request.environ.get("werkzeug.server.shutdown")
        if func is None:
            raise RuntimeError("Not running with the Werkzeug Server")
        func()
        return "OpenMM Setup has stopped running.  You can close this window."


@app.route("/")
def showSelectFileType():
    return render_template("selectFileType.html")


@app.route("/selectFiles")
def selectFiles():
    session["fileType"] = request.args.get("type", "")  # get the value of `type` from the url
    return showConfigureFiles()


def showConfigureFiles():
    try:
        fileType = session["fileType"]
        if fileType == "pdb":
            return render_template("configurePdbFile.html")
        elif fileType == "amber":
            return render_template("configureAmberFiles.html")
    except Exception:
        app.logger.error("Error displaying configure files page", exc_info=True)
    # The file type is invalid, so send them back to the select file type page.
    return showSelectFileType()


#################################################################################################
@app.route("/configureFiles", methods=["POST"])
def configureFiles():
    fileType = session["fileType"]
    if fileType == "pdb":
        if "file" not in request.files or request.files["file"].filename == "":
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        saveUploadedFiles()
        session["forcefield"] = request.form.get("forcefield", "")
        session["ml_forcefield"] = request.form.get("ml_forcefield", "")
        session["waterModel"] = request.form.get("waterModel", "")
        session["smallMoleculeForceField"] = request.form.get("smallMoleculeForceField", "")
        session["ligandMinimization"] = request.form.get("ligandMinimization", "")
        session["ligandSanitization"] = request.form.get("ligandSanitization", "")
        session["sdfFile"] = uploadedFiles["sdfFile"][0][1]
        configureDefaultOptions()
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
                # if the user doesn't select prmtop or incprd file.  Send them back.
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
            saveUploadedFiles()
        elif has_files == "no":
            configureDefaultAmberOptions()
            return showAmberOptions()

    configureDefaultOptions()
    return showSimulationOptions()


@app.route("/showAmberOptions")
def showAmberOptions():
    return render_template("AmberOptions.html")


@app.route("/setAmberOptions", methods=["POST"])
def setAmberOptions():
    for key in request.form:
        session[key] = request.form[key]
    ######## Receptor ########
    session["rcpType"] = request.form.get(
        "rcpType", ""
    )  # store the value of rcpType in session, e.g. protRcp, dnaRcp, rnaRcp, carboRcp
    session["prot_ff"] = request.form.get("prot_ff", "")
    session["other_prot_ff_input"] = request.form.get("other_prot_ff_input", "")
    session["dna_ff"] = request.form.get("dna_ff", "")
    session["other_dna_ff_input"] = request.form.get("other_dna_ff_input", "")
    session["rna_ff"] = request.form.get("rna_ff", "")
    session["other_rna_ff_input"] = request.form.get("other_rna_ff_input", "")
    session["carbo_ff"] = request.form.get("carbo_ff", "")
    session["other_carbo_ff_input"] = request.form.get("other_carbo_ff_input", "")
    # save uploaded pdb file for receptor
    rcpType = session["rcpType"]
    if rcpType == "protRcp":
        if "protFile" not in request.files or request.files["protFile"].filename == "":
            showAmberOptions()
        saveUploadedFiles()
    elif rcpType == "dnaRcp":
        if "dnaFile" not in request.files or request.files["dnaFile"].filename == "":
            showAmberOptions()
        saveUploadedFiles()
    elif rcpType == "rnaRcp":
        if "rnaFile" not in request.files or request.files["rnaFile"].filename == "":
            showAmberOptions()
        saveUploadedFiles()
    elif rcpType == "carboRcp":
        if "carboFile" not in request.files or request.files["carboFile"].filename == "":
            showAmberOptions()
        saveUploadedFiles()

    ######## Ligand ########
    session["nmLig"] = "nmLig" in request.form  # store whether the nmLig checkbox is checked, e.g. True or False
    session["spLig"] = "spLig" in request.form
    # save uploaded pdb or sdf file for ligand
    ## for normal ligand
    if session["nmLig"]:
        if "nmLigFile" not in request.files or request.files["nmLigFile"].filename == "":
            showAmberOptions()
        saveUploadedFiles()

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
        saveUploadedFiles()

    ######## Add Water/Membrane ########
    session["addType"] = request.form.get("addType", "")
    session["boxType"] = request.form.get("boxType", "")
    session["dist"] = request.form.get("dist", "")
    session["lipid_tp"] = request.form.get("lipid_tp", "")
    session["other_lipid_tp_input"] = request.form.get("other_lipid_tp_input", "")
    session["lipid_ratio"] = request.form.get("lipid_ratio", "")
    session["lipid_ff"] = request.form.get("lipid_ff", "")
    session["dist2Border"] = request.form.get("dist2Border", "")
    session["padDist"] = request.form.get("padDist", "")
    session["water_ff"] = request.form.get("water_ff", "")
    session["pos_ion"] = request.form.get("pos_ion", "")
    session["neg_ion"] = request.form.get("neg_ion", "")
    session["ionConc"] = request.form.get("ionConc", "")

    return createAmberBashScript()


@app.route("/downloadAmberBashScript")
def downloadAmberBashScript():
    response = make_response(createAmberBashScript())
    response.headers["Content-Disposition"] = 'attachment; filename="run_ambertools.sh"'
    return response


def configureDefaultAmberOptions():
    """Select default options based on the file format and force field."""
    # Ligand
    session["nmLig"] = ""
    session["spLig"] = ""
    session["lig_ff"] = "gaff2"
    session["charge_value"] = "0"
    session["charge_method"] = "bcc"

    # Receptor
    session["prot_ff"] = "ff19SB"
    session["dna_ff"] = "OL15"
    session["rna_ff"] = "OL3"
    session["carbo_ff"] = "GLYCAM_06j"

    # AddWaterMembrane
    session["addType"] = "addWater"
    session["boxType"] = "cube"
    session["dist"] = "10"

    session["lipid_tp"] = "POPC"
    session["other_lipid_tp_input"] = "POPC:TOPC"
    session["lipid_ratio"] = "1:1"
    session["lipid_ff"] = "lipid21"
    session["dist2Border"] = "15"
    session["padDist"] = "17"

    session["water_ff"] = "opc"
    session["pos_ion"] = "Na+"
    session["neg_ion"] = "Cl-"
    session["ionConc"] = "0.15"


def createAmberBashScript():
    a_script = []
    a_script.append("# This script was generated by OpenMMDL-Setup on %s.\n" % datetime.date.today())
    a_script.append(
        """
        #       ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      
        #     .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      
        #    / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      
        #   ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    
        #   |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    
        #   : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    
        #    \ `"/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  
        #     '. \_/``".'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ 
        #       '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` 
                                                                                                      
                                                                                                      
    """
    )

    a_script.append("#!/bin/bash\n")

    # Receptor
    a_script.append("################################## Receptor ######################################")
    rcpType = session["rcpType"]
    if rcpType == "protRcp":
        protFile = uploadedFiles["protFile"][0][1]
        protFile = protFile[:-4]
        a_script.append("rcp_nm=%s # the file name of ligand without suffix `pdb`" % protFile)

        prot_ff = session["prot_ff"]
        if prot_ff != "other_prot_ff":
            a_script.append("rcp_ff=%s" % session["prot_ff"])
        elif prot_ff == "other_prot_ff":
            a_script.append(
                "rcp_ff=%s  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                % session["other_prot_ff_input"]
            )
        a_script.append("\n")

    elif rcpType == "dnaRcp":
        dnaFile = uploadedFiles["dnaFile"][0][1]
        dnaFile = dnaFile[:-4]
        a_script.append("rcp_nm=%s # the file name of ligand without suffix `pdb`" % dnaFile)

        dna_ff = session["dna_ff"]
        if dna_ff != "other_dna_ff":
            a_script.append("rcp_ff=%s" % session["dna_ff"])
        elif dna_ff == "other_dna_ff":
            a_script.append(
                "rcp_ff=%s  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                % session["other_dna_ff_input"]
            )
        a_script.append("\n")

    elif rcpType == "rnaRcp":
        rnaFile = uploadedFiles["rnaFile"][0][1]
        rnaFile = rnaFile[:-4]
        a_script.append("rcp_nm=%s # the file name of ligand without suffix `pdb`" % rnaFile)

        rna_ff = session["rna_ff"]
        if rna_ff != "other_rna_ff":
            a_script.append("rcp_ff=%s" % session["rna_ff"])
        elif rna_ff == "other_rna_ff":
            a_script.append(
                "rcp_ff=%s  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                % session["other_rna_ff_input"]
            )
        a_script.append("\n")

    elif rcpType == "carboRcp":
        carboFile = uploadedFiles["carboFile"][0][1]
        carboFile = carboFile[:-4]
        a_script.append("rcp_nm=%s # the file name of ligand without suffix `pdb`" % carboFile)

        carbo_ff = session["carbo_ff"]
        if carbo_ff != "other_carbo_ff":
            a_script.append("rcp_ff=%s" % session["carbo_ff"])
        elif carbo_ff == "other_carbo_ff":
            a_script.append(
                "rcp_ff=%s  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                % session["other_carbo_ff_input"]
            )
        a_script.append("\n")

    a_script.append("## Clean the PDB file by pdb4amber")
    a_script.append(("pdb4amber -i ${rcp_nm}.pdb -o ${rcp_nm}_amber.pdb"))
    a_script.append(
        """
## `tleap` requires that all residues and atoms have appropriate types to ensure compatibility with the specified force field.
## To avoid `tleap` failing, we delete non-essential atoms, such as hydrogens, but preserve important atoms like carbon and nitrogen within the caps residues.
## Don' worry about the missing atoms as tleap has the capability to reconstruct them automatically. """
    )
    a_script.append(
        """awk '! ($2 ~ "(CH3|HH31|HH32|HH33)" || $3 ~ "(CH3|HH31|HH32|HH33)" )' ${rcp_nm}_amber.pdb > ${rcp_nm}_amber_f.pdb """
    )
    a_script.append("grep -v '^CONECT' ${rcp_nm}_amber_f.pdb > ${rcp_nm}_cnt_rmv.pdb\n")

    # Ligand
    if session["nmLig"] or session["spLig"]:
        a_script.append("################################## Ligand ######################################")
    if session["nmLig"]:
        a_script.append("# Normal Ligand that is compatible with GAFF force field")
        nmLigFile = uploadedFiles["nmLigFile"][0][1]
        a_script.append("nmLigFile=%s # the file name of ligand without suffix `.pdb` or `.sdf`" % nmLigFile[:-4])
        # depending on the uploaded file format,convert it to pdb or sdf file.
        if nmLigFile[-4:] == ".sdf":
            a_script.append(
                "obabel ${nmLigFile}.sdf -O ${nmLigFile}.pdb -p # convert to pdb file for tleap, -p: add hydrogens appropriate for pH7.4"
            )
        elif nmLigFile[-4:] == ".pdb":
            a_script.append(
                "obabel ${nmLigFile}.pdb -O ${nmLigFile}.sdf -p # convert to sdf file for openmmdl_analysis, -p: add hydrogens appropriate for pH7.4"
            )

        a_script.append(
            "charge_method=%s # refers to the charge method that antechamber will adopt" % session["charge_method"]
        )
        a_script.append(
            "charge_value=%s # Enter the net molecular charge of the ligand as integer (e.g. 1 or -2)"
            % session["charge_value"]
        )
        a_script.append("lig_ff=%s # Ligand force field \n" % session["lig_ff"])

        a_script.append("## Clean the PDB file by pdb4amber")
        a_script.append(("pdb4amber -i ${nmLigFile}.pdb -o ${nmLigFile}_amber.pdb\n"))

        a_script.append("## Generate a prepc file and an additional frcmod file by `antechamber`")
        a_script.append(
            "antechamber -fi pdb -fo prepc -i ${nmLigFile}_amber.pdb -o ${nmLigFile}.prepc -c ${charge_method} -at ${lig_ff} -nc ${charge_value} -pf y"
        )
        a_script.append("parmchk2 -f prepc -i ${nmLigFile}.prepc -o ${nmLigFile}.frcmod\n")
        a_script.append("## Rename ligand pdb")
        a_script.append("antechamber -i ${nmLigFile}.prepc -fi prepc -o rename_${nmLigFile}.pdb -fo pdb\n")

    if session["spLig"]:
        a_script.append("# Special Ligand that is incompatible with GAFF force field")
        spLigFile = uploadedFiles["spLigFile"][0][1]
        a_script.append("spLigFile=%s # the file name of ligand without suffix `.pdb`" % spLigFile[:-4])
        prepcFile = uploadedFiles["prepcFile"][0][1]
        prepcFile = prepcFile[:-6]
        a_script.append("prepc=%s # the file name without suffix `prepc`" % prepcFile)

        frcmodFile = uploadedFiles["frcmodFile"][0][1]
        frcmodFile = frcmodFile[:-7]
        a_script.append("frcmod=%s # the file name without suffix `frcmod`\n" % frcmodFile)

        a_script.append("## Clean the PDB file by pdb4amber")
        a_script.append(("pdb4amber -i ${spLigFile}.pdb -o ${spLigFile}_amber.pdb\n"))

        # get the name of ligand in the pdb file: it is the fourth column of the first line
        a_script.append("spLigName=$(awk 'NR==1 {print $4}' ${spLigFile}_amber.pdb)\n")

    # Combine all components to be modelled.
    if session["nmLig"] or session["spLig"]:
        a_script.append("######################  Combine All Components to Be Modelled ####################")
        a_script.append("cat > tleap.combine.in <<EOF\n")
        a_script.append("source ${rcp_ff}")
        a_script.append("source leaprc.${lig_ff}")
        ## load the prepc and frcmod file for either normal or special ligand
        if session["nmLig"]:
            a_script.append("\nloadamberprep ${nmLigFile}.prepc")
            a_script.append("loadamberparams ${nmLigFile}.frcmod\n")
        if session["spLig"]:
            a_script.append("loadamberprep ${prepc}.prepc")
            a_script.append("loadamberparams ${frcmod}.frcmod\n")
        ## load both receptor and ligand pdb file
        a_script.append("rcp = loadpdb ${rcp_nm}_cnt_rmv.pdb")
        if session["nmLig"] and session["spLig"]:
            a_script.append("nmLig = loadpdb rename_${nmLigFile}.pdb ")
            a_script.append("spLig = loadpdb ${spLigFile}_amber.pdb ")
            a_script.append("comp = combine{rcp nmLig spLig}")
        elif session["nmLig"]:
            a_script.append("nmLig = loadpdb rename_${nmLigFile}.pdb ")
            a_script.append("comp = combine{rcp nmLig}")
        elif session["spLig"]:
            a_script.append("spLig = loadpdb ${spLigFile}_amber.pdb")
            a_script.append("comp = combine {rcp spLig}")

        a_script.append("savepdb comp comp.pdb")
        a_script.append("\nquit")
        a_script.append("\nEOF\n")
        a_script.append("tleap -s -f tleap.combine.in > tleap.combine.out")
        ## remove 'CONECT' line in the pdb file
        a_script.append("grep -v '^CONECT' comp.pdb > comp_cnt_rmv.pdb\n")

    # Add Water/Membrane
    a_script.append("################################ Add Water/Membrane ##############################")

    ## Box setting
    addType = session["addType"]
    if addType == "addWater":
        boxType = session["boxType"]
        if boxType == "cube":
            a_script.append("boxType=solvatebox # `solvatebox`, a command in tleap, creates a cubic box ")
            a_script.append(
                "dist=%s # the minimum distance between any atom originally present in solute and the edge of the periodic box."
                % session["dist"]
            )
        elif boxType == "octahedron":
            a_script.append(
                "boxType=solvateoct # `solvateoct`, a command in tleap, creates a truncated octahedron box."
            )
            a_script.append(
                "dist=%s # the minimum distance between any atom originally present in solute and the edge of the periodic box"
                % session["dist"]
            )
        elif boxType == "cap":
            a_script.append(
                "boxType=solvatecap # `solvatecap`, a command in tleap, creates a solvent cap around solute. In development!"
            )
            a_script.append("radius=%s # the radius of the sphere" % session["dist"])
        elif boxType == "shell":
            a_script.append(
                "boxType=solvateshell # `solvatecap`, a command in tleap, adds a solent shell to solute, which reflect the contours of the original solute molecule. "
            )
            a_script.append("thickness=%s # the thickness of the shell" % session["dist"])
    elif addType == "addMembrane":
        lipid_tp = session["lipid_tp"]
        if lipid_tp != "other_lipid_tp":
            a_script.append("lipid_tp=%s" % session["lipid_tp"])
            a_script.append("lipid_ratio=1")
        elif lipid_tp == "other_lipid_tp":
            a_script.append(
                "lipid_tp=%s  # The command to check supported lipids: packmol-memgen --available_lipids"
                % session["other_lipid_tp_input"]
            )
            a_script.append(("lipid_ratio=%s # Set to 1 if only one lipid required" % session["lipid_ratio"]))

        lipid_ff = session["lipid_ff"]
        if lipid_ff != "other_lipid_ff":
            a_script.append("lipid_ff=%s" % session["lipid_ff"])
        elif lipid_ff == "other_lipid_ff":
            a_script.append(
                "lipid_ff=%s  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                % session["other_lipid_ff_input"]
            )

        a_script.append(
            "dist2Border=%s  # The minimum distance between the maxmin values for x y and z to the box boundaries. Flag --dist"
            % session["dist2Border"]
        )
        a_script.append(
            "padDist=%s  # The width of the water layer over the membrane or protein in the z axis. Flag --dist_wat"
            % session["padDist"]
        )

    ## Water Setting
    water_ff = session["water_ff"]
    if water_ff != "other_water_ff":
        a_script.append("water_ff=%s" % session["water_ff"])
        if addType == "addWater":
            water_ff = session["water_ff"]
            if water_ff == "tip3p":
                a_script.append("solvent=%sBOX  # set the water box" % session["water_ff"].upper())
            elif water_ff == "fb3":
                a_script.append("solvent=TIP3PFBOX # set the water box")
            elif water_ff == "spce":
                a_script.append("solvent=SPCBOX # set the water box")
            elif water_ff == "tip4pew":
                a_script.append("solvent=TIP4PEWBOX # set the water box")
            elif water_ff == "fb4":
                a_script.append("solvent=TIP4PBOX # set the water box")
            elif water_ff == "opc":
                a_script.append("solvent=OPCBOX # set the water box")
            elif water_ff == "opc3":
                a_script.append("solvent=OPC3BOX # set the water box")
    elif water_ff == "other_water_ff":
        a_script.append(
            "water_ff=%s  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
            % session["other_water_ff_input"]
        )
        if addType == "addWater":
            a_script.append("solvent=%sBOX  # set the water box" % session["other_water_ff_input"].upper())

    ## Ion Setting
    pos_ion = session["pos_ion"]
    if pos_ion != "other_pos_ion":
        a_script.append("pos_ion=%s" % session["pos_ion"])
    elif pos_ion == "other_pos_ion":
        a_script.append("pos_ion=%s  # In development!" % session["other_pos_ion_input"])

    neg_ion = session["neg_ion"]
    if neg_ion != "other_neg_ion":
        a_script.append("neg_ion=%s" % session["neg_ion"])
    elif neg_ion == "other_neg_ion":
        a_script.append("neg_ion=%s  # In development!" % session["other_neg_ion_input"])

    if addType == "addWater":
        a_script.append(
            "numIon=0 # `numIon` is the flag for `addions` in tleap. When set to 0, the system will be neutralized"
        )
        a_script.append("\n")
    elif addType == "addMembrane":
        a_script.append("ionConc=%s" % session["ionConc"])
        a_script.append("\n")

    ## Build the membrane
    if addType == "addMembrane":
        a_script.append("## Build the membrane")
        if not session["nmLig"] and not session["spLig"]:
            a_script.append(
                "packmol-memgen --pdb ${rcp_nm}_cnt_rmv.pdb --lipids ${lipid_tp} --ratio ${lipid_ratio} --preoriented --dist ${dist2Border} --dist_wat ${padDist} --salt --salt_c ${pos_ion} --saltcon ${ionConc} --nottrim --overwrite --notprotonate\n"
            )
            a_script.append("## Clean the complex pdb by `pdb4amber` for further `tleap` process")
            a_script.append("pdb4amber -i bilayer_${rcp_nm}_cnt_rmv.pdb -o clean_bilayer_${rcp_nm}.pdb")
            ## remove 'CONECT' line in the pdb file
            a_script.append("grep -v '^CONECT' clean_bilayer_${rcp_nm}.pdb > clean_bilayer_${rcp_nm}_cnt_rmv.pdb")
            a_script.append("\n")
        if session["nmLig"] or session["spLig"]:
            a_script.append(
                "packmol-memgen --pdb comp.pdb --lipids ${lipid_tp} --ratio ${lipid_ratio} --preoriented --dist ${dist2Border} --dist_wat ${padDist} --salt --salt_c ${pos_ion} --saltcon ${ionConc} --nottrim --overwrite --notprotonate\n"
            )
            a_script.append("## Clean the complex pdb by `pdb4amber` for further `tleap` process")
            a_script.append("pdb4amber -i bilayer_comp.pdb -o clean_bilayer_comp.pdb")
            ## remove 'CONECT' line in the pdb file
            a_script.append("grep -v '^CONECT' clean_bilayer_comp.pdb > clean_bilayer_comp_cnt_rmv.pdb")
            a_script.append("\n")

    # Generate the prmtop and frcmod file for the complex.
    a_script.append("##################### Generate Prmtop and Frcmod File for the Complex ###################### ")
    a_script.append("cat > tleap.in <<EOF\n")
    ## source the force field
    a_script.append("source ${rcp_ff}")
    a_script.append("source leaprc.water.${water_ff}")
    if session["nmLig"] or session["spLig"]:
        a_script.append("source leaprc.${lig_ff}")
    if addType == "addMembrane":
        a_script.append("source leaprc.${lipid_ff}")
    ## load the prepc and frcmod file
    if session["nmLig"]:
        a_script.append("\nloadamberprep ${nmLigFile}.prepc")
        a_script.append("loadamberparams ${nmLigFile}.frcmod\n")
    if session["spLig"]:
        a_script.append("loadamberprep ${prepc}.prepc")
        a_script.append("loadamberparams ${frcmod}.frcmod\n")
    ## load the complex pdb which comtains all the components to be modelled.
    if addType == "addWater":
        if not session["nmLig"] and not session["spLig"]:
            a_script.append("\nsystem = loadpdb ${rcp_nm}_cnt_rmv.pdb\n ")
        else:
            a_script.append("system = loadpdb comp_cnt_rmv.pdb\n")
    elif addType == "addMembrane":
        if not session["nmLig"] and not session["spLig"]:
            a_script.append("\nsystem = loadpdb clean_bilayer_${rcp_nm}_cnt_rmv.pdb\n")
        if session["nmLig"] or session["spLig"]:
            a_script.append("system = loadpdb clean_bilayer_comp_cnt_rmv.pdb\n")
    ## add box
    if addType == "addWater":
        if boxType == "cube":
            a_script.append("solvatebox system ${solvent} ${dist} ")
        elif boxType == "octahedron":
            a_script.append("solvateoct system ${solvent} ${dist}")
        elif boxType == "cap":
            a_script.append("solvatecap system ${solvent} ${radius}")
        elif boxType == "shell":
            a_script.append("solvateshell system ${solvent} ${thickness}")
        a_script.append("addions2 system ${neg_ion} ${numIon}")
        a_script.append("addions2 system ${pos_ion} ${numIon}")
    elif addType == "addMembrane":
        a_script.append('setBox system "vdw"')
    a_script.append("check system ")
    a_script.append("charge system\n")
    ## save pdb, prmtop and inpcrd file.
    a_script.append("savepdb system system.${water_ff}.pdb")
    a_script.append("saveamberparm system system.${water_ff}.prmtop system.${water_ff}.inpcrd")
    a_script.append("\nquit")
    a_script.append("\nEOF")
    a_script.append("\ntleap -s -f tleap.in > tleap.out")

    return "\n".join(a_script)


def extractLigName(LigFileName):
    """
    Extract the ligand name from the pdb file as a string, which is the fourth column of the first line in the pdb file.
    This string can be used for openmmdl_analysis in later function `createScript`.

    Params:
    -------
    LigFile: the tuple that stores both the buffered file and its name, uploadedFiles['nmLigFile'][0][1]
    """

    if LigFileName[-4:] == ".sdf":
        LigName = "UNL"
    elif LigFileName[-4:] == ".pdb":
        LigName = LigFileName[:-4]

    return LigName


########################################################################################################################


@app.route("/getCurrentStructure")
def getCurrentStructure():
    pdb = StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, pdb)
    return pdb.getvalue()


def showSelectChains():
    chains = []
    hasHeterogen = False
    for chain in fixer.topology.chains():
        residues = list(r.name for r in chain.residues())
        if any(r in proteinResidues for r in residues):
            content = "Protein"
        elif any(r in rnaResidues for r in residues):
            content = "RNA"
        elif any(r in dnaResidues for r in residues):
            content = "DNA"
        else:
            content = ", ".join(set(residues))
            hasHeterogen = True
        chains.append((chain.id, len(residues), content))
    if len(chains) < 2 and not hasHeterogen:
        session["heterogens"] = "all"
        return showAddResidues()
    return render_template("selectChains.html", chains=chains)


@app.route("/selectChains", methods=["POST"])
def selectChains():
    session["heterogens"] = request.form.get("heterogens", "")
    numChains = len(list(fixer.topology.chains()))
    request.form.getlist("include")
    deleteIndices = [i for i in range(numChains) if str(i) not in request.form.getlist("include")]
    fixer.removeChains(deleteIndices)
    return showAddResidues()


def showAddResidues():
    spans = []
    chains = list(fixer.topology.chains())
    fixer.findMissingResidues()
    if len(fixer.missingResidues) == 0:
        return showConvertResidues()
    for i, key in enumerate(sorted(fixer.missingResidues)):
        residues = fixer.missingResidues[key]
        chain = chains[key[0]]
        chainResidues = list(chain.residues())
        if key[1] < len(chainResidues):
            offset = int(chainResidues[key[1]].id) - len(residues) - 1
        else:
            offset = int(chainResidues[-1].id)
        spans.append((chain.id, offset + 1, offset + len(residues), ", ".join(residues)))
    return render_template("addResidues.html", spans=spans)


@app.route("/addResidues", methods=["POST"])
def addResidues():
    keys = [key for key in sorted(fixer.missingResidues)]
    for i, key in enumerate(keys):
        if str(i) not in request.form.getlist("add"):
            del fixer.missingResidues[key]
    return showConvertResidues()


def showConvertResidues():
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) == 0:
        return showAddHeavyAtoms()
    residues = []
    nucleotides = ["DA", "DC", "DG", "DT", "A", "C", "G", "T"]
    for i in range(len(fixer.nonstandardResidues)):
        residue, replaceWith = fixer.nonstandardResidues[i]
        if replaceWith in proteinResidues:
            replacements = proteinResidues
        else:
            replacements = nucleotides
        residues.append((residue.chain.id, residue.name, residue.id, replacements, replaceWith))
    return render_template("convertResidues.html", residues=residues)


@app.route("/convertResidues", methods=["POST"])
def convertResidues():
    for i in range(len(fixer.nonstandardResidues)):
        if str(i) in request.form.getlist("convert"):
            fixer.nonstandardResidues[i] = (
                fixer.nonstandardResidues[i][0],
                request.form["residue" + str(i)],
            )
    fixer.replaceNonstandardResidues()
    return showAddHeavyAtoms()


def showAddHeavyAtoms():
    if session["heterogens"] == "none":
        fixer.removeHeterogens(False)
    elif session["heterogens"] == "water":
        fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    allResidues = list(set(fixer.missingAtoms.keys()).union(fixer.missingTerminals.keys()))
    allResidues.sort(key=lambda x: x.index)
    if len(allResidues) == 0:
        return addHeavyAtoms()
    residues = []
    for residue in allResidues:
        atoms = []
        if residue in fixer.missingAtoms:
            atoms.extend(atom.name for atom in fixer.missingAtoms[residue])
        if residue in fixer.missingTerminals:
            atoms.extend(atom for atom in fixer.missingTerminals[residue])
        residues.append((residue.chain.id, residue.name, residue.id, ", ".join(atoms)))
    return render_template("addHeavyAtoms.html", residues=residues)


@app.route("/addHeavyAtoms", methods=["POST"])
def addHeavyAtoms():
    fixer.addMissingAtoms()
    return showAddHydrogens()


def showAddHydrogens():
    unitCell = fixer.topology.getUnitCellDimensions()
    if unitCell is not None:
        unitCell = unitCell.value_in_unit(unit.nanometer)
    boundingBox = tuple(
        (max((pos[i] for pos in fixer.positions)) - min((pos[i] for pos in fixer.positions))).value_in_unit(
            unit.nanometer
        )
        for i in range(3)
    )
    return render_template("addHydrogens.html", unitCell=unitCell, boundingBox=boundingBox)


@app.route("/addHydrogens", methods=["POST"])
def addHydrogens():
    session["solvent"] = False
    if "addHydrogens" in request.form:
        pH = float(request.form.get("ph", "7"))
        fixer.addMissingHydrogens(pH)
    if "addWater" in request.form:
        session["solvent"] = True
        session["add_membrane"] = False
        if request.form["boxType"] == "geometry":
            session["water_padding"] = True
            session["water_padding_distance"] = float(request.form["geomPadding"])
            session["water_boxShape"] = request.form["geometryDropdown"]
        else:
            session["water_padding"] = False
            session["box_x"] = float(request.form["boxx"])
            session["box_y"] = float(request.form["boxy"])
            session["box_z"] = float(request.form["boxz"])
        session["water_ionicstrength"] = float(request.form["ionicstrength"])
        session["water_positive"] = request.form["positiveion"] + "+"
        session["water_negative"] = request.form["negativeion"] + "-"
    elif "addMembrane" in request.form:
        session["solvent"] = True
        session["add_membrane"] = True
        session["lipidType"] = request.form["lipidType"]
        session["membrane_padding"] = float(request.form["membranePadding"])
        session["membrane_ionicstrength"] = float(request.form["ionicstrength"])
        session["membrane_positive"] = request.form["positiveion"] + "+"
        session["membrane_negative"] = request.form["negativeion"] + "-"

    # Save the new PDB file.

    uploadedFiles["originalFile"] = uploadedFiles["file"]
    pdb = StringIO()
    if session["pdbType"] == "pdb":
        try:
            PDBFile.writeFile(fixer.topology, fixer.positions, pdb, True)
        except Exception:
            # This can happen if the ids are too large to fit in the allowed space.
            pdb = StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, pdb, False)
    else:
        PDBxFile.writeFile(fixer.topology, fixer.positions, pdb, True)
    temp = tempfile.TemporaryFile()
    temp.write(pdb.getvalue().encode("utf-8"))
    name = uploadedFiles["file"][0][1]
    dotIndex = name.rfind(".")
    if dotIndex == -1:
        prefix = name
        suffix = ""
    else:
        prefix = name[:dotIndex]
        suffix = name[dotIndex:]
    uploadedFiles["file"] = [(temp, prefix + "-processed_openMMDL" + suffix)]
    return showSimulationOptions()


@app.route("/showSimulationOptions")
def showSimulationOptions():
    file_type = session.get("fileType", "")

    # render buttons based on the fileType
    if file_type == "pdb":
        return render_template(
            "simulationOptions.html",
            display_save_script=True,
            display_processed_pdb=True,
            display_save_all_files=True,
        )
    elif file_type == "amber":
        return render_template(
            "simulationOptions.html",
            display_save_script=True,
            display_processed_pdb=False,
            display_save_all_files=False,
        )


@app.route("/setSimulationOptions", methods=["POST"])
def setSimulationOptions():
    for key in request.form:
        session[key] = request.form[key]
    session["ligand"] = (
        "ligand" in request.form
    )  # store whether the ligand is present, so the retruned value can be 'True' or 'False'
    session["writeDCD"] = "writeDCD" in request.form
    session["writeData"] = "writeData" in request.form
    session["writeCheckpoint"] = "writeCheckpoint" in request.form
    session["dataFields"] = request.form.getlist("dataFields")
    session["hmr"] = "hmr" in request.form
    session["writeSimulationXml"] = "writeSimulationXml" in request.form
    session["writeFinalState"] = "writeFinalState" in request.form
    return createScript()


@app.route("/downloadScript")
def downloadScript():
    response = make_response(createScript())
    response.headers["Content-Disposition"] = 'attachment; filename="OpenMMDL_Simulation.py"'
    return response


@app.route("/downloadStructuralfiles")
def downloadStructuralfiles():
    file, name = uploadedFiles["file"][0]
    file.seek(0, 0)
    response = make_response(file.read())
    response.headers["Content-Disposition"] = 'attachment; filename="%s"' % name
    return response


@app.route("/downloadPackage")
def downloadPackage():
    temp = tempfile.NamedTemporaryFile()
    with zipfile.ZipFile(temp, "w", zipfile.ZIP_DEFLATED) as zip:
        zip.writestr("openmmdl_simulation/OpenMMDL_Simulation.py", createScript())
        for key in uploadedFiles:
            for file, name in uploadedFiles[key]:
                file.seek(0, 0)
                zip.writestr("openmmdl_simulation/%s" % name, file.read())
    temp.seek(0, 0)
    return send_file(temp, "application/zip", True, "openmmdl_simulation.zip", max_age=0)


def configureDefaultOptions():
    """Select default options based on the file format and force field."""
    implicitWater = False
    session["restart_checkpoint"] = "no"
    session["mdtraj_output"] = "mdtraj_pdb_dcd"
    session["mdtraj_removal"] = "False"
    session["mda_output"] = "mda_pdb_dcd"
    session["mda_selection"] = "mda_prot_lig_all"
    session["openmmdl_analysis"] = "No"
    session["analysis_selection"] = "analysis_all"
    session["binding_mode"] = "40"
    session["min_transition"] = "1"
    session["rmsd_diff"] = "No"
    session["pml_generation"] = "True"
    session["stable_water"] = "Yes"
    session["wc_distance"] = "1.0"
    if session["fileType"] == "pdb" and session["waterModel"] == "implicit":
        implicitWater = True
    session["ensemble"] = "nvt" if implicitWater else "npt"
    session["platform"] = "CUDA"
    session["precision"] = "mixed"
    session["cutoff"] = "2.0" if implicitWater else "1.0"
    session["ewaldTol"] = "0.0005"
    session["constraintTol"] = "0.000001"
    session["hmr"] = True
    session["hmrMass"] = "1.5"
    session["dt"] = "0.002"
    session["sim_length"] = "50"
    session["equilibration_length"] = "0.5"
    session["temperature"] = "300"
    session["friction"] = "1.0"
    session["pressure"] = "1.0"
    session["barostatInterval"] = "25"
    session["nonbondedMethod"] = "CutoffNonPeriodic" if implicitWater else "PME"
    session["writeDCD"] = True
    session["dcdFilename"] = "trajectory.dcd"
    session["dcdFrames"] = "1000"
    session["pdbInterval_ns"] = "10"
    session["writeData"] = True
    session["dataFilename"] = "log.txt"
    session["dataInterval"] = "1000"
    session["dataFields"] = [
        "step",
        "speed",
        "progress",
        "potentialEnergy",
        "temperature",
    ]
    session["writeCheckpoint"] = True
    session["checkpointFilename"] = "checkpoint.chk"
    session["checkpointInterval_ns"] = "0.02"
    session["writeSimulationXml"] = False
    session["systemXmlFilename"] = "system.xml"
    session["integratorXmlFilename"] = "integrator.xml"
    session["writeFinalState"] = False
    session["finalStateFileType"] = "stateXML"
    session["finalStateFilename"] = "final_state.xml"
    session["constraints"] = "hbonds"
    session["rmsd"] = "True"
    session["md_postprocessing"] = "True"


def createScript(isInternal=False):
    script = []

    # If we are creating this script for internal use to run a simulation directly, add extra code at the top
    # to set the working directory and redirect stdout to the pipe.

    if isInternal:
        script.append(
            """
import os
import sys
import time

class PipeOutput(object):
    def write(self, string):
        output.send(string)

sys.stdout = PipeOutput()
sys.stderr = PipeOutput()
os.chdir(outputDir)"""
        )

    # Header

    script.append("# This script was generated by OpenMM-MDL Setup on %s.\n" % datetime.date.today())
    script.append(
        """
#       ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      
#     .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      
#    / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      
#   ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    
#   |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    
#   : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    
#    \ `"/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  
#     '. \_/``".'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ 
#       '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` 
                                                                                                      
                                                                                                      
"""
    )
    script.append(
        "from openmmdl.openmmdl_simulation.scripts.forcefield_water import ff_selection, water_forcefield_selection, water_model_selection, generate_forcefield, generate_transitional_forcefield"
    )
    script.append(
        "from openmmdl.openmmdl_simulation.scripts.protein_ligand_prep import prepare_ligand, rdkit_to_openmm, merge_protein_and_ligand, water_padding_solvent_builder, water_absolute_solvent_builder, membrane_builder, water_conversion"
    )
    script.append(
        "from openmmdl.openmmdl_simulation.scripts.post_md_conversions import mdtraj_conversion, MDanalysis_conversion"
    )
    script.append(
        "from openmmdl.openmmdl_simulation.scripts.cleaning_procedures import cleanup, create_directory_if_not_exists, copy_file, organize_files, post_md_file_movement \n"
    )

    script.append("import simtk.openmm.app as app")
    script.append(
        "from simtk.openmm.app import PDBFile, Modeller, PDBReporter, StateDataReporter, DCDReporter, CheckpointReporter, AmberPrmtopFile, AmberInpcrdFile"
    )
    script.append("from simtk.openmm import unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator")
    script.append("from openmm.openmm import XmlSerializer")
    script.append("from simtk.openmm import Vec3")
    script.append("import simtk.openmm as mm")
    script.append("import pdbfixer")
    script.append("import sys")
    script.append("import os")
    script.append("import shutil")
    if session["openmmdl_analysis"] == "Yes":
        script.append("import subprocess")

    # Input files
    script.append("\n# Input Files")
    fileType = session["fileType"]
    if fileType == "pdb":
        script.append("""############# Ligand and Protein Data ###################""")
        script.append(
            """########   Add the Ligand SDF File and Protein PDB File in the Folder with the Script  ######### \n"""
        )
        pdbType = session["pdbType"]
        if pdbType == "pdb":
            script.append('protein = "%s"' % uploadedFiles["file"][0][1])
            if session["sdfFile"] != "":
                script.append("ligand = '%s'" % session["sdfFile"])
                script.append('ligand_name = "UNK"')
                script.append("minimization = %s" % session["ligandMinimization"])
                script.append("smallMoleculeForceField = '%s'" % session["smallMoleculeForceField"])
                script.append("sanitization = %s" % session["ligandSanitization"])
            water = session["waterModel"]
    elif fileType == "amber":
        script.append("""####### Add the Amber Files in the Folder with this Script ####### \n""")
        # amber_files related variables
        if session["has_files"] == "yes":
            script.append("prmtop_file = '%s'" % uploadedFiles["prmtopFile"][0][1])
            script.append('inpcrd_file = "%s"' % uploadedFiles["inpcrdFile"][0][1])

            # ligand related variables
            if session["nmLig"]:
                nmLigName = session["nmLigName"]  # e.g. 'UNL'
            else:
                nmLigName = None

            if session["spLig"]:  # success
                spLigName = session["spLigName"]  # e.g. 'HEME'
            else:
                spLigName = None

        elif session["has_files"] == "no":
            script.append("prmtop_file = 'system.%s.prmtop'" % session["water_ff"])
            script.append("inpcrd_file = 'system.%s.inpcrd' " % session["water_ff"])

            # ligand related variables
            if session["nmLig"]:
                nmLigFileName = uploadedFiles["nmLigFile"][0][1]  # e.g. '8QY.pdb'
                nmLigName = extractLigName(nmLigFileName)  # e.g '8QY' or 'UNL' # resname in topology
            else:
                nmLigFileName = None
                nmLigName = None

            if session["spLig"]:  # success
                spLigFileName = uploadedFiles["spLigFile"][0][1]
                spLigName = extractLigName(spLigFileName)
            else:
                spLigFileName = None
                spLigName = None

        ## Feed prmtop_file and inpcrd_file to OpenMM Reader
        script.append("prmtop = AmberPrmtopFile(prmtop_file)")
        script.append("inpcrd = AmberInpcrdFile(inpcrd_file)")

        ##################### For debugging #################################
        # print all key-value pairs in session
        # print(f"session is {session}")
        print(f"fileType is {fileType}")
        print(f"session['has_files'] is {session['has_files']}")
        print(f"session['nmLig'] is {session['nmLig']}")
        print(f"session['spLig'] is {session['spLig']}")
        print(f"nmLigName is {nmLigName}")
        print(f"spLigName is {spLigName}")
        # print(f"nmLigFileName is {nmLigFileName}")
        # print(f"nmLigName is {nmLigName}")
        # print(f"spLigFileName is {spLigFileName}")
        # print(f"spLigName is {spLigName}")
        ######################################################################

    if fileType == "pdb":
        script.append("""\n############# Forcefield, Water and Membrane Model Selection ###################\n""")
        script.append("ff = '%s'" % session["forcefield"])
        if water != "None":
            script.append("water = '%s'" % water)
        else:
            script.append("water = %s" % water)

    ################################## IF CLEANING WAS PERFORMED ##############################################
    ###########################################################################################################
    ###########################################################################################################
    if fileType == "pdb":
        if session["solvent"]:
            if session["add_membrane"]:
                script.append("""\n############# Membrane Settings ###################\n""")
                script.append("add_membrane = %s" % session["add_membrane"])
                script.append("membrane_lipid_type = '%s'" % session["lipidType"])
                script.append("membrane_padding = %s" % session["membrane_padding"])
                script.append("membrane_ionicstrength = %s" % session["membrane_ionicstrength"])
                script.append("membrane_positive_ion = '%s'" % session["membrane_positive"])
                script.append("membrane_negative_ion = '%s'" % session["membrane_negative"])
            elif not session["add_membrane"]:
                script.append("""\n############# Water Box Settings ###################\n""")
                script.append("add_membrane = %s" % session["add_membrane"])
                if session["water_padding"]:
                    script.append('Water_Box = "Buffer"')
                    script.append("water_padding_distance = %s" % session["water_padding_distance"])
                    script.append("water_boxShape = '%s'" % session["water_boxShape"])
                else:
                    script.append('Water_Box = "Absolute"')
                    script.append("water_box_x = %s" % session["box_x"])
                    script.append("water_box_y = %s" % session["box_y"])
                    script.append("water_box_z = %s" % session["box_z"])
                script.append("water_ionicstrength = %s" % session["water_ionicstrength"])
                script.append("water_positive_ion = '%s'" % session["water_positive"])
                script.append("water_negative_ion = '%s'" % session["water_negative"])
        else:
            if not session["solvent"]:
                script.append("Solvent = %s" % session["solvent"])

    ################################## IF CLEANING WAS NOT PERFORMED ##########################################
    ###########################################################################################################
    ###########################################################################################################

    # System configuration
    script.append("\n# System Configuration\n")
    nonbondedMethod = session["nonbondedMethod"]
    script.append("nonbondedMethod = app.%s" % nonbondedMethod)
    if nonbondedMethod != "NoCutoff":
        script.append("nonbondedCutoff = %s*unit.nanometers" % session["cutoff"])
    if nonbondedMethod == "PME":
        script.append("ewaldErrorTolerance = %s" % session["ewaldTol"])
    constraints = session["constraints"]
    constraintMethods = {
        "none": "None",
        "water": "None",
        "hbonds": "HBonds",
        "allbonds": "AllBonds",
    }
    if constraints != "none" and constraints != "water":
        script.append("constraints = app.%s" % constraintMethods[constraints])
    if constraints == "none":
        script.append("constraints = %s" % constraintMethods[constraints])
    script.append("rigidWater = %s" % ("False" if constraints == "none" else "True"))
    if constraints != "none":
        script.append("constraintTolerance = %s" % session["constraintTol"])
    if session["hmr"]:
        script.append("hydrogenMass = %s*unit.amu" % session["hmrMass"])

    # Integration options

    script.append("\n# Integration Options\n")
    script.append("step_time = %s" % session["dt"])
    script.append("dt = %s*unit.picoseconds" % session["dt"])
    script.append("temperature = %s*unit.kelvin" % session["temperature"])
    script.append("friction = %s/unit.picosecond" % session["friction"])
    ensemble = session["ensemble"]
    if ensemble == "npt":
        script.append("pressure = %s*unit.atmospheres" % session["pressure"])
        script.append("barostatInterval = %s" % session["barostatInterval"])

    # Simulation options

    script.append("\n# Simulation Options\n")
    script.append("sim_length = %s" % session["sim_length"])
    script.append("steps = int(sim_length / step_time * 1000)")
    script.extend(
        [
            "dcdFrames = %s" % session["dcdFrames"],
            "dcdInterval = int(steps / dcdFrames)",
        ]
    )
    script.extend(
        [
            "pdbInterval_ns = %s" % session["pdbInterval_ns"],
            "pdbInterval = int(steps * (pdbInterval_ns / sim_length))",
        ]
    )
    if session["restart_checkpoint"] == "yes":
        script.append("restart_step = %s" % session["restart_step"])
    script.append("equilibration_length = %s" % session["equilibration_length"])
    script.append("equilibrationSteps = int(equilibration_length / step_time * 1000)")
    script.append("platform = Platform.getPlatformByName('%s')" % session["platform"])
    if session["platform"] in ("CUDA", "OpenCL"):
        script.append("platformProperties = {'Precision': '%s'}" % session["precision"])
    if session["writeDCD"]:
        if session["restart_checkpoint"] == "yes":
            script.append(
                "dcdReporter = DCDReporter('%s_%s', dcdInterval)" % (session["restart_step"], session["dcdFilename"])
            )
        else:
            script.append("dcdReporter = DCDReporter('%s', dcdInterval)" % (session["dcdFilename"]))
    if session["writeData"]:
        args = ", ".join("%s=True" % field for field in session["dataFields"])
        if session["restart_checkpoint"] == "yes":
            script.append(
                "dataReporter = StateDataReporter('%s_%s', %s, totalSteps=steps,"
                % (
                    session["restart_step"],
                    session["dataFilename"],
                    session["dataInterval"],
                )
            )
        else:
            script.append(
                "dataReporter = StateDataReporter('%s', %s, totalSteps=steps,"
                % (session["dataFilename"], session["dataInterval"])
            )
        script.append("    %s, separator='\\t')" % args)
        if isInternal:
            # Create a second reporting sending to stdout so we can display it in the browser.
            script.append(
                "consoleReporter = StateDataReporter(sys.stdout, %s, totalSteps=steps, %s, separator='\\t')"
                % (session["dataInterval"], args)
            )
    if session["writeCheckpoint"]:
        script.append("checkpointInterval = int(1000 * %s / %s)" % (session["checkpointInterval_ns"], session["dt"]))
        if session["restart_checkpoint"] == "yes":
            script.append(
                "checkpointReporter = CheckpointReporter('%s_%s', checkpointInterval)"
                % (session["restart_step"], session["checkpointFilename"])
            )
            script.append(
                "checkpointReporter10 = CheckpointReporter('10x_%s__%s', checkpointInterval *10)"
                % (session["restart_step"], session["checkpointFilename"])
            )
            script.append(
                "checkpointReporter100 = CheckpointReporter('100x_%s_%s', checkpointInterval *100)"
                % (session["restart_step"], session["checkpointFilename"])
            )
        else:
            script.append(
                "checkpointReporter = CheckpointReporter('%s', checkpointInterval)" % (session["checkpointFilename"])
            )
            script.append(
                "checkpointReporter10 = CheckpointReporter('10x_%s', checkpointInterval* 10)"
                % (session["checkpointFilename"])
            )
            script.append(
                "checkpointReporter100 = CheckpointReporter('100x_%s', checkpointInterval* 100)"
                % (session["checkpointFilename"])
            )

    # Prepare the simulation

    if fileType == "pdb":
        if session["sdfFile"] != "":
            script.append(
                """
print("Preparing MD Simulation with ligand")
ligand_prepared = prepare_ligand(ligand,minimize_molecule=minimization)
omm_ligand = rdkit_to_openmm(ligand_prepared, ligand_name)
protein_pdb = pdbfixer.PDBFixer(str(protein))
forcefield_selected = ff_selection(ff)
water_selected = water_forcefield_selection(water=water,forcefield_selection=ff_selection(ff))
model_water = water_model_selection(water=water,forcefield_selection=ff_selection(ff))
print("Forcefield and Water Model Selected")
if add_membrane:
    transitional_forcefield = generate_transitional_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, smallMoleculeForceField=smallMoleculeForceField, rdkit_mol=ligand_prepared)
forcefield = generate_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, smallMoleculeForceField=smallMoleculeForceField, rdkit_mol=ligand_prepared)
complex_topology, complex_positions = merge_protein_and_ligand(protein_pdb, omm_ligand)
print("Complex topology has", complex_topology.getNumAtoms(), "atoms.")     """
            )
        elif session["sdfFile"] == "":
            script.append(
                """
protein_pdb = PDBFile(protein)     
forcefield_selected = ff_selection(ff)
water_selected = water_forcefield_selection(water=water,forcefield_selection=ff_selection(ff))
model_water = water_model_selection(water=water,forcefield_selection=ff_selection(ff))
print("Forcefield and Water Model Selected")
if water_selected != None:
    forcefield = generate_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, smallMoleculeForceField=smallMoleculeForceField, rdkit_mol=None) 
else:
    forcefield = app.ForceField(forcefield_selected)    
if add_membrane:
        transitional_forcefield = generate_transitional_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, smallMoleculeForceField=smallMoleculeForceField, rdkit_mol=None)     """
            )
        if session["sdfFile"] == "":
            script.append(
                """
forcefield = generate_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, smallMoleculeForceField=smallMoleculeForceField, rdkit_mol=None)        
modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
if add_membrane:
    membrane_builder(ff, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, membrane_lipid_type, membrane_padding, membrane_positive_ion, membrane_negative_ion, membrane_ionicstrength, protein)
elif not add_membrane:
    if Water_Box == "Buffer":
        water_padding_solvent_builder(model_water, forcefield, water_padding_distance, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
    elif Water_Box == "Absolute":
        water_absolute_solvent_builder(model_water, forcefield, water_box_x, water_box_y, water_box_z, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
if add_membrane:
    if model_water == 'tip4pew' or model_water == 'tip5p':
        water_conversion(model_water, modeller, protein)
topology = modeller.topology
positions = modeller.positions """
            )
        elif session["sdfFile"] != "":
            script.append(
                """
modeller = app.Modeller(complex_topology, complex_positions)
if add_membrane:
    membrane_builder(ff, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, membrane_lipid_type, membrane_padding, membrane_positive_ion, membrane_negative_ion, membrane_ionicstrength, protein)
elif not add_membrane:
    if Water_Box == "Buffer":
        water_padding_solvent_builder(model_water, forcefield, water_padding_distance, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
    elif Water_Box == "Absolute":
        water_absolute_solvent_builder(model_water, forcefield, water_box_x, water_box_y, water_box_z, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)  
if add_membrane:
    if model_water == 'tip4pew' or model_water == 'tip5p':
        water_conversion(model_water, modeller, protein)
topology = modeller.topology
positions = modeller.positions  """
            )
    elif fileType == "amber":
        script.append("topology = prmtop.topology")
        script.append("positions = inpcrd.positions")

    script.append("\n# Prepare the Simulation\n")
    script.append("print('Building system...')")
    hmrOptions = ", hydrogenMass=hydrogenMass" if session["hmr"] else ""
    if fileType == "pdb":
        script.append(
            "system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod,%s"
            % (" nonbondedCutoff=nonbondedCutoff," if nonbondedMethod != "NoCutoff" else "")
        )
        script.append(
            "    constraints=constraints, rigidWater=rigidWater%s%s)"
            % (
                (", ewaldErrorTolerance=ewaldErrorTolerance" if nonbondedMethod == "PME" else ""),
                hmrOptions,
            )
        )
    elif fileType == "amber":
        script.append(
            "system = prmtop.createSystem(nonbondedMethod=nonbondedMethod,%s"
            % (" nonbondedCutoff=nonbondedCutoff," if nonbondedMethod != "NoCutoff" else "")
        )
        script.append(
            "    constraints=constraints, rigidWater=rigidWater%s%s)"
            % (
                (", ewaldErrorTolerance=ewaldErrorTolerance" if nonbondedMethod == "PME" else ""),
                hmrOptions,
            )
        )
    if ensemble == "npt":
        script.append("system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))")
    script.append("integrator = LangevinMiddleIntegrator(temperature, friction, dt)")
    if constraints != "none":
        script.append("integrator.setConstraintTolerance(constraintTolerance)")
    script.append(
        "simulation = app.Simulation(topology, system, integrator, platform%s)"
        % (", platformProperties" if session["platform"] in ("CUDA", "OpenCL") else "")
    )
    script.append("simulation.context.setPositions(positions)")
    if fileType == "amber":
        script.append("if inpcrd.boxVectors is not None:")
        script.append("    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)")
    # Output XML files for system and integrator

    if session["writeSimulationXml"]:

        def _xml_script_segment(to_serialize, target_file):
            if target_file == "":
                # if filename is blank, we cannot create the file
                return []
            return [
                f'with open("{target_file}", mode="w") as file:',
                f"    file.write(XmlSerializer.serialize({to_serialize}))",
            ]

        script.append("\n# Write XML serialized objects\n")
        script.extend(_xml_script_segment("system", session["systemXmlFilename"]))
        script.extend(_xml_script_segment("integrator", session["integratorXmlFilename"]))

    # Minimize and equilibrate

    script.append("\n# Minimize and Equilibrate\n")
    script.append("print('Performing energy minimization...')")
    script.append("simulation.minimizeEnergy()")
    if fileType == "pdb":
        script.append(
            """
with open(f'Energyminimization_{protein}', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    """
        )
    elif fileType == "amber":
        script.append(
            """
with open(f'Energyminimization_{prmtop_file[:-7]}.pdb', 'w') as outfile:
    PDBFile.writeFile(prmtop.topology, inpcrd.positions, outfile)
    """
        )
    script.append("print('Equilibrating...')")
    script.append("simulation.context.setVelocitiesToTemperature(temperature)")
    script.append("simulation.step(equilibrationSteps)")
    if fileType == "pdb":
        script.append(
            """
with open(f'Equilibration_{protein}', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    """
        )
    elif fileType == "amber":
        script.append(
            """
with open(f'Equilibration_{prmtop_file[:-7]}.pdb', 'w') as outfile:
    PDBFile.writeFile(prmtop.topology, inpcrd.positions, outfile)
    """
        )
    if session["restart_checkpoint"] == "yes":
        script.append("simulation.loadCheckpoint('%s')" % session["checkpointFilename"])

    # Simulate

    script.append("\n# Simulate\n")
    script.append("print('Simulating...')")
    if session["restart_checkpoint"] == "yes":
        if fileType == "pdb":
            script.append("simulation.reporters.append(PDBReporter(f'restart_output_{protein}', pdbInterval))")
        elif fileType == "amber":
            script.append("simulation.reporters.append(PDBReporter(f'restart_output_{prmtop_file}', pdbInterval))")
    else:
        if fileType == "pdb":
            script.append("simulation.reporters.append(PDBReporter(f'output_{protein}', pdbInterval))")
        elif fileType == "amber":
            script.append("simulation.reporters.append(PDBReporter(f'output_{prmtop_file[:-7]}.pdb', pdbInterval))")
    if session["writeDCD"]:
        script.append("simulation.reporters.append(dcdReporter)")
    if session["writeData"]:
        script.append("simulation.reporters.append(dataReporter)")
        if isInternal:
            script.append("simulation.reporters.append(consoleReporter)")
    if session["writeCheckpoint"]:
        script.append("simulation.reporters.append(checkpointReporter)")
        script.append("simulation.reporters.append(checkpointReporter10)")
        script.append("simulation.reporters.append(checkpointReporter100)")
    script.append(
        "simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))"
    )
    if session["restart_checkpoint"] == "yes":
        script.append("simulation.currentStep = %s" % session["restart_step"])
    else:
        script.append("simulation.currentStep = 0")
    script.append("simulation.step(steps)")

    # Output final simulation state
    if session["writeFinalState"]:
        script.append("\n# Write file with final simulation state\n")
        state_script = {
            "checkpoint": ['simulation.saveCheckpoint("{filename}")'],
            "stateXML": ['simulation.saveState("{filename}")'],
            "pdbx": [
                "state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())",
                'with open("{filename}", mode="w") as file:',
                "    PDBxFile.writeFile(simulation.topology, state.getPositions(), file)",
            ],
        }[session["finalStateFileType"]]
        lines = [line.format(filename=session["finalStateFilename"]) for line in state_script]
        script.extend(lines)

    # session[md_postprocessing]
    if session["md_postprocessing"] == "True":
        # mdtraj_conversion() and MDanalysis_conversion()
        if fileType == "pdb":
            script.append("mdtraj_conversion(f'Equilibration_{protein}', '%s')" % session["mdtraj_output"])
            if session["sdfFile"]:
                if session["mdtraj_output"] != "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.pdb', 'centered_old_coordinates.dcd', mda_output='%s', output_selection='%s', ligand_name='UNK')"
                        % (session["mda_output"], session["mda_selection"])
                    )
                elif session["mdtraj_output"] == "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.gro', 'centered_old_coordinates.xtc', mda_output='%s', output_selection='%s', ligand_name='UNK')"
                        % (session["mda_output"], session["mda_selection"])
                    )
            elif session["sdfFile"] == "":
                if session["mdtraj_output"] != "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.pdb', 'centered_old_coordinates.dcd', mda_output='%s', output_selection='%s')"
                        % (session["mda_output"], session["mda_selection"])
                    )
                elif session["mdtraj_output"] == "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.gro', 'centered_old_coordinates.xtc', mda_output='%s', output_selection='%s')"
                        % (session["mda_output"], session["mda_selection"])
                    )
        elif fileType == "amber":
            script.append("mdtraj_conversion(prmtop_file, '%s')" % session["mdtraj_output"])
            if not session["nmLig"] and not session["spLig"]:
                if session["mdtraj_output"] != "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.pdb', 'centered_old_coordinates.dcd', mda_output='%s', output_selection='%s')"
                        % (session["mda_output"], session["mda_selection"])
                    )
                elif session["mdtraj_output"] == "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.gro', 'centered_old_coordinates.xtc', mda_output='%s', output_selection='%s')"
                        % (session["mda_output"], session["mda_selection"])
                    )
            elif session["nmLig"] and not session["spLig"]:
                if session["mdtraj_output"] != "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.pdb', 'centered_old_coordinates.dcd', mda_output='%s', output_selection='%s', ligand_name=%s)"
                        % (
                            session["mda_output"],
                            session["mda_selection"],
                            f"'{nmLigName}'",
                        )
                    )
                elif session["mdtraj_output"] == "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.gro', 'centered_old_coordinates.xtc', mda_output='%s', output_selection='%s', ligand_name=%s)"
                        % (
                            session["mda_output"],
                            session["mda_selection"],
                            f"'{nmLigName}'",
                        )
                    )
            elif session["nmLig"] and session["spLig"]:
                if session["mdtraj_output"] != "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.pdb', 'centered_old_coordinates.dcd', mda_output='%s', output_selection='%s', ligand_name=%s, special_ligname=%s)"
                        % (
                            session["mda_output"],
                            session["mda_selection"],
                            f"'{nmLigName}'",
                            f"'{spLigName}'",
                        )
                    )
                elif session["mdtraj_output"] == "mdtraj_gro_xtc":
                    script.append(
                        "MDanalysis_conversion('centered_old_coordinates_top.gro', 'centered_old_coordinates.xtc', mda_output='%s', output_selection='%s', ligand_name=%s, special_ligname=%s)"
                        % (
                            session["mda_output"],
                            session["mda_selection"],
                            f"'{nmLigName}'",
                            f"'{spLigName}'",
                        )
                    )
        # cleanup()
        if session["mdtraj_removal"] == "True":
            if fileType == "pdb":
                script.append("cleanup(f'{protein}')")
            elif fileType == "amber":
                script.append("cleanup(f'{prmtop_file}')")

    # post_md_file_movement()
    if fileType == "pdb":
        if session["sdfFile"]:
            script.append("post_md_file_movement(protein,ligands=[ligand])")
        elif session["sdfFile"] == "":
            script.append("post_md_file_movement(protein)")
    elif fileType == "amber":
        if (
            session["has_files"] == "yes"
        ):  # In this case, no ligand file will be uploaded, thus not neccessary to assign value to argument `ligands`
            script.append(
                "post_md_file_movement(protein_name=f'{prmtop_file[:-7]}.pdb', prmtop=prmtop_file, inpcrd=inpcrd_file)"
            )
        elif session["has_files"] == "no":
            if not session["nmLig"] and not session["spLig"]:
                script.append(
                    "post_md_file_movement(protein_name=f'{prmtop_file[:-7]}.pdb', prmtop=prmtop_file, inpcrd=inpcrd_file)"
                )
            elif session["nmLig"] and not session["spLig"]:
                script.append(
                    "post_md_file_movement(protein_name=f'{prmtop_file[:-7]}.pdb', prmtop=prmtop_file, inpcrd=inpcrd_file, ligands=['%s'])"
                    % nmLigFileName
                )
            elif session["nmLig"] and session["spLig"]:
                script.append(
                    "post_md_file_movement(protein_name=f'{prmtop_file[:-7]}.pdb', prmtop=prmtop_file, inpcrd=inpcrd_file, ligands=['%s', '%s'])"
                    % (nmLigFileName, spLigFileName)
                )

    # session[openmmdl_analysis]
    if session["openmmdl_analysis"] == "Yes":
        if session["mdtraj_output"] != "mdtraj_gro_xtc":
            top_ext = ".pdb"
            traj_ext = ".dcd"
        elif session["mdtraj_output"] == "mdtraj_gro_xtc":
            top_ext = ".gro"
            traj_ext = ".xtc"
        # session[analysis_selection] == 'analysis_all'
        if session["analysis_selection"] == "analysis_all":
            script.append("os.chdir('Final_Output/All_Atoms')")
            if fileType == "pdb":
                if session["sdfFile"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -l %s -n UNK -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["sdfFile"],
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["sdfFile"] == "":
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
            elif fileType == "amber":
                if not session["nmLig"] and not session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["nmLig"] and not session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -n %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["nmLig"] and session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -n %s -s %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            spLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
        # session[analysis_selection] == 'analysis_prot'
        elif session["analysis_selection"] == "analysis_prot_lig":
            script.append("os.chdir('Final_Output/Prot_Lig')")
            if fileType == "pdb":
                if session["sdfFile"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -l %s -n UNK -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["sdfFile"],
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["sdfFile"] == "":
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
            elif fileType == "amber":
                if not session["nmLig"] and not session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["nmLig"] and not session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -n %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["nmLig"] and session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -n %s -s %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            spLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
        # session[analysis_selection] == 'analysis_all_prot'
        elif session["analysis_selection"] == "analysis_all_prot_lig":
            if fileType == "pdb":
                script.append("os.chdir('Final_Output/All_Atoms')")
                if session["sdfFile"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -l %s -n UNK -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["sdfFile"],
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                    script.append("os.chdir('../Prot_Lig')")
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -l %s -n UNK -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["sdfFile"],
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["sdfFile"] == "":
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                    script.append("os.chdir('../Prot_Lig')")
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
            elif fileType == "amber":
                script.append("os.chdir('Final_Output/All_Atoms')")
                if not session["nmLig"] and not session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                    script.append("os.chdir('../Prot_Lig')")
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["nmLig"] and not session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -n %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                    script.append("os.chdir('../Prot_Lig')")
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -n %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                elif session["nmLig"] and session["spLig"]:
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t centered_top%s -d centered_traj%s -n %s -s %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            spLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )
                    script.append("os.chdir('../Prot_Lig')")
                    script.append(
                        "analysis_run_command = 'openmmdl_analysis -t prot_lig_top%s -d prot_lig_traj%s -n %s -s %s -b %s -m %s -r %s -p %s -w %s --watereps %s' "
                        % (
                            top_ext,
                            traj_ext,
                            nmLigName,
                            spLigName,
                            session["binding_mode"],
                            session["min_transition"],
                            session["rmsd_diff"],
                            session["pml_generation"],
                            session["stable_water"],
                            session["wc_distance"],
                        )
                    )

    return "\n".join(script)


def main():
    def open_browser():
        # Give the server a moment to start before opening the browser.
        time.sleep(1)
        url = "http://127.0.0.1:5000"
        webbrowser.open(url)

    threading.Thread(target=open_browser).start()
    app.run(debug=False)


if __name__ == "__main__":
    main()

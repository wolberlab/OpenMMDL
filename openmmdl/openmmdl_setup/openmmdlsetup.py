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

from openmmdl.openmmdl_setup.setup_options import SetupOptionsConfigurator
from openmmdl.openmmdl_setup.configfile_creator import ConfigCreator
from openmmdl.openmmdl_setup.amberscript_creator import AmberScriptGenerator


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
    session["fileType"] = request.args.get(
        "type", ""
    )  # get the value of `type` from the url
    return showConfigureFiles()


def showConfigureFiles():
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
        if (
            "carboFile" not in request.files
            or request.files["carboFile"].filename == ""
        ):
            showAmberOptions()
        saveUploadedFiles()

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
    
    simulation_configurator = SetupOptionsConfigurator(session)
    
    simulation_configurator.configureDefaultAmberOptions()


def createAmberBashScript():
    a_script = []
    a_script.append(
        "# This script was generated by OpenMMDL-Setup on %s.\n" % datetime.date.today()
    )
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
    boxType = session.get("boxType")

    a_script.append("#!/bin/bash\n")
    
    amber_script = AmberScriptGenerator(session, uploadedFiles)
    
    

    # Receptor
    amber_script.add_receptor_type(a_script)


    amber_script.add_clean_pdb_commands(a_script)
    # Ligand
    amber_script.add_ligand_commands(a_script)


    # Combine all components to be modelled.
    amber_script.add_combine_components_commands(a_script)

    # Add Water/Membrane
    a_script.append(
        "################################ Add Water/Membrane ##############################"
    )

    ## Box setting
    amber_script.add_solvation_commands(a_script)
    amber_script.add_membrane_commands(a_script)
    addType = session["addType"]

    ## Water Setting
    amber_script.add_water_ff_commands(a_script)
    water_ff = session["water_ff"]

    ## Ion Setting
    pos_ion = session["pos_ion"]
    if pos_ion != "other_pos_ion":
        a_script.append("pos_ion=%s" % session["pos_ion"])
    elif pos_ion == "other_pos_ion":
        a_script.append(
            "pos_ion=%s  # In development!" % session["other_pos_ion_input"]
        )

    neg_ion = session["neg_ion"]
    if neg_ion != "other_neg_ion":
        a_script.append("neg_ion=%s" % session["neg_ion"])
    elif neg_ion == "other_neg_ion":
        a_script.append(
            "neg_ion=%s  # In development!" % session["other_neg_ion_input"]
        )

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
        if session["nmLig"] == False and session["spLig"] == False:
            a_script.append(
                "packmol-memgen --pdb ${rcp_nm}_cnt_rmv.pdb --lipids ${lipid_tp} --ratio ${lipid_ratio} --preoriented --dist ${dist2Border} --dist_wat ${padDist} --salt --salt_c ${pos_ion} --saltcon ${ionConc} --nottrim --overwrite --notprotonate\n"
            )
            a_script.append(
                "## Clean the complex pdb by `pdb4amber` for further `tleap` process"
            )
            a_script.append(
                "pdb4amber -i bilayer_${rcp_nm}_cnt_rmv.pdb -o clean_bilayer_${rcp_nm}.pdb"
            )
            ## remove 'CONECT' line in the pdb file
            a_script.append(
                "grep -v '^CONECT' clean_bilayer_${rcp_nm}.pdb > clean_bilayer_${rcp_nm}_cnt_rmv.pdb"
            )
            a_script.append("\n")
        if session["nmLig"] or session["spLig"]:
            a_script.append(
                "packmol-memgen --pdb comp.pdb --lipids ${lipid_tp} --ratio ${lipid_ratio} --preoriented --dist ${dist2Border} --dist_wat ${padDist} --salt --salt_c ${pos_ion} --saltcon ${ionConc} --nottrim --overwrite --notprotonate\n"
            )
            a_script.append(
                "## Clean the complex pdb by `pdb4amber` for further `tleap` process"
            )
            a_script.append("pdb4amber -i bilayer_comp.pdb -o clean_bilayer_comp.pdb")
            ## remove 'CONECT' line in the pdb file
            a_script.append(
                "grep -v '^CONECT' clean_bilayer_comp.pdb > clean_bilayer_comp_cnt_rmv.pdb"
            )
            a_script.append("\n")

    # Generate the prmtop and frcmod file for the complex.
    a_script.append(
        "##################### Generate Prmtop and Frcmod File for the Complex ###################### "
    )
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
        if session["nmLig"] == False and session["spLig"] == False:
            a_script.append("\nsystem = loadpdb ${rcp_nm}_cnt_rmv.pdb\n ")
        else:
            a_script.append("system = loadpdb comp_cnt_rmv.pdb\n")
    elif addType == "addMembrane":
        if session["nmLig"] == False and session["spLig"] == False:
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
    a_script.append(
        "saveamberparm system system.${water_ff}.prmtop system.${water_ff}.inpcrd"
    )
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
    deleteIndices = [
        i for i in range(numChains) if str(i) not in request.form.getlist("include")
    ]
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
        spans.append(
            (chain.id, offset + 1, offset + len(residues), ", ".join(residues))
        )
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
        residues.append(
            (residue.chain.id, residue.name, residue.id, replacements, replaceWith)
        )
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
    allResidues = list(
        set(fixer.missingAtoms.keys()).union(fixer.missingTerminals.keys())
    )
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
        (
            max((pos[i] for pos in fixer.positions))
            - min((pos[i] for pos in fixer.positions))
        ).value_in_unit(unit.nanometer)
        for i in range(3)
    )
    return render_template(
        "addHydrogens.html", unitCell=unitCell, boundingBox=boundingBox
    )


@app.route("/addHydrogens", methods=["POST"])
def addHydrogens():
    session["solvent"] = False
    if "addHydrogens" in request.form:
        pH = float(request.form.get("ph", "7"))
        fixer.addMissingHydrogens(pH)
    if "addWater" in request.form:
        session["solvent"] = True
        session["add_membrane"] = False
        padding, boxSize, boxShape = None, None, None
        if request.form["boxType"] == "geometry":
            session["water_padding"] = True
            session["water_padding_distance"] = float(request.form["geomPadding"])
            session["water_boxShape"] = request.form["geometryDropdown"]
        else:
            session["water_padding"] = False
            session["box_x"] = float(request.form["boxx"])
            session["box_y"] = float(request.form["boxy"])
            session["box_z"] = float(request.form["boxz"])
            boxSize = (
                float(request.form["boxx"]),
                float(request.form["boxy"]),
                float(request.form["boxz"]),
            )
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
        except:
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
    response.headers["Content-Disposition"] = (
        'attachment; filename="OpenMMDL_Simulation.conf"'
    )
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
    return send_file(
        temp, "application/zip", True, "openmmdl_simulation.zip", max_age=0
    )


def configureDefaultOptions():
    """Select default options based on the file format and force field."""
    
    # Initialize the configurator with the session
    simulation_configurator = SetupOptionsConfigurator(session)

    simulation_configurator.configure_default_options()

def createScript():
    script = []
    config_creator = ConfigCreator(session, uploadedFiles)

    # Input files
    
    config_creator.add_openmmdl_ascii_art_logo(script)
    
    config_creator.add_ascii_config_art(script)

    config_creator.add_pdb_input_files_configuration(script)
    
    config_creator.add_amber_file_configuration(script)
    
    config_creator.add_forcefield_and_water_model_configuration(script)

    
    config_creator.add_solvent_configuration(script)
    

    # System configuration
    config_creator.add_system_configuration(script)

    # Integration options

    config_creator.add_integration_configuration(script)

    # Simulation options


    config_creator.add_simulation_time_and_steps_configuration(script)

    config_creator.add_equilibration_configuration(script)

    config_creator.add_simulation_configuration(script)
    
    config_creator.add_checkpoint_configuration(script)


    config_creator.add_xml_serialization_configuration(script)

    config_creator.add_postprocessing_configuration(script)

    config_creator.add_openmmdl_analysis_configuration(script)





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

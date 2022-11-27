import openmm as mm
import openmm.unit as unit
from openmm.app import PDBFile, PDBxFile
from pdbfixer.pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat
from flask import Flask, request, session, g, render_template, make_response, send_file, url_for
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

if sys.version_info >= (3,0):
    from io import StringIO
else:
    from cStringIO import StringIO


app = Flask(__name__)
app.config.from_object(__name__)
app.config.update({'SECRET_KEY':'development key'})
app.jinja_env.globals['mm'] = mm

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

@app.route('/headerControls')
def headerControls():
    if 'startOver' in request.args:
        return showSelectFileType()
    if 'quit' in request.args:
        func = request.environ.get('werkzeug.server.shutdown')
        if func is None:
            raise RuntimeError('Not running with the Werkzeug Server')
        func()
        return "OpenMM Setup has stopped running.  You can close this window."

@app.route('/')
def showSelectFileType():
    return render_template('selectFileType.html')

@app.route('/selectFiles')
def selectFiles():
    session['fileType'] = request.args.get('type', '')
    return showConfigureFiles()

def showConfigureFiles():
    try:
        fileType = session['fileType']
        if fileType == 'pdb':
            return render_template('configurePdbFile.html')
        elif fileType == 'amber':
            return render_template('configureAmberFiles.html')
    except:
        app.logger.error('Error displaying configure files page', exc_info=True)
    # The file type is invalid, so send them back to the select file type page.
    return showSelectFileType()

@app.route('/configureFiles', methods=['POST'])
def configureFiles():
    fileType = session['fileType']
    if fileType == 'pdb':
        if 'file' not in request.files or request.files['file'].filename == '':
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        saveUploadedFiles()
        session['forcefield'] = request.form.get('forcefield', '')
        session['waterModel'] = request.form.get('waterModel', '')
        session['cleanup'] = request.form.get('cleanup', '')
        configureDefaultOptions()
        file, name = uploadedFiles['file'][0]
        file.seek(0, 0)
        session['pdbType'] = _guessFileFormat(file, name)
        if session['cleanup'] == 'yes':
            global fixer
            if session['pdbType'] == 'pdb':
                fixer = PDBFixer(pdbfile=file)
            else:
                fixer = PDBFixer(pdbxfile=StringIO(file.read().decode()))
            return showSelectChains()
    elif fileType == 'amber':
        if 'prmtopFile' not in request.files or request.files['prmtopFile'].filename == '' or 'inpcrdFile' not in request.files or request.files['inpcrdFile'].filename == '':
            # They didn't select a file.  Send them back.
            return showConfigureFiles()
        saveUploadedFiles()
    configureDefaultOptions()
    return showSimulationOptions()

@app.route('/getCurrentStructure')
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
            content = ', '.join(set(residues))
            hasHeterogen = True
        chains.append((chain.id, len(residues), content))
    if len(chains) < 2 and not hasHeterogen:
        session['heterogens'] = 'all'
        return showAddResidues()
    return render_template('selectChains.html', chains=chains)

@app.route('/selectChains', methods=['POST'])
def selectChains():
    session['heterogens'] = request.form.get('heterogens', '')
    numChains = len(list(fixer.topology.chains()))
    request.form.getlist('include')
    deleteIndices = [i for i in range(numChains) if str(i) not in request.form.getlist('include')]
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
            offset = int(chainResidues[key[1]].id)-len(residues)-1
        else:
            offset = int(chainResidues[-1].id)
        spans.append((chain.id, offset+1, offset+len(residues), ', '.join(residues)))
    return render_template('addResidues.html', spans=spans)

@app.route('/addResidues', methods=['POST'])
def addResidues():
    keys = [key for key in sorted(fixer.missingResidues)]
    for i, key in enumerate(keys):
        if str(i) not in request.form.getlist('add'):
            del fixer.missingResidues[key]
    return showConvertResidues()

def showConvertResidues():
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) == 0:
        return showAddHeavyAtoms()
    residues = []
    nucleotides = ['DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'T']
    for i in range(len(fixer.nonstandardResidues)):
        residue, replaceWith = fixer.nonstandardResidues[i]
        if replaceWith in proteinResidues:
            replacements = proteinResidues
        else:
            replacements = nucleotides
        residues.append((residue.chain.id, residue.name, residue.id, replacements, replaceWith))
    return render_template('convertResidues.html', residues=residues)

@app.route('/convertResidues', methods=['POST'])
def convertResidues():
    for i in range(len(fixer.nonstandardResidues)):
        if str(i) in request.form.getlist('convert'):
            fixer.nonstandardResidues[i] = (fixer.nonstandardResidues[i][0], request.form['residue'+str(i)])
    fixer.replaceNonstandardResidues()
    return showAddHeavyAtoms()

def showAddHeavyAtoms():
    if session['heterogens'] == 'none':
        fixer.removeHeterogens(False)
    elif session['heterogens'] == 'water':
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
        residues.append((residue.chain.id, residue.name, residue.id, ', '.join(atoms)))
    return render_template('addHeavyAtoms.html', residues=residues)

@app.route('/addHeavyAtoms', methods=['POST'])
def addHeavyAtoms():
    fixer.addMissingAtoms()
    return showAddHydrogens()

def showAddHydrogens():
    unitCell = fixer.topology.getUnitCellDimensions()
    if unitCell is not None:
        unitCell = unitCell.value_in_unit(unit.nanometer)
    boundingBox = tuple((max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions))).value_in_unit(unit.nanometer) for i in range(3))
    return render_template('addHydrogens.html', unitCell=unitCell, boundingBox=boundingBox)

@app.route('/addHydrogens', methods=['POST'])
def addHydrogens():
    session['solvent'] = False
    if 'addHydrogens' in request.form:
        pH = float(request.form.get('ph', '7'))
        fixer.addMissingHydrogens(pH)
    if 'addWater' in request.form:
        session['solvent'] = True
        session['add_membrane'] = False
        padding, boxSize, boxShape = None, None, None
        if request.form['boxType'] == 'geometry':
            session['water_padding'] = True
            session['water_padding_distance'] = float(request.form['geomPadding'])
            session['water_boxShape'] = request.form['geometryDropdown']
        else:
            session['water_padding'] = False
            session['box_x'] = (float(request.form['boxx']))
            session['box_y'] = (float(request.form['boxy']))
            session['box_z'] = (float(request.form['boxz']))
            boxSize = (float(request.form['boxx']), float(request.form['boxy']), float(request.form['boxz']))
        session['water_ionicstrength'] = float(request.form['ionicstrength'])
        session['water_positive'] = request.form['positiveion']+'+'
        session['water_negative'] = request.form['negativeion']+'-'
    elif 'addMembrane' in request.form:
        session['solvent'] = True
        session['add_membrane'] = True
        session['lipidType'] = request.form['lipidType']
        session['membrane_padding'] = float(request.form['membranePadding'])
        session['membrane_ionicstrength'] = float(request.form['ionicstrength'])
        session['membrane_positive'] = request.form['positiveion']+'+'
        session['membrane_negative'] = request.form['negativeion']+'-'
    
    # Save the new PDB file.
    
    uploadedFiles['originalFile'] = uploadedFiles['file']
    pdb = StringIO()
    if session['pdbType'] == 'pdb':
        try:
            PDBFile.writeFile(fixer.topology, fixer.positions, pdb, True)
        except:
            # This can happen if the ids are too large to fit in the allowed space.
            pdb = StringIO()
            PDBFile.writeFile(fixer.topology, fixer.positions, pdb, False)
    else:
        PDBxFile.writeFile(fixer.topology, fixer.positions, pdb, True)
    temp = tempfile.TemporaryFile()
    temp.write(pdb.getvalue().encode('utf-8'))
    name = uploadedFiles['file'][0][1]
    dotIndex = name.rfind('.')
    if dotIndex == -1:
        prefix = name
        suffix = ''
    else:
        prefix = name[:dotIndex]
        suffix = name[dotIndex:]
    uploadedFiles['file'] = [(temp, prefix+'-processed_openMMDL'+suffix)]
    return showSimulationOptions()

def showSimulationOptions():
    return render_template('simulationOptions.html')

@app.route('/setSimulationOptions', methods=['POST'])
def setSimulationOptions():
    for key in request.form:
        session[key] = request.form[key]
    session['ligand'] = 'ligand' in request.form
    session['writeDCD'] = 'writeDCD' in request.form
    session['writeData'] = 'writeData' in request.form
    session['writeCheckpoint'] = 'writeCheckpoint' in request.form
    session['dataFields'] = request.form.getlist('dataFields')
    session['hmr'] = 'hmr' in request.form
    session['writeSimulationXml'] = 'writeSimulationXml' in request.form
    session['writeFinalState'] = 'writeFinalState' in request.form
    return createScript()

@app.route('/downloadScript')
def downloadScript():
    response = make_response(createScript())
    response.headers['Content-Disposition'] = 'attachment; filename="OpenMMDL_Simulation.py"'
    return response

@app.route('/downloadPDB')
def downloadPDB():
    file, name = uploadedFiles['file'][0]
    file.seek(0, 0)
    response = make_response(file.read())
    response.headers['Content-Disposition'] = 'attachment; filename="%s"' % name
    return response

@app.route('/downloadPackage')
def downloadPackage():
    temp = tempfile.NamedTemporaryFile()
    with zipfile.ZipFile(temp, 'w', zipfile.ZIP_DEFLATED) as zip:
        zip.writestr('openmm_simulation/OpenMMDL_Simulation.py', createScript())
        for key in uploadedFiles:
            for file, name in uploadedFiles[key]:
                file.seek(0, 0)
                zip.writestr('openmm_simulation/%s' % name, file.read())
    temp.seek(0, 0)
    return send_file(temp, 'application/zip', True, 'openmm_simulation.zip', max_age=0)

@app.route('/showRunSimulation')
def showRunSimulation():
    homeDir = os.path.expanduser('~')
    defaultDir = os.path.join(homeDir, 'openmm_simulation')
    return render_template('runSimulation.html', defaultDir=defaultDir)

@app.route('/startSimulation', methods=['POST'])
def startSimulation():
    global scriptOutput, simulationProcess
    conn1, conn2 = Pipe()
    scriptOutput = conn1
    # Create the simulation directory and copy files.
    try:
        outputDir = request.form['directory']
        if not os.path.isdir(outputDir):
            os.makedirs(outputDir)
    except:
        conn2.send('An error occurred while creating the simulation directory: %s' % sys.exc_info()[1])
        conn2.send(None)
        return ""
    try:
        for key in uploadedFiles:
            for file, name in uploadedFiles[key]:
                file.seek(0, 0)
                with open(os.path.join(outputDir, name), 'wb') as outputFile:
                    shutil.copyfileobj(file, outputFile)
        with open(os.path.join(outputDir, 'OpenMMDL_Simulation.py'), 'w') as outputFile:
            outputFile.write(createScript())
    except:
        conn2.send('An error occurred while copying the input files: %s' % sys.exc_info()[1])
        conn2.send(None)
        return ""
    # Run the simulation in a subprocess.
    simulationProcess = Process(target=simulate, args=(conn2, outputDir, createScript(True)))
    simulationProcess.start()
    return ""

@app.route('/stopSimulation', methods=['POST'])
def stopSimulation():
    global scriptOutput, simulationProcess
    os.kill(simulationProcess.pid, signal.SIGKILL)
    scriptOutput = None
    return ""

@app.route('/getSimulationOutput')
def getSimulationOutput():
    global scriptOutput
    if scriptOutput is None:
        return "", 404
    output = []
    try:
        while scriptOutput.poll():
            data = scriptOutput.recv()
            if data is None:
                scriptOutput = None
                break
            else:
                output.append(data)
    except EOFError:
        scriptOutput = None
    return "".join(output)

def simulate(output, outputDir, script):
    try:
        exec(script, {"output":output, "outputDir":outputDir})
    except Exception as e:
        output.send('\nThe simulation failed with the following error:\n\n')
        output.send(str(e))
        output.send('\n\nDetails:\n\n')
        output.send(traceback.format_exc())
    output.send(None)

def configureDefaultOptions():
    """Select default options based on the file format and force field."""
    session['ligand_select'] = 'no'
    session['minimization'] = 'False'
    implicitWater = False
    session['restart_checkpoint'] = 'no'
    session['manual_sys'] = 'no'
    session['md_postprocessing'] = 'True'
    session['mdtraj_remove'] = 'False'
    session['box_choice'] = 'no_solvent'
    session['wat_padding'] = '1.0'
    session['box_vector_x'] = '10.0'
    session['box_vector_y'] = '10.0'
    session['box_vector_z'] = '10.0'
    session['membrane_type_system'] = 'POPC'
    session['mem_padding'] = '1.0'
    session['ionic_str'] = '0.15'
    session['positive_ion'] = 'Na+'
    session['negative_ion'] = 'Cl-'
    if session['fileType'] == 'pdb' and session['waterModel'] == 'implicit':
        implicitWater = True
    session['ensemble'] = 'nvt' if implicitWater else 'npt'
    session['platform'] = 'CUDA'
    session['precision'] = 'single'
    session['cutoff'] = '2.0' if implicitWater else '1.0'
    session['ewaldTol'] = '0.0005'
    session['constraintTol'] = '0.000001'
    session['hmr'] = True
    session['hmrMass'] = '1.5'
    session['dt'] = '0.004'
    session['steps'] = '1000000'
    session['equilibrationSteps'] = '1000'
    session['temperature'] = '300'
    session['friction'] = '1.0'
    session['pressure'] = '1.0'
    session['barostatInterval'] = '25'
    session['nonbondedMethod'] = 'CutoffNonPeriodic' if implicitWater else 'PME'
    session['writeDCD'] = True
    session['dcdFilename'] = 'trajectory.dcd'
    session['dcdInterval'] = '10000'
    session['pdbinterval'] = '10000'
    session['writeData'] = True
    session['dataFilename'] = 'log.txt'
    session['dataInterval'] = '1000'
    session['dataFields'] = ['step', 'speed' ,'progress', 'potentialEnergy', 'temperature']
    session['writeCheckpoint'] = True
    session['checkpointFilename'] = 'checkpoint.chk'
    session['checkpointInterval'] = '10000'
    session['writeSimulationXml'] = False
    session['systemXmlFilename'] = 'system.xml'
    session['integratorXmlFilename'] = 'integrator.xml'
    session['writeFinalState'] = False
    session['finalStateFileType'] = 'stateXML'
    session['finalStateFilename'] = "final_state.xml"
    session['constraints'] = 'hbonds'
    session['rmsd'] = 'True'
    session['interaction_analysis'] = 'False'

def createScript(isInternal=False):
    script = []

    # If we are creating this script for internal use to run a simulation directly, add extra code at the top
    # to set the working directory and redirect stdout to the pipe.

    if isInternal:
        script.append("""
import os
import sys
import time

class PipeOutput(object):
    def write(self, string):
        output.send(string)

sys.stdout = PipeOutput()
sys.stderr = PipeOutput()
os.chdir(outputDir)""")

    # Header
    
    script.append('# This script was generated by OpenMM-MDL Setup on %s.\n' % datetime.date.today())
    script.append('''
#       ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      
#     .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      
#    / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      
#   ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    
#   |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    
#   : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    
#    \ `"/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  
#     '. \_/``".'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ 
#       '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` 
                                                                                                      
                                                                                                      
''')
    script.append('from scripts.forcefield_water import ff_selection, water_forecfield_selection, water_model_selection, generate_forcefield, generate_transitional_forcefield')
    script.append('from scripts.protein_ligand_prep import protein_choice, prepare_ligand, rdkit_to_openmm, merge_protein_and_ligand, water_padding_solvent_builder, water_absolute_solvent_builder, membrane_builder, water_conversion')
    script.append('from scripts.post_md_conversions import mdtraj_conversion, MDanalysis_conversion, rmsd_for_atomgroups, RMSD_dist_frames, atomic_distance')
    script.append('from scripts.cleaning_procedures import cleanup, post_md_file_movement \n')
    
    script.append('import simtk.openmm.app as app')
    script.append('from simtk.openmm.app import PDBFile, Modeller, PDBReporter, StateDataReporter, DCDReporter, CheckpointReporter')
    script.append('from simtk.openmm import unit, Platform, Platform_getPlatformByName, MonteCarloBarostat, LangevinMiddleIntegrator')
    script.append('from simtk.openmm import Vec3')
    script.append('import simtk.openmm as mm')
    script.append('import sys')
    script.append('import os')
    script.append('import shutil')
    
   
    # Input files
    
    script.append('\n# Input Files')
    script.append('''############# Ligand and Protein Data ###################''')
    script.append('''########   Add the Ligand SDf File and Protein PDB File in the Folder with the Script  ######### \n''')
    if session['ligand_select'] == 'yes':
        script.append('ligand_select = "%s"' % session['ligand_select'])
        script.append('ligand_name = "UNK"')
        script.append('ligand_sdf = "%s"' % session['sdf_file'])
        script.append('\nminimize = %s '% session['minimization'])
    fileType = session['fileType']
    if fileType == 'pdb':
        pdbType = session['pdbType']
        if pdbType == 'pdb':
            script.append('protein = "%s"' % uploadedFiles['file'][0][1])
            forcefield = session['forcefield']
            water_model = session['waterModel']
            water = session['waterModel']
    elif fileType == 'amber':
        script.append("prmtop = AmberPrmtopFile('%s')" % uploadedFiles['prmtopFile'][0][1])
        script.append("inpcrd = AmberInpcrdFile('%s')" % uploadedFiles['inpcrdFile'][0][1])

    script.append('''\n############# Ligand and Protein Preparation ###################\n''')
    script.append('protein_prepared = "Yes"')
    
    script.append('''\n############# Forcefield, Water and Membrane Model Selection ###################\n''')
    if fileType == 'pdb':
        script.append("ff = '%s'" % session['forcefield'])
        if water != 'None':
            script.append("water = '%s'" % water)
        else:
            script.append("water = %s" % water)



################################## IF CLEANING WAS PERFORMED ##############################################
###########################################################################################################
###########################################################################################################
    
    
    
    if fileType == 'pdb':
        if session['cleanup'] == 'yes':
            if session['solvent'] == True:
                if session['add_membrane'] == True and session['manual_sys'] == 'yes' and session['box_choice'] == 'membrane':
                    script.append('''\n############# Membrane Settings ###################\n''')
                    script.append("add_membrane = %s" % session['add_membrane'])
                    script.append("membrane_lipid_type = '%s'" % session['membrane_type_system'])
                    script.append("membrane_padding = %s" % session['mem_padding'])
                    script.append("membrane_ionicstrength = %s" % session['ionic_str'])
                    script.append("membrane_positive_ion = '%s'" % session['positive_ion'])
                    script.append("membrane_negative_ion = '%s'" % session['negative_ion'])
                elif session['add_membrane'] == True and session['manual_sys'] == 'no':
                    script.append('''\n############# Membrane Settings ###################\n''')
                    script.append("add_membrane = %s" % session['add_membrane'])
                    script.append("membrane_lipid_type = '%s'" % session['lipidType'])
                    script.append("membrane_padding = %s" % session['membrane_padding'])
                    script.append("membrane_ionicstrength = %s" % session['membrane_ionicstrength'])
                    script.append("membrane_positive_ion = '%s'" % session['membrane_positive'])
                    script.append("membrane_negative_ion = '%s'" % session['membrane_negative'])
                elif session['add_membrane'] == False and session['manual_sys'] == 'yes' and session['box_choice'] == 'membrane':
                    script.append('''\n############# Membrane Settings ###################\n''')
                    script.append("add_membrane = %s" % session['add_membrane'])
                    script.append("membrane_lipid_type = '%s'" % session['membrane_type_system'])
                    script.append("membrane_padding = %s" % session['mem_padding'])
                    script.append("membrane_ionicstrength = %s" % session['ionic_str'])
                    script.append("membrane_positive_ion = '%s'" % session['positive_ion'])
                    script.append("membrane_negative_ion = '%s'" % session['negative_ion'])

                elif session['add_membrane'] == False and session['manual_sys'] == 'no':
                    script.append('''\n############# Water Box Settings ###################\n''')
                    script.append("add_membrane = %s" % session['add_membrane'])
                    if session['water_padding'] == True:
                        script.append('Water_Box = "Buffer"')
                        script.append("water_padding_distance = %s" % session['water_padding_distance'])
                        script.append("water_boxShape = '%s'" % session['water_boxShape'])
                    else:
                        script.append('Water_Box = "Absolute"')
                        script.append("water_box_x = %s" % session['box_x'])
                        script.append("water_box_y = %s" % session['box_y']) 
                        script.append("water_box_z = %s" % session['box_z'])   
                    script.append("water_ionicstrength = %s" % session['water_ionicstrength'])
                    script.append("water_positive_ion = '%s'" % session['water_positive'])
                    script.append("water_negative_ion = '%s'" % session['water_negative'])
            
                elif session['add_membrane'] == False and session['manual_sys'] == 'yes' and session['box_choice'] == 'water':
                    script.append('''\n############# Water Box Settings ###################\n''')
                    script.append("add_membrane = False" % session['add_membrane'])
                    if session['water_type_system'] == 'Buffer':
                        script.append('Water_Box = "Buffer"')
                        script.append("water_padding_distance = %s" % session['wat_padding'])
                        script.append("water_boxShape = 'cube'")        
                    else:
                        script.append('Water_Box = "Absolute"')
                        script.append("water_box_x = %s" % session['box_vector_x'])
                        script.append("water_box_y = %s" % session['box_vector_y']) 
                        script.append("water_box_z = %s" % session['box_vector_z'])  
                    script.append("water_ionicstrength = %s" % session['ionic_str'])
                    script.append("water_positive_ion = '%s'" % session['positive_ion'])
                    script.append("water_negative_ion = '%s'" % session['negative_ion'])
                
                elif session['add_membrane'] == True and session['manual_sys'] == 'yes' and session['box_choice'] == 'water':
                    script.append('''\n############# Water Box Settings ###################\n''')
                    script.append("add_membrane = False")
                    if session['water_type_system'] == 'Buffer':
                        script.append('Water_Box = "Buffer"')
                        script.append("water_padding_distance = %s" % session['wat_padding'])
                        script.append("water_boxShape = 'cube'")        
                    else:
                        script.append('Water_Box = "Absolute"')
                        script.append("water_box_x = %s" % session['box_vector_x'])
                        script.append("water_box_y = %s" % session['box_vector_y']) 
                        script.append("water_box_z = %s" % session['box_vector_z'])  
                    script.append("water_ionicstrength = %s" % session['ionic_str'])
                    script.append("water_positive_ion = '%s'" % session['positive_ion'])
                    script.append("water_negative_ion = '%s'" % session['negative_ion'])   
                    
            else:
                if session['solvent'] == False:
                    if session['manual_sys'] == 'no' or session['box_choice'] == 'no_solvent':
                        script.append("Solvent = %s" % session['solvent'])
                    elif session['manual_sys'] == 'yes' and session['box_choice'] == 'membrane':
                        script.append('''\n############# Membrane Settings ###################\n''')
                        script.append("add_membrane = %s" % session['add_membrane'])
                        script.append("membrane_lipid_type = '%s'" % session['membrane_type_system'])
                        script.append("membrane_padding = %s" % session['mem_padding'])
                        script.append("membrane_ionicstrength = %s" % session['ionic_str'])
                        script.append("membrane_positive_ion = '%s'" % session['positive_ion'])
                        script.append("membrane_negative_ion = '%s'" % session['negative_ion'])
                    elif session['manual_sys'] == 'yes' and session['box_choice'] == 'water':
                        if session['water_type_system'] == 'Buffer':
                            script.append('Water_Box = "Buffer"')
                            script.append("water_padding_distance = %s" % session['wat_padding'])
                            script.append("water_boxShape = 'cube'")        
                        else:
                            script.append('Water_Box = "Absolute"')
                            script.append("water_box_x = %s" % session['box_vector_x'])
                            script.append("water_box_y = %s" % session['box_vector_y']) 
                            script.append("water_box_z = %s" % session['box_vector_z'])  
                        script.append("water_ionicstrength = %s" % session['ionic_str'])
                        script.append("water_positive_ion = '%s'" % session['positive_ion'])
                        script.append("water_negative_ion = '%s'" % session['negative_ion'])  
    

################################## IF CLEANING WAS NOT PERFORMED ##########################################
###########################################################################################################
###########################################################################################################

        elif session['cleanup'] == 'no':
            if session['manual_sys'] == 'yes' and session['box_choice'] == 'membrane':
                script.append('''\n############# Membrane Settings ###################\n''')
                script.append("add_membrane = %s" % session['add_membrane'])
                script.append("membrane_lipid_type = '%s'" % session['membrane_type_system'])
                script.append("membrane_padding = %s" % session['mem_padding'])
                script.append("membrane_ionicstrength = %s" % session['ionic_str'])
                script.append("membrane_positive_ion = '%s'" % session['positive_ion'])
                script.append("membrane_negative_ion = '%s'" % session['negative_ion'])
 
            elif session['manual_sys'] == 'yes' and session['box_choice'] == 'water':
                script.append('''\n############# Water Box Settings ###################\n''')
                script.append("add_membrane = False")
                if session['water_type_system'] == 'Buffer':
                    script.append('Water_Box = "Buffer"')
                    script.append("water_padding_distance = %s" % session['wat_padding'])
                    script.append("water_boxShape = 'cube'")        
                else:
                    script.append('Water_Box = "Absolute"')
                    script.append("water_box_x = %s" % session['box_vector_x'])
                    script.append("water_box_y = %s" % session['box_vector_y']) 
                    script.append("water_box_z = %s" % session['box_vector_z'])  
                script.append("water_ionicstrength = %s" % session['ionic_str'])
                script.append("water_positive_ion = '%s'" % session['positive_ion'])
                script.append("water_negative_ion = '%s'" % session['negative_ion'])  
    
    script.append('''\n############# Post MD Processing ###################\n''')
    script.append('MDAnalysis_Postprocessing = %s' % session['md_postprocessing'])
    script.append('MDTraj_Cleanup = %s' % session['md_postprocessing'])

    # System configuration
    script.append('\n# System Configuration\n')
    nonbondedMethod = session['nonbondedMethod']
    script.append('nonbondedMethod = app.%s' % nonbondedMethod)
    if nonbondedMethod != 'NoCutoff':
        script.append('nonbondedCutoff = %s*unit.nanometers' % session['cutoff'])
    if nonbondedMethod == 'PME':
        script.append('ewaldErrorTolerance = %s' % session['ewaldTol'])
    constraints = session['constraints']
    constraintMethods = {'none': 'None',
                         'water': 'None',
                         'hbonds': 'HBonds',
                         'allbonds': 'AllBonds'}
    if constraints != 'none' and constraints != 'water':
        script.append('constraints = app.%s' % constraintMethods[constraints])
    if constraints == 'none':
        script.append('constraints = %s' % constraintMethods[constraints])
    script.append('rigidWater = %s' % ('False' if constraints == 'none' else 'True'))
    if constraints != 'none':
        script.append('constraintTolerance = %s' % session['constraintTol'])
    if session['hmr']:
        script.append('hydrogenMass = %s*unit.amu' % session['hmrMass'])

    # Integration options

    script.append('\n# Integration Options\n')
    script.append('dt = %s*unit.picoseconds' % session['dt'])
    script.append('temperature = %s*unit.kelvin' % session['temperature'])
    script.append('friction = %s/unit.picosecond' % session['friction'])
    ensemble = session['ensemble']
    if ensemble == 'npt':
            script.append('pressure = %s*unit.atmospheres' % session['pressure'])
            script.append('barostatInterval = %s' % session['barostatInterval'])

    # Simulation options

    script.append('\n# Simulation Options\n')
    script.append('steps = %s' % session['steps'])
    if session['restart_checkpoint'] == 'yes':
        script.append('restart_step = %s' % session['restart_step'])
    script.append('equilibrationSteps = %s' % session['equilibrationSteps'])
    script.append("platform = Platform.getPlatformByName('%s')" % session['platform'])
    if session['platform'] in ('CUDA', 'OpenCL'):
        script.append("platformProperties = {'Precision': '%s'}" % session['precision'])
    if session['writeDCD']:
        if session['restart_checkpoint'] == 'yes':
            script.append("dcdReporter = DCDReporter('%s_%s', %s)" % (session['restart_step'], session['dcdFilename'], session['dcdInterval']))
        else:
            script.append("dcdReporter = DCDReporter('%s', %s)" % (session['dcdFilename'], session['dcdInterval']))
    if session['writeData']:
        args = ', '.join('%s=True' % field for field in session['dataFields'])
        if session['restart_checkpoint'] == 'yes':
            script.append("dataReporter = StateDataReporter('%s_%s', %s, totalSteps=steps," % (session['restart_step'], session['dataFilename'], session['dataInterval']))
        else:
            script.append("dataReporter = StateDataReporter('%s', %s, totalSteps=steps," % (session['dataFilename'], session['dataInterval']))
        script.append("    %s, separator='\\t')" % args)
        if isInternal:
            # Create a second reporting sending to stdout so we can display it in the browser.
            script.append("consoleReporter = StateDataReporter(sys.stdout, %s, totalSteps=steps, %s, separator='\\t')" % (session['dataInterval'], args))
    if session['writeCheckpoint']:
        if session['restart_checkpoint'] == 'yes':
            script.append("checkpointReporter = CheckpointReporter('%s_%s', %s)" % (session['restart_step'], session['checkpointFilename'], session['checkpointInterval']))
            script.append("checkpointReporter10 = CheckpointReporter('10x_%s__%s', %s0)" % (session['restart_step'], session['checkpointFilename'], session['checkpointInterval']))
            script.append("checkpointReporter100 = CheckpointReporter('100x_%s_%s', %s00)" % (session['restart_step'], session['checkpointFilename'], session['checkpointInterval']))
        else:
            script.append("checkpointReporter = CheckpointReporter('%s', %s)" % (session['checkpointFilename'], session['checkpointInterval']))
            script.append("checkpointReporter10 = CheckpointReporter('10x_%s', %s0)" % (session['checkpointFilename'], session['checkpointInterval']))
            script.append("checkpointReporter100 = CheckpointReporter('100x_%s', %s00)" % (session['checkpointFilename'], session['checkpointInterval']))
    
    # Prepare the simulation
    
    if fileType  == 'pdb':
        if session['ligand_select'] == 'yes':
            script.append('''
if ligand_select == 'yes':
    
    print("Preparing MD Simulation with ligand")
    
    ligand_prepared = prepare_ligand(ligand_sdf,minimize_molecule=minimize)
     
    omm_ligand = rdkit_to_openmm(ligand_prepared, ligand_name) ''')
        script.append('''
protein_pdb = protein_choice(protein_is_prepared=protein_prepared,protein=protein)
forcefield_selected = ff_selection(ff)
water_selected = water_forecfield_selection(water=water,forcefield_selection=ff_selection(ff))
model_water = water_model_selection(water=water,forcefield_selection=ff_selection(ff))
print("Forcefield and Water Model Selected")   ''')

        if session['ligand_select'] == 'yes':
            script.append('''
if ligand_select == 'yes':
    
    if add_membrane == True:
        transitional_forcefield = generate_transitional_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, rdkit_mol=ligand_prepared)
    
    forcefield = generate_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, rdkit_mol=ligand_prepared)
    
    complex_topology, complex_positions = merge_protein_and_ligand(protein_pdb, omm_ligand)
    
    print("Complex topology has", complex_topology.getNumAtoms(), "atoms.")     ''')
        elif session['ligand_select'] == 'no':
            script.append('''      
if water_selected != None:
    forcefield = generate_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, rdkit_mol=None) 
else:
    forcefield = app.ForceField(forcefield_selected)    
    
if add_membrane == True:
        transitional_forcefield = generate_transitional_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, rdkit_mol=None)     ''')


        
        if session['manual_sys'] == 'yes' and session['box_choice'] == 'no_solvent':
            script.append('''
topology = protein_pdb.topology
positions = protein_pdb.positions                       ''') 
    
    
        elif session['manual_sys'] == 'yes' and session['box_choice'] == 'membrane': 
            script.append('''
modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)

membrane_builder(ff, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, membrane_lipid_type, membrane_padding, membrane_positive_ion, membrane_negative_ion, membrane_ionicstrength, protein)

if model_water == 'tip4pew' or model_water == 'tip5p':
    water_conversion(model_water, modeller, protein)

topology = modeller.topology
positions = modeller.positions
           ''')
       
        elif session['manual_sys'] == 'yes' and session['box_choice'] == 'water':  
            script.append('''
modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)       
if Water_Box == "Buffer":
    water_padding_solvent_builder(model_water, forcefield, water_padding_distance, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)    	        
elif Water_Box == "Absolute":
    water_absolute_solvent_builder(model_water, forcefield, water_box_x, water_box_y, water_box_z, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)

topology = modeller.topology
positions = modeller.positions
       ''')
       
        if session['cleanup'] == 'yes' and session['ligand_select'] == 'no':
            script.append('''
forcefield = generate_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, rdkit_mol=None)        

modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
 
if add_membrane == True:
    membrane_builder(ff, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, membrane_lipid_type, membrane_padding, membrane_positive_ion, membrane_negative_ion, membrane_ionicstrength, protein)

elif add_membrane == False:
    if Water_Box == "Buffer":
        water_padding_solvent_builder(model_water, forcefield, water_padding_distance, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
    elif Water_Box == "Absolute":
        water_absolute_solvent_builder(model_water, forcefield, water_box_x, water_box_y, water_box_z, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
    
if add_membrane == True:
    if model_water == 'tip4pew' or model_water == 'tip5p':
        water_conversion(model_water, modeller, protein)
    
topology = modeller.topology

positions = modeller.positions
        ''')
    
        elif session['cleanup'] == 'yes' and session['ligand_select'] == 'yes':
            script.append('''
modeller = app.Modeller(complex_topology, complex_positions)
 
if add_membrane == True:
    membrane_builder(ff, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, membrane_lipid_type, membrane_padding, membrane_positive_ion, membrane_negative_ion, membrane_ionicstrength, protein)

elif add_membrane == False:
    if Water_Box == "Buffer":
        water_padding_solvent_builder(model_water, forcefield, water_padding_distance, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
    elif Water_Box == "Absolute":
        water_absolute_solvent_builder(model_water, forcefield, water_box_x, water_box_y, water_box_z, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
    
if add_membrane == True:
    if model_water == 'tip4pew' or model_water == 'tip5p':
        water_conversion(model_water, modeller, protein)

topology = modeller.topology

positions = modeller.positions

        ''')
    
    
    
        elif session['cleanup'] == 'no' and session['ligand_select'] == 'yes':
            script.append('''
forcefield_selected = ff_selection(ff)
    
water_selected =water_selection(water=water,force_selection=ff_selection(ff))
    
forcefield = generate_forcefield(protein_ff=forcefield_selected, solvent_ff=water_selected, add_membrane=add_membrane, rdkit_mol=ligand_prepared)

modeller = app.Modeller(complex_topology, complex_positions)

if add_membrane == True:
    membrane_builder(ff, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, membrane_lipid_type, membrane_padding, membrane_positive_ion, membrane_negative_ion, membrane_ionicstrength, protein)
elif add_membrane == False:
    if Water_Box == "Buffer":
        water_padding_solvent_builder(model_water, forcefield, water_padding_distance, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
    elif Water_Box == "Absolute":
        water_absolute_solvent_builder(model_water, forcefield, water_box_x, water_box_y, water_box_z, protein_pdb, modeller, water_positive_ion, water_negative_ion, water_ionicstrength, protein)
            ''')
            
    elif fileType == 'amber':
        script.append('topology = prmtop.topology')
        script.append('positions = inpcrd.positions')

    script.append('\n# Prepare the Simulation\n')
    script.append("print('Building system...')")
    hmrOptions = ', hydrogenMass=hydrogenMass' if session['hmr'] else ''
    if fileType  == 'pdb':
        script.append('system = forcefield.createSystem(topology, nonbondedMethod=nonbondedMethod,%s' % (' nonbondedCutoff=nonbondedCutoff,' if nonbondedMethod != 'NoCutoff' else ''))
        script.append('    constraints=constraints, rigidWater=rigidWater%s%s)' % (', ewaldErrorTolerance=ewaldErrorTolerance' if nonbondedMethod == 'PME' else '', hmrOptions))
    elif fileType == 'amber':
        script.append('system = prmtop.createSystem(nonbondedMethod=nonbondedMethod,%s' % (' nonbondedCutoff=nonbondedCutoff,' if nonbondedMethod != 'NoCutoff' else ''))
        script.append('    constraints=constraints, rigidWater=rigidWater%s%s)' % (', ewaldErrorTolerance=ewaldErrorTolerance' if nonbondedMethod == 'PME' else '', hmrOptions))
    if ensemble == 'npt':
            script.append('system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))')
    script.append('integrator = LangevinMiddleIntegrator(temperature, friction, dt)')
    if constraints != 'none':
        script.append('integrator.setConstraintTolerance(constraintTolerance)')
    script.append('simulation = app.Simulation(topology, system, integrator, platform%s)' % (', platformProperties' if session['platform'] in ('CUDA', 'OpenCL') else ''))
    script.append('simulation.context.setPositions(positions)')
    if fileType == 'amber':
        script.append('if inpcrd.boxVectors is not None:')
        script.append('    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)')
    # Output XML files for system and integrator

    if session['writeSimulationXml']:
        def _xml_script_segment(to_serialize, target_file):
            if target_file == "":
                # if filename is blank, we cannot create the file
                return []
            return [
                f'with open("{target_file}", mode="w") as file:',
                f'    file.write(XmlSerializer.serialize({to_serialize}))'
            ]

        script.append("\n# Write XML serialized objects\n")
        script.extend(_xml_script_segment('system', session['systemXmlFilename']))
        script.extend(_xml_script_segment('integrator', session['integratorXmlFilename']))

    # Minimize and equilibrate
    
    script.append('\n# Minimize and Equilibrate\n')
    script.append("print('Performing energy minimization...')")
    script.append('simulation.minimizeEnergy()')
    if fileType  == 'pdb':
        script.append("""
with open(f'Energyminimization_{protein}', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    """)
    script.append("print('Equilibrating...')")
    script.append('simulation.context.setVelocitiesToTemperature(temperature)')
    script.append('simulation.step(equilibrationSteps)')
    if fileType  == 'pdb':
        script.append("""
with open(f'Equilibration_{protein}', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)
    """)    
    if session['restart_checkpoint'] == 'yes':
        script.append("simulation.loadCheckpoint('%s')" % session['checkpointFilename'])
    
    # Simulate
    
    script.append('\n# Simulate\n')
    script.append("print('Simulating...')")
    if session['restart_checkpoint'] == 'yes':
        script.append("simulation.reporters.append(PDBReporter(f'restart_output_{protein}', %s))" % session['pdbinterval'])
    else:
        script.append("simulation.reporters.append(PDBReporter(f'output_{protein}', %s))" % session['pdbinterval'])
    if session['writeDCD']:
        script.append('simulation.reporters.append(dcdReporter)')
    if session['writeData']:
        script.append('simulation.reporters.append(dataReporter)')
        if isInternal:
            script.append('simulation.reporters.append(consoleReporter)')
    if session['writeCheckpoint']:
        script.append('simulation.reporters.append(checkpointReporter)')
        script.append('simulation.reporters.append(checkpointReporter10)')
        script.append('simulation.reporters.append(checkpointReporter100)')
    script.append('simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))')
    if session['restart_checkpoint'] == 'yes':
        script.append('simulation.currentStep = %s' % session ['restart_step'])
    else:
        script.append('simulation.currentStep = 0')
    script.append('simulation.step(steps)')

    # Output final simulation state
    if session['writeFinalState']:
        script.append("\n# Write file with final simulation state\n")
        state_script = {
            'checkpoint': ['simulation.saveCheckpoint("{filename}")'],
            'stateXML': ['simulation.saveState("{filename}")'],
            'pdbx': ['state = simulation.context.getState(getPositions=True, enforcePeriodicBox=system.usesPeriodicBoundaryConditions())',
                     'with open("{filename}", mode="w") as file:',
                     '    PDBxFile.writeFile(simulation.topology, state.getPositions(), file)'],
        }[session['finalStateFileType']]
        lines = [line.format(filename=session['finalStateFilename']) for line in state_script]
        script.extend(lines)
    
    if session ['md_postprocessing'] == 'True':
            script.append("mdtraj_conversion(f'Equilibration_{protein}')")
            script.append("MDanalysis_conversion(f'centered_old_coordinates.pdb', f'centered_old_coordinates.dcd', ligand_name='UNK')")
    
    if session['rmsd'] == 'True':
            script.append("rmsd_for_atomgroups(f'prot_lig_top.pdb', f'prot_lig_traj.dcd', selection1='backbone', selection2=['protein', 'resname UNK'])")
            script.append("RMSD_dist_frames(f'prot_lig_top.pdb', f'prot_lig_traj.dcd', lig='UNK')")
    
    script.append('post_md_file_movement(protein,ligand_sdf)')
    return "\n".join(script)


def main():

    def open_browser():
        # Give the server a moment to start before opening the browser.
        time.sleep(1)
        url = 'http://127.0.0.1:5000'
        webbrowser.open(url)

    threading.Thread(target=open_browser).start()
    app.run(debug=False)

if __name__ == '__main__':
    main()

import pytest
from openmmdl.openmmdl_setup.configfile_creator import ConfigCreator, ConfigWriter  # Adjust import as necessary

def test_add_openmmdl_ascii_art_logo():
    """Test if add_openmmdl_ascii_art_logo correctly adds the ASCII art logo."""

    # Create a ConfigCreator instance
    session = {}
    uploadedFiles = {}
    config_creator = ConfigCreator(session, uploadedFiles)

    # Initialize the script list
    script = []
    config_creator.add_openmmdl_ascii_art_logo(script)

    # Expected ASCII art
    expected_logo = """
       ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      
     .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      
    / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      
   ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    
   |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    
   : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    
    \ `"/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  
     '. \_/``".'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ 
       '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------`                                                                                                  
            """
    
    # Compare script content with expected output
    assert script[0].strip() == expected_logo.strip()
    
def test_add_ascii_config_art():
    """Test if add_ascii_config_art correctly adds the configuration file ASCII art header."""

    # Create a ConfigCreator instance
    session = {}
    uploadedFiles = {}
    config_creator = ConfigCreator(session, uploadedFiles)

    # Initialize the script list
    script = []
    config_creator.add_ascii_config_art(script)

    # Expected ASCII art
    expected_art = """
                            __   __        ___    __      ___         ___ 
                           /  ` /  \ |\ | |__  | / _`    |__  | |    |__  
                           \__, \__/ | \| |    | \__>    |    | |___ |___                                                                                                                              
            """
    
    # Compare script content with expected output
    assert script[0].strip() == expected_art.strip()
    
@pytest.fixture
def config_creator_pdb():
    session = {
        "fileType": "pdb",
        "pdbType": "pdb",
        "sdfFile": "ligand.sdf",
        "ligandMinimization": "minimization_method",
        "smallMoleculeForceField": "force_field",
        "ligandSanitization": "sanitization_method",
        "waterModel": "spce"
    }
    uploadedFiles = {
        "file": [("protein.pdb", "path/to/protein.pdb")]
    }
    return ConfigCreator(session, uploadedFiles)

def test_add_pdb_input_files_configuration(config_creator_pdb):
    """Test if add_pdb_input_files_configuration correctly configures PDB input files."""

    # Initialize the script list
    script = []
    config_creator_pdb.add_pdb_input_files_configuration(script)

    # Expected script content
    expected_lines = [
        "\n# Input Files",
        "############# Ligand and Protein Data ###################",
        "input_file_type = pdb",
        "protein = path/to/protein.pdb",
        "ligand = ligand.sdf",
        "ligand_name = UNK",
        "minimization = minimization_method",
        "smallMoleculeForceField = force_field",
        "sanitization = sanitization_method"
    ]

    # Compare script content with expected output
    assert script == expected_lines

@pytest.fixture
def config_creator_amber():
    session = {
        "fileType": "amber",
        "has_files": "yes",
        "nmLig": True,
        "spLig": False,
        "nmLigName": "nmLig",
        "spLigName": None,
        "water_ff": "tip3p"
    }
    uploadedFiles = {
        "prmtopFile": [("prmtop.prmtop", "path/to/prmtop.prmtop")],
        "inpcrdFile": [("inpcrd.inpcrd", "path/to/inpcrd.inpcrd")],
        "nmLigFile": [("nmLig.sdf", "path/to/nmLig.sdf")]
    }
    return ConfigCreator(session, uploadedFiles)

def test_add_amber_file_configuration(config_creator_amber):
    """Test if add_amber_file_configuration correctly configures Amber files."""

    # Initialize the script list
    script = []
    config_creator_amber.add_amber_file_configuration(script)

    # Expected script content
    expected_lines = [
        """####### Add the Amber Files in the Folder with this Script ####### \n""",
        "input_file_type = amber",
        "prmtop_file = path/to/prmtop.prmtop",
        "inpcrd_file = path/to/inpcrd.inpcrd",
        "prmtop = AmberPrmtopFile(prmtop_file)",
        "inpcrd = AmberInpcrdFile(inpcrd_file)"
    ]

    # Compare script content with expected output
    assert script == expected_lines
    
@pytest.fixture
def config_creator_ff():
    session = {
        "fileType": "pdb",
        "forcefield": "ff99SB",
        "waterModel": "spce"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_forcefield_and_water_model_configuration(config_creator_ff):
    """Test if add_forcefield_and_water_model_configuration correctly configures forcefield and water model."""

    # Initialize the script list
    script = []
    config_creator_ff.add_forcefield_and_water_model_configuration(script)

    # Expected script content
    expected_lines = [
        "\n############# Forcefield, Water and Membrane Model Selection ###################\n",
        "forcefield = ff99SB",
        "water = spce"
    ]

    # Compare script content with expected output
    assert script == expected_lines
    
@pytest.fixture
def config_creator_solvent():
    session = {
        "fileType": "pdb",
        "solvent": True,
        "add_membrane": False,
        "water_padding": True,
        "water_padding_distance": "10.0",
        "water_boxShape": "cubic",
        "box_x": None,
        "box_y": None,
        "box_z": None,
        "water_ionicstrength": "0.15",
        "water_positive": "Na+",
        "water_negative": "Cl-"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_solvent_configuration(config_creator_solvent):
    """Test if add_solvent_configuration correctly configures solvent or membrane settings."""

    # Initialize the script list
    script = []
    config_creator_solvent.add_solvent_configuration(script)

    # Expected script content
    expected_lines = [
        "\n############# Water Box Settings ###################\n",
        "add_membrane = False",
        "Water_Box = Buffer",
        "water_padding_distance = 10.0",
        "water_boxShape = cubic",
        "water_ionicstrength = 0.15",
        "water_positive_ion = Na+",
        "water_negative_ion = Cl-"
    ]

    # Compare script content with expected output
    assert script == expected_lines
    
@pytest.fixture
def config_creator_system():
    session = {
        "nonbondedMethod": "PME",
        "cutoff": "1.0",
        "ewaldTol": "0.0005",
        "hmr": True,
        "hmrMass": "1.008",
        "constraints": "hbonds",
        "constraintTol": "0.0001"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_system_configuration(config_creator_system):
    """Test if add_system_configuration correctly configures system settings."""

    # Initialize the script list
    script = []
    config_creator_system.add_system_configuration(script)
    print(script)

    # Expected script content
    expected_lines = [
        "\n# System Configuration\n",
        "nonbondedMethod = app.PME",
        "nonbondedCutoff = 1.0*unit.nanometers",
        "ewaldErrorTolerance = 0.0005",
        "hmrOptions = ', hydrogenMass=hydrogenMass'",
        "constraints = app.HBonds",
        "rigidWater = True",
        "constraintTolerance = 0.0001",
        "hydrogenMass = 1.008*unit.amu"
    ]

    # Compare script content with expected output
    assert script == expected_lines
    
@pytest.fixture
def config_creator_integrator():
    session = {
        "dt": "0.002",
        "temperature": "300",
        "friction": "1.0",
        "ensemble": "npt",
        "pressure": "1.0",
        "barostatInterval": "100"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_integration_configuration(config_creator_integrator):
    """Test if add_integration_configuration correctly configures integration settings."""

    # Initialize the script list
    script = []
    config_creator_integrator.add_integration_configuration(script)

    # Expected script content
    expected_lines = [
        "\n# Integration Configuration\n",
        "step_time = 0.002",
        "dt = 0.002*unit.picoseconds",
        "temperature = 300*unit.kelvin",
        "friction = 1.0/unit.picosecond",
        "pressure = 1.0*unit.atmospheres",
        "barostatInterval = 100"
    ]

    # Compare script content with expected output
    assert script == expected_lines
    
@pytest.fixture
def config_creator_sim_time():
    session = {
        "sim_length": "1000",  # Simulation length in ns
        "dt": "0.002",         # Time step in ps
        "dcdFrames": "5000",   # Number of frames for DCD output
        "pdbInterval_ns": "10" # Interval for PDB output in ns
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_simulation_time_and_steps_configuration(config_creator_sim_time):
    """Test if add_simulation_time_and_steps_configuration correctly configures simulation time and steps."""

    # Initialize the script list
    script = []
    config_creator_sim_time.add_simulation_time_and_steps_configuration(script)

    # Calculate expected values
    steps = int(float("1000") / float("0.002") * 1000)  # Total steps
    dcdinterval = int(steps / int("5000"))  # DCD interval
    pdbInterval = int(
        steps * (float("10") / float("1000"))  # PDB interval
    )

    # Expected script content
    expected_lines = [
        "\n# Simulation Time and Steps Configuration\n",
        "sim_length = 1000",
        "steps = %s" % steps,
        "\n# Frames and Interval Configuration\n",
        "dcdFrames = 5000",
        "dcdInterval = %s" % dcdinterval,
        "pdbInterval_ns = 10",
        "pdbInterval = %s" % pdbInterval
    ]

    # Compare script content with expected output
    assert script == expected_lines
    
@pytest.fixture
def config_creator_equil():
    session = {
        "equilibration": "equilibration"  # Change to "minimization" for testing that case
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_equilibration_configuration(config_creator_equil):
    """Test if add_equilibration_configuration correctly configures equilibration or minimization settings."""

    # Initialize the script list
    script = []
    config_creator_equil.add_equilibration_configuration(script)

    # Expected script content
    expected_lines = [
        "\n# Equilibration & Minimization Configuration\n",
        "preparation_type = equilibration"
    ]

    # Compare script content with expected output
    assert script == expected_lines

    # Change the session parameter to test "minimization"
    config_creator_equil.session["equilibration"] = "minimization"
    script = []
    config_creator_equil.add_equilibration_configuration(script)

    # Expected script content for minimization
    expected_lines = [
        "\n# Equilibration & Minimization Configuration\n",
        "preparation_type = minimization"
    ]

    # Compare script content with expected output
    assert script == expected_lines
    
@pytest.fixture
def config_creator_simulation():
    session = {
        "platform": "CUDA",
        "precision": "mixed",
        "writeDCD": True,
        "dcdFilename": "simulation",
        "writeData": True,
        "dataFields": ["energy", "temperature"],
        "dataInterval": "1000",
        "restart_checkpoint": "yes",
        "restart_step": "step1",
        "dataFilename": "data_reporter"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_simulation_configuration(config_creator_simulation):
    """Test if add_simulation_configuration correctly configures simulation platform, precision, and file outputs."""

    # Initialize the script list
    script = []
    config_creator_simulation.add_simulation_configuration(script)

    # Expected script content
    expected_lines = [
        "\n# Simulation Configuration\n",
        "platform = CUDA",
        "platformProperties = {'Precision': 'mixed'}",
        "dcd_name = step1_simulation",
        "dataReporter = StateDataReporter('step1_data_reporter', 1000, totalSteps=steps,",
        "    energy=True, temperature=True, separator='\\t')"
    ]

    # Compare script content with expected output
    assert script == expected_lines

    # Test case with no data reporting
    config_creator_simulation.session["writeData"] = False
    script = []
    config_creator_simulation.add_simulation_configuration(script)
    
    # Expected script content without dataReporter
    expected_lines = [
        "\n# Simulation Configuration\n",
        "platform = CUDA",
        "platformProperties = {'Precision': 'mixed'}",
        "dcd_name = step1_simulation"
    ]
    
    assert script == expected_lines
    

@pytest.fixture
def config_creator_checkpoint():
    session = {
        "writeCheckpoint": True,
        "checkpointInterval_ns": "500",  # Checkpoint interval in ns
        "dt": "0.002",                   # Time step in ps
        "checkpointFilename": "checkpoint",
        "restart_checkpoint": "yes",
        "restart_step": "step1"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_checkpoint_configuration(config_creator_checkpoint):
    """Test if add_checkpoint_configuration correctly configures checkpoint and restart settings."""

    # Initialize the script list
    script = []
    config_creator_checkpoint.add_checkpoint_configuration(script)

    # Calculate expected checkpoint interval
    checkpointInterval = int(
        1000 * float("500") / float("0.002")
    )

    # Expected script content
    expected_lines = [
        "\n# Checkpoint and Restart Configuration\n",
        "checkpointInterval = %s" % checkpointInterval,
        "checkpoint_name = checkpoint",
        "restart_step = step1"
    ]

    # Compare script content with expected output
    assert script == expected_lines

    # Test case with no checkpoint
    config_creator_checkpoint.session["writeCheckpoint"] = False
    script = []
    config_creator_checkpoint.add_checkpoint_configuration(script)
    
    # Expected script content without checkpoint configuration
    expected_lines = []

    assert script == expected_lines
    
@pytest.fixture
def config_creator_xml():
    session = {
        "writeSimulationXml": True,
        "systemXmlFilename": "system.xml",
        "integratorXmlFilename": "integrator.xml"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_xml_serialization_configuration(config_creator_xml):
    """Test if add_xml_serialization_configuration correctly configures XML serialization settings."""

    # Initialize the script list
    script = []
    config_creator_xml.add_xml_serialization_configuration(script)

    # Expected script content
    expected_lines = [
        "\n# Write XML Serialized Objects\n",
        "xmlsystem_filename = system.xml",
        "xmlintegrator_filename = integrator.xml"
    ]

    # Compare script content with expected output
    assert script == expected_lines

    # Test case with no XML serialization
    config_creator_xml.session["writeSimulationXml"] = False
    script = []
    config_creator_xml.add_xml_serialization_configuration(script)
    
    # Expected script content without XML serialization configuration
    expected_lines = []

    assert script == expected_lines
    
@pytest.fixture
def config_creator_postprocessing():
    session = {
        "md_postprocessing": "enabled",
        "mdtraj_output": "mdtraj_output.pdb",
        "mdtraj_removal": "mdtraj_removal.pdb",
        "mda_output": "mda_output.h5",
        "mda_selection": "resname LIG"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_postprocessing_configuration(config_creator_postprocessing):
    """Test if add_postprocessing_configuration correctly configures MD post-processing settings."""

    # Initialize the script list
    script = []
    config_creator_postprocessing.add_postprocessing_configuration(script)

    # Expected script content
    expected_lines = [
        "\n# Post-Processing Configuration\n",
        "postprocessing = enabled",
        "old_output = mdtraj_output.pdb",
        "old_removal = mdtraj_removal.pdb",
        "mda_output = mda_output.h5",
        "mda_selection = resname LIG"
    ]

    # Compare script content with expected output
    assert script == expected_lines
    

@pytest.fixture
def config_creator_analysis():
    session = {
        "openmmdl_analysis": "Yes",
        "analysis_selection": "all",
        "binding_mode": "flexible",
        "min_transition": "5.0",
        "rmsd_diff": "0.1",
        "pml_generation": "enabled"
    }
    uploadedFiles = {}
    return ConfigCreator(session, uploadedFiles)

def test_add_openmmdl_analysis_configuration(config_creator_analysis):
    """Test if add_openmmdl_analysis_configuration correctly configures OpenMMDL Analysis settings."""

    # Initialize the script list
    script = []
    config_creator_analysis.add_openmmdl_analysis_configuration(script)

    # Expected script content
    expected_lines = [
        "\n# OpenMMDL Analysis Configuration\n",
        "openmmdl_analysis = Yes",
        "analysis_selection = all",
        "binding_mode = flexible",
        "min_transition = 5.0",
        "rmsd_diff = 0.1",
        "pml_generation = enabled"
    ]

    # Compare script content with expected output
    assert script == expected_lines

    # Test case with OpenMMDL Analysis disabled
    config_creator_analysis.session["openmmdl_analysis"] = "No"
    script = []
    config_creator_analysis.add_openmmdl_analysis_configuration(script)
    
    # Expected script content without OpenMMDL Analysis configuration
    expected_lines = [
        "\n# OpenMMDL Analysis Configuration\n",
        "openmmdl_analysis = No"
    ]

    assert script == expected_lines
    


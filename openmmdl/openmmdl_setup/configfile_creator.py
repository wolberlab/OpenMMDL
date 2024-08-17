import datetime
from openmmdl.openmmdl_setup.file_operator import LigandExtractor


class ConfigCreator:
    """
    ConfigCreator is responsible for generating different sections of a simulation
    configuration script based on the session parameters and uploaded files.

    Attributes:
        session (dict): A dictionary containing user-defined session parameters.
        uploadedFiles (dict): A dictionary containing information about the files uploaded by the user.
    """

    def __init__(self, session, uploadedFiles):
        """
        Initializes the ConfigCreator with session parameters and uploaded files from OpenMMDL Setup.

        Args:
            session (dict): A dictionary containing user-defined session parameters.
            uploadedFiles (dict): A dictionary containing information about the files uploaded by the user during the usage of OpenMMDL Setup.
        """
        self.session = session
        self.uploadedFiles = uploadedFiles

    def add_openmmdl_ascii_art_logo(self, script):
        """
        Adds the OpenMMDL logo to the configuration file.

        Args:
            script (list): A list of strings representing the configuration file being built.
        """
        script.append(
            """
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
        )

    def add_ascii_config_art(self, script):
        """
        Adds a Config File ASCII art header for the configuration file.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append(
            """
                            __   __        ___    __      ___         ___ 
                           /  ` /  \ |\ | |__  | / _`    |__  | |    |__  
                           \__, \__/ | \| |    | \__>    |    | |___ |___                                                                                                                              
            """
        )

    def add_pdb_input_files_configuration(self, script):
        """
        Adds configuration settings for PDB input files, including protein and ligand data.
        Handles both protein-only and protein-ligand systems, setting up the necessary variables
        for simulation input.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append("\n# Input Files")
        if self.session["fileType"] == "pdb":
            script.append(
                """############# Ligand and Protein Data ###################"""
            )
            if self.session["pdbType"] == "pdb":
                script.append('input_file_type = "pdb"')
                script.append('protein = "%s"' % self.uploadedFiles["file"][0][1])
                if self.session["sdfFile"] != "":
                    script.append("ligand = '%s'" % self.session["sdfFile"])
                    script.append('ligand_name = "UNK"')
                    script.append(
                        "minimization = %s" % self.session["ligandMinimization"]
                    )
                    script.append(
                        "smallMoleculeForceField = '%s'"
                        % self.session["smallMoleculeForceField"]
                    )
                    script.append(
                        "sanitization = %s" % self.session["ligandSanitization"]
                    )
            water = self.session["waterModel"]

    def add_amber_file_configuration(self, script):
        """
        Adds configuration settings for Amber input files, including topology and coordinate files.
        This method handles both scenarios where Amber files are pre-existing or need to be generated.
        It also extracts ligand names if ligands are present.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        if self.session["fileType"] == "amber":
            script.append(
                """####### Add the Amber Files in the Folder with this Script ####### \n"""
            )

            # amber_files related variables
            if self.session["has_files"] == "yes":
                script.append('input_file_type = "amber"')
                script.append(
                    "prmtop_file = '%s'" % self.uploadedFiles["prmtopFile"][0][1]
                )
                script.append(
                    'inpcrd_file = "%s"' % self.uploadedFiles["inpcrdFile"][0][1]
                )

                # ligand related variables
                nmLigName = self.session["nmLigName"] if self.session["nmLig"] else None
                spLigName = self.session["spLigName"] if self.session["spLig"] else None

            elif self.session["has_files"] == "no":
                script.append(
                    "prmtop_file = 'system.%s.prmtop'" % self.session["water_ff"]
                )
                script.append(
                    "inpcrd_file = 'system.%s.inpcrd' " % self.session["water_ff"]
                )

                # ligand related variables
                if self.session["nmLig"]:
                    nmLigFileName = self.uploadedFiles["nmLigFile"][0][1]
                    nmLigName = LigandExtractor.extract_ligand_name(nmLigFileName)
                else:
                    nmLigFileName = None
                    nmLigName = None

                if self.session["spLig"]:
                    spLigFileName = self.uploadedFiles["spLigFile"][0][1]
                    spLigName = LigandExtractor.extract_ligand_name(spLigFileName)
                else:
                    spLigFileName = None
                    spLigName = None

            ## Feed prmtop_file and inpcrd_file to OpenMM Reader
            script.append("prmtop = AmberPrmtopFile(prmtop_file)")
            script.append("inpcrd = AmberInpcrdFile(inpcrd_file)")

    def add_forcefield_and_water_model_configuration(self, script):
        """
        Adds forcefield and water model configuration settings to the script based on the file type.
        This includes selecting the appropriate forcefield and water model if applicable.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        if self.session["fileType"] == "pdb":
            script.append(
                """\n############# Forcefield, Water and Membrane Model Selection ###################\n"""
            )
            script.append("forcefield = '%s'" % self.session["forcefield"])
            if self.session["waterModel"] != "None":
                script.append("water = '%s'" % self.session["waterModel"])
            else:
                script.append("water = %s" % self.session["waterModel"])

    def add_solvent_configuration(self, script):
        """
        Adds solvent or membrane configuration settings to the script, depending on whether
        a membrane or water box is being used, including parameters like padding, ionic strength,
        and ion types.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        if self.session["fileType"] == "pdb":
            if self.session["solvent"]:
                if self.session["add_membrane"]:
                    script.append(
                        "\n############# Membrane Settings ###################\n"
                    )
                    script.append("add_membrane = %s" % self.session["add_membrane"])
                    script.append(
                        "membrane_lipid_type = '%s'" % self.session["lipidType"]
                    )
                    script.append(
                        "membrane_padding = %s" % self.session["membrane_padding"]
                    )
                    script.append(
                        "membrane_ionicstrength = %s"
                        % self.session["membrane_ionicstrength"]
                    )
                    script.append(
                        "membrane_positive_ion = '%s'"
                        % self.session["membrane_positive"]
                    )
                    script.append(
                        "membrane_negative_ion = '%s'"
                        % self.session["membrane_negative"]
                    )
                else:
                    script.append(
                        "\n############# Water Box Settings ###################\n"
                    )
                    script.append("add_membrane = %s" % self.session["add_membrane"])
                    if self.session["water_padding"]:
                        script.append('Water_Box = "Buffer"')
                        script.append(
                            "water_padding_distance = %s"
                            % self.session["water_padding_distance"]
                        )
                        script.append(
                            "water_boxShape = '%s'" % self.session["water_boxShape"]
                        )
                    else:
                        script.append('Water_Box = "Absolute"')
                        script.append("water_box_x = %s" % self.session["box_x"])
                        script.append("water_box_y = %s" % self.session["box_y"])
                        script.append("water_box_z = %s" % self.session["box_z"])
                    script.append(
                        "water_ionicstrength = %s" % self.session["water_ionicstrength"]
                    )
                    script.append(
                        "water_positive_ion = '%s'" % self.session["water_positive"]
                    )
                    script.append(
                        "water_negative_ion = '%s'" % self.session["water_negative"]
                    )
            else:
                script.append("Solvent = %s" % self.session["solvent"])

    def add_system_configuration(self, script):
        """
        Adds system configuration settings to the script, including nonbonded method, cutoff distances,
        ewald error tolerance, hydrogen mass repartitioning (HMR), and constraints settings.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append("\n# System Configuration\n")
        nonbondedMethod = self.session["nonbondedMethod"]
        script.append("nonbondedMethod = app.%s" % nonbondedMethod)

        if nonbondedMethod != "NoCutoff":
            script.append(
                "nonbondedCutoff = %s*unit.nanometers" % self.session["cutoff"]
            )

        if nonbondedMethod == "PME":
            script.append("ewaldErrorTolerance = %s" % self.session["ewaldTol"])

        hmrOptions = ", hydrogenMass=hydrogenMass" if self.session["hmr"] else ""
        script.append("hmrOptions = '%s'" % hmrOptions)

        constraints = self.session["constraints"]
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

        script.append(
            "rigidWater = %s" % ("False" if constraints == "none" else "True")
        )

        if constraints != "none":
            script.append("constraintTolerance = %s" % self.session["constraintTol"])

        if self.session["hmr"]:
            script.append("hydrogenMass = %s*unit.amu" % self.session["hmrMass"])

    def add_integration_configuration(self, script):
        """
        Adds the integration configuration to the script, specifying time step, temperature,
        friction, and pressure (for NPT ensemble).

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append("\n# Integration Configuration\n")
        script.append("step_time = %s" % self.session["dt"])
        script.append("dt = %s*unit.picoseconds" % self.session["dt"])
        script.append("temperature = %s*unit.kelvin" % self.session["temperature"])
        script.append("friction = %s/unit.picosecond" % self.session["friction"])

        ensemble = self.session["ensemble"]
        if ensemble == "npt":
            script.append("pressure = %s*unit.atmospheres" % self.session["pressure"])
            script.append("barostatInterval = %s" % self.session["barostatInterval"])

    def add_simulation_time_and_steps_configuration(self, script):
        """
        Adds the simulation time and step configuration to the script, including
        the total simulation length, number of steps, and intervals for frames and PDB output.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append("\n# Simulation Time and Steps Configuration\n")
        script.append("sim_length = %s" % self.session["sim_length"])
        steps = int(
            float(self.session["sim_length"]) / float(self.session["dt"]) * 1000
        )
        script.append("steps = %s" % steps)

        script.append("\n# Frames and Interval Configuration\n")
        script.append("dcdFrames = %s" % self.session["dcdFrames"])
        dcdinterval = int(steps / int(self.session["dcdFrames"]))
        script.append("dcdInterval = %s" % dcdinterval)
        script.append("pdbInterval_ns = %s" % self.session["pdbInterval_ns"])
        pdbInterval = int(
            steps
            * (
                float(self.session["pdbInterval_ns"])
                / float(self.session["sim_length"])
            )
        )
        script.append("pdbInterval = %s" % pdbInterval)

    def add_equilibration_configuration(self, script):
        """
        Adds equilibration or minimization configuration to the script. It defines different
        equilibration stages based on the selected method (equilibration, minimization, or none).

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        equilibration_type = self.session["equilibration"]
        script.append("\n# Equilibration & Minimization Configuration\n")

        if equilibration_type == "equilibration":
            script.append(
                """equilibration_stages = [
                {'EOM': 'minimize', 'n_steps': 10000, 'temperature': 300*unit.kelvin, 'ensemble': None, 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD_interpolate', 'n_steps': 100000, 'temperature': 100*unit.kelvin, 'temperature_end': 300*unit.kelvin, 'ensemble': 'NVT', 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 10/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300, 'ensemble': 'NPT', 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 10/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 250000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and not type H', 'force_constant': 10*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'minimize', 'n_steps': 10000, 'temperature': 300*unit.kelvin, 'ensemble': None, 'restraint_selection': 'protein and backbone', 'force_constant': 10*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and backbone', 'force_constant': 10*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and backbone', 'force_constant': 1*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and backbone', 'force_constant': 0.1*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 2500000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': None, 'force_constant': 0*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 2*unit.femtoseconds},
                ]
            """
            )
        elif equilibration_type == "minimization":
            script.append(
                """equilibration_stages = [
                {'EOM': 'minimize', 'n_steps': 10000, 'temperature': 300*unit.kelvin, 'ensemble': None, 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                ]
            """
            )
        elif equilibration_type == "None":
            script.append("equilibration_stages = None")

    def add_simulation_configuration(self, script):
        """
        Adds the simulation configuration to the script, including platform settings,
        precision, and file outputs like DCD and data reporters.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append("\n# Simulation Configuration\n")

        script.append(
            "platform = Platform.getPlatformByName('%s')" % self.session["platform"]
        )
        if self.session["platform"] in ("CUDA", "OpenCL"):
            script.append(
                "platformProperties = {'Precision': '%s'}" % self.session["precision"]
            )
        if self.session["writeDCD"]:
            if self.session["restart_step"] == "":
                script.append(
                    "dcdReporter = DCDReporter('%s_%s', dcdInterval)"
                    % (self.session["restart_step"], self.session["dcdFilename"])
                )
            else:
                script.append(
                    "dcdReporter = DCDReporter('%s', dcdInterval)"
                    % (self.session["dcdFilename"])
                )
        if self.session["writeData"]:
            args = ", ".join("%s=True" % field for field in self.session["dataFields"])
            if self.session["restart_checkpoint"] == "yes":
                script.append(
                    "dataReporter = StateDataReporter('%s_%s', %s, totalSteps=steps,"
                    % (
                        self.session["restart_step"],
                        self.session["dataFilename"],
                        self.session["dataInterval"],
                    )
                )
            else:
                script.append(
                    "dataReporter = StateDataReporter('%s', %s, totalSteps=steps,"
                    % (self.session["dataFilename"], self.session["dataInterval"])
                )
            script.append("    %s, separator='\\t')" % args)

    def add_checkpoint_configuration(self, script):
        """
        Adds checkpoint and restart configuration settings to the script.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        if self.session["writeCheckpoint"]:
            checkpointInterval = int(
                1000
                * float(self.session["checkpointInterval_ns"])
                / float(self.session["dt"])
            )

            script.append("\n# Checkpoint and Restart Configuration\n")
            script.append("checkpointInterval = %s" % checkpointInterval)
            script.append("checkpoint_name = %s" % self.session["checkpointFilename"])

            if self.session["restart_checkpoint"]:
                script.append("restart_step = %s" % self.session["restart_step"])

    def add_xml_serialization_configuration(self, script):
        """
        Adds XML serialization settings for the system and integrator objects.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        if self.session["writeSimulationXml"]:
            script.append("\n# Write XML Serialized Objects\n")
            script.append("xmlsystem_filename = %s" % self.session["systemXmlFilename"])
            script.append(
                "xmlintegrator_filename = %s" % self.session["integratorXmlFilename"]
            )

    def add_postprocessing_configuration(self, script):
        """
        Adds post-processing configuration, including settings for MD post-processing,
        with MDTraj and MDAnalysis outputs.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append("\n# Post-Processing Configuration\n")
        script.append("postprocessing = %s" % self.session["md_postprocessing"])
        script.append("old_output = %s" % self.session["mdtraj_output"])
        script.append("old_removal = %s" % self.session["mdtraj_removal"])
        script.append("mda_output = %s" % self.session["mda_output"])
        script.append("mda_selection = %s" % self.session["mda_selection"])

    def add_openmmdl_analysis_configuration(self, script):
        """
        Adds OpenMMDL Analysis configuration to the script, including various analysis settings
        if OpenMMDL Analysis is enabled.

        Args:
            script (list): A list of strings representing the configuration script being built.
        """
        script.append("\n# OpenMMDL Analysis Configuration\n")
        script.append("openmmdl_analysis = %s" % self.session["openmmdl_analysis"])

        if self.session["openmmdl_analysis"] == "Yes":
            script.append(
                "analysis_selection = %s" % self.session["analysis_selection"]
            )
            script.append("binding_mode = %s" % self.session["binding_mode"])
            script.append("min_transition = %s" % self.session["min_transition"])
            script.append("rmsd_diff = %s" % self.session["rmsd_diff"])
            script.append("pml_generation = %s" % self.session["pml_generation"])


class ConfigWriter:
    """
    This class is responsible for generating a configuration script based on the provided session
    and uploaded files. It acts as a wrapper around the ConfigCreator class, orchestrating the
    various configuration sections that need to be included in the final script.

    Attributes:
        script (list): A list of strings that make up the configuration script.
        config_creator (ConfigCreator): An instance of the ConfigCreator class used to add different
                                        configuration sections to the script.
    """

    def __init__(self, session, uploaded_files):
        """
        Initializes the ConfigWriter with a session and uploaded files.

        Args:
            session (dict): A dictionary containing the session settings and parameters.
            uploaded_files (dict): A dictionary containing the uploaded file paths and names.
        """
        self.script = []
        self.config_creator = ConfigCreator(session, uploaded_files)

    def create_config_script(self):
        """
        Generates the complete configuration script by sequentially adding various configuration
        sections. The method integrates several configurations such as PDB input files, Amber input
        files, forcefield, solvent, system settings, integration options, simulation parameters,
        equilibration steps, postprocessing and analysis options

        Returns:
            str: The complete configuration script as a single string.
        """
        # OpenMMDL Logo
        self.config_creator.add_openmmdl_ascii_art_logo(self.script)

        # Config Logo
        self.config_creator.add_ascii_config_art(self.script)

        # PDB Input Files
        self.config_creator.add_pdb_input_files_configuration(self.script)

        # Amber Input Files
        self.config_creator.add_amber_file_configuration(self.script)

        # Forcefield and Water Model Selection
        self.config_creator.add_forcefield_and_water_model_configuration(self.script)

        # Add Solvent or Membrane Configuration
        self.config_creator.add_solvent_configuration(self.script)

        # System Configuration
        self.config_creator.add_system_configuration(self.script)

        # Integration Options
        self.config_creator.add_integration_configuration(self.script)

        # Simulation Time & Step Configuration
        self.config_creator.add_simulation_time_and_steps_configuration(self.script)

        # Equilibartion Options
        self.config_creator.add_equilibration_configuration(self.script)

        # Simulation Configuration
        self.config_creator.add_simulation_configuration(self.script)

        # Checkpoint Options
        self.config_creator.add_checkpoint_configuration(self.script)

        # XML Options
        self.config_creator.add_xml_serialization_configuration(self.script)

        # Postprocessing Options
        self.config_creator.add_postprocessing_configuration(self.script)

        # OpenMMDL Analysis Options
        self.config_creator.add_openmmdl_analysis_configuration(self.script)

        return "\n".join(self.script)


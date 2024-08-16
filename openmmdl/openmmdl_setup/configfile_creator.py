import datetime
from openmmdl.openmmdl_setup.file_operator import LigandExtractor

from typing import List, Dict, Optional
from typing import Dict, List


class ConfigCreator:
    """
    ConfigCreator is responsible for generating different sections of a simulation
    configuration script based on the session parameters and uploaded files.

    Attributes:
        session (Dict[str, str]): A dictionary containing user-defined session parameters.
        uploadedFiles (Dict[str, List[tuple[str, str]]]): A dictionary containing information about the files uploaded by the user.
    """

    def __init__(
        self, session: Dict[str, str], uploadedFiles: Dict[str, List[tuple[str, str]]]
    ):
        """
        Initializes the ConfigCreator with session parameters and uploaded files from OpenMMDL Setup.

        Args:
            session (Dict[str, str]): A dictionary containing user-defined session parameters.
            uploadedFiles (Dict[str, List[tuple[str, str]]]): A dictionary containing information about the files uploaded by the user during the usage of OpenMMDL Setup.
        """
        self.session = session
        self.uploadedFiles = uploadedFiles

    def add_openmmdl_ascii_art_logo(self, script: List[str]) -> None:
        """
        Adds the OpenMMDL logo to the configuration file.

        Args:
            script (List[str]): A list of strings representing the configuration file being built.
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

    def add_ascii_config_art(self, script: List[str]) -> None:
        """
        Adds a Config File ASCII art header for the configuration file.

        Args:
            script (List[str]): A list of strings representing the configuration script being built.
        """
        script.append(
            """
                            __   __        ___    __      ___         ___ 
                           /  ` /  \ |\ | |__  | / _`    |__  | |    |__  
                           \__, \__/ | \| |    | \__>    |    | |___ |___                                                                                                                              
            """
        )

    def add_pdb_input_files_configuration(self, script: List[str]) -> None:
        """
        Adds configuration settings for PDB input files, including protein and ligand data.
        Handles both protein-only and protein-ligand systems, setting up the necessary variables
        for simulation input.

        Args:
            script (List[str]): A list of strings representing the configuration script being built.
        """
        script.append("\n# Input Files")
        if self.session["fileType"] == "pdb":
            script.append(
                """############# Ligand and Protein Data ###################"""
            )
            if self.session["pdbType"] == "pdb":
                script.append('input_file_type = "pdb"')
                script.append('protein = "%s"' % self.uploadedFiles["file"][0][1])
                if self.session["sdfFile"]:
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

    def add_amber_file_configuration(self, script: List[str]) -> None:
        """
        Adds configuration settings for Amber input files, including topology and coordinate files.
        This method handles both scenarios where Amber files are pre-existing or need to be generated.
        It also extracts ligand names if ligands are present.

        Args:
            script (List[str]): A list of strings representing the configuration script being built.
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
                nmLigName: Optional[str] = (
                    self.session["nmLigName"] if self.session["nmLig"] else None
                )
                spLigName: Optional[str] = (
                    self.session["spLigName"] if self.session["spLig"] else None
                )

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

            # Feed prmtop_file and inpcrd_file to OpenMM Reader
            script.append("prmtop = AmberPrmtopFile(prmtop_file)")
            script.append("inpcrd = AmberInpcrdFile(inpcrd_file)")

    def add_forcefield_and_water_model_configuration(self, script: List[str]) -> None:
        """
        Adds forcefield and water model configuration settings to the script based on the file type.
        This includes selecting the appropriate forcefield and water model if applicable.

        Args:
            script (List[str]): A list of strings representing the configuration script being built.
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

    def add_solvent_configuration(self, script: List[str]) -> None:
        """
        Adds solvent or membrane configuration settings to the script, depending on whether
        a membrane or water box is being used, including parameters like padding, ionic strength,
        and ion types.

        Args:
            script (List[str]): A list of strings representing the configuration script being built.
        """
        if self.session["fileType"] == "pdb":
            if self.session["solvent"]:
                if self.session["add_membrane"]:
                    script.append(
                        "\n############# Membrane Settings ###################\n"
                    )
                    script.append("add_membrane = %s" % self.session["add_membrane"])
                    script.append(
                        "membrane_lipid_type = '%s'"
                        % self.session["membrane_lipid_type"]
                    )
                else:
                    script.append(
                        "\n############# Solvent Settings ###################\n"
                    )
                    script.append("padding = %s" % self.session["padding"])
                    script.append(
                        "ionic_strength = %s" % self.session["ionic_strength"]
                    )
                    script.append("cation = %s" % self.session["cation"])
                    script.append("anion = %s" % self.session["anion"])


class ConfigWriter:
    """
    This class is responsible for generating a configuration script based on the provided session
    and uploaded files. It acts as a wrapper around the ConfigCreator class, orchestrating the
    various configuration sections that need to be included in the final script.

    Attributes:
        script (List[str]): A list of strings that make up the configuration script.
        config_creator (ConfigCreator): An instance of the ConfigCreator class used to add different
                                        configuration sections to the script.
    """

    def __init__(self, session: Dict[str, any], uploaded_files: Dict[str, str]) -> None:
        """
        Initializes the ConfigWriter with a session and uploaded files.

        Args:
            session (Dict[str, any]): A dictionary containing the session settings and parameters.
            uploaded_files (Dict[str, str]): A dictionary containing the uploaded file paths and names.
        """
        self.script: List[str] = []
        self.config_creator = ConfigCreator(session, uploaded_files)

    def create_config_script(self) -> str:
        """
        Generates the complete configuration script by sequentially adding various configuration
        sections. The method integrates several configurations such as PDB input files, Amber input
        files, forcefield, solvent, system settings, integration options, simulation parameters,
        equilibration steps, postprocessing and analysis options.

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

        # Equilibration Options
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

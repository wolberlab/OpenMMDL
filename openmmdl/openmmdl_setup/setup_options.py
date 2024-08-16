from flask import session, request


from typing import TypedDict, List, Union

class SessionDict(TypedDict, total=False):
    restart_checkpoint: bool
    mdtraj_output: str
    mdtraj_removal: str
    mda_output: str
    mda_selection: str
    openmmdl_analysis: str
    analysis_selection: str
    binding_mode: str
    min_transition: str
    rmsd_diff: str
    pml_generation: str
    stable_water: str
    wc_distance: str
    fileType: str
    waterModel: str
    ensemble: str
    platform: str
    precision: str
    cutoff: str
    ewaldTol: str
    constraintTol: str
    hmr: bool
    hmrMass: str
    dt: str
    sim_length: str
    equilibration: str
    temperature: str
    friction: str
    pressure: str
    barostatInterval: str
    nonbondedMethod: str
    writeDCD: bool
    dcdFilename: str
    dcdFrames: str
    pdbInterval_ns: str
    writeData: bool
    dataFilename: str
    dataInterval: str
    dataFields: List[str]
    writeCheckpoint: bool
    checkpointFilename: str
    checkpointInterval_ns: str
    writeSimulationXml: bool
    systemXmlFilename: str
    integratorXmlFilename: str
    writeFinalState: bool
    finalStateFileType: str
    finalStateFilename: str
    constraints: str
    rmsd: str
    md_postprocessing: str
    nmLig: str
    spLig: str
    lig_ff: str
    charge_value: str
    charge_method: str
    prot_ff: str
    dna_ff: str
    rna_ff: str
    carbo_ff: str
    addType: str
    boxType: str
    dist: str
    lipid_tp: str
    other_lipid_tp_input: str
    lipid_ratio: str
    lipid_ff: str
    dist2Border: str
    padDist: str
    water_ff: str
    pos_ion: str
    neg_ion: str
    ionConc: str

class SetupOptionsConfigurator:
    """
    Configures and initializes default simulation options based on Flask session data.

    This class is responsible for setting up various simulation parameters such as
    ensemble types, platform, precision, cutoffs, force fields, and other essential
    options. It primarily deals with setting default values for these parameters,
    which can then be used throughout the simulation process.

    Attributes:
        session (SessionDict): The session object used to store configuration data
                               for the simulation setup.
    """

    def __init__(self, session: SessionDict) -> None:
        """
        Initializes the SetupOptionsConfigurator with the provided session object.

        Args:
            session (SessionDict): The session object to store simulation configuration data.
        """
        self.session: SessionDict = session

    def configure_default_options(self) -> None:
        """
        Sets default simulation options based on session data.

        Updates parameters for ensemble type, platform, cutoff, tolerances,
        simulation length, and output files.
        """
        implicitWater: bool = False
        self.session["restart_checkpoint"] = False
        self.session["mdtraj_output"] = "mdtraj_pdb_dcd"
        self.session["mdtraj_removal"] = "False"
        self.session["mda_output"] = "mda_pdb_dcd"
        self.session["mda_selection"] = "mda_prot_lig_all"
        self.session["openmmdl_analysis"] = "No"
        self.session["analysis_selection"] = "analysis_all"
        self.session["binding_mode"] = "40"
        self.session["min_transition"] = "1"
        self.session["rmsd_diff"] = "No"
        self.session["pml_generation"] = "True"
        self.session["stable_water"] = "Yes"
        self.session["wc_distance"] = "1.0"

        if (
            self.session["fileType"] == "pdb"
            and self.session["waterModel"] == "implicit"
        ):
            implicitWater = True

        self.session["ensemble"] = "nvt" if implicitWater else "npt"
        self.session["platform"] = "CUDA"
        self.session["precision"] = "mixed"
        self.session["cutoff"] = "2.0" if implicitWater else "1.0"
        self.session["ewaldTol"] = "0.0005"
        self.session["constraintTol"] = "0.000001"
        self.session["hmr"] = True
        self.session["hmrMass"] = "1.5"
        self.session["dt"] = "0.002"
        self.session["sim_length"] = "50"
        self.session["equilibration"] = "equilibration"
        self.session["temperature"] = "300"
        self.session["friction"] = "1.0"
        self.session["pressure"] = "1.0"
        self.session["barostatInterval"] = "25"
        self.session["nonbondedMethod"] = (
            "CutoffNonPeriodic" if implicitWater else "PME"
        )
        self.session["writeDCD"] = True
        self.session["dcdFilename"] = "trajectory.dcd"
        self.session["dcdFrames"] = "1000"
        self.session["pdbInterval_ns"] = "10"
        self.session["writeData"] = True
        self.session["dataFilename"] = "log.txt"
        self.session["dataInterval"] = "1000"
        self.session["dataFields"] = [
            "step",
            "speed",
            "progress",
            "potentialEnergy",
            "temperature",
        ]
        self.session["writeCheckpoint"] = True
        self.session["checkpointFilename"] = "checkpoint.chk"
        self.session["checkpointInterval_ns"] = "0.02"
        self.session["writeSimulationXml"] = False
        self.session["systemXmlFilename"] = "system.xml"
        self.session["integratorXmlFilename"] = "integrator.xml"
        self.session["writeFinalState"] = False
        self.session["finalStateFileType"] = "stateXML"
        self.session["finalStateFilename"] = "final_state.xml"
        self.session["constraints"] = "hbonds"
        self.session["rmsd"] = "True"
        self.session["md_postprocessing"] = "True"

    def configureDefaultAmberOptions(self) -> None:
        """
        Sets default options for ligand, receptor, and solvation in Amber.

        Configures force fields, ion types, water model, and lipid properties.
        """
        # Ligand
        self.session["nmLig"] = ""
        self.session["spLig"] = ""
        self.session["lig_ff"] = "gaff2"
        self.session["charge_value"] = "0"
        self.session["charge_method"] = "bcc"

        # Receptor
        self.session["prot_ff"] = "ff19SB"
        self.session["dna_ff"] = "OL15"
        self.session["rna_ff"] = "OL3"
        self.session["carbo_ff"] = "GLYCAM_06j"

        # AddWaterMembrane
        self.session["addType"] = "addWater"
        self.session["boxType"] = "cube"
        self.session["dist"] = "10"

        self.session["lipid_tp"] = "POPC"
        self.session["other_lipid_tp_input"] = "POPC:TOPC"
        self.session["lipid_ratio"] = "1:1"
        self.session["lipid_ff"] = "lipid21"
        self.session["dist2Border"] = "15"
        self.session["padDist"] = "17"

        self.session["water_ff"] = "opc"
        self.session["pos_ion"] = "Na+"
        self.session["neg_ion"] = "Cl-"
        self.session["ionConc"] = "0.15"


from typing import Any, Dict, Optional
from werkzeug.datastructures import ImmutableMultiDict

class RequestSessionManager:
    """
    Manages the configuration of session variables based on form data received in HTTP requests.

    This class is responsible for setting up diverse options related to receptors, water, membranes,
    and other general simulation settings. The configurations are stored in the session object,
    which tracks the state of the simulation setup across different user interactions.

    Attributes:
        form (ImmutableMultiDict[str, str]): The form data from an HTTP request.
    """

    def __init__(self, form: ImmutableMultiDict[str, str]) -> None:
        """
        Initializes the RequestSessionManager with form data.

        Args:
            form (ImmutableMultiDict[str, str]): The form data from an HTTP request.
        """
        self.form: ImmutableMultiDict[str, str] = form

    def setAmberOptions_rcp_session(self) -> None:
        """
        Sets receptor-related Amber options from form data.
        """
        session["rcpType"] = self.form.get("rcpType", "")
        session["prot_ff"] = self.form.get("prot_ff", "")
        session["other_prot_ff_input"] = self.form.get("other_prot_ff_input", "")
        session["dna_ff"] = self.form.get("dna_ff", "")
        session["other_dna_ff_input"] = self.form.get("other_dna_ff_input", "")
        session["rna_ff"] = self.form.get("rna_ff", "")
        session["other_rna_ff_input"] = self.form.get("other_rna_ff_input", "")
        session["carbo_ff"] = self.form.get("carbo_ff", "")
        session["other_carbo_ff_input"] = self.form.get("other_carbo_ff_input", "")

    def setAmberOptions_water_membrane_session(self) -> None:
        """
        Sets water and membrane-related Amber options from form data.
        """
        session["addType"] = self.form.get("addType", "")
        session["boxType"] = self.form.get("boxType", "")
        session["dist"] = self.form.get("dist", "")
        session["lipid_tp"] = self.form.get("lipid_tp", "")
        session["other_lipid_tp_input"] = self.form.get("other_lipid_tp_input", "")
        session["lipid_ratio"] = self.form.get("lipid_ratio", "")
        session["lipid_ff"] = self.form.get("lipid_ff", "")
        session["dist2Border"] = self.form.get("dist2Border", "")
        session["padDist"] = self.form.get("padDist", "")
        session["water_ff"] = self.form.get("water_ff", "")
        session["pos_ion"] = self.form.get("pos_ion", "")
        session["neg_ion"] = self.form.get("neg_ion", "")
        session["ionConc"] = self.form.get("ionConc", "")

    def addhydrogens_add_water_or_membrane_session(self) -> None:
        """
        Configures water or membrane settings based on form data.
        """
        if "addWater" in self.form:
            session["solvent"] = True
            session["add_membrane"] = False

            if self.form["boxType"] == "geometry":
                session["water_padding"] = True
                session["water_padding_distance"] = self._parse_float(self.form.get("geomPadding"))
                session["water_boxShape"] = self.form["geometryDropdown"]
            else:
                session["water_padding"] = False
                session["box_x"] = self._parse_float(self.form.get("boxx"))
                session["box_y"] = self._parse_float(self.form.get("boxy"))
                session["box_z"] = self._parse_float(self.form.get("boxz"))
                
            session["water_ionicstrength"] = self._parse_float(self.form.get("ionicstrength"))
            session["water_positive"] = self.form.get("positiveion", "") + "+"
            session["water_negative"] = self.form.get("negativeion", "") + "-"

        elif "addMembrane" in self.form:
            session["solvent"] = True
            session["add_membrane"] = True
            session["lipidType"] = self.form["lipidType"]
            session["membrane_padding"] = self._parse_float(self.form.get("membranePadding"))
            session["membrane_ionicstrength"] = self._parse_float(self.form.get("ionicstrength"))
            session["membrane_positive"] = self.form.get("positiveion", "") + "+"
            session["membrane_negative"] = self.form.get("negativeion", "") + "-"

    def simulationoptions_add_general_settings(self) -> None:
        """
        Adds general simulation settings from form data to the session.
        """
        for key in self.form:
            session[key] = self.form[key]
        session["ligand"] = "ligand" in self.form
        session["writeDCD"] = "writeDCD" in self.form
        session["writeData"] = "writeData" in self.form
        session["writeCheckpoint"] = "writeCheckpoint" in self.form
        session["dataFields"] = self.form.getlist("dataFields")
        session["hmr"] = "hmr" in self.form
        session["writeSimulationXml"] = "writeSimulationXml" in self.form
        session["writeFinalState"] = "writeFinalState" in self.form

    def configureFiles_add_forcefield_ligand_settings(self) -> None:
        """
        Adds forcefield and ligand-related settings from form data to the session.
        """
        session["forcefield"] = self.form.get("forcefield", "")
        session["ml_forcefield"] = self.form.get("ml_forcefield", "")
        session["waterModel"] = self.form.get("waterModel", "")
        session["smallMoleculeForceField"] = self.form.get("smallMoleculeForceField", "")
        session["ligandMinimization"] = self.form.get("ligandMinimization", "")
        session["ligandSanitization"] = self.form.get("ligandSanitization", "")

    def _parse_float(self, value: Optional[str]) -> Optional[float]:
        """
        Safely parses a string to a float, returning None if parsing fails.

        Args:
            value (Optional[str]): The string value to parse.

        Returns:
            Optional[float]: The parsed float, or None if parsing fails.
        """
        try:
            return float(value) if value is not None else None
        except (ValueError, TypeError):
            return None

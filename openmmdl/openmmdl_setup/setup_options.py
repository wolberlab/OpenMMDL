from flask import session, request

class SetupOptionsConfigurator:
    def __init__(self, session):
        self.session = session

    def configure_default_options(self):
        """
        Sets default simulation options based on session data.
        
        Updates parameters for ensemble type, platform, cutoff, tolerances,
        simulation length, and output files.
        """
        implicitWater = False
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
        
        if self.session["fileType"] == "pdb" and self.session["waterModel"] == "implicit":
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
        self.session["nonbondedMethod"] = "CutoffNonPeriodic" if implicitWater else "PME"
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

    def configureDefaultAmberOptions(self):
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
        
        
class RequestSessionManager:
    def __init__(self, form):
        self.form = form

    def setAmberOptions_rcp_session(self):
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

    def setAmberOptions_water_membrane_session(self):
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
        
    def addhydrogens_add_water_or_membrane_session(self):
        """
        Configures water or membrane settings based on form data.
        """
        if "addWater" in self.form:
            session["solvent"] = True
            session["add_membrane"] = False
            padding, boxSize, boxShape = None, None, None

            if self.form["boxType"] == "geometry":
                session["water_padding"] = True
                session["water_padding_distance"] = float(self.form["geomPadding"])
                session["water_boxShape"] = self.form["geometryDropdown"]
            else:
                session["water_padding"] = False
                session["box_x"] = float(self.form["boxx"])
                session["box_y"] = float(self.form["boxy"])
                session["box_z"] = float(self.form["boxz"])
                boxSize = (
                    float(self.form["boxx"]),
                    float(self.form["boxy"]),
                    float(self.form["boxz"]),
                )
            session["water_ionicstrength"] = float(self.form["ionicstrength"])
            session["water_positive"] = self.form["positiveion"] + "+"
            session["water_negative"] = self.form["negativeion"] + "-"
        
        elif "addMembrane" in self.form:
            session["solvent"] = True
            session["add_membrane"] = True
            session["lipidType"] = self.form["lipidType"]
            session["membrane_padding"] = float(self.form["membranePadding"])
            session["membrane_ionicstrength"] = float(self.form["ionicstrength"])
            session["membrane_positive"] = self.form["positiveion"] + "+"
            session["membrane_negative"] = self.form["negativeion"] + "-"
            
    def simulationoptions_add_general_settings(self):
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

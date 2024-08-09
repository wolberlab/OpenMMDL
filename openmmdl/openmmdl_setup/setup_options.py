

class SetupOptionsConfigurator:
    def __init__(self, session):
        self.session = session

    def configure_default_options(self):
        """Select default options based on the file format and force field."""
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
        self.session["sim_length"] = "30"
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
        """Select default options based on the file format and force field."""
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

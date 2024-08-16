class AmberScriptGenerator:
    """
    Generates the bash script for setting up molecular dynamics simulations using Amber tools.

    Attributes:
        session (dict): A dictionary containing session-specific parameters and settings.
        uploadedFiles (dict): A dictionary containing information about uploaded files, organized by file type.
    """

    def __init__(self, session, uploadedFiles):
        """
        Initializes the AmberScriptGenerator with session parameters and uploaded files from OpenMMDL Setup.

        Args:
            session (dict): Session-specific parameters and settings.
            uploadedFiles (dict): Uploaded files organized by their respective types.
        """

        self.session = session
        self.uploadedFiles = uploadedFiles

    def add_openmmdl_logo(self, amber_script):
        """
        Appends the OpenMMDL logo to the Amber bash script.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        logo = """
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
        amber_script.append(logo)

    def add_openmmdl_amber_logo(self, amber_script):
        """
        Appends the OpenMMDL Amber name to the Amber script.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        logo = """
#				      _              _               
#				     / \   _ __ ___ | |__   ___ _ __ 
#				    / _ \ | '_ ` _ \| '_ \ / _ \ '__|
#				   / ___ \| | | | | | |_) |  __/ |   
#				  /_/   \_\_| |_| |_|_.__/ \___|_|    
        """
        amber_script.append(logo)

    def add_receptor_type(self, amber_script):
        """
        Determines the receptor type and appends corresponding Ambertools commands to the Amber bash script.

        This method adds commands to the Amber bash script based on the type of receptor specified.
        It handles various receptor types including protein, DNA, RNA, and carbohydrate, and sets up
        the appropriate force field for each type. It also updates the receptor file names and force fields
        based on the user's input and session settings.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        amber_script.append("#!/bin/bash\n")
        amber_script.append(
            "################################## Receptor ######################################\n"
        )

        rcpType = self.session["rcpType"]

        if rcpType == "protRcp":
            protFile = self.uploadedFiles["protFile"][0][1][:-4]
            amber_script.append(
                f"rcp_nm={protFile} # the file name of ligand without suffix `pdb`"
            )

            prot_ff = self.session["prot_ff"]
            if prot_ff != "other_prot_ff":
                amber_script.append(f"rcp_ff={prot_ff}")
            else:
                amber_script.append(
                    f"rcp_ff={self.session['other_prot_ff_input']}  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                )
            amber_script.append("\n")

        elif rcpType == "dnaRcp":
            dnaFile = self.uploadedFiles["dnaFile"][0][1][:-4]
            amber_script.append(
                f"rcp_nm={dnaFile} # the file name of ligand without suffix `pdb`"
            )

            dna_ff = self.session["dna_ff"]
            if dna_ff != "other_dna_ff":
                amber_script.append(f"rcp_ff={dna_ff}")
            else:
                amber_script.append(
                    f"rcp_ff={self.session['other_dna_ff_input']}  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                )
            amber_script.append("\n")

        elif rcpType == "rnaRcp":
            rnaFile = self.uploadedFiles["rnaFile"][0][1][:-4]
            amber_script.append(
                f"rcp_nm={rnaFile} # the file name of ligand without suffix `pdb`"
            )

            rna_ff = self.session["rna_ff"]
            if rna_ff != "other_rna_ff":
                amber_script.append(f"rcp_ff={rna_ff}")
            else:
                amber_script.append(
                    f"rcp_ff={self.session['other_rna_ff_input']}  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                )
            amber_script.append("\n")

        elif rcpType == "carboRcp":
            carboFile = self.uploadedFiles["carboFile"][0][1][:-4]
            amber_script.append(
                f"rcp_nm={carboFile} # the file name of ligand without suffix `pdb`"
            )

            carbo_ff = self.session["carbo_ff"]
            if carbo_ff != "other_carbo_ff":
                amber_script.append(f"rcp_ff={carbo_ff}")
            else:
                amber_script.append(
                    f"rcp_ff={self.session['other_carbo_ff_input']}  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                )
            amber_script.append("\n")

    def add_clean_pdb_commands(self, amber_script):
        """
        Appends commands to clean the PDB file using pdb4amber to the Amber bash script.

        This method adds commands to clean up the receptor PDB file using `pdb4amber`, including:
            - Generating an amber-compatible PDB file by removing unnecessary atoms.
            - Ensuring that the PDB file is cleaned of `CONECT` records and non-essential atoms
              to ensure compatibility with TLEAP and proper processing in the Amber simulation.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        amber_script.append("## Clean the PDB file by pdb4amber")
        amber_script.append("pdb4amber -i ${rcp_nm}.pdb -o ${rcp_nm}_amber.pdb")
        amber_script.append(
            """
## `tleap` requires that all residues and atoms have appropriate types to ensure compatibility with the specified force field.
## To avoid `tleap` failing, we delete non-essential atoms, such as hydrogens, but preserve important atoms like carbon and nitrogen within the caps residues.
## Don' worry about the missing atoms as tleap has the capability to reconstruct them automatically."""
        )
        amber_script.append(
            """awk '! ($2 ~ "(CH3|HH31|HH32|HH33)" || $3 ~ "(CH3|HH31|HH32|HH33)" )' ${rcp_nm}_amber.pdb > ${rcp_nm}_amber_f.pdb"""
        )
        amber_script.append(
            "grep -v '^CONECT' ${rcp_nm}_amber_f.pdb > ${rcp_nm}_cnt_rmv.pdb\n"
        )

    def add_ligand_commands(self, amber_script):
        """
        Appends commands related to ligand processing to the Amber bash script, handling both normal and special ligands.

        This method adds commands for processing ligands in the Amber simulation setup. It handles:
            - Normal ligands compatible with the GAFF force field by converting file formats and generating necessary parameter files.
            - Special ligands that are incompatible with GAFF, including file cleaning and processing steps.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        if self.session.get("nmLig") or self.session.get("spLig"):
            amber_script.append(
                "################################## Ligand ######################################"
            )

        if self.session.get("nmLig"):
            amber_script.append(
                "# Normal Ligand that is compatible with GAFF force field"
            )
            nmLigFile = self.uploadedFiles["nmLigFile"][0][1]
            amber_script.append(
                f"nmLigFile={nmLigFile[:-4]} # the file name of ligand without suffix `.pdb` or `.sdf`"
            )

            # Depending on the uploaded file format, convert it to pdb or sdf file.
            if nmLigFile.endswith(".sdf"):
                amber_script.append(
                    "obabel ${nmLigFile}.sdf -O ${nmLigFile}.pdb -p # convert to pdb file for tleap, -p: add hydrogens appropriate for pH7.4"
                )
            elif nmLigFile.endswith(".pdb"):
                amber_script.append(
                    "obabel ${nmLigFile}.pdb -O ${nmLigFile}.sdf -p # convert to sdf file for openmmdl_analysis, -p: add hydrogens appropriate for pH7.4"
                )

            amber_script.append(
                f"charge_method={self.session['charge_method']} # refers to the charge method that antechamber will adopt"
            )
            amber_script.append(
                f"charge_value={self.session['charge_value']} # Enter the net molecular charge of the ligand as integer (e.g. 1 or -2)"
            )
            amber_script.append(
                f"lig_ff={self.session['lig_ff']} # Ligand force field\n"
            )

            amber_script.append("## Clean the PDB file by pdb4amber")
            amber_script.append(
                "pdb4amber -i ${nmLigFile}.pdb -o ${nmLigFile}_amber.pdb\n"
            )

            amber_script.append(
                "## Generate a prepc file and an additional frcmod file by `antechamber`"
            )
            amber_script.append(
                "antechamber -fi pdb -fo prepc -i ${nmLigFile}_amber.pdb -o ${nmLigFile}.prepc -c ${charge_method} -at ${lig_ff} -nc ${charge_value} -pf y"
            )
            amber_script.append(
                "parmchk2 -f prepc -i ${nmLigFile}.prepc -o ${nmLigFile}.frcmod\n"
            )
            amber_script.append("## Rename ligand pdb")
            amber_script.append(
                "antechamber -i ${nmLigFile}.prepc -fi prepc -o rename_${nmLigFile}.pdb -fo pdb\n"
            )

        if self.session.get("spLig"):
            amber_script.append(
                "# Special Ligand that is incompatible with GAFF force field"
            )
            spLigFile = self.uploadedFiles["spLigFile"][0][1]
            amber_script.append(
                f"spLigFile={spLigFile[:-4]} # the file name of ligand without suffix `.pdb`"
            )

            prepcFile = self.uploadedFiles["prepcFile"][0][1][:-6]
            amber_script.append(
                f"prepc={prepcFile} # the file name without suffix `prepc`"
            )

            frcmodFile = self.uploadedFiles["frcmodFile"][0][1][:-7]
            amber_script.append(
                f"frcmod={frcmodFile} # the file name without suffix `frcmod`\n"
            )

            amber_script.append("## Clean the PDB file by pdb4amber")
            amber_script.append(
                "pdb4amber -i ${spLigFile}.pdb -o ${spLigFile}_amber.pdb\n"
            )

            # Get the name of ligand in the pdb file: it is the fourth column of the first line
            amber_script.append(
                "spLigName=$(awk 'NR==1 {print $4}' ${spLigFile}_amber.pdb)\n"
            )

    def add_combine_components_commands(self, amber_script):
        """
        Appends commands to combine all components (receptor and ligands) into a single complex.

        This method generates commands to combine various components (receptor and ligands) into a single
        complex in the Amber simulation. It handles loading of the receptor and ligands (normal and/or special),
        and combines them using TLEAP. The resulting complex is saved in PDB format and cleaned by removing
        CONECT lines.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        if self.session.get("nmLig") or self.session.get("spLig"):
            amber_script.append(
                "######################  Combine All Components to Be Modelled ####################"
            )
            amber_script.append("cat > tleap.combine.in <<EOF\n")
            amber_script.append("source ${rcp_ff}")
            amber_script.append("source leaprc.${lig_ff}")

            # Load the prepc and frcmod file for either normal or special ligand
            if self.session.get("nmLig"):
                amber_script.append("loadamberprep ${nmLigFile}.prepc")
                amber_script.append("loadamberparams ${nmLigFile}.frcmod\n")
            if self.session.get("spLig"):
                amber_script.append("loadamberprep ${prepc}.prepc")
                amber_script.append("loadamberparams ${frcmod}.frcmod\n")

            # Load both receptor and ligand pdb file
            amber_script.append("rcp = loadpdb ${rcp_nm}_cnt_rmv.pdb")
            if self.session.get("nmLig") and self.session.get("spLig"):
                amber_script.append("nmLig = loadpdb rename_${nmLigFile}.pdb ")
                amber_script.append("spLig = loadpdb ${spLigFile}_amber.pdb ")
                amber_script.append("comp = combine{rcp nmLig spLig}")
            elif self.session.get("nmLig"):
                amber_script.append("nmLig = loadpdb rename_${nmLigFile}.pdb ")
                amber_script.append("comp = combine{rcp nmLig}")
            elif self.session.get("spLig"):
                amber_script.append("spLig = loadpdb ${spLigFile}_amber.pdb")
                amber_script.append("comp = combine{rcp spLig}")

            amber_script.append("savepdb comp comp.pdb")
            amber_script.append("quit\nEOF\n")
            amber_script.append("tleap -s -f tleap.combine.in > tleap.combine.out")
            # Remove 'CONECT' line in the pdb file
            amber_script.append("grep -v '^CONECT' comp.pdb > comp_cnt_rmv.pdb\n")

    def add_solvation_commands(self, amber_script):
        """
        Appends commands related to solvation settings, including box type and dimensions for the water box.

        This method generates the necessary bash commands to set up the solvation box around the solute
        in an Amber simulation. It supports various types of solvent boxes such as cubic, truncated octahedral,
        spherical cap, and shell, and includes parameters for the box dimensions or shell thickness.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        addType = self.session.get("addType")
        if addType == "addWater":
            boxType = self.session.get("boxType")
            if boxType == "cube":
                amber_script.append(
                    "boxType=solvatebox # `solvatebox`, a command in tleap, creates a cubic box "
                )
                amber_script.append(
                    f"dist={self.session['dist']} # the minimum distance between any atom originally present in solute and the edge of the periodic box."
                )
            elif boxType == "octahedron":
                amber_script.append(
                    "boxType=solvateoct # `solvateoct`, a command in tleap, creates a truncated octahedron box."
                )
                amber_script.append(
                    f"dist={self.session['dist']} # the minimum distance between any atom originally present in solute and the edge of the periodic box"
                )
            elif boxType == "cap":
                amber_script.append(
                    "boxType=solvatecap # `solvatecap`, a command in tleap, creates a solvent cap around solute. In development!"
                )
                amber_script.append(
                    f"radius={self.session['dist']} # the radius of the sphere"
                )
            elif boxType == "shell":
                amber_script.append(
                    "boxType=solvateshell # `solvatecap`, a command in tleap, adds a solvent shell to solute, which reflects the contours of the original solute molecule."
                )
                amber_script.append(
                    f"thickness={self.session['dist']} # the thickness of the shell"
                )

    def add_membrane_commands(self, amber_script):
        """
        Appends commands related to membrane settings, including lipid types and force fields.

        This method generates bash commands required to configure the membrane system in an Amber simulation.
        It specifies lipid types, lipid force fields, and related parameters such as lipid ratios and distances
        for the membrane setup. Custom lipid types and force fields can also be specified.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        addType = self.session.get("addType")
        if addType == "addMembrane":
            lipid_tp = self.session.get("lipid_tp")
            if lipid_tp != "other_lipid_tp":
                amber_script.append(f"lipid_tp={lipid_tp}")
                amber_script.append("lipid_ratio=1")
            elif lipid_tp == "other_lipid_tp":
                amber_script.append(
                    f"lipid_tp={self.session['other_lipid_tp_input']}  # The command to check supported lipids: packmol-memgen --available_lipids"
                )
                amber_script.append(
                    f"lipid_ratio={self.session['lipid_ratio']} # Set to 1 if only one lipid required"
                )

            lipid_ff = self.session.get("lipid_ff")
            if lipid_ff != "other_lipid_ff":
                amber_script.append(f"lipid_ff={lipid_ff}")
            elif lipid_ff == "other_lipid_ff":
                amber_script.append(
                    f"lipid_ff={self.session['other_lipid_ff_input']}  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
                )

            amber_script.append(
                f"dist2Border={self.session['dist2Border']}  # The minimum distance between the maxmin values for x y and z to the box boundaries. Flag --dist"
            )
            amber_script.append(
                f"padDist={self.session['padDist']}  # The width of the water layer over the membrane or protein in the z axis. Flag --dist_wat"
            )

    def add_water_ff_commands(self, amber_script):
        """
        Appends commands related to water force field settings.

        This method generates the bash commands needed to configure the water force
        field (FF) in an Amber file preparation. Depending on the selected water model,
        it sets up the appropriate solvent box. It also accommodates custom water force
        fields specified by the user.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        water_ff = self.session.get("water_ff")
        addType = self.session.get("addType")

        if water_ff != "other_water_ff":
            amber_script.append(f"water_ff={water_ff}")
            if addType == "addWater":
                if water_ff == "tip3p":
                    amber_script.append(
                        f"solvent={water_ff.upper()}BOX  # set the water box"
                    )
                elif water_ff == "fb3":
                    amber_script.append("solvent=TIP3PFBOX # set the water box")
                elif water_ff == "spce":
                    amber_script.append("solvent=SPCBOX # set the water box")
                elif water_ff == "tip4pew":
                    amber_script.append("solvent=TIP4PEWBOX # set the water box")
                elif water_ff == "fb4":
                    amber_script.append("solvent=TIP4PBOX # set the water box")
                elif water_ff == "opc":
                    amber_script.append("solvent=OPCBOX # set the water box")
                elif water_ff == "opc3":
                    amber_script.append("solvent=OPC3BOX # set the water box")
        elif water_ff == "other_water_ff":
            amber_script.append(
                f"water_ff={self.session['other_water_ff_input']}  # See the supported force fields in the original file at `$AMBERHOME/dat/leap/cmd/`"
            )
            if addType == "addWater":
                amber_script.append(
                    f"solvent={self.session['other_water_ff_input'].upper()}BOX  # set the water box"
                )

    def add_ion_commands(self, amber_script):
        """
        Appends commands related to ion settings, including positive and negative ions.

        This method constructs the commands needed to set up ion types and their
        concentrations in the Amber simulation script. It handles standard ions
        (positive and negative) and allows for custom ion types. Additionally, it
        configures ion concentration based on the system setup, whether it includes
        water or a membrane.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        pos_ion = self.session.get("pos_ion")
        if pos_ion != "other_pos_ion":
            amber_script.append(f"pos_ion={pos_ion}")
        elif pos_ion == "other_pos_ion":
            amber_script.append(
                f"pos_ion={self.session['other_pos_ion_input']}  # In development!"
            )

        neg_ion = self.session.get("neg_ion")
        if neg_ion != "other_neg_ion":
            amber_script.append(f"neg_ion={neg_ion}")
        elif neg_ion == "other_neg_ion":
            amber_script.append(
                f"neg_ion={self.session['other_neg_ion_input']}  # In development!"
            )

        addType = self.session.get("addType")
        if addType == "addWater":
            amber_script.append(
                "numIon=0 # `numIon` is the flag for `addions` in tleap. When set to 0, the system will be neutralized"
            )
        elif addType == "addMembrane":
            amber_script.append(f"ionConc={self.session['ionConc']}")

        amber_script.append("\n")

    def add_membrane_building_commands(self, amber_script):
        """
        Appends commands to build the membrane using packmol-memgen.

        This method generates the necessary bash commands to build a membrane system in
        an Amber simulation using the packmol-memgen tool. It configures the membrane
        building process based on the session's ligand and system settings.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        addType = self.session.get("addType")
        if addType == "addMembrane":
            amber_script.append("## Build the membrane")
            if not self.session["nmLig"] and not self.session["spLig"]:
                amber_script.append(
                    "packmol-memgen --pdb ${rcp_nm}_cnt_rmv.pdb --lipids ${lipid_tp} --ratio ${lipid_ratio} --preoriented --dist ${dist2Border} --dist_wat ${padDist} --salt --salt_c ${pos_ion} --saltcon ${ionConc} --nottrim --overwrite --notprotonate\n"
                )
                amber_script.append(
                    "## Clean the complex pdb by `pdb4amber` for further `tleap` process"
                )
                amber_script.append(
                    "pdb4amber -i bilayer_${rcp_nm}_cnt_rmv.pdb -o clean_bilayer_${rcp_nm}.pdb"
                )
                # Remove 'CONECT' line in the pdb file
                amber_script.append(
                    "grep -v '^CONECT' clean_bilayer_${rcp_nm}.pdb > clean_bilayer_${rcp_nm}_cnt_rmv.pdb"
                )
                amber_script.append("\n")
            if self.session["nmLig"] or self.session["spLig"]:
                amber_script.append(
                    "packmol-memgen --pdb comp.pdb --lipids ${lipid_tp} --ratio ${lipid_ratio} --preoriented --dist ${dist2Border} --dist_wat ${padDist} --salt --salt_c ${pos_ion} --saltcon ${ionConc} --nottrim --overwrite --notprotonate\n"
                )
                amber_script.append(
                    "## Clean the complex pdb by `pdb4amber` for further `tleap` process"
                )
                amber_script.append(
                    "pdb4amber -i bilayer_comp.pdb -o clean_bilayer_comp.pdb"
                )
                # Remove 'CONECT' line in the pdb file
                amber_script.append(
                    "grep -v '^CONECT' clean_bilayer_comp.pdb > clean_bilayer_comp_cnt_rmv.pdb"
                )
                amber_script.append("\n")

    def generate_tleap_commands(self, amber_script):
        """
        Appends commands to generate TLEAP input files and execute TLEAP for system preparation.

        This method constructs the necessary TLEAP commands to prepare the system for simulation
        in Amber. It sources the appropriate force fields, loads the relevant ligand and system files,
        configures the box type, and saves the system's topology and coordinate files. The method
        also accounts for different system setups such as water-only or membrane-included simulations.

        Args:
            amber_script (list): The list representing lines of the Amber bash script.
        """
        amber_script.append("cat > tleap.in <<EOF\n")

        # Source the force field
        amber_script.append("source ${rcp_ff}")
        amber_script.append("source leaprc.water.${water_ff}")
        if self.session.get("nmLig") or self.session.get("spLig"):
            amber_script.append("source leaprc.${lig_ff}")
        if self.session.get("addType") == "addMembrane":
            amber_script.append("source leaprc.${lipid_ff}")

        # Load the prepc and frcmod files
        if self.session.get("nmLig"):
            amber_script.append("\nloadamberprep ${nmLigFile}.prepc")
            amber_script.append("loadamberparams ${nmLigFile}.frcmod\n")
        if self.session.get("spLig"):
            amber_script.append("loadamberprep ${prepc}.prepc")
            amber_script.append("loadamberparams ${frcmod}.frcmod\n")

        # Load the complex pdb which contains all the components to be modeled
        if self.session.get("addType") == "addWater":
            if not self.session.get("nmLig") and not self.session.get("spLig"):
                amber_script.append("\nsystem = loadpdb ${rcp_nm}_cnt_rmv.pdb\n")
            else:
                amber_script.append("system = loadpdb comp_cnt_rmv.pdb\n")
        elif self.session.get("addType") == "addMembrane":
            if not self.session.get("nmLig") and not self.session.get("spLig"):
                amber_script.append(
                    "\nsystem = loadpdb clean_bilayer_${rcp_nm}_cnt_rmv.pdb\n"
                )
            if self.session.get("nmLig") or self.session.get("spLig"):
                amber_script.append("system = loadpdb clean_bilayer_comp_cnt_rmv.pdb\n")

        # Add box commands based on the type of system
        if self.session.get("addType") == "addWater":
            boxType = self.session.get("boxType")
            if boxType == "cube":
                amber_script.append("solvatebox system ${solvent} ${dist} ")
            elif boxType == "octahedron":
                amber_script.append("solvateoct system ${solvent} ${dist}")
            elif boxType == "cap":
                amber_script.append("solvatecap system ${solvent} ${radius}")
            elif boxType == "shell":
                amber_script.append("solvateshell system ${solvent} ${thickness}")
            amber_script.append("addions2 system ${neg_ion} ${numIon}")
            amber_script.append("addions2 system ${pos_ion} ${numIon}")
        elif self.session.get("addType") == "addMembrane":
            amber_script.append('setBox system "vdw"')

        # Final TLEAP commands
        amber_script.append("check system")
        amber_script.append("charge system\n")

        # Save the PDB, prmtop, and inpcrd files
        amber_script.append("savepdb system system.${water_ff}.pdb")
        amber_script.append(
            "saveamberparm system system.${water_ff}.prmtop system.${water_ff}.inpcrd"
        )

        # Quit TLEAP and EOF
        amber_script.append("\nquit")
        amber_script.append("\nEOF")
        amber_script.append("\ntleap -s -f tleap.in > tleap.out")

 # Load the prepc and frcmod files
        if self.session.get("nmLig"):
            amber_script.append("\nloadamberprep ${nmLigFile}.prepc")
            amber_script.append("loadamberparams ${nmLigFile}.frcmod")

        if self.session.get("spLig"):
            amber_script.append("\nloadamberprep ${prepc}.prepc")
            amber_script.append("loadamberparams ${frcmod}.frcmod")

        # Load the receptor and/or ligand PDB files
        amber_script.append(f"\nloadpdb {self.session['rcp_nm']}_cnt_rmv.pdb")
        if self.session.get("nmLig"):
            amber_script.append(f"loadpdb rename_${{nmLigFile}}.pdb")
        if self.session.get("spLig"):
            amber_script.append(f"loadpdb ${spLigFile}_amber.pdb")

        # Combine the components
        if self.session.get("nmLig") and self.session.get("spLig"):
            amber_script.append("combine {rcp nmLig spLig}")
        elif self.session.get("nmLig"):
            amber_script.append("combine {rcp nmLig}")
        elif self.session.get("spLig"):
            amber_script.append("combine {rcp spLig}")

        amber_script.append("savepdb comp comp.pdb")

        # Set up the water or membrane
        if self.session.get("addType") == "addWater":
            amber_script.append(f"\nsolvate {self.session['boxType']} {self.session['dist']}")
            amber_script.append(f"savepdb comp_solvated.pdb")
        elif self.session.get("addType") == "addMembrane":
            amber_script.append(f"\nsource leaprc.{self.session['lipid_ff']}")
            amber_script.append(f"savepdb {self.session['rcp_nm']}_with_membrane.pdb")

        amber_script.append("quit\nEOF\n")
        amber_script.append("tleap -f tleap.in > tleap.out\n")


from typing import Dict, List
import logging

class AmberScriptGenerator:
    # Assuming the methods used in AmberScriptGenerator
    def __init__(self, session: Dict[str, any], uploaded_files: Dict[str, any]):
        pass
    
    def add_openmmdl_logo(self, script: List[str]) -> None:
        pass

    def add_openmmdl_amber_logo(self, script: List[str]) -> None:
        pass

    def add_receptor_type(self, script: List[str]) -> None:
        pass

    def add_clean_pdb_commands(self, script: List[str]) -> None:
        pass

    def add_ligand_commands(self, script: List[str]) -> None:
        pass

    def add_combine_components_commands(self, script: List[str]) -> None:
        pass

    def add_solvation_commands(self, script: List[str]) -> None:
        pass

    def add_membrane_commands(self, script: List[str]) -> None:
        pass

    def add_water_ff_commands(self, script: List[str]) -> None:
        pass

    def add_ion_commands(self, script: List[str]) -> None:
        pass

    def add_membrane_building_commands(self, script: List[str]) -> None:
        pass

    def generate_tleap_commands(self, script: List[str]) -> None:
        pass


class AmberBashScriptWriter:
    """
    Writes the complete Amber bash script by aggregating various components generated by AmberScriptGenerator.

    Attributes:
        script (List[str]): A list representing lines of the Amber bash script.
        amber_script (AmberScriptGenerator): An instance of AmberScriptGenerator to generate script components.
    """

    def __init__(self, session: Dict[str, any], uploaded_files: Dict[str, any]):
        """
        Initializes the AmberBashScriptWriter with session parameters and uploaded files from OpenMMDL Setup.

        Args:
            session (Dict[str, any]): Session-specific parameters and settings.
            uploaded_files (Dict[str, any]): Uploaded files organized by their respective types.
        """
        self.script: List[str] = []
        self.amber_script: AmberScriptGenerator = AmberScriptGenerator(session, uploaded_files)

    def write_amber_script(self) -> str:
        """
        Aggregates all components of the Amber bash script and returns the complete script as a string.

        Returns:
            str: The complete Amber bash script.
        """
        # Add the OpenMM/Amber logo
        self.amber_script.add_openmmdl_logo(self.script)

        # Add the Amber name
        self.amber_script.add_openmmdl_amber_logo(self.script)

        # Receptor
        self.amber_script.add_receptor_type(self.script)
        self.amber_script.add_clean_pdb_commands(self.script)

        # Ligand
        self.amber_script.add_ligand_commands(self.script)

        # Combine all components to be modeled
        self.amber_script.add_combine_components_commands(self.script)

        # Waterbox & Membrane setting
        self.amber_script.add_solvation_commands(self.script)
        self.amber_script.add_membrane_commands(self.script)

        # Water setting
        self.amber_script.add_water_ff_commands(self.script)

        # Ion setting
        self.amber_script.add_ion_commands(self.script)

        # Build the membrane
        self.amber_script.add_membrane_building_commands(self.script)

        # Generate the prmtop and frcmod file for the complex
        self.amber_script.generate_tleap_commands(self.script)

        return "\n".join(self.script)

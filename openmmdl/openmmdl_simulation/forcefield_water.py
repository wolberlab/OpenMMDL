import simtk.openmm.app as app
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import (
    GAFFTemplateGenerator,
    SMIRNOFFTemplateGenerator,
)


class ForcefieldSetup:
    def __init__(self, forcefield_name, water_model, ligand_file=None, minimization=False):
        """
        Initialize the SimulationRunner with the given forcefield, water model, and optional ligand parameters.

        Parameters:
        - forcefield_name (str): The name of the forcefield to use.
        - water_model (str): The water model to use.
        - ligand_file (str, optional): Path to the ligand file (e.g., SDF). Default is None.
        - minimization (bool): Whether to perform minimization on the ligand. Default is False.
        """
        self.forcefield_name = forcefield_name
        self.water_model = water_model
        self.ligand_file = ligand_file
        self.minimization = minimization

    def setup_forcefield(self):
        """
        Set up the forcefield and water model.

        Returns:
        - tuple: (forcefield, water_forcefield, model_water)
        """
        # Initialize ForcefieldPreparation
        prep = ForcefieldPreparation(
            forcefield_name=self.forcefield_name,
            water_model=self.water_model
        )

        # Select forcefield and water model
        forcefield = prep.select_forcefield()
        water_forcefield, model_water = prep.select_water_model()

        return forcefield, water_forcefield, model_water


class ForcefieldPreparation:
    def __init__(self, forcefield_name, water_model):
        """
        Initialize the setup for simulation.

        Parameters:
        - forcefield_name (str): The name of the forcefield to use.
        - water_model (str): The water model to use.
        """
        self.forcefield_name = forcefield_name
        self.water_model = water_model

        # Initialize the objects for forcefield and water model selection
        self.forcefield_selector = ForcefieldSelector()

    def select_forcefield(self):
        """Select the forcefield based on the provided name."""
        forcefield = self.forcefield_selector.ff_selection(self.forcefield_name)
        return forcefield

    def select_water_model(self):
        """Select the water model and forcefield."""
        print("Selecting forcefield...")
        forcefield = self.select_forcefield()
        print(f"Selected forcefield: {forcefield}")
        print("Selecting water forcefield and model...")
        water_forcefield = self.forcefield_selector.water_forcefield_selection(
            water=self.water_model, 
            forcefield_selection=forcefield
        )
        model_water = self.forcefield_selector.water_model_selection(
            water=self.water_model, 
            forcefield_selection=forcefield
        )
        print(f"Selected water forcefield: {water_forcefield}")
        print(f"Selected water model: {model_water}")
        return water_forcefield, model_water


class ForcefieldConfigurator:
    def __init__(self, config_parser, forcefield_selected, water_forcefield, water_selected, smallMoleculeForceField=None, prepared_ligand=None):
        """
        Initialize the ForcefieldConfigurator with configuration and generation parameters.

        Parameters:
        - config_parser (ConfigParser): The configuration parser instance.
        - forcefield_generator (ForcefieldGenerator): The forcefield generator instance.
        - forcefield_selected (str): The selected protein forcefield.
        - water_selected (str): The selected water forcefield.
        - smallMoleculeForceField (str): The forcefield for small molecules.
        - prepared_ligand (object, optional): The prepared ligand in OpenMM format. Default is None.
        """
        self.config_parser = config_parser
        self.forcefield_generator = ForcefieldGenerator()
        self.forcefield_selected = forcefield_selected
        self.water_selected = water_selected
        self.water_forcefield = water_forcefield
        self.smallMoleculeForceField = smallMoleculeForceField or self.config_parser.smallMoleculeForceField
        self.prepared_ligand = prepared_ligand
        self.add_membrane = self.config_parser.add_membrane

    def check_and_generate_forcefield(self):
        """
        Check if add_membrane is True and generate the transitional forcefield if so.

        Returns:
        - transitional_forcefield (object) or None
        """

        if self.add_membrane == "True":
            print("Generating transitional forcefield with membrane...")
            transitional_forcefield = self.forcefield_generator.generate_transitional_forcefield(
                protein_ff=self.forcefield_selected,
                solvent_ff=self.water_forcefield,
                add_membrane=self.add_membrane,
                smallMoleculeForceField=self.smallMoleculeForceField,
                rdkit_mol=self.prepared_ligand
            )
            print("Transitional forcefield generation complete.")
            return transitional_forcefield
        else:
            print("add_membrane is False; no transitional forcefield generated.")
            return None
        
    def create_forcefield(self):
        """
        Create the forcefield based on configuration.

        Returns:
        - forcefield (object): The generated or selected forcefield.
        """
        ligand = self.config_parser.ligand

        if ligand is not None:
            print("Generating forcefield with ligand...")
            forcefield = self.forcefield_generator.generate_forcefield(
                protein_ff=self.forcefield_selected,
                solvent_ff=self.water_forcefield,
                add_membrane=self.add_membrane,
                smallMoleculeForceField=self.smallMoleculeForceField,
                rdkit_mol=self.prepared_ligand
            )
            print("Forcefield generation with ligand complete.")
        else:
            if self.water_model is not None:
                print("Generating forcefield without ligand but with water model...")
                forcefield = self.forcefield_generator.generate_forcefield(
                    protein_ff=self.forcefield_selected,
                    solvent_ff=self.water_forcefield,
                    add_membrane=self.add_membrane,
                    smallMoleculeForceField=self.smallMoleculeForceField,
                    rdkit_mol=None
                )
                print("Forcefield generation without ligand but with water model complete.")
            else:
                print("Using default forcefield.")
                forcefield = app.ForceField(self.forcefield_selected)

        return forcefield


class ForcefieldSelector:
    """
    A class to handle the selection of forcefields and water models for simulations.

    Attributes:
        forcefield_dict (dict): Dictionary mapping forcefield names to their XML files.
        old_amber (set): Set of old AMBER forcefield XML file names.
        water_forcefield_mapping (dict): Dictionary mapping water model names to their XML files.
        water_model_mapping (dict): Dictionary mapping water model names to short model names.
        water_forcefields (dict): Dictionary mapping AMBER14 and CHARMM forcefield XML files to their respective water models.
        charmm_water_mapping (dict): Dictionary mapping CHARMM water models to their short model names.
    """

    def __init__(self):
        """
        Initializes the ForcefieldSelector class with predefined forcefield and water model mappings.
        """
        self.forcefield_dict = {
            "AMBER14": "amber14-all.xml",
            "AMBER99SB": "amber99sb.xml",
            "AMBER99SB-ILDN": "amber99sbildn.xml",
            "AMBER03": "amber03.xml",
            "AMBER10": "amber10.xml",
            "CHARMM36": "charmm36.xml",
        }
        self.old_amber = {
            "amber99sb.xml",
            "amber99sbildn.xml",
            "amber03.xml",
            "amber10.xml",
        }
        self.water_forcefield_mapping = {
            "TIP3P": "tip3p.xml",
            "TIP3P-FB": "tip3pfb.xml",
            "SPC/E": "spce.xml",
            "TIP4P-Ew": "tip4pew.xml",
            "TIP4P-FB": "tip4pfb.xml",
            "TIP5P": "tip5p.xml",
        }
        self.water_model_mapping = {
            "TIP3P": "tip3p",
            "TIP3P-FB": "tip3pfb",
            "SPC/E": "spce",
            "TIP4P-Ew": "tip4pew",
            "TIP4P-FB": "tip4pfb",
        }
        self.water_forcefields = {
            "amber14-all.xml": {
                "TIP3P": "amber14/tip3p.xml",
                "TIP3P-FB": "amber14/tip3pfb.xml",
                "SPC/E": "amber14/spce.xml",
                "TIP4P-Ew": "amber14/tip4pew.xml",
                "TIP4P-FB": "amber14/tip4pfb.xml",
            },
            "charmm36.xml": {
                "CHARMM default": "charmm36/water.xml",
                "TIP3P-PME-B": "charmm36/tip3p-pme-b.xml",
                "TIP3P-PME-F": "charmm36/tip3p-pme-f.xml",
                "SPC/E": "charmm36/spce.xml",
                "TIP4P-Ew": "charmm36/tip4pew.xml",
                "TIP4P-2005": "charmm36/tip4p2005.xml",
                "TIP5P": "charmm36/tip5p.xml",
                "TIP5P-Ew": "charmm36/tip5pew.xml",
            },
        }
        self.charmm_water_mapping = {
            "CHARMM default": "charmm",
            "TIP3P-PME-B": "charmm",
            "TIP3P-PME-F": "charmm",
            "SPC/E": "charmm",
            "TIP4P-Ew": "tip4pew",
            "TIP4P-2005": "tip4pew",
            "TIP5P": "tip5p",
            "TIP5P-Ew": "tip5p",
        }

    def ff_selection(self, ff):
        """
        Selects the required XML forcefield file.

        Args:
            ff (str): Input forcefield.

        Returns:
            str: Selected XML forcefield file.
        """
        return self.forcefield_dict.get(ff, None)

    def water_forcefield_selection(self, water, forcefield_selection):
        """
        Selects the XML filename for water force field parameters based on the chosen force field and water model.

        Args:
            water (str): The chosen water model.
            forcefield_selection (str): The selected force field.

        Returns:
            str: The XML filename of the water forcefield.
        """
        if forcefield_selection in self.old_amber:
            water_model = self.water_forcefield_mapping.get(water, None)
        else:
            water_model = self.water_forcefields.get(forcefield_selection, {}).get(
                water, None
            )
        return water_model

    def water_model_selection(self, water, forcefield_selection):
        """
        Selects the required water model forcefield XML file according to water selection and previous force field selection.

        Args:
            water (str): Water model input.
            forcefield_selection (str): Input of selected forcefield XML file.

        Returns:
            str: Water model forcefield XML file.
        """
        if forcefield_selection in self.old_amber:
            water_model = self.water_model_mapping.get(water)
        elif forcefield_selection == "amber14-all.xml":
            if water == "TIP5P":
                return None  # 'TIP5P' is not available in 'amber14-all.xml'
            water_model = self.water_model_mapping.get(water)
        elif forcefield_selection == "charmm36.xml":
            water_model = self.charmm_water_mapping.get(water)
        else:
            return None
        return water_model


class ForcefieldGenerator:
    """
    A class to generate OpenMM forcefields and register small molecules.

    Attributes:
        old_amber (set): Set of old AMBER forcefield XML file names.
    """

    def __init__(self):
        """
        Initializes the ForcefieldGenerator class with predefined old AMBER forcefield XML file names.
        """
        self.old_amber = {
            "amber99sb.xml",
            "amber99sbildn.xml",
            "amber03.xml",
            "amber10.xml",
        }

    def generate_forcefield(
        self,
        protein_ff,
        solvent_ff,
        add_membrane,
        smallMoleculeForceField,
        rdkit_mol=None,
    ):
        """
        Generate an OpenMM Forcefield object and register a small molecule.

        Args:
            protein_ff (str): Input of selected forcefield XML File.
            solvent_ff (str): Input of selected water model forcefield XML File.
            add_membrane (bool): Selection if the system should be built with a membrane.
            rdkit_mol (rdkit.Chem.rdchem.Mol): Small molecule to register in the force field.

        Returns:
            simtk.openmm.app.Forcefield: Forcefield with a registered small molecule.
        """
        if add_membrane:
            if protein_ff in self.old_amber:
                forcefield = app.ForceField(
                    protein_ff, solvent_ff, "amber14/lipid17.xml"
                )
            else:
                forcefield = app.ForceField(protein_ff, solvent_ff)
        else:
            forcefield = app.ForceField(protein_ff, solvent_ff)

        if rdkit_mol is not None:
            if smallMoleculeForceField == "gaff":
                gaff = GAFFTemplateGenerator(
                    molecules=Molecule.from_rdkit(
                        rdkit_mol, allow_undefined_stereo=True
                    ),
                    forcefield="gaff-2.11",
                )
                forcefield.registerTemplateGenerator(gaff.generator)
            elif smallMoleculeForceField == "smirnoff":
                smirnoff = SMIRNOFFTemplateGenerator(
                    molecules=Molecule.from_rdkit(
                        rdkit_mol, allow_undefined_stereo=True
                    ),
                    forcefield="openff-2.2.0",
                )
                forcefield.registerTemplateGenerator(smirnoff.generator)

        return forcefield

    def generate_transitional_forcefield(
        self,
        protein_ff,
        solvent_ff,
        add_membrane,
        smallMoleculeForceField,
        rdkit_mol=None,
    ):
        """
        Generate an OpenMM transitional forcefield object with TIP3P water model for membrane building and register a small molecule.

        Args:
            protein_ff (str): Name of the force field in XML format.
            solvent_ff (str): Name of the water model force field in XML format.
            add_membrane (bool): Selection if the system should be built with a membrane.
            rdkit_mol (rdkit.Chem.rdchem.Mol): Small molecule to register in the force field.

        Returns:
            simtk.openmm.app.Forcefield: A transitional forcefield with TIP3P water and a registered small molecule.
        """
        if add_membrane:
            if protein_ff in self.old_amber:
                transitional_forcefield = app.ForceField(
                    protein_ff, "tip3p.xml", "amber14/lipid17.xml"
                )
            else:
                transitional_forcefield = app.ForceField(
                    protein_ff, "amber14/tip3p.xml"
                )
        else:
            transitional_forcefield = app.ForceField(protein_ff, solvent_ff)

        if rdkit_mol is not None:
            if smallMoleculeForceField == "gaff":
                gaff = GAFFTemplateGenerator(
                    molecules=Molecule.from_rdkit(
                        rdkit_mol, allow_undefined_stereo=True
                    ),
                    forcefield="gaff-2.11",
                )
                transitional_forcefield.registerTemplateGenerator(gaff.generator)
            elif smallMoleculeForceField == "smirnoff":
                smirnoff = SMIRNOFFTemplateGenerator(
                    molecules=Molecule.from_rdkit(
                        rdkit_mol, allow_undefined_stereo=True
                    ),
                    forcefield="openff-2.2.0",
                )
                transitional_forcefield.registerTemplateGenerator(smirnoff.generator)

        return transitional_forcefield


    def prepare_forcefield_and_water(self):
        """
        Prepare the forcefield and water model for MD simulation.

        Returns:
            forcefield_selected: The selected forcefield.
            water_selected: The selected water forcefield.
            model_water: The selected water model.
        """
        forcefield_selected = self.forcefield_selector.ff_selection(self.ff)
        water_selected = self.forcefield_selector.water_forcefield_selection(
            water=self.water, forcefield_selection=forcefield_selected
        )
        model_water = self.forcefield_selector.water_model_selection(
            water=self.water, forcefield_selection=forcefield_selected
        )
        print("Forcefield and Water Model Selected")
        return forcefield_selected, water_selected, model_water

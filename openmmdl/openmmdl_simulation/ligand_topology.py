from openmmdl.openmmdl_simulation.forcefield_water import ForcefieldSelector, ForcefieldGenerator, ForcefieldPreparation
from openmmdl.openmmdl_simulation.protein_ligand_prep import MoleculePreparation
import numpy as np
from simtk.openmm import unit


class SimulationRunner:
    def __init__(self, config_parser):
        """
        Initialize the SimulationRunner with the given forcefield, water model, and optional ligand parameters.

        Parameters:
        - forcefield_name (str): The name of the forcefield to use.
        - water_model (str): The water model to use.
        - ligand_file (str, optional): Path to the ligand file (e.g., SDF). Default is None.
        - minimization (bool): Whether to perform minimization on the ligand. Default is False.
        """
        self.config_parser = config_parser
        self.forcefield_name = config_parser.forcefield
        self.water_model = config_parser.water
        self.ligand_file = config_parser.ligand
        self.minimization = config_parser.minimization

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

    def prepare_ligand(self):
        """
        Prepare the ligand if a ligand file is provided.

        Returns:
        - tuple: (prepared_ligand, omm_ligand) or (None, None)
        """
        if self.ligand_file is not None:
            print("Preparing ligand...")
            molecule_preparation = MoleculePreparation(self.ligand_file, self.minimization)
            prepared_ligand = molecule_preparation.prepare_ligand()
            omm_ligand = molecule_preparation.rdkit_to_openmm(prepared_ligand, name="UNK")
            print("Ligand preparation complete.")
            return prepared_ligand, omm_ligand
        else:
            print("No ligand file provided; skipping ligand preparation.")
            return None, None

class TopologyAndPositionsSelector:
    def __init__(self, config_parser, modeller=None, prmtop=None, inpcrd=None):
        """
        Initialize the TopologyAndPositionsSelector with the necessary parameters.

        Parameters:
        - config_parser (ConfigParser): The configuration parser instance.
        - modeller (app.Modeller, optional): OpenMM Modeller object.
        - prmtop (app.AmberPrmtopFile, optional): Amber prmtop file.
        - inpcrd (app.AmberInpcrdFile, optional): Amber inpcrd file.
        """
        self.config_parser = config_parser
        self.modeller = modeller
        self.prmtop = prmtop
        self.inpcrd = inpcrd

    def select_topology_and_positions(self):
        """
        Select the appropriate topology and positions based on the input file type.

        Returns:
        - topology (app.Topology): The selected topology.
        - positions (unit.Quantity): The selected positions.
        """
        input_file_type = self.config_parser.get_variable('input_file_type', None)
        
        if input_file_type == "pdb":
            if self.modeller is None:
                raise ValueError("Modeller object must be provided for PDB input filetype.")
            topology = self.modeller.topology
            positions = self.modeller.positions
            positions = np.array([[pos[i].value_in_unit(unit.nanometer) for i in range(3)] for pos in self.modeller.positions])
        elif input_file_type == "amber":
            if self.prmtop is None or self.inpcrd is None:
                raise ValueError("Prmtop and Inpcrd files must be provided for Amber input filetype.")
            topology = self.prmtop.topology
            positions = self.inpcrd.positions
        else:
            raise ValueError("Unsupported input filetype. Supported types are 'pdb' and 'amber'.")

        print("Topology and positions obtained")
        return topology, positions

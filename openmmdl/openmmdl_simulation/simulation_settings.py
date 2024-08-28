import simtk.openmm.app as app
from simtk.openmm import unit, Platform, MonteCarloBarostat, LangevinMiddleIntegrator, XmlSerializer
from openmmdl.openmmdl_simulation.minimization_equilibration import MinimizerEquilibrator
import ast
import sys
import simtk.openmm.app as app
from simtk.openmm.app import PDBFile, Modeller, PDBReporter, StateDataReporter, DCDReporter, CheckpointReporter, AmberPrmtopFile, AmberInpcrdFile


class SystemCreator:
    def __init__(self, config_parser, forcefield=None, prmtop=None):
        """
        Initialize the SystemCreator with the necessary parameters.

        Parameters:
        - config_parser (ConfigParser): The configuration parser instance.
        - forcefield (app.ForceField, optional): OpenMM ForceField object.
        - prmtop (app.AmberPrmtopFile, optional): Amber prmtop file.
        """
        self.config_parser = config_parser
        self.forcefield = forcefield
        self.prmtop = prmtop
        self.nonbondedMethod_app, self.nonbondedMethod_name = config_parser.nonbondedMethod.split(".")
        self.nonbondedMethod = getattr(app, self.nonbondedMethod_name)

        self.platformProperties = ast.literal_eval(config_parser.platformProperties)
        # do it like this self.pme_check2 = getattr(app, string)

    def create_system(self, nonbondedMethod, topology=None, positions=None, nonbondedCutoff=None):
        """
        Create the system based on the input file type.

        Parameters:
        - topology (app.Topology): The topology object.
        - nonbondedMethod (str): The nonbonded method for the system.
        - nonbondedCutoff (Quantity, optional): The nonbonded cutoff distance.

        Returns:
        - system (openmm.System): The created system.
        """
        input_filetype = self.config_parser.get_variable('input_file_type', None)

        # Define a mapping to convert constraint strings to actual objects
        constraint_mapping = {
            "app.HBonds": app.HBonds,
            "app.AllBonds": app.AllBonds,
            "None": None
        }

        # Fetch and convert the constraint directly when getting the variable
        constraints_str = self.config_parser.get_variable('constraints', None)
        constraints = constraint_mapping.get(constraints_str, None)  # Apply mapping

        # Define parameters for system creation
        system_params = {
            'nonbondedMethod': self.config_parser.nonbondedMethod,
            'rigidWater': self.config_parser.get_variable('rigidWater', False),
            'hydrogenMass': self.config_parser.hydrogenMass
        }
        hydrogenMass_value = float(system_params['hydrogenMass'].split('*')[0])


        if nonbondedMethod != "NoCutoff":
            system_params['nonbondedCutoff'] = nonbondedCutoff
            nonbondedCutoff_value = float(system_params['nonbondedCutoff'].split('*')[0])
        
        if nonbondedMethod == "app.PME":
            system_params['ewaldErrorTolerance'] = float(self.config_parser.ewaldErrorTolerance)
        
        if input_filetype == "pdb":
            if self.forcefield is None:
                raise ValueError("Forcefield object must be provided for PDB input filetype.")
            if nonbondedMethod == "app.PME":
                ### Cover all cases
                system = self.forcefield.createSystem(
                    topology,
                    nonbondedMethod = self.nonbondedMethod,
                    nonbondedCutoff = nonbondedCutoff_value*unit.nanometers,
                    constraints = constraints,
                    rigidWater = system_params['rigidWater'],
                    ewaldErrorTolerance = system_params['ewaldErrorTolerance'],
                    hydrogenMass = hydrogenMass_value*unit.amu,
                )
        elif input_filetype == "amber":
            ### Cover all cases
            if self.prmtop is None:
                raise ValueError("Prmtop file must be provided for Amber input filetype.")
            if nonbondedMethod == "app.PME":
                system = self.prmtop.createSystem(
                    nonbondedMethod = self.nonbondedMethod,
                    nonbondedCutoff = nonbondedCutoff_value*unit.nanometers,
                    constraints = constraints,
                    rigidWater = system_params['rigidWater'],
                    ewaldErrorTolerance = system_params['ewaldErrorTolerance'],
                    hydrogenMass = hydrogenMass_value*unit.amu,
                )
        else:
            raise ValueError("Unsupported input filetype. Supported types are 'pdb' and 'amber'.")

        # Minimization and Equilibration
        minimizer_equilibrator = MinimizerEquilibrator(topology, positions, system, self.config_parser)
        system, positions = minimizer_equilibrator.minimize_and_equilibrate()

        # Add Values
        temperature_value = float(self.config_parser.temperature.split('*')[0])
        friction_value = float(self.config_parser.friction.split('/')[0])
        dt_value = float(self.config_parser.dt.split('*')[0])
        
        # Add MonteCarloBarostat if ensemble is 'npt'
        if self.config_parser.ensemble == "npt":
            system.addForce(MonteCarloBarostat(self.config_parser.pressure, temperature_value*unit.kelvin, self.config_parser.barostatInterval))

        # Create integrator
 
        integrator = LangevinMiddleIntegrator(temperature_value*unit.kelvin, friction_value/unit.picosecond, dt_value*unit.picoseconds)
        if self.config_parser.constraints != None:
            integrator.setConstraintTolerance(float(self.config_parser.constraintTolerance))

        # Create simulation
        #Cover Platforms
        simulation = app.Simulation(
            topology, system, integrator, Platform.getPlatformByName(self.config_parser.platform), self.platformProperties
        )
        simulation.context.setPositions(positions)
        
        return system, integrator, simulation

    def run_simulation(self, simulation):
        """
        Run the simulation, adding necessary reporters and executing the simulation steps.

        Parameters:
        - simulation (app.Simulation): The simulation object created by create_simulation.
        - protein (str): The name of the protein, used for naming the output PDB file.
        - steps (int): The number of simulation steps to run.
        - pdbInterval (int): Interval at which PDB frames should be written.
        - dcdReporter (DCDReporter): The DCD reporter object.
        - dataReporter (StateDataReporter): The data reporter object.
        - checkpointReporter (CheckpointReporter): The primary checkpoint reporter object.
        - checkpointReporter10 (CheckpointReporter): The secondary checkpoint reporter for 10 steps interval.
        - checkpointReporter100 (CheckpointReporter): The tertiary checkpoint reporter for 100 steps interval.
        """
        print('Simulating...')
        
        # Append reporters to the simulation
        simulation.reporters.append(PDBReporter(f'output_{self.config_parser.protein}.pdb', int(self.config_parser.pdbInterval)))
        simulation.reporters.append(DCDReporter(self.config_parser.dcd_name, int(self.config_parser.dcdInterval)))
        simulation.reporters.append(CheckpointReporter(self.config_parser.checkpoint_name, int(self.config_parser.checkpointInterval)))
        simulation.reporters.append(StateDataReporter(sys.stdout, 1000, step=True, potentialEnergy=True, temperature=True))

        simulation.currentStep = 0
        print("run")
        print(int(self.config_parser.steps))
        simulation.step(int(self.config_parser.steps))


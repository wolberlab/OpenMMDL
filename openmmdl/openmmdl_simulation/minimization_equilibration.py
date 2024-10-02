from simtk.openmm.app import PDBFile
import numpy as np
from openmm import unit, OpenMMException
import openmm.unit as unit
import logging
import openmm
import copy
import re
import time
import mdtraj as md
from simtk.openmm import unit, Platform

class MinimizerEquilibrator:
    def __init__(self, topology, positions, system, config_parser):
        """
        Initialize the MinimizerEquilibrator with the necessary parameters.

        Parameters:
        - simulation (openmm.app.Simulation): The simulation object.
        - modeller (openmm.app.Modeller): The modeller object with the topology and positions.
        - protein_name (str): The name of the protein (used for output file naming).
        - temperature (Quantity): The temperature for equilibration.
        - equilibration_steps (int): The number of steps for equilibration.
        - file_type (str): The input file type, either 'pdb' or 'amber'.
        - prmtop (AmberPrmtopFile, optional): The AMBER prmtop object, required if file_type is 'amber'.
        - inpcrd (AmberInpcrdFile, optional): The AMBER inpcrd object, required if file_type is 'amber'.
        """
        self.topology = topology
        self.positions = positions
        self.system= system
        self.config_parser = config_parser
        
    def minimize_and_equilibrate(self):
        """
        Run gentle equilibration. Serialize results in an external PDB or CIF file.
        Utility functions useful for equilibration.
        Adapted from zhang-ivy(GitHub) and openmmtools(GitHub)
        Add gentle equilibration protocol as utility function (#669)

        Parameters
        ----------
        topology : openmm.app.Topology
            topology
        positions : np.array in unit.nanometer
            positions
        system : openmm.System
            system
        stages : list of dicts
            each dict corresponds to a stage of equilibration and contains the equilibration parameters for that stage

            equilibration parameters:
                EOM : str
                    'minimize' or 'MD' or 'MD_interpolate' (the last one will allow interpolation between 'temperature' and 'temperature_end')
                n_steps : int
                    number of steps of MD
                temperature : openmm.unit.kelvin
                    temperature (kelvin)
                temperature_end : openmm.unit.kelvin, optional
                    the temperature (kelvin) at which to finish interpolation, if 'EOM' is 'MD_interpolate'
                ensemble : str or None
                    'NPT' or 'NVT'
                restraint_selection : str or None
                    to be used by mdtraj to select atoms for which to apply restraints
                force_constant : openmm.unit.kilocalories_per_mole/openmm.unit.angstrom**2
                    force constant (kcal/molA^2)
                collision_rate : 1/openmm.unit.picoseconds
                    collision rate (1/picoseconds)
                timestep : openmm.unit.femtoseconds
                    timestep (femtoseconds)
        filename : str
            path to save the equilibrated structure
        platform_name : str, default 'CUDA'
            name of platform to be used by OpenMM. If not specified, OpenMM will select the fastest available platform
        save_box_vectors : bool
            Whether to save the box vectors in a box_vectors.npy file in the working directory, after execution.
            Defaults to True.

        """
        _logger = logging.getLogger(__name__)
        equilibration_manager = EquilibrationManager(self.config_parser)
        stages = equilibration_manager.get_equilibration_stages()
        filename=f"Equilibration_{self.config_parser.protein}"
        platform_name = self.config_parser.platform
        platform = Platform.getPlatformByName(platform_name)
        save_box_vectors=bool(self.config_parser.save_box_vectors_equilibration)
        positions = self.positions
        if stages:
            print("Minimizing/Equilibrating...")
            for stage_index, parameters in enumerate(stages):

                initial_time = time.time()
                print(f"Executing stage {stage_index + 1} of equilibration...")

                # Make a copy of the system
                system_copy = copy.deepcopy(self.system)

                # Add restraint
                if parameters['restraint_selection'] is not None:
                    traj = md.Trajectory(positions, md.Topology.from_openmm(self.topology))
                    selection_indices = traj.topology.select(parameters['restraint_selection'])

                    custom_cv_force = openmm.CustomCVForce('(K_RMSD/2)*(RMSD)^2')
                    custom_cv_force.addGlobalParameter('K_RMSD', parameters['force_constant'] * 2)
                    rmsd_force = openmm.RMSDForce(positions, selection_indices)
                    custom_cv_force.addCollectiveVariable('RMSD', rmsd_force)
                    system_copy.addForce(custom_cv_force)

                # Set barostat update interval to 0 (for NVT)
                if parameters['ensemble'] == 'NVT':
                    force_dict = {force.__class__.__name__: index for index, force in enumerate(system_copy.getForces())}
                    try:
                        system_copy.getForce(force_dict['MonteCarloBarostat']).setFrequency(0)  # This requires openmm 8
                    except KeyError:
                        # No MonteCarloBarostat found
                        _logger.debug("No MonteCarloBarostat found in forces. Continuing.")
                    except OpenMMException:
                        # Must be Openmm<8
                        system_copy.removeForce(force_dict['MonteCarloBarostat'])

                elif parameters['ensemble'] == 'NPT' or parameters['ensemble'] is None:
                    pass

                else:
                    raise ValueError("Invalid parameter supplied for 'ensemble'")

                # Set up integrator
                temperature = parameters['temperature']
                collision_rate = parameters['collision_rate']
                timestep = parameters['timestep']

                integrator = openmm.LangevinMiddleIntegrator(temperature, collision_rate, timestep)

                # Set up context
                print(platform_name)
                if platform_name in ['CUDA', 'OpenCL']:
                    platform.setPropertyDefaultValue('Precision', 'mixed')
                if platform_name in ['CUDA']:
                    platform.setPropertyDefaultValue('DeterministicForces', 'true')

                context = openmm.Context(system_copy, integrator, platform)
                context.setPeriodicBoxVectors(*system_copy.getDefaultPeriodicBoxVectors())
                context.setPositions(positions)
                context.setVelocitiesToTemperature(temperature)

                # Run minimization or MD
                n_steps = parameters['n_steps']
                n_steps_per_iteration = 100

                if parameters['EOM'] == 'minimize':
                    openmm.LocalEnergyMinimizer.minimize(context, maxIterations=n_steps)

                elif parameters['EOM'] == 'MD':
                    for _ in range(int(n_steps / n_steps_per_iteration)):
                        integrator.step(n_steps_per_iteration)

                elif parameters['EOM'] == 'MD_interpolate':
                    temperature_end = parameters['temperature_end']
                    temperature_unit = unit.kelvin
                    temperatures = np.linspace(temperature / temperature_unit, temperature_end / temperature_unit,
                                            int(n_steps / n_steps_per_iteration)) * temperature_unit
                    for temperature in temperatures:
                        integrator.setTemperature(temperature)
                        integrator.step(n_steps_per_iteration)

                else:
                    raise ValueError("Invalid parameter supplied for 'EOM'")

                # Retrieve positions after this stage of equil
                state = context.getState(getPositions=True)
                positions = state.getPositions(asNumpy=True)

                # Update default box vectors for next iteration
                box_vectors = state.getPeriodicBoxVectors()
                self.system.setDefaultPeriodicBoxVectors(*box_vectors)

                # Delete context and integrator
                del context, integrator, system_copy

                elapsed_time = time.time() - initial_time
                print(f"\tStage {stage_index + 1} took {elapsed_time} seconds")

            # Save the final equilibrated positions
            
            if filename.endswith('pdb'):
                openmm.app.PDBFile.writeFile(self.topology, positions, open(filename, "w"), keepIds=True)
            elif filename.endswith('cif'):
                openmm.app.PDBxFile.writeFile(self.topology, positions, open(filename, "w"), keepIds=True)

            # Save the box vectors
            if save_box_vectors:
                with open(filename[:-4] + '_box_vectors.npy', 'wb') as f:
                    np.save(f, box_vectors)
            print("Equilibration/Minimization complete.")
            return self.system, positions
        else:
            print("No minimization or equilibration stages provided; skipping.")
            return self.system, positions


class EquilibrationManager:
    def __init__(self, config_parser):
        self.config_parser = config_parser
        self.equilibration_stages = []

        # Initialize equilibration stages based on equilibration type
        self.setup_equilibration_stages()

    def setup_equilibration_stages(self):
        equilibration_type = self.config_parser.preparation_type

        if equilibration_type == "equilibration":
            self.equilibration_stages = [
                {'EOM': 'minimize', 'n_steps': 10000, 'temperature': 300*unit.kelvin, 'ensemble': None, 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD_interpolate', 'n_steps': 100000, 'temperature': 100*unit.kelvin, 'temperature_end': 300*unit.kelvin, 'ensemble': 'NVT', 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 10/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 10/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 250000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and not type H', 'force_constant': 10*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'minimize', 'n_steps': 10000, 'temperature': 300*unit.kelvin, 'ensemble': None, 'restraint_selection': 'protein and backbone', 'force_constant': 10*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and backbone', 'force_constant': 10*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and backbone', 'force_constant': 1*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 100000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': 'protein and backbone', 'force_constant': 0.1*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
                {'EOM': 'MD', 'n_steps': 2500000, 'temperature': 300*unit.kelvin, 'ensemble': 'NPT', 'restraint_selection': None, 'force_constant': 0*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 2*unit.femtoseconds},
            ]
        
        elif equilibration_type == "minimization":
            self.equilibration_stages = [
                {'EOM': 'minimize', 'n_steps': 10000, 'temperature': 300*unit.kelvin, 'ensemble': None, 'restraint_selection': 'protein and not type H', 'force_constant': 100*unit.kilocalories_per_mole/unit.angstrom**2, 'collision_rate': 2/unit.picoseconds, 'timestep': 1*unit.femtoseconds},
            ]
        
        elif equilibration_type == "None":
            self.equilibration_stages = None

        else:
            raise ValueError(f"Unknown equilibration_type: {equilibration_type}")

    def get_equilibration_stages(self):
        return self.equilibration_stages
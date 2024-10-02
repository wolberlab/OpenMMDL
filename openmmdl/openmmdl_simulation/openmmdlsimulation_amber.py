"""
mmdl_simulation.py
Perform Simulations of Protein-ligand complexes with OpenMM
"""

import argparse
import sys
import os
import shutil
import argparse

from simtk.openmm.app import PDBFile, Modeller, PDBReporter, StateDataReporter, DCDReporter, CheckpointReporter, AmberPrmtopFile, AmberInpcrdFile
from openmmdl.openmmdl_simulation.parser import ConfigParser
from openmmdl.openmmdl_simulation.runner import SimulationRunner, TopologyAndPositionsSelector
from openmmdl.openmmdl_simulation.minimization_equilibration import MinimizerEquilibrator
from openmmdl.openmmdl_simulation.simulation_settings import SystemCreator
from openmmdl.openmmdl_simulation.forcefield_water import ForcefieldSelector, ForcefieldGenerator, ForcefieldPreparation, ForcefieldConfigurator
from openmmdl.openmmdl_simulation.protein_ligand_prep import MoleculePreparation, WaterMembraneComplexBuilder, PDBFormatter, ComplexBuilder, ModellerCreator, EnvironmentBuilder


parser = argparse.ArgumentParser()


logo = "\n".join(
    [
        "     ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      ",
        "   .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      ",
        "  / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      ",
        " ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    ",
        " |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    ",
        " : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    ",
        "  \ `_/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  ",
        "   '. \_/``'.'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ ",
        "     '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` ",
        "              Prepare and Perform OpenMM Protein-Ligand MD Simulations                                 ",
        "                                     Alpha Version                                                     ",
    ]
)


def main():
    parser = argparse.ArgumentParser(
        prog="openmmdl_simulation",
        description=logo,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-f",
        dest="folder",
        type=str,
        help="Folder Name for MD Simulation",
        required=True,
    )
    parser.add_argument(
        "-s",
        dest="config",
        type=str,
        help="MD Simulation config file",
        required=True,
    )
    parser.add_argument(
        "-t", dest="topology", help="Protein Topology PDB/Amber File", required=True
    )
    parser.add_argument("-l", dest="ligand", help="SDF File of Ligand", default=None)
    parser.add_argument(
        "-c", dest="coordinate", help="Amber coordinates file", default=None
    )
    input_formats = [".conf", ".pdb", ".sdf", ".mol", ".prmtop", ".inpcrd"]

    args = parser.parse_args()

    # Create or replace the folder
    if not os.path.exists(args.folder):
        os.mkdir(args.folder)
    else:
        shutil.rmtree(args.folder)
        os.mkdir(args.folder)

    script_dir = os.path.abspath(os.path.dirname(__file__))

    # Copy config file
    if input_formats[0] not in args.config:
        print("Wrong Format for config file. Expected .conf")
        exit(1)

    if not os.path.exists(args.config):
        print("Wrong python config path, try the absolute path")
        exit(1)

    shutil.copy(args.config, args.folder)

    # Copy topology file
    if input_formats[1] not in args.topology and input_formats[4] not in args.topology:
        print("Wrong Format for topology file. Expected .pdb or .prmtop")
        exit(1)

    if not os.path.exists(args.topology):
        print("Wrong topology file path, try the absolute path")
        exit(1)

    shutil.copy(args.topology, args.folder)

    # Copy ligand file if provided
    if args.ligand:
        if input_formats[2] not in args.ligand and input_formats[3] not in args.ligand:
            print("Wrong Format for ligand file. Expected .sdf or .mol")
            exit(1)

        if not os.path.exists(args.ligand):
            print("Wrong ligand file path, try the absolute path")
            exit(1)

        shutil.copy(args.ligand, args.folder)

    # Copy coordinate file if provided
    if args.coordinate:
        if input_formats[5] not in args.coordinate:
            print("Wrong Format for coordinate file. Expected .inpcrd")
            exit(1)

        if not os.path.exists(args.coordinate):
            print("Wrong coordinates file path, try the absolute path")
            exit(1)

        shutil.copy(args.coordinate, args.folder)

    os.chdir(args.folder)
    print("hello")
    # Obtain Data from Configuration File
    config_file_path = args.config
    config_parser = ConfigParser(config_file_path)
    config_parser.print_variables()

    print(args.topology)

    prmtop = AmberPrmtopFile(config_parser.prmtop_file)
    inpcrd = AmberInpcrdFile(config_parser.inpcrd_file)

    print(prmtop)
    print(inpcrd)


    topology = prmtop.topology
    positions = inpcrd.positions


    # Create the system creator instance
    system_creator = SystemCreator(config_parser, prmtop=prmtop)

    # Create the system
    system, integrator, simulation = system_creator.create_system(config_parser.nonbondedMethod, topology, positions, config_parser.nonbondedCutoff)


    # For PDB input
    temperature_value = float(config_parser.temperature.split('*')[0])
    minimizer_equilibrator = MinimizerEquilibrator(simulation, modeller, config_parser.protein, temperature_value*unit.kelvin, int(config_parser.equilibrationSteps), file_type=config_parser.input_filetype)
    minimizer_equilibrator.minimize_and_equilibrate()

if __name__ == "__main__":
    main()

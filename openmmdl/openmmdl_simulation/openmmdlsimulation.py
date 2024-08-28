"""
mmdl_simulation.py
Perform Simulations of Protein-ligand complexes with OpenMM
"""

import argparse
import sys
import os
import shutil
import argparse

from openmmdl.openmmdl_simulation.parser import ConfigParser
from openmmdl.openmmdl_simulation.ligand_topology import SimulationRunner, TopologyAndPositionsSelector
from openmmdl.openmmdl_simulation.simulation_settings import SystemCreator
from openmmdl.openmmdl_simulation.post_md_conversions import MDPostProcessingHandler
from openmmdl.openmmdl_simulation.cleaning_procedures import PostMDProcessor, Cleanup
from openmmdl.openmmdl_simulation.minimization_equilibration import MinimizerEquilibrator
from openmmdl.openmmdl_simulation.forcefield_water import ForcefieldSelector, ForcefieldGenerator, ForcefieldPreparation, ForcefieldConfigurator
from openmmdl.openmmdl_simulation.protein_ligand_prep import MoleculePreparation, WaterMembraneComplexBuilder, PDBFormatter, ComplexBuilder, ModellerCreator, EnvironmentBuilder
import numpy as np
from openmm import unit, OpenMMException

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
    if not os.path.exists(args.folder):
        os.mkdir(args.folder)
    else:
        shutil.rmtree(args.folder)
        os.mkdir(args.folder)
    script_dir = os.path.abspath(os.path.dirname(__file__))
    if os.path.exists(args.folder):
        if input_formats[0] in args.config:
            if os.path.exists(args.config):
                shutil.copy(args.config, args.folder)
            else:
                print("Wrong python config path, try the absolute path")
        if input_formats[1] in args.topology:
            if os.path.exists(args.topology):
                shutil.copy(args.topology, args.folder)
            else:
                print("Wrong topology file path, try the absolute path")
        elif input_formats[4] in args.topology:
            if os.path.exists(args.topology):
                shutil.copy(args.topology, args.folder)
            else:
                print("Wrong topology file path, try the absolute path")
        else:
            print("Wrong Format, don't forget the .pdb/.prmtop of the file")
        if args.ligand != None:
            if input_formats[2] in args.ligand or input_formats[3] in args.ligand:
                if os.path.exists(args.ligand):
                    shutil.copy(args.ligand, args.folder)
                else:
                    print("Wrong ligand file path, try the absolute path")
            else:
                print("Wrong Format, don't forget the .sdf of the ligand file")
        if args.coordinate != None:
            if input_formats[5] in args.coordinate:
                if os.path.exists(args.coordinate):
                    shutil.copy(args.coordinate, args.folder)
                else:
                    print("Wrong coordinates file path, try the absolute path")
            else:
                print("Wrong Format, don't forget the .inpcrd of the coordinate file")
        os.chdir(args.folder)
        # Obtain Data from Configuration File
        config_file_path = args.config
        config_parser = ConfigParser(config_file_path)

        forcefield_name = config_parser.forcefield

        pdb_formatter = PDBFormatter(config_parser)
        protein_pdb = pdb_formatter.get_protein_pdb()

        runner = SimulationRunner(
        config_parser=config_parser
        )


        forcefield, water_forcefield, model_water = runner.setup_forcefield()
        prepared_ligand, omm_ligand = runner.prepare_ligand()

        forcefield_configurator = ForcefieldConfigurator(
            config_parser=config_parser,
            water_selected=water_forcefield,
            water_forcefield=water_forcefield,
            forcefield_selected = forcefield,
            smallMoleculeForceField=config_parser.smallMoleculeForceField,
            prepared_ligand=prepared_ligand
        )  

        transitional_forcefield = forcefield_configurator.check_and_generate_forcefield()

        forcefield = forcefield_configurator.create_forcefield()
        
        # Create the ComplexBuilder instance
        complex_builder = ComplexBuilder(
            config_parser=config_parser,
            protein_pdb=protein_pdb,
            omm_ligand=omm_ligand
        )

        # Build the complex
        complex_topology, complex_positions = complex_builder.build_complex()

        modeller_creator = ModellerCreator(
            config_parser=config_parser,
            protein_pdb=protein_pdb,
            complex_topology=complex_topology,  # Optional
            complex_positions=complex_positions  # Optional
        )

        modeller = modeller_creator.create_modeller()

        # Create the environment builder instance
        env_builder = EnvironmentBuilder(
            config_parser=config_parser,
            forcefield_name=forcefield_name,
            model_water=model_water,
            forcefield=forcefield,
            transitional_forcefield=transitional_forcefield,
            protein_pdb=protein_pdb,
            protein=config_parser.protein,
            modeller=modeller,
            ligand=config_parser.ligand,
        )

        env_builder.build_environment()

        modeller_complex = env_builder.modeller

        # Create the topology and positions selector
        selector = TopologyAndPositionsSelector(config_parser, modeller_complex, config_parser.prmtop, config_parser.inpcrd)

        topology, positions = selector.select_topology_and_positions()
        # Select the appropriate topology and positions

        # Create the system creator instance
        system_creator = SystemCreator(config_parser, forcefield=forcefield, prmtop=config_parser.prmtop)

        # Create the system
        system, integrator, simulation = system_creator.create_system( config_parser.nonbondedMethod, topology, positions, config_parser.nonbondedCutoff)

        system_creator.run_simulation(simulation)

        # Initialize and run the MDPostProcessingHandler
        md_postprocessing_handler = MDPostProcessingHandler(
            config_parser=config_parser
        )

        # Perform the MD post-processing
        md_postprocessing_handler.process()

        PostMDProcessor.post_md_file_movement(config_parser.protein, config_parser.prmtop, config_parser.inpcrd, config_parser.ligand)

if __name__ == "__main__":
    main()

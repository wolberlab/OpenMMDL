<img src="https://openmm.org/images/logo.svg" height="200">


# OpenMMDL Simulation


Interface to OpenMM for the easy setup of molecular dynamic simulations of protein-ligand complexes

## Install

#### Clone this repository

Open a new terminal and clone this repository
    cd ~
    git clone 

#### Install all required dependencies in a separate environment


OpenMMDL is written in python 3.10 and uses diverse packages, which can be easily installed on a separate environment using conda:

    cd openmmdl
    conda create -n openmmdl --file requirements.txt

#### Create alias for setup and simulation for your bash

    echo 'alias openmmdl-setup="python3 ~/openmmdl/openmmdl_setup/openmmdlsetup.py"' >> ~/.bashrc
    echo 'alias openmmdl-simulation="python3 ~/openmmdl/openmmdl_simulation/openmmdlsimulation.py"' >> ~/.bashrc

## Using OpenMMDL-Setup

Activate conda environment:

    conda activate openmmdl


Start the Setup by executing the command:

    openmmdl-setup

The setup interface is displayed through a web browser, but it is still a single-user desktop application, not a web application. It should automatically open a web browser displaying the user interface. If for any reason that does not happen, open a browser yourself and point it to the address displayed in the console window (usually http://127.0.0.1:5000).

Download the processed PDB File and python script, which serve as input for the simulation script.

## Using OpenMMDL-Simulation

Activate conda environment if not activated:

    conda activate openmmdl

Start the Simulation with the following Inputs:
Mandatory:
-f = Name of the folder, where the MD simulation files will be stored
-p = PDB file of your protein from the Setup including the path
-s = Python script from the Setup including the path
Optional:
-l = SDF File of the ligand, if the ligand was selected during the Setup. sdf file name should be consistent with the input in setup

Command line example without ligand

    openmmdl-simulation -f {path/to/folder_name} -p {path/to/protein.pdb} -s {path/to/script.py}

Command line example with ligand

    openmmdl-simulation -f {path/to/folder_name} -p {path/to/protein.pdb} -s {path/to/script.py} -l {path/to/ligand.sdf}

## Running OpenMMDL-Simulation test simulations

There are two Systems prepared for the testing of the simulation.

 1: A 10ns simulation of the 6b73 Protein-ligand complex with POPC membrane and TIP3P-FB water. To run the testing of 6b73 enter the following command line:

    openmmdl-simulation -f 6b73_testing_simulation -p ~/openmmdl/openmmdl_simulation/testing_sytems/6b73_membrane/6b73-moe-processed_openMMDL.pdb -s ~/openmmdl/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_simulation.py -l  ~/openmmdl/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_lig.sdf

 1: A 10ns simulation of the 5wyz Protein-ligand complex with TIP3P water. To run the testing of 5wyz enter the following command line:

    openmmdl-simulation -f 5wyz_testing_simulation -p ~/openmmdl/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz-moe-processed_openMMDL.pdb -s ~/openmmdl/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz_simulation.py -l  ~/openmmdl/openmmdl_simulation/testing_sytems/5wyz_solvent/5VF.sdf


## Copyright
Copyright (c) 2022, Valerij Talagayev & Yu Chen

#### Acknowledgements

The Script is based upon the OpenMM Toolkit [OpenMM](https://github.com/openmm)
Many thanks to all of the creators and contributors of the OpenMM Toolkit, especially a big thank you to [jchodera](https://github.com/jchodera), [peastman](https://github.com/peastman), [mikemhenry](https://github.com/mikemhenry) and the whole [choderalab](https://github.com/choderalab) 

The Simulation script is heavily inspired by the CADD T019 Talkatorial
(https://github.com/volkamerlab/teachopencadd/tree/master/teachopencadd/talktorials/T019_md_simulation)
Thanks to the members of the [Volkamer Lab](https://volkamerlab.org/),
especially [dominiquesydow](https://github.com/dominiquesydow/), [schallerdavid](https://github.com/schallerdavid) and [AndreaVolkamer](https://github.com/andreavolkamer).

The Setup script is a modified version of [openmm-setup](https://github.com/openmm/openmm-setup).
 
Project is based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.


[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<img src="https://github.com/pacificore/OpenMMDL/blob/main/openmmdl_setup/static/OpenMMDL_logo.svg" height="250">


# OpenMMDL @ Molecular Design Lab

Interface to OpenMM for easy setup of molecular dynamic simulations of
protein-ligand complexes

## Install

#### Clone this repository

Open a new terminal and clone this repository

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL

#### Install all required dependencies in a separate environment

OpenMMDL is written in python 3.10 and uses several packages, which can
be easily installed on a separate environment using conda (we recommend
using miniconda):

    cd OpenMMDL
    conda create -n openmmdl --file requirements.txt

After installation, activate the conda environment:

    conda activate openmmdl

#### Create alias for setup and simulation for your bash

    echo 'alias openmmdl-setup="python3 ~/OpenMMDL/openmmdl_setup/openmmdlsetup.py"' >> ~/.bashrc
    echo 'alias openmmdl-simulation="python3 ~/OpenMMDL/openmmdl_simulation/openmmdlsimulation.py"' >> ~/.bashrc
    echo 'alias openmmdl-analysis="python3 ~/OpenMMDL/openmmdl_analysis/openmmdlanalysis.py"' >> ~/.bashrc

## OpenMMDL-Setup

OpenMMDL-Setup will setup an MD simulation for OpenMM (i.e. create all
the needed files for running the simulation)

### Usage

Start the setup process by executing the command:

    openmmdl-setup

The setup interface is displayed through a web browser, but it is still
a single-user desktop application, not a web application. It should
automatically open a web browser displaying the user interface. If for
any reason that does not happen, open a browser yourself and point it to
the address displayed in the console window (e.g. http://127.0.0.1:5000).

Download the processed PDB File and python script, which serve as input
for the simulation script.

## OpenMMDL-Simulation

OpenMMDL-Simulation starts the MD simulation with the inputs acquired
from OpenMMDL-Setup

### Usage

Start the simulation with the following Inputs:

Mandatory:
-f = name of the folder, where the MD simulation files are stored
-t = Topology file of your protein from the setup including the path
-s = Python script from the setup including the path
Optional:
-l = SDF File of the ligand, if the ligand was selected during the
Setup. sdf file name should be consistent with the input in setup
-c = Cooridnates file of Amber


Command line example with ligand

    openmmdl-simulation -f {path/to/folder_name} -t {path/to/topology} -s {path/to/script} -l {path/to/ligand}

## OpenMMDL-Analysis

OpenMMDL-Analysis performs an analysis of the MD simulation obtained from OpenMMDL-Simulation.
It analyzes the protein-ligand complex interactions throughout the whole MD trajectory, delivering the list of
all possible interactions. In addition, it generates interaction fingerprints, provides the most occurring binding modes
and displays the transitions between the separate binding modes.

### Usage

Start the analysis with the following Inputs:

Mandatory:
-t = Topology file
-d = Trajectory file
-l = SDF File of the ligand
-n = Name of the ligand in the pdf file
Optional:
-b = Binding mode treshold, which is used to remove interactions under a certain procentual occurence treshold from the binding mode generation, default is 40 (values: 0-100)
-df = Dataframe File, if all of the interaction already were calculated. The default name of this file, which is obtained after calculating the interactions is interactions_gathered.csv
-m = Minimal Transition, which is a treshold applied for the display of the binding mode transitions via a markov chains network figure, the defaul value is 1.

Command line example with default values

    openmmdl-analysis -t {path/to/topology} -d {path/to/trajectory} -l {path/to/sdf_ligand} -n {Ligand_name}


## Running OpenMMDL-Simulation test simulations

There are two Systems prepared for the testing of the simulation.

 1: A 10ns simulation of the 6b73 Protein-ligand complex with POPC membrane and TIP3P-FB water. To run the testing of 6b73 enter the following command line:

    openmmdl-simulation -f 6b73_testing_simulation -t ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73-moe-processed_openMMDL.pdb -s ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_simulation.py -l  ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_lig.sdf

 2: A 10ns simulation of the 5wyz Protein-ligand complex with TIP3P water. To run the testing of 5wyz enter the following command line:

    openmmdl-simulation -f 5wyz_testing_simulation -t ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz-moe-processed_openMMDL.pdb -s ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz_simulation.py -l  ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5VF.sdf

## Running OpenMMDL-Simulations using slurm (multiple replicas are possible)

Two scripts are needed to run simulations via slurm. Start using the runOpenMM_slurm.sh bash script when being in the repository folder. It has several inputs. For help just type: 

    bash runOpenMM_slurm.sh
    
without any flags, it will list the flags needed. The simplest way to run a simulation is to use the "-i" flag, which takes an input directory including the simulationscript, the topology and optionally the ligand file and it will create the outputs folder within the given directory. NOTE: ake sure only one topology is present in the input folder, so that it finds it automatically.

The script calls a second script (you don't need to do that) that is used for slurms "sbatch" command to run multiple replicas. The second script can be left where it is and named how it is (SlurmWrap.sh). 

One example line of how to start a five replicas on cn-gpus would be: 

    bash runOpenMM.sh -i /mdspace/leon-moveWIP/Ligand-search-stability-check/ALk2R206H-D207Q-backmut/simulation -n 5 -c

## Copyright
Copyright (c) 2022, Valerij Talagayev & Yu Chen (Wolber lab)

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

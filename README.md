OpenMMDL @ Molecular Design Lab
==============================
[//]: # (Badges)
[![CI-CD](https://github.com/wolberlab/OpenMMDL/actions/workflows/CI-CD.yml/badge.svg)](https://github.com/wolberlab/OpenMMDL/actions/workflows/CI-CD.yml)
[![codecov](https://codecov.io/gh/talagayev/OpenMMDL/graph/badge.svg?token=950HZ93CKS)](https://codecov.io/gh/talagayev/OpenMMDL)
[![Documentation Status](https://readthedocs.org/projects/openmmdl/badge/?version=latest)](https://openmmdl.readthedocs.io)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![CodeQL](https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml/badge.svg)](https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml)


<img src="https://github.com/wolberlab/OpenMMDL/blob/main/openmmdl/openmmdl_setup/static/OpenMMDL_logo_2.png" height="250">

Interface to OpenMM for easy setup of molecular dynamic simulations of
protein-ligand complexes

http://openmmdl.readthedocs.io/

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
    conda create -n openmmdl -c conda-forge --file requirements.txt

After installation, activate the conda environment:

    conda activate openmmdl

#### Install the openmmdl package with pip

    pip install .

## OpenMMDL-Setup

OpenMMDL-Setup will setup an MD simulation for OpenMM (i.e. create all
the needed files for running the simulation)

### Usage

Start the setup process by executing the command:

    openmmdl_setup

The setup interface is displayed through a web browser, but it is still
a single-user desktop application, not a web application. It should
automatically open a web browser displaying the user interface. If that does not happen for any reason, open a browser and point it to
the address displayed in the console window (e.g. http://127.0.0.1:5000).

Download the processed PDB File and python script, which serve as input
for the simulation script.

## OpenMMDL-Simulation

OpenMMDL-Simulation starts the MD simulation with the inputs acquired
from OpenMMDL-Setup

### Usage

Start the simulation with the following Inputs:

#### Mandatory:
-f = name of the folder, where the MD simulation files are stored

-t = Topology file of your protein from the setup including the path

-s = Python script from the setup including the path

#### Optional:
-l = SDF File of the ligand, if the ligand was selected during the
Setup. The SDF file name should be consistent with the input in the setup

-c = Coordinates file of Amber


#### Command line example with ligand

    openmmdl_simulation -f {path/to/folder_name} -t {path/to/topology} -s {path/to/script} -l {path/to/ligand}

## OpenMMDL-Analysis

OpenMMDL-Analysis performs an analysis of the MD simulation obtained from OpenMMDL-Simulation.
It analyzes the protein-ligand complex interactions throughout the MD trajectory, delivering the list of
all possible interactions. In addition, it generates interaction fingerprints, provides the most occurring binding modes
, and displays the transitions between the separate binding modes.

If there is no ligand given, OpenMMDL-Analysis will instead analyze the trajectory on stable watermolecules
and cluster those at positions where in at least 75% of the MD a watermolecule is present. Outputs include a PDB with representative waters
and a CSV of nearby protein Residuenumbers and chains as well as PDBs of each water cluster.


### Usage

Start the analysis with the following Inputs:

#### Mandatory:
-t = Topology file (PDB)

-d = Trajectory file (DCD)

-l = SDF File of the ligand

-n = Name of the ligand in the PDB file

#### Optional:
-b = Binding mode treshold, which is used to remove interactions under a certain procentual occurence treshold from the binding mode generation, default is 40 (values: 0-100)

-df = Dataframe File, if all of the interaction already were calculated. The default name of this file, which is obtained after calculating the interactions is interactions_gathered.csv

-m = Minimal Transition, which is a treshold applied for the display of the binding mode transitions via a markov chains network figure, the defaul value is 1.

-c = CPU count, which will be used to calculate interactions. The default value is half of the CPU count.

#### Command line example with default values

    openmmdl_analysis -t {path/to/topology} -d {path/to/trajectory} -l {path/to/sdf_ligand} -n {Ligand_name}


#### Visualization
Most of the analysis outputs are JEPG images and do not need any further preparation to be viewed.

For the visualization of your trajectory with interaction pointclouds you can use the jupyter notebook prepared in the OpenMMDL repository.

Or use this command:
```
openmmdl_visualization
```
Then edit the notebook to include the output of your analyis.
## Copyright
Copyright (c) 2022, Valerij Talagayev, Yu Chen,  Niklas Piet Doering & Leon Obendorf (Wolber lab)

#### Acknowledgements

The Script is based upon the OpenMM Toolkit [OpenMM](https://github.com/openmm)
Many thanks to all of the creators and contributors of the OpenMM Toolkit, especially a big thank you to [jchodera](https://github.com/jchodera), [peastman](https://github.com/peastman), [mikemhenry](https://github.com/mikemhenry) and the whole [choderalab](https://github.com/choderalab) 

The Simulation script is heavily inspired by the CADD T019 Talkatorial
(https://github.com/volkamerlab/teachopencadd/tree/master/teachopencadd/talktorials/T019_md_simulation)
Thanks to the members of the [Volkamer Lab](https://volkamerlab.org/),
especially [dominiquesydow](https://github.com/dominiquesydow/), [schallerdavid](https://github.com/schallerdavid) and [AndreaVolkamer](https://github.com/andreavolkamer).

The Setup script is a modified version of [openmm-setup](https://github.com/openmm/openmm-setup).
 
This Project is based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.


[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Documentation Status](https://readthedocs.org/projects/openmmdl/badge/?version=latest)](https://openmmdl.readthedocs.io/en/latest/?badge=latest)
[![CodeQL](https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml/badge.svg)](https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml)

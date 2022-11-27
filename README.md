<img src="https://openmm.org/images/logo.svg" height="200">


# OpenMMDL Simulation


Interface to OpenMM for the easy setup of molecular dymaic simulations of protein-ligand complexes

## Install

#### Clone this repository

Open a new terminal and clone this repository
```bash   
    cd ~
    git clone 
```
#### Install all required dependencies in an separate enviorement


OpenMMDL is written in python 3.10 and uses diverse packages, which can be easily be installed on a separate enviroment using conda:
```bash
    cd openmmdl
    conda create -n openmmdl --file requirements.txt
```

#### Create alias for setup and simulation for your bash
```bash
    echo 'alias openmmdl-setup="python3 ~/openmmdl/openmmdl_setup/openmmdlsetup.py"' >> ~/.bashrc
    echo 'alias openmmdl-simulation="python3 ~/openmmdl/openmmdl_simulation/openmmdlsimulation.py"' >> ~/.bashrc
```

## Using OpenMMDL-Setup

Activate conda enviorement:

```bash
    conda activate openmmdl
```

Start the Setup by excecuting the command:

```bash
    openmmdl-setup
```

The setup interface is displayed through a web browser, but it is still a single user desktop application, not a web application. It should automatically open a web browser displaying the user interface. If for any reason that does not happen, open a browser yourself and point it to the address displayed in the console window (usually http://127.0.0.1:5000).

Download the processed PDB File and python script, which serve as input for the simulation script.

## Using OpenMMDL-Simulation

Activate conda enviorement if not activated:

```bash
    conda activate openmmdl
```

Start the Simulation with the following Inputs:
Mandatory:
-f = Name of folder, where the MD simulation files will be stored
-p = PDB file of your protein from the Setup including path
-s = Python script from the Setup including path
Optional:
-l = SDF File of the ligand, if the ligand was selected during the Setup. sdf file name should be consistent with the input in setup

Command line example without ligand
```bash
    openmmdl-simulation -f {path/to/folder_name} -p {path/to/protein.pdb} -s {path/to/script.py}
```

Command line example with ligand
```bash
    openmmdl-simulation -f {path/to/folder_name} -p {path/to/protein.pdb} -s {path/to/script.py} -l {path/to/ligand.sdf}
```

## Running OpenMMDL-Simulation test simulations

There are two Systems prepared for the testing of the simulation.




#### Acknowledgements

The Script is based upon the OpenMM Toolkit [OpenMM](https://github.com/openmm)
Many thanks to all of the creators and contributers of the OpenMM Toolkit, especially a big thank you to [jchodera](https://github.com/jchodera), [peastman](https://github.com/peastman), [mikemhenry](https://github.com/mikemhenry) and the whole [choderalab](https://github.com/choderalab) 

The Simulation script is heavly inspired by the CADD T019 Talkatorial
(https://github.com/volkamerlab/teachopencadd/tree/master/teachopencadd/talktorials/T019_md_simulation)
Thanks to the members of the [Volkamer Lab](https://volkamerlab.org/),
especially [dominiquesydow](https://github.com/dominiquesydow/), [schallerdavid](https://github.com/schallerdavid) and [AndreaVolkamer](https://github.com/andreavolkamer).

The Setup script is a modified version of [openmm-setup](https://github.com/openmm/openmm-setup).
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.


[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
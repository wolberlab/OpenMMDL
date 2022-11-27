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
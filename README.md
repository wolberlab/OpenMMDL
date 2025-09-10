OpenMMDL @ Molecular Design Lab
==============================

<img src="https://github.com/wolberlab/OpenMMDL/blob/main/openmmdl/openmmdl_setup/static/OpenMMDL_logo_2.png" height="250">

[//]: # (Badges)

| **Release**        | [![Last release tag][badge_release]][url_latest_release] ![GitHub commits since latest release (by date) for a branch][badge_commits_since]|
| :----------------- | :------- |
| **Availability**   | [![Condaforge][badge_conda_forge]][url_conda_forge] ![Docker][badge_docker]|
| **Workflows**      | ![CI_CD][badge_ci_cd] ![CodeQL][badge_codeql]|
| **Documentation**  | [![docs][badge_docs]][url_docs] [![codecov][badge_codecov]][url_codecov]|
| **Dependencies**   | [![rdkit][badge_rdkit]][url_rdkit] [![mdanalysis][badge_mda]][url_mda] [![Code style: black][badge_black]][url_black]|
| **License**        | [![License: GPL v2][badge_license]][url_license]|


Interface for easy setup and analysis of molecular dynamic (MD) simulations of
protein-ligand complexes with **OpenMM**

http://openmmdl.readthedocs.io/

## Installation via conda-forge

**OpenMMDL** is implemented in conda-forge and can be installed for linux-based system with the following command line

    conda install -c conda-forge openmmdl

## Installation via docker

**OpenMMDL** is mainly supported for Linux distribution systems, thus for Windows and MacOS the installation with docker may be preferred,
due to docker creating an image with ubuntu:22.04

For this first clone this repository

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL

change the path to the OpenMMDL folder and build the docker image from there

    cd OpenMMDL
    docker build -t openmmdl_env .

This will build the OpenMMDL image with docker. Now that it is build you can access it through an interactive terminal

    docker run -it --name openmmdl_test openmmdl_env

From there you can access all the OpenMMDL entry points. Currently due to OpenMMDL Setup using flask it can be difficult
to access it through the docker image.

## Installation via repository

#### Clone this repository

Open a new terminal and clone this repository

    cd ~
    git clone https://github.com/wolberlab/OpenMMDL

#### Install all required dependencies in a separate environment

**OpenMMDL** is written in Python 3.10 and uses several packages, which can
be easily installed on a separate environment using conda (we recommend
using miniconda):

    cd OpenMMDL
    conda create -n openmmdl -c conda-forge --file requirements.txt

After installation, activate the conda environment:

    conda activate openmmdl

#### Install the openmmdl package with pip

    pip install .

## OpenMMDL Setup

**OpenMMDL Setup** will allow you to prepare the files needed to perform a  protein-ligand complex MD simulation with **OpenMM**.

### Usage

Start **OpenMMDL Setup** by executing the command:

    openmmdl_setup

The **OpenMMDL Setup** interface is displayed through a web browser, but it is still
a single-user desktop application, not a web application. It should
automatically open a web browser displaying the user interface. If that does not happen for any reason, open a browser and point it to
the address displayed in the console window (e.g. http://127.0.0.1:5000).

Download the processed PDB file and Python script, which will serve as input
for the **OpenMMDL Simulation** script.

## OpenMMDL Simulation

**OpenMMDL Simulation** starts the MD simulation with the inputs acquired
from **OpenMMDL Setup**.

### Usage

Start the simulation with the following Inputs:

#### Mandatory:
-f = Name of the folder, where the MD simulation files are stored

-t = Topology file of your protein from the setup including the path

-s = Python script from the setup including the path

#### Optional:
-l = SDF file of the ligand, if the ligand was selected during OpenMMDL
Setup. The SDF file name should be consistent with the input in the setup

-c = Coordinates file of Amber


#### Command line example with ligand

    openmmdl_simulation -f {path/to/folder_name} -t {path/to/topology} -s {path/to/script} -l {path/to/ligand}

## OpenMMDL Analysis

**OpenMMDL Analysis** performs an analysis of the MD simulation obtained from **OpenMMDL Simulation**.
It analyzes the protein-ligand complex interactions throughout the MD trajectory, delivering the list of
all possible interactions. In addition, it generates interaction fingerprints, provides the most occurring so-called Binding Modes
and displays the transitions between the separate binding modes.

If there is no ligand given, **OpenMMDL Analysis** will instead analyze the trajectory on stable water molecules
and cluster those at positions where in at least 75% frames of the MD a water molecule is present. Outputs include a PDB with representative waters
and a CSV of nearby protein residue numbers and chains as well as PDBs of each water cluster.


### Usage

Start the analysis with the following Inputs:

#### Mandatory:
-t = topology file of the simulation (in .pdb format)

-d = trajectory file of the simulation (in .dcd format)

-n = Ligand name (3 letter code in PDB)

#### Optional:

-l = Ligand in SDF format

-b = binding mode threshold. Is used to remove interactions under the defined procentual occurence from the binding mode generation. The default is 40% (accepted values: 0-100)

-df = Dataframe (use if the interactions were already calculated, default name would be "interactions_gathered.csv")

-m = minimal transition threshold. Is used for the display of the binding mode transitions in the Markov state chains network figure. The default value is 1

-c = CPU count, specify how many CPUs should be used, default is half of the CPU count.

-p = Generate .pml files for pharmacophore visualization. The default is False (accepted values: True/False)

-s = special ligand name to calculate interactions with special ligands.

-nuc = Treat nucleic acids as receptor

-pep = Calculate interactions with peptides. Give the peptide chain name as input. Defaults to None

-ref = Add a reference PDB to renumber the residue numbers. Defaults to None (accepted values: str of PDB)

-r = Calculate the RMSD difference between frames. The default is False (accepted values: True/False)

-w = stable-water-analysis. Defines if the analysis of stable water molecules should be performed. The default is False (accepted values: True/False)

--watereps = the EPS of the clustering part during the water analysis. will only result in something if "-w True" is added. Accepts float (in Angstrom).

--figure = File type for the figures, default is png. Can be changed to all file types supported by matplotlib.

#### Command line example with default values

    openmmdl_analysis -t {path/to/topology} -d {path/to/trajectory} -n {Ligand_name}


#### Visualization
Most of the analysis outputs are JPEG images and do not need any further preparation to be viewed.

For the visualization of your complex with interaction pointclouds you can use NGLView with the jupyter notebook prepared in the **OpenMMDL** repository or visualize the pointclouds in PyMOL.

### Usage
```
openmmdl_visualization
```
#### Optional:
--type = Software you wish to visualize openmmdl interactions with. Options: nglview, pymol. Default: nglview
#### NGLView
After running the start comand a jupyter notebook will automatically open.
Edit the notebook to include the output files of your analysis.
Then run all cells.
#### PyMOL
After running the start comand a python skript will apear in your directory.
Open up PyMOL then run these two comands in the PyMOL console:
```
run visualization_pymol.py
openmdl_visualization PATH_TO_interacting_waters.pdb, modulePATH_TO_clouds.json
```
## Copyright
Copyright (c) 2022, Valerij Talagayev, Yu Chen,  Niklas Piet Doering & Leon Obendorf (Wolber lab)

### Citation
If you use OpenMMDL in your research, please cite the following [paper](https://doi.org/10.1021/acs.jcim.4c02158).

#### Acknowledgments

The Script is based upon the OpenMM Toolkit [OpenMM](https://github.com/openmm)
Many thanks to all of the creators and contributors of the OpenMM Toolkit, especially a big thank you to [jchodera](https://github.com/jchodera), [peastman](https://github.com/peastman), [mikemhenry](https://github.com/mikemhenry) and the whole [choderalab](https://github.com/choderalab) 

The Simulation script is heavily inspired by the CADD T019 Talkatorial
(https://github.com/volkamerlab/teachopencadd/tree/master/teachopencadd/talktorials/T019_md_simulation)
Thanks to the members of the [Volkamer Lab](https://volkamerlab.org/),
especially [dominiquesydow](https://github.com/dominiquesydow/), [schallerdavid](https://github.com/schallerdavid) and [AndreaVolkamer](https://github.com/andreavolkamer).

The Setup script is a modified version of [openmm-setup](https://github.com/openmm/openmm-setup).
 
This Project is based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.


[badge_release]: https://img.shields.io/github/release-pre/wolberlab/openmmdl.svg
[badge_commits_since]: https://img.shields.io/github/commits-since/wolberlab/openmmdl/latest
[badge_ci_cd]: https://github.com/wolberlab/OpenMMDL/actions/workflows/CI-CD.yml/badge.svg
[badge_codeql]: https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml/badge.svg
[badge_docs]: https://readthedocs.org/projects/openmmdl/badge/?version=latest
[badge_codecov]: https://codecov.io/gh/talagayev/OpenMMDL/graph/badge.svg?token=950HZ93CKS
[badge_conda_forge]: https://img.shields.io/conda/vn/conda-forge/openmmdl.svg
[badge_rdkit]: https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC
[badge_mda]: https://img.shields.io/badge/powered%20by-MDAnalysis-orange.svg?logoWidth=16&logo=data:image/x-icon;base64,AAABAAEAEBAAAAEAIAAoBAAAFgAAACgAAAAQAAAAIAAAAAEAIAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJD+XwCY/fEAkf3uAJf97wGT/a+HfHaoiIWE7n9/f+6Hh4fvgICAjwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACT/yYAlP//AJ///wCg//8JjvOchXly1oaGhv+Ghob/j4+P/39/f3IAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAJH8aQCY/8wAkv2kfY+elJ6al/yVlZX7iIiI8H9/f7h/f38UAAAAAAAAAAAAAAAAAAAAAAAAAAB/f38egYF/noqAebF8gYaagnx3oFpUUtZpaWr/WFhY8zo6OmT///8BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAgICAn46Ojv+Hh4b/jouJ/4iGhfcAAADnAAAA/wAAAP8AAADIAAAAAwCj/zIAnf2VAJD/PAAAAAAAAAAAAAAAAICAgNGHh4f/gICA/4SEhP+Xl5f/AwMD/wAAAP8AAAD/AAAA/wAAAB8Aov9/ALr//wCS/Z0AAAAAAAAAAAAAAACBgYGOjo6O/4mJif+Pj4//iYmJ/wAAAOAAAAD+AAAA/wAAAP8AAABhAP7+FgCi/38Axf4fAAAAAAAAAAAAAAAAiIiID4GBgYKCgoKogoB+fYSEgZhgYGDZXl5e/m9vb/9ISEjpEBAQxw8AAFQAAAAAAAAANQAAADcAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAjo6Mb5iYmP+cnJz/jY2N95CQkO4pKSn/AAAA7gAAAP0AAAD7AAAAhgAAAAEAAAAAAAAAAACL/gsAkv2uAJX/QQAAAAB9fX3egoKC/4CAgP+NjY3/c3Nz+wAAAP8AAAD/AAAA/wAAAPUAAAAcAAAAAAAAAAAAnP4NAJL9rgCR/0YAAAAAfX19w4ODg/98fHz/i4uL/4qKivwAAAD/AAAA/wAAAP8AAAD1AAAAGwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAALGxsVyqqqr/mpqa/6mpqf9KSUn/AAAA5QAAAPkAAAD5AAAAhQAAAAEAAAAAAAAAAAAAAAAAAAAAAAAAAAAAADkUFBSuZ2dn/3V1df8uLi7bAAAATgBGfyQAAAA2AAAAMwAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAB0AAADoAAAA/wAAAP8AAAD/AAAAWgC3/2AAnv3eAJ/+dgAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA9AAAA/wAAAP8AAAD/AAAA/wAKDzEAnP3WAKn//wCS/OgAf/8MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAIQAAANwAAADtAAAA7QAAAMAAABUMAJn9gwCe/e0Aj/2LAP//AQAAAAAAAAAA
[badge_black]: https://img.shields.io/badge/code%20style-black-000000.svg
[badge_license]: https://img.shields.io/badge/License-MIT-blue.svg
[badge_docker]: https://img.shields.io/badge/Docker-2496ED?logo=docker&logoColor=fff

[url_latest_release]: https://github.com/wolberlab/openmmdl/releases
[url_ci_cd]: https://github.com/wolberlab/OpenMMDL/actions/workflows/CI-CD.yml
[url_codeql]: https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml
[url_docs]: https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml
[url_codecov]: https://github.com/wolberlab/OpenMMDL/actions/workflows/codeql.yml
[url_conda_forge]: https://anaconda.org/conda-forge/openmmdl
[url_rdkit]: https://www.rdkit.org/
[url_mda]: https://www.mdanalysis.org
[url_black]: https://github.com/psf/black
[url_license]: https://opensource.org/licenses/MIT

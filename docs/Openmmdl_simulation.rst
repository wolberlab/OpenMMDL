**Running OpenMMDL-Simulation**
===============

OpenMMDL-Simulation starts the MD simulation with the inputs acquired from OpenMMDL-Setup.

This page details the variables that are required for starting the script and showcases the application of OpenMMDL-Simulation. 

Variables
------------------------------
OpenMMDL-Simulation consists of mandatory and optional variables. The following are listed down below:


Mandatory:

.. code-block:: text

    -f = name of the folder, where the MD simulation files are stored
    -t = Topology file of your protein from the setup including the path
    -s = Python script from the setup including the path

Optional:

.. code-block:: text

    -l = SDF File of the ligand. SDF file name should be consistent with the input in setup
    -c = Cooridnates file of Ambe

Application
------------------------------

An example of how a command line of OpenMMDL-Simulation should look is:

.. code-block:: text

    openmmdl-simulation -f {path/to/folder_name} -t {path/to/topology} -s {path/to/script} -l {path/to/ligand}


For help during usage of OpenMMDL-Simulation use the following command line:

.. code-block:: text

    openmmdl-simulation -h all

Running OpenMMDL-Simulation test simulations
------------------------------
There are two Systems prepared for the testing of the simulation.

1: A 10ns simulation of the 6b73 Protein-ligand complex with POPC membrane and TIP3P-FB water. To run the testing of 6b73 enter the following command line:

.. code-block:: text

    openmmdl-simulation -f 6b73_testing_simulation -t ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73-moe-processed_openMMDL.pdb -s ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_simulation.py -l  ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_lig.sdf

2: A 10ns simulation of the 5wyz Protein-ligand complex with TIP3P water. To run the testing of 5wyz enter the following command line:

.. code-block:: text

    openmmdl-simulation -f 5wyz_testing_simulation -t ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz-moe-processed_openMMDL.pdb -s ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz_simulation.py -l  ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5VF.sdf

Each of the command lines should generate a folder, where the script and the input data will be moved and further perform a MD simulation and postprocessing of the systems.

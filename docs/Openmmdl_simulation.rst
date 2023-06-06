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

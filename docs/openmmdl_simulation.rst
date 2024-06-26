**Running OpenMMDL Simulation**
===============
This page details the variables required to run the simulation and showcases the application of **OpenMMDL Simulation**.

.. figure:: /_static/images/OpenMMDL_Simulation_logo.png
    :figwidth: 600px
    :align: center

|  
**OpenMMDL Simulation** starts the MD simulation with the inputs acquired from **OpenMMDL Setup**.

This page details the variables that are required for starting the script and showcases the application of **OpenMMDL Simulation**. 

Variables
------------------------------
**OpenMMDL Simulation** consists of mandatory and optional variables. The following are listed down below:


Mandatory:

.. code-block:: text

    -f = name of the folder, where the MD simulation files are stored
    -t = Topology file of your protein from the setup including the path
    -s = Python script from the setup including the path

Optional:

.. code-block:: text

    -l = SDF file of the ligand. The SDF file name should be consistent with the input in OpenMMDL Setup
    -c = Coordinates file of Amber

Application
------------------------------

An example of how a command line of **OpenMMDL Simulation** should look is:

.. code-block:: text

    openmmdl_simulation -f {path/to/folder_name} -t {path/to/topology} -s {path/to/script} -l {path/to/ligand}


For help during the usage of **OpenMMDL Simulation** use the following command line:

.. code-block:: text

    openmmdl_simulation -h all

Running OpenMMDL Simulation test simulations
------------------------------
There are two Systems prepared for the testing of the simulation.

1: A 10ns simulation of the 6B73 Protein-ligand complex with POPC membrane and TIP3P-FB water. To run the testing of 6b73 enter the following command line:

.. code-block:: text

    openmmdl_simulation -f 6b73_testing_simulation -t ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73-moe-processed_openMMDL.pdb -s ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_simulation.py -l  ~/OpenMMDL/openmmdl_simulation/testing_sytems/6b73_membrane/6b73_lig.sdf

2: A 10ns simulation of the 5WYZ Protein-ligand complex with TIP3P water. To run the testing of 5WYZ enter the following command line:

.. code-block:: text

    openmmdl_simulation -f 5wyz_testing_simulation -t ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz-moe-processed_openMMDL.pdb -s ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5wyz_simulation.py -l  ~/OpenMMDL/openmmdl_simulation/testing_sytems/5wyz_solvent/5VF.sdf

Each of the command lines should generate a folder, where the script and the input data will be moved and further perform a MD simulation and postprocessing of the systems.

Running OpenMMDL Simulation using slurm
------------------------------
Two scripts are needed to run simulations via slurm. Start using the `runOpenMM_slurm.sh` bash script when being in the repository folder. It has several inputs. For help just type:

.. code-block:: text

    bash runOpenMM_slurm.sh
    
without any flags, it will list the flags needed. The simplest way to run a simulation is to use the `-i` flag, which takes an input directory including the simulation script, the topology and optionally the ligand file and it will create the outputs folder within the given directory. NOTE: make sure only one topology is present in the input folder so that it finds it automatically.

The script calls a second script (you don't need to do that) that is used for slurms `sbatch` command to run multiple replicas. The second script can be left where it is and named how it is (SlurmWrap.sh).

One example line of how to start five replicas on cn-gpus would be:

.. code-block:: text

    bash runOpenMM.sh -i /path/to/Input/Folder -n 5 -c

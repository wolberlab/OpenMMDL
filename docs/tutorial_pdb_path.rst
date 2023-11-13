**PDB Path**
==============

This is a tutorial how to use the PDB Path
The system that we will use for this tutorial the Toll-like receptor 8.

Starting OpenMMDL-Setup
------------------------------

We start the tutorial by creating an folder where we will copy all the input files and later run the MD Simulation.

To create a new folder we enter the following command lines:


.. code-block:: text

    mkdir openmmdl_pdb_tutorial


activating our openmmdl environment to start the OpenMMDL-Setup.
To start the OpenMMDL-Setup we need to activate the openmmdl environment. to do this we have to enter the following command lines:

.. code-block:: text

    conda activate openmmdl

Now that we have activated the openmmdl environment we can start OpenMMDL-Setup. To do this you need to type the following:

.. code-block:: text

    openmmdl_setup


Selecting Input Files
------------------------------

This will open the OpenMMDL-Setup, which you can use for the creation of the input files for OpenMMDL-Simulation.



.. figure:: /_static/images/tutorials/PDB_Path/Selectfiletype.png
   :figwidth: 725px
   :align: center


Here We select the PDB File Option and press continue.


This leads us to the next page, where we can see the selection of the PDB and SDF File of the protein and ligand respectively.


.. figure:: /_static/images/tutorials/PDB_Path/Selectfiles1.png
   :figwidth: 725px
   :align: center


Go to the PDB File Browse button and select 5wyz.pdb in the folder. 

Go to the SDF File Browse button and select the 5VF.sdf file in the folder.

The page should look like the page below.

   
.. figure:: /_static/images/tutorials/PDB_Path/Selectfiles2.png
   :figwidth: 725px
   :align: center

Now that we selected the PDB and SDF File we select continue and go to the next page.

Selecting Chains
------------------------------

The next page shows us the protein with all of the chains that are in the pdb file.


.. figure:: /_static/images/tutorials/PDB_Path/Modifychains1.png
   :figwidth: 725px
   :align: center


During this simulation we only require the chains A and B, thus we deselect the remaining chains.


.. figure:: /_static/images/tutorials/PDB_Path/Modifychains2.png
   :figwidth: 725px
   :align: center

Now that we have only the chains A and B selected we select Continue and go to the next page.

Optional: Adding Residues
------------------------------

The next page shows that the PDB is missing residues, which can be added.

.. figure:: /_static/images/tutorials/PDB_Path/Addresidues1.png
   :figwidth: 725px
   :align: center

We add the missing residues from 101 to 111 and from 435/438 to 460 in the chains A and B.

.. figure:: /_static/images/tutorials/PDB_Path/Addresidues2.png
   :figwidth: 725px
   :align: center

Now that we selected the residues that we want to add we continue to the next page.

Optional: Adding Heavy Atoms
------------------------------

The next page shows the heavy atoms that are missing and need to be added.

.. figure:: /_static/images/tutorials/PDB_Path/Addheavyatoms.png
   :figwidth: 725px
   :align: center

We press continue and go to the next page.

Adding Hydrogens and Water Box
------------------------------

The next page allows us to add the missing hydrogens and add a water box or an membrane.

.. figure:: /_static/images/tutorials/PDB_Path/Addwater1.png
   :figwidth: 725px
   :align: center

We need to add hydrogens for the MD Simulation at a pH at 7.4 so we change the number from 7.0 to 7.4.

Additionally we want to add a Water Box, so we select the Water Box Option.

In the Water Box we also change the Ionic strenght to 0.15.

.. figure:: /_static/images/tutorials/PDB_Path/Addwater2.png
   :figwidth: 725px
   :align: center

Now that we selected the following options we continue to the next page.

Adjusting MD Simulation Options
------------------------------

The final page displays the MD simulation options.

.. figure:: /_static/images/tutorials/PDB_Path/Simulationoptions1.png
   :figwidth: 725px
   :align: center

We will simulate the protein for 10 ns.

For this we go to the Simulation Tab and Change Simulation lenght (steps) to 250000.

If you want to have a longer or shorter simulation you can increase or reduce the amount of steps.

.. figure:: /_static/images/tutorials/PDB_Path/Simulationoptions2.png
   :figwidth: 725px
   :align: center

Now that we changed the amount of steps we select the Save Script button to save the script for the simulation.

Select the Save Processed PDF File to save the PDB File that will be the input for the MD simulation.

Running Tutorial OpenMMDL-Simulation
------------------------------

Create a separate folder and copy the Simulation script, Processed PDB File and the Ligand SDF File into the folder.

The SDF File should be the same that was used as an input for the Openmm-Setup.

.. figure:: /_static/images/tutorials/PDB_Path/Inputfiles.png
   :figwidth: 725px
   :align: center

Now that we have the files in one folder we can start the MD simulation.

For this we start by activating the environment

.. code-block:: text

    conda activate openmmdl

Now that activated the environment we start the simulation.

For this enter the following command

.. code-block:: text

    openmmdl-simulation -f tutorial_simulation -s OpenMMDL_Simulation.py -t 5wyz-processed_openMMDL.pdb -l 5VF.sdf

By entering the command we create a folder called tutorial_simulation, where the Output of the MD simulation will appear.

As the Input for the MD simulation we used the -t to select 5wyz-processed_openMMDL.pdb as the topology file for the simulation, -l to select the ligand 5VF.sdf and -s to specify the OpenMMDL_Simulation.py script that will run the MD simulation.

.. figure:: /_static/images/tutorials/PDB_Path/Outputfiles1.png
   :figwidth: 725px
   :align: center

During and after simulation you can open the folder to see the progress.

After the simulation is finished the tutorial_simulation should look like the picture below.

.. figure:: /_static/images/tutorials/PDB_Path/Outputfiles2.png
   :figwidth: 725px
   :align: center

If there are files or folders missing, repeat the MD simulation.

This concludes the Tutorial for the OpenMMDL PDB Path simulations.

To see what the separate files in the Output represent follow this page:

* :doc:`MD Simulation Output </simulation_output>`

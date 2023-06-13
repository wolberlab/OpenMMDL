**PDB Path**
==============

This is a tutorial how to use the PDB Path
The system that we will use for this tutorial the Toll-like receptor 8.


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

    openmmdl-setup

This will open the OpenMMDL-Setup, which you can use for the creation of the input files for OpenMMDL-Simulation.



.. figure:: /_static/images/tutorials/PDB_Path/Selectfiletype.png
   :figwidth: 725px
   :align: center


Here We select the PDB File Option and press continue.


This leads us to the next page, where we can see the selection of the PDB and SDF File of the protein and ligand respectively.

.. figure:: /_static/images/tutorials/PDB_Path/Selectfiles1.png
   :figwidth: 725px
   :align: center

Go to the PDB File Browse button and select 5wyz.pdb in the folder. Afterwards go to the SDF File Browse button and select the 5VF.sdf file in the folder.

The page should look like the page below.


   
.. figure:: /_static/images/tutorials/PDB_Path/Selectfiles2.png
   :figwidth: 725px
   :align: center

Now that we selected the PDB and SDF File we select continue and go to the next page.

The next page shows us the protein with all of the chains that are in the pdb file.
.. figure:: /_static/images/tutorials/PDB_Path/Modifychains1.png
   :figwidth: 725px
   :align: center

During this simulation we only require the chains A and B, thus we deselect the remaining chains.


.. figure:: /_static/images/tutorials/PDB_Path/Modifychains2.png
   :figwidth: 725px
   :align: center

sus5


.. figure:: /_static/images/tutorials/PDB_Path/Addresidues1.png
   :figwidth: 725px
   :align: center

sus6

.. figure:: /_static/images/tutorials/PDB_Path/Addresidues2.png
   :figwidth: 725px
   :align: center

sus7

.. figure:: /_static/images/tutorials/PDB_Path/Addheavyatoms.png
   :figwidth: 725px
   :align: center

sus8

.. figure:: /_static/images/tutorials/PDB_Path/Addwater1.png
   :figwidth: 725px
   :align: center

sus9

.. figure:: /_static/images/tutorials/PDB_Path/Addwater2.png
   :figwidth: 725px
   :align: center

sus10

.. figure:: /_static/images/tutorials/PDB_Path/Simulationoptions1.png
   :figwidth: 725px
   :align: center

sus11

.. figure:: /_static/images/tutorials/PDB_Path/Simulationoptions2.png
   :figwidth: 725px
   :align: center

sus12

.. figure:: /_static/images/tutorials/PDB_Path/Inputfiles.png
   :figwidth: 725px
   :align: center

sus13

.. figure:: /_static/images/tutorials/PDB_Path/Outputfiles1.png
   :figwidth: 725px
   :align: center

sus14

.. figure:: /_static/images/tutorials/PDB_Path/Outputfiles2.png
   :figwidth: 725px
   :align: center

sus15

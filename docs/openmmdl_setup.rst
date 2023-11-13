**Running OpenMMDL-Setup**
=============================

This is a detailed explanation on how to run the OpenMMDL Setup.

To start the OpenMMDL-Setup we need to activate the openmmdl environment. to do this we have to enter the following command lines:

.. code-block:: text

    conda activate openmmdl

Now that we have activated the openmmdl environment we can start OpenMMDL-Setup. To do this you need to type the following:

.. code-block:: text

    openmmdl_setup

This will open the OpenMMDL-Setup, which you can use for the creation of the input files for OpenMMDL-Simulation.

There are two possible options to create the input files for OpenMMDL-Simulation:

1. The PDB Path, where a PDB of the protein is used as an input for the preparation and simulation.
The tutorial for the PDB Path can be found :doc:`here </tutorial_pdb_path>`.

Here is the table of the currently avaiable forcefields and watermodels for the PDB path: 

.. figure:: /_static/images/Forcefield_watermodels.png
   :figwidth: 725px
   :align: center

2. The Amber Path, where `prmtop` and `inpcrd` files are used the preparation and simulation. This path allows to either use already prepared `prmtop` and `inpcrd` as an input or create the `prmtop` and `inpcrd` from PDB files of the receptor and ligand.
The tutorial for the Amber path can be found :doc:`here </tutorial_amber_path>`.

.. figure:: /_static/images/amber_ff.png
   :figwidth: 725px
   :align: center


In the table, the first row is the default setting, and the term `other` allows users to type their desired forcefields from those accessible in AmberTools 22.0 into the designated textbox.

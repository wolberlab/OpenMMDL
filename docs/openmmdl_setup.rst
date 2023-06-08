**Running OpenMMDL-Setup**
=============================

This is a detailed explanation on how to run the OpenMMDL Setup.

To start the OpenMMDL-Setup we need to activate the openmmdl environment. to do this we have to enter the following command lines:

.. code-block:: text

    conda activate openmmdl

Now that we have activated the openmmdl environment we can start OpenMMDL-Setup. To do this you need to type the following:

.. code-block:: text

    openmmdl-setup

This will open the OpenMMDL-Setup, which you can use for the creation of the input files for OpenMMDL-Simulation.

There are two possible options to create the input files for OpenMMDL-Simulation:

1. The PDB Path, where an PDB of the protein is used as an Input for the preparation and simulation.
The tutorial for the PDB Path can be found here:

* :doc:`PDB Path Tutorial </tutorial_pdb_path>`

2. The Amber Path, where Amber prmtop and inpcrd are used the preparation and simulation. This path allows to either use already prepared prmtop and inpcrd as an input or create the prmtop and inpcrd from an PDB File

The tutorial for the Amber Path can be found here:

* :doc:`Amber Path Tutorial </tutorial_amber_path>`

Here is the table of the currently avaiable Forcefields and Watermodels for the PDB Path: 

.. figure:: /_static/images/Forcefield_watermodels.png
   :figwidth: 725px
   :align: center


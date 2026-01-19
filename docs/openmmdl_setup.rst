**Running OpenMMDL Setup**
=============================

This page displays the preparation paths and forcefields available in **OpenMMDL** and showcases the application of **OpenMMDL Setup**.

.. figure:: /_static/images/OpenMMDL_Setup.png
    :figwidth: 600px
    :height: 100px
    :align: center
|  
To start the **OpenMMDL Setup** we need to activate the `openmmdl` environment. To do this we have to enter the following command lines:

.. code-block:: text

    conda activate openmmdl

Now that we have activated the `openmmdl` environment we can start **OpenMMDL Setup**. To do this you need to type the following:

.. code-block:: text

    openmmdl setup

This will open the **OpenMMDL Setup**, which you can use for the creation of the input files for **OpenMMDL Simulation**.

There are two possible options to create the input files for **OpenMMDL Simulation**:

1. The PDBFixer path, where a `pdb` file of the protein is used as an input for the preparation and simulation.
The tutorial for the PDBFixer path can be found :doc:`here </tutorial_pdb_path>`.

Here is the table of the currently available forcefields and water models for the PDBFixer path: 

.. figure:: /_static/images/Forcefield_watermodels.png
   :figwidth: 725px
   :align: center

.. list-table:: Supported water models by protein force field (PDBFixer path)
   :header-rows: 1
   :widths: 25 75

   * - Protein force field
     - Water models selectable in OpenMMDL Setup
   * - AMBER19
     - TIP3P; TIP3P-FB; SPC/E; TIP4P-Ew; TIP4P-FB; OPC3; OPC
   * - AMBER14
     - TIP3P; TIP3P-FB; SPC/E; TIP4P-Ew; TIP4P-FB; OPC3; OPC
   * - AMBER99SB
     - TIP3P; TIP3P-FB; SPC/E; TIP4P-Ew; TIP5P; OPC3; OPC
   * - AMBER99SB-ILDN
     - TIP3P; TIP3P-FB; SPC/E; TIP4P-Ew; TIP5P; OPC3; OPC
   * - AMBER03
     - TIP3P; TIP3P-FB; SPC/E; TIP4P-Ew; TIP5P; OPC3; OPC
   * - AMBER10
     - TIP3P; TIP3P-FB; SPC/E; TIP4P-Ew; TIP5P; OPC3; OPC
   * - CHARMM36
     - CHARMM default; TIP3P-PME-B; TIP3P-PME-F; SPC/E; TIP4P-Ew; TIP4P-2005; TIP5P; TIP5P-Ew
   * - CHARMM36 2024
     - CHARMM default; TIP3P-PME-B; TIP3P-PME-F; SPC/E; TIP4P-Ew; TIP4P-2005; TIP5P; TIP5P-Ew


2. The Amber path, where `prmtop` and `inpcrd` files are used the preparation and simulation. This path allows us to either use already prepared `prmtop` and `inpcrd` as an input or create the `prmtop` and `inpcrd` from PDB files of the receptor and ligand.
The tutorial for the Amber path can be found :doc:`here </tutorial_amber_path>`.

.. figure:: /_static/images/amber_ff.png
   :figwidth: 725px
   :align: center

In the table, the first row is the default setting, and the term `other` allows users to type their desired forcefields from those accessible in AmberTools 22.0 into the designated textbox.

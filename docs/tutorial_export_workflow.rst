**Export Workflow**
===================

This tutorial shows how to use the **Export Workflow** in **OpenMMDL Setup** to prepare a protein or protein-ligand complex and write the prepared files in formats that are suitable for other molecular simulation workflows.

The export workflow is intended for users who want to carry out the structure preparation in **OpenMMDL Setup**, but continue the simulation in other engines such as **Amber**, **GROMACS**, or **CHARMM/NAMD**.

Starting OpenMMDL Setup
-----------------------

To start the export workflow, first activate the ``openmmdl`` environment:

.. code-block:: text

    conda activate openmmdl

Then launch **OpenMMDL Setup**:

.. code-block:: text

    openmmdl setup

Selecting the Input Type
------------------------

After starting the setup interface, select the **PDB File** option.

.. figure:: /_static/images/tutorials/PDB_Path/Selectfiletype.png
   :figwidth: 725px
   :align: center

This opens the PDB path.

Selecting the Workflow
----------------------

After choosing the PDB path, the next page allows the user to choose between the standard **OpenMMDL Simulation** workflow and the **Export Workflow**.

.. figure:: /_static/images/tutorials/Export_Workflow/Workflow_Selection.png
   :figwidth: 725px
   :align: center

To prepare files for use in another molecular dynamics workflow, select the **Export Workflow** option.

Preparing the Input Structure
-----------------------------

In the next step, choose the protein PDB file and, if required, the ligand SDF file. The preferred force field, water model, and small molecule force field can also be selected here.

.. figure:: /_static/images/tutorials/Export_Workflow/File_Selection.png
   :figwidth: 725px
   :align: center

At this stage, the same preparation steps used for the normal **OpenMMDL Simulation** route are available. This ensures that the structure cleanup and system preparation are performed consistently before exporting the files.

.. figure:: /_static/images/tutorials/Export_Workflow/Protein_Selection.png
   :figwidth: 725px
   :align: center

Selecting the Environment
-------------------------

The export workflow uses the same environment preparation steps as the standard simulation workflow. The user can therefore prepare the system with:

- a water box,
- or a membrane environment.

.. figure:: /_static/images/tutorials/Export_Workflow/Environment_Selection.png
   :figwidth: 725px
   :align: center

Exporting the Prepared System
-----------------------------

After the structure has been prepared, the final page of the export workflow allows the user to select the desired output format.

.. figure:: /_static/images/tutorials/Export_Workflow/Export_Options.png
   :figwidth: 725px
   :align: center

The export page provides two main sections:

- **Export format selection** on the left
- **Preparation summary** on the right

The summary reports the selected input structure, ligand, force field, solvent model, ionization, and preparation settings.

Available Export Formats
------------------------

The current export workflow supports the following output formats:

- **Processed PDB / mmCIF**
- **Amber topology and coordinates** (``.prmtop`` and ``.inpcrd``)
- **GROMACS topology and coordinates** (``.top`` and ``.gro``)
- **PSF and companion coordinates** for CHARMM/NAMD-style workflows
- **OpenMM XML** for reuse in OpenMM-based workflows

These files are generated from the prepared OpenMM system after structure cleanup, ligand merging, solvent or membrane construction, and force-field assignment.

Preparing and Saving Files
--------------------------

The export process is performed in two steps:

1. **Prepare Files**
2. **Save Files**

First, click **Prepare Files**. This runs the preparation and file writing in the background. Once the export has completed successfully, the interface reports that the prepared files are ready.

After preparation is complete, click **Save Files** to download the generated files as a compressed archive.

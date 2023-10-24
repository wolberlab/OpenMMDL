**Amber Path**
==============

Introduction
------------------

Amber files is hard to get for those don't familar with Amber language and Bash script. OpenMMDL-Setup provides a user-friendly interface to help users to generate Amber files for later MD simulation.
This tutorial provides a step-by-step guide on setting up an MD simulation using the Amber files using OpenMMDL-Setup.

1. [Starting OpenMMDL-Setup](#starting-openmmdl-setup)
2. [Selecting the Amber Path](#selecting-the-amber-path)
3. [Selecting the input files](#selecting-the-input-files)
4. [Amber Configuration](#amber-configuration)
5. [Simulation Setup](#simulation-setup)
6. [Run Simulation](#run-simulation)

Starting OpenMMDL-Setup
------------------------------
We start the tutorial by creating an folder where we will copy all the input files and later run the MD Simulation.

To create a new folder we enter the following command lines:

.. code-block:: text
    
    mkdir openmmdl_amber_tutorial

To copy or move the pdb files for receptor and ligand to the folder we just created:

.. code-block:: text
    
    copy receptor.pdb openmmdl_amber_tutorial
    copy ligand.pdb openmmdl_amber_tutorial
Replace the receptor.pdb and ligand.pdb with the pdb files you want to simulate.

To start the OpenMMDL-Setup we need to activate the openmmdl environment. To do this we have to enter the following command lines:
.. code-block:: text

    conda activate openmmdl

Now that we have activated the openmmdl environment we can launch OpenMMDL-Setup. To do this you need to type the following:

.. code-block:: text

    openmmdl-setup

Selecting the Amber Path
------------------------------
The first step is to select the Amber path. To do this we have to click on the "Amber" button.

Selecting the input files
------------------------------
The second step is to select the input files. There are two options available:

.. figure:: /_static/images/tutorials/Amber_Path/selectAmberFiles.png
   :figwidth: 700px
   :align: center

1. **Yes**, my files are already to simulate.
When this button is selected, the user can select the Prmtop and Inpcrd files from the file browser. Once the files are selected, the user can click on the "Continue" button, which will take the user to the 'Simulation  Setup' page.
2. **No**, I want to prepare them here.
When this button is selected, and the user clicks on the "Continue" button, the user will be taken to the 'Amber Configuration' page. 

Amber Configuration
------------------------------
The user can configure the Amber setup and generate the bash script to generate the Amber files for later-on simulation.

.. figure:: /_static/images/tutorials/Amber_Path/AmberOption.png
   :figwidth: 700px
   :align: center

There are three tabs in this page: 
1. **Receptor**
Depends on the type of macromolecule the user wants to simulate, there are four types are available:

.. figure:: /_static/images/tutorials/Amber_Path/receptor.png
   :figwidth: 700px
   :align: center

- Protein. Click on the "Browse..." button to select the PDB file, and then select the force field from the drop-down menu of "Force Field". If the force filed the user wants is not listed in the menu, the user can click on the 'other', and the 'Other Force Field' textbox will appear. 

.. figure:: /_static/images/tutorials/Amber_Path/receptor_otherFF.png
   :figwidth: 700px
   :align: center

The user can type in the force field name in the textbox. Note: only the force field provided by the AmberTools can be used.``
- DNA;
- RNA;
- Carbohydrate.

The user can select one of the above type of macromolecule at one time.

2. **Ligand**
Depends on the type of ligand the user wants to simulate, there are two types are available:

.. figure:: /_static/images/tutorials/Amber_Path/ligand.png
   :figwidth: 700px
   :align: center

- Normal Ligand. It refers to the small molecule that is made up of C, N, O, S, P, H, F, Cl, Br and I. Here `Antechamber` and `Parmchk2`` will be used to generate the Amber files for the ligand.

.. figure:: /_static/images/tutorials/Amber_Path/normalLigand.png
   :figwidth: 700px
   :align: center

After check the 'Normal Ligand' option, the parameter setting for the ligand will be unfold. 
Firstly, click on the "Browse..." button to select the ligand pdb file. 
And then fill up the charge value of the ligand in the textbox of 'Charge Value'. This value should be a integer, e.g., -1 or 2. and can be calculating easily by openning the ligand pdb file in one text editor, and then suming up the value of the last column of the pdb file.
For 'Ligand Force Field', two force fields are available: General Amber Force Field (GAFF) and GAFF2. The user can select either one of them.
Lastly, pick up 'Charge Method' in the drop-down menu.

- Special Ligand. When it comes to the ligand that `Antechamber` is not able to processed, e.g. cofactor heme in CYP450 enzymes, the user can check the 'Special Ligand' option. 
AMBER parameter database is a good source to find the Amber files for these special ligands. Follow the guidance listed in this application and finally setup for generating the Amber files for the special ligand.

.. figure:: /_static/images/tutorials/Amber_Path/specialLigand.png
   :figwidth: 700px
   :align: center

The user can select either one of or both of the above type of ligand at one time.

3. **Add Water/Membrane**
Depends on the environment the biosystem is in, the user should consider to add water or membrane to the system. 
The user can select either 'Add water Box' or 'Add Membrane and Water' in this tab.
- When 'Add water Box' is selected, the user can further select the 'Box Type' in the drop-down list, and then type the 'Distance(Å)' value in the textbox.

.. figure:: /_static/images/tutorials/Amber_Path/addWater.png
   :figwidth: 700px
   :align: center

- When 'Add Membrane and Water' is selected, the user can further select the 'Lipid Type' and 'Lipid Force Field' in the drop-down list. 
If the listed lipid type is not the one the user wants, the user can click on the 'Other Type or Mixture' option, and then type in the lipid type in the pop-up textbox of 'Other Types or Mixture' and 'Lipid Ratio'. 
For example, 'POPC:TOPC' in 'Other Types or Mixture' and '1:1' in the 'Lipid Ratio' means the membrane is made up of 1 POPC and 1 TOPC. When only one type of lipid is selected, set the 'Lipid Ratio' to 1.

.. figure:: /_static/images/tutorials/Amber_Path/adMembrane.png
   :figwidth: 700px
   :align: center

Warning: The input structure, encompassing both the receptor and ligand, must be aligned with its respective PDB structure available in the OPM database. This alignment is essential for adding the membrane properly using this application.

- 'Water and Ions Setting' is a must for both 'Add water Box' and 'Add Membrane and Water' options. The Ions will be added to neutralize the model. The user can select the 'Water Force Field', 'Positive Ion' and 'Negative Ion' in the drop-down list, and then type the 'Ion Concentration (molar)' value in the textbox.

.. figure:: /_static/images/tutorials/Amber_Path/water_ion_setting.png
   :figwidth: 700px
   :align: center

`tleap` in AmberTools is used to create water boxes around solute. For more information, see AmberTools22 .
`PACKMOL-Memgen` is used to build the membrane in this application. For more information, see the literature.

4. **Save Script**
Last, click on 'Save Script' to download the generated bash script based on the configuration the user set up in this page and save it in the folder the user created at the beginning of this tutorial.
Click 'Continue' to go to the 'Simulation  Setup' page.

Simulation Setup
------------------------------
Firstly, configure the simulation options including the five tabs available in this page: System, Integrator, Simulation, Output, MDAnalysis.
And then click on 'Save Script' to download the generated python script based on the configuration the user set up in this page and save it in the folder the user created at the beginning of this tutorial.

Run Simulation
------------------------------
In the terminal, go to the folder 'openmmdl_amber_tutorial', type the following command lines to run the bash script to generate the Amber files:

.. code-block:: text

    bash run_ambertools.sh

Once the 'Prmtop' and 'Inpcrd' files are generated, the user can run the MD simulation by typing the following command lines:

.. code-block:: text

    python3 OpenMMDL_Simulation.py

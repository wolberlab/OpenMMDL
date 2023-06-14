**OpenMMDL-Simulation Output Files**
=================================

The OpenMMDL Simulation script creates output folders and files during and after the simulation.

The following list contains an detailed overview of the folders and files:


**Analysis**: A folder that contains RMSD Analysis of the MD Simulation performed with MDAnaylsis:


.. list-table::
   :header-rows: 1
   :widths: 25 50

   * - Name
     - Description
   * - RMSD_between_the_frames.png
     - Displays the RMSD between the frames of the protein and ligand.
   * - RMSD_over_time.png
     - Displays the RMSD of the protein, backbone and ligand during the simulation.
   * - RMSD_over_time.csv
     -  Displays the Data used in RMSD_over_time.png.



**Checkpoints**: A folder that contains Checkpoints, which can be used to restart the MD Simulation.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - checkpoint.chk
     - Checkpoint saved every 10000 steps.
   * - 10x_checkpoint.chk
     - Checkpoint saved every 100000 steps.
   * - 100x_checkpoint.chk
     - Checkpoint saved every 1000000 steps.


**Input_Files**: A folder that contains the PDF and SDF files, which served as Input Files for the MD Simulation.



.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - Protein_name.pdb
     - Protein PDB File that served as the input for the MD Simulation.
   * - ligand_name.sdf
     - Ligand SDF File that served as the input for the MD Simulation.



**MD_Files**: A folder that contains the files that were generated during the MD Simulation. This folder contains 3 subfolders.

*Pre_MD*: Folder that contains the files that were prepared by the script before the MD Simulation.

*Minimization_Equilibration*: Folder that contains topology files after the minimization and equilibration of OpenMM.

*MD_Output*: Folder that contains the Output trajectory files generated during the MD Simulation.




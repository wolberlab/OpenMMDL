**OpenMMDL-Simulation Output Files**
=================================

The OpenMMDL Simulation script creates output folders and files during and after the simulation.

The following list contains an detailed overview of the folders and files:
The Optional files are highlighted with an *asteriks**


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
   * - 100x_checkpoint.chk*
     - Checkpoint saved every 1000000 steps.


**Input_Files**: A folder that contains the PDF and SDF files, which served as Input Files for the MD Simulation.



.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - Protein_name.pdb*
     - Protein PDB File that served as the input for the MD Simulation.
   * - ligand_name.sdf*
     - Ligand SDF File that served as the input for the MD Simulation.



**MD_Files**: A folder that contains the files that were generated during the MD Simulation. This folder contains 3 subfolders.

*Pre_MD*: Folder that contains the files that were prepared by the script before the MD Simulation.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - prepared_no_solvent_protein_name.pdb
     - Prepared PDB File without solvent or membrane.
   * - solvent_padding_protein_name.pdb*
     - Prepared PDB File with padding solvent.
   * - solvent_absoulte_protein_name.pdb*
     - Prepared PDB File with absolute solvent.
   * - membrane_protein_name.pdb*
     - Prepared PDB File with membrane.

*Minimization_Equilibration*: Folder that contains topology files after the minimization and equilibration of OpenMM.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - Energyminimization_protein_name.pdb
     - Prepared PDB File after OpenMM energy minimization.
   * - Equilibration_protein_name.pdb
     - Prepared PDB File after OpenMM energy minimization and equilibration.

*MD_Output*: Folder that contains the Output trajectory files generated during the MD Simulation.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - output_protein_name.pdb
     - PDB Topology File of the first frame of the simulation.
   * - trajectory.dcd
     - Trajectory of the OpenMM Simulation.


**MD_Postprocessing**: A folder that contains the postprocessing files after the MD Simulation. This folder contains 2 subfolders.

*1_MDTraj*: Folder that contains the centered topology files with old coordinates.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - centered_old_coordinates_top.pdb
     - Topology File of the centered protein in PDB Format.
   * - centered_old_coordinates_top.gro
     - Topology File of the centered protein in Gromacs GRO Format.
   * - centered_old_coordinates.dcd
     - Trajectory File of the centered protein in DCD Format.
   * - centered_old_coordinates.xtc
     - Trajectory File of the centered protein in XTC Format.
     
     
*2_MDTraj*: Folder that contains the centered topology files with new coordinates according to the center of mass.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - centered_top.pdb
     - Topology File of the centered protein in PDB Format with all atoms and new coordinates.
   * - centered_top.gro
     - Topology File of the centered protein in Gromacs GRO Format with all atoms and new coordinates.
   * - centered_traj.dcd
     - Trajectory File of the centered protein in DCD Format with all atoms and new coordinates.
   * - centered_traj.xtc
     - Trajectory File of the centered protein in XTC Format with all atoms and new coordinates.
   * - prot_lig_top.pdb
     - Topology File of the centered protein in PDB Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_top.gro
     - Topology File of the centered protein in Gromacs GRO Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.dcd
     - Trajectory File of the centered protein in DCD Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.xtc
     - Trajectory File of the centered protein in XTC Format with only protein and ligand atoms and new coordinates.

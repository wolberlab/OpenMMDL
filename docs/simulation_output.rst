**OpenMMDL Simulation Output Files**
=================================

The **OpenMMDL Simulation** script creates output folders and files during and after the simulation.

The following list contains a detailed overview of the folders and files:
The Optional files are highlighted with an *asterisk**

Input Files
------------------------------
**Input Files**: A folder that contains the PDB and SDF files, which served as input files for the MD simulation.



.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - Protein_name.pdb*
     - Protein PDB file that served as the input for the MD simulation.
   * - ligand_name.sdf*
     - Ligand SDF file that served as the input for the MD simulation.

Checkpoints
------------------------------

**Checkpoints**: A folder that contains checkpoints, which can be used to restart the MD Simulation.


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

MD Files
------------------------------
**MD Files**: A folder that contains the files that were generated during the MD simulation. This folder contains 3 subfolders.

*Pre MD*: Folder that contains the files that were prepared by the script before the MD simulation.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - prepared_no_solvent_protein_name.pdb
     - Prepared PDB file without solvent or membrane.
   * - solvent_padding_protein_name.pdb*
     - Prepared PDB file with padding solvent.
   * - solvent_absoulte_protein_name.pdb*
     - Prepared PDB file with absolute solvent.
   * - membrane_protein_name.pdb*
     - Prepared PDB file with membrane.

*Minimization Equilibration*: Folder that contains topology files after the minimization and equilibration with **OpenMM**.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - Energyminimization_protein_name.pdb
     - Prepared PDB file after OpenMM energy minimization.
   * - Equilibration_protein_name.pdb
     - Prepared PDB file after OpenMM energy minimization and equilibration.

*MD Output*: Folder that contains the output trajectory files generated during the MD simulation.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - output_protein_name.pdb
     - PDB topology file of the first frame of the simulation.
   * - trajectory.dcd
     - Trajectory of the OpenMM Simulation.


MD Postprocessing
------------------------------
**MD Postprocessing**: A folder that contains the postprocessing files after the MD simulation.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - centered_old_coordinates_top.pdb
     - Topology file of the centered protein in PDB format.
   * - centered_old_coordinates_top.gro
     - Topology file of the centered protein in Gromacs GRO format.
   * - centered_old_coordinates.dcd
     - Trajectory file of the centered protein in DCD format.
   * - centered_old_coordinates.xtc
     - Trajectory file of the centered protein in XTC format.
   * - centered_traj_unaligned.dcd
     - Trajectory file of the unaligned centered protein in DCD format with all atoms and new coordinates.
   * - centered_traj_unaligned.xtc
     - Trajectory file of the unaligned centered protein in XTC format with all atoms and new coordinates.
   * - prot_lig_traj_unaligned.dcd*
     - Trajectory file of the unaligned centered protein in DCD format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj_unaligned.xtc*
     - Trajectory file of the unaligned centered protein in XTC format with only protein and ligand atoms and new coordinates.
     
Final Output
------------------------------
**Final Output**: A folder that contains the final files after the MD simulation, ready to be analyzed. This folder contains 2 subfolders.

1. *All Atoms*: Folder that contains the centered topology files of all atoms with new coordinates according to the center of mass.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - centered_top.pdb*
     - Topology file of the centered protein in PDB format with all atoms and new coordinates.
   * - centered_top.gro*
     - Topology file of the centered protein in Gromacs GRO format with all atoms and new coordinates.
   * - centered_traj.dcd*
     - Trajectory file of the aligned centered protein in DCD format with all atoms and new coordinates.
   * - centered_traj.xtc*
     - Trajectory file of the aligned centered protein in XTC format with all atoms and new coordinates.



2. *Prot Lig*: Folder that contains the centered topology files of only the protein and ligand atoms with new coordinates according to the center of mass.



.. list-table::
   :header-rows: 1
   :widths: 25 75


   * - Name
     - Description
   * - prot_lig_top.pdb*
     - Topology file of the centered protein in PDB format with only protein and ligand atoms and new coordinates.
   * - prot_lig_top.gro*
     - Topology file of the centered protein in Gromacs GRO format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.dcd*
     - Trajectory file of the centered protein in DCD format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.xtc*
     - Trajectory file of the centered protein in XTC format with only protein and ligand atoms and new coordinates.

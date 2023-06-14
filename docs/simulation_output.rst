**OpenMMDL-Simulation Output Files**
=================================

The OpenMMDL Simulation script creates output folders and files during and after the simulation.

The following list contains an detailed overview of the folders and files:

**Analysis**: A folder that contains RMSD Analysis of the MD Simulation performed with MDAnaylsis.

  1.1 **RMSD_between_the_frames.png**: Displays the RMSD between the frames of the protein and ligand
  
  1.2 **RMSD_over_time.png**: Displays the RMSD of the protein, backbone and ligand during the simulation.
  
  1.3 **RMSD_over_time.csv**: Displays the Data used in RMSD_over_time.png.

**Checkpoints**: A folder that contains Checkpoints, which can be used to restart the MD Simulation.

  1.1 **checkpoint.chk**: Checkpoint saved every 10000 steps.
  
  1.2 **10x_checkpoint.chk**: Checkpoint saved every 100000 steps.
  
  1.2 **100x_checkpoint.chk**: Checkpoint saved every 1000000 steps.

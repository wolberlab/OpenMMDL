**OpenMMDL Analysis Output Files**
=================================

The **OpenMMDL Analysis** script creates output folders and files during and after the analysis.

The following list contains an detailed overview of the folders and files:
The Optional files are highlighted with an *asteriks**

Barcodes
------------------------------
**Barcodes**: A folder that contains the figures of the barcodes for the interactions.



.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - acceptor_interactions.png*
     - File displaying the barcodes for all hydrogen bond acceptor interactions.
   * - donor_interactions.png*
     - File displaying the barcodes for all hydrogen bond donor interactions.
   * - halogen_interactions.png*
     - File displaying the barcodes for all halogen interactions.
   * - hydrophobic_interactions.png*
     - File displaying the barcodes for all hydrophobic interactions.
   * - metal_interactions.png*
     - File displaying the barcodes for all metal interactions.
   * - saltbridge_interactions.png*
     - File displaying the barcodes for all saltbridge interactions.
   * - pication_interactions.png*
     - File displaying the barcodes for all pication interactions.
   * - pistacking_interactions.png*
     - File displaying the barcodes for all pistacking interactions.
   * - waterbridge_interactions.png*
     - File displaying the barcodes for all waterbridge interactions.

Binding_Modes_Markov_States
------------------------------

**Binding_Modes_Markov_States**: A folder that contains the figures of the markov state figures and 2D depiction figure of the top 10 occuring binding modes in addition to PDB and PML files of the representative frames of the binding modes.


.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - all_binding_modes_arranged.png
     - Figure of the 2D depictions of the top 10 binding modes with the interacting protein partners.
   * - Binding_mode_x.pdb
     - PDB File of the representative frame of the binding mode.
   * - Binding_mode_x.pml*
     - PML File of the representative frame of the binding mode.
   * - markov_chain_plot_1.png
     - Markov Chain Figure displaying Binding modes with a minimum transition of 1% of the frames between the binding modes.
   * - markov_chain_plot_2.png
     - Markov Chain Figure displaying Binding modes with a minimum transition of 2% of the frames between the binding modes.
   * - markov_chain_plot_5.png
     - Markov Chain Figure displaying Binding modes with a minimum transition of 5% of the frames between the binding modes.
   * - markov_chain_plot_10.png
     - Markov Chain Figure displaying Binding modes with a minimum transition of 10% of the frames between the binding modes.

RMSD
------------------------------
**RMSD**: A folder that contains the files for the RMSD calculation during the simulation. The RMSD over time is default, while the RMSD between the frames is optional

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - RMSD_over_time.csv
     - CSV File with the values of the RMSD  of the protein/backbone/ligand during the simulation time.
   * - RMSD_over_time.png
     - Figure of the RMSD of the protein/backbone/ligand during the simulation time.
   * - RMSD_between_the_frames.png*
     - Matrix figure displaying the RMSD between each consecutive frame.

Binding Modes Markov States
------------------------------
**Binding Modes Markov States**: A folder that contains the Figures for the binding modes, the pdb files and pml files of the representative frames of the binding modes and the markov chain figure with the transitions between the binding modes.

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
   * - centered_traj_unaligned.dcd
     - Trajectory File of the unaligned centered protein in DCD Format with all atoms and new coordinates.
   * - centered_traj_unaligned.xtc
     - Trajectory File of the unaligned centered protein in XTC Format with all atoms and new coordinates.
   * - prot_lig_traj_unaligned.dcd*
     - Trajectory File of the unaligned centered protein in DCD Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj_unaligned.xtc*
     - Trajectory File of the unaligned centered protein in XTC Format with only protein and ligand atoms and new coordinates.
     
Final Output
------------------------------
**Final Output**: A folder that contains the final files after the MD Simulation, ready to be analyzed. This folder contains 2 subfolders.

1. *All Atoms*: Folder that contains the centered topology files of all atoms with new coordinates according to the center of mass.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - centered_top.pdb*
     - Topology File of the centered protein in PDB Format with all atoms and new coordinates.
   * - centered_top.gro*
     - Topology File of the centered protein in Gromacs GRO Format with all atoms and new coordinates.
   * - centered_traj.dcd*
     - Trajectory File of the aligned centered protein in DCD Format with all atoms and new coordinates.
   * - centered_traj.xtc*
     - Trajectory File of the aligned centered protein in XTC Format with all atoms and new coordinates.



2. *Prot Lig*: Folder that contains the centered topology files of only the protein and ligand atoms with new coordinates according to the center of mass.



.. list-table::
   :header-rows: 1
   :widths: 25 75


   * - Name
     - Description
   * - prot_lig_top.pdb*
     - Topology File of the centered protein in PDB Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_top.gro*
     - Topology File of the centered protein in Gromacs GRO Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.dcd*
     - Trajectory File of the centered protein in DCD Format with only protein and ligand atoms and new coordinates.
   * - prot_lig_traj.xtc*
     - Trajectory File of the centered protein in XTC Format with only protein and ligand atoms and new coordinates.

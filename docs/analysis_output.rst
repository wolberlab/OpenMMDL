**OpenMMDL Analysis Output Files**
=================================

The **OpenMMDL Analysis** script creates output folders and files during and after the analysis.

The following list contains a detailed overview of the folders and files:
The Optional files are highlighted with an *asterisk**


The following Files are generated and stored directly in the working folder where you start the **OpenMMDL Analysis**:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - complex.pdb
     - PDB file of the protein, ligand, waters 10Å around the ligand and special ligand if selected).
   * - df_all.csv
     - CSV file of the calculated interaction together with the generated fingerprints.
   * - interactions_gathered.csv
     - CSV File of the calculated interaction together without the generated fingerprints.
   * - lig.pdb
     - PDB file of the ligand with hydrogens.
   * - lig.smi
     - Smiles file of the ligand.
   * - ligand_numbering.png*
     - File that displays the ligand in 2D depiction with the correct assignment of atom numbers, which makes it easier to understand the interacting atoms.
   * - ligand_special.pdb*
     - PDB file of the ligand with the special ligand.
   * - lig_no_h.pdb
     - PDB file of the ligand without hydrogens.
   * - missing_frames_filled.csv
     - CSV file of the calculated interactions with added frames for frames without interactions
   * - top_10_binding_modes.csv
     - CSV file of the top 10 most occurring binding modes.


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

**Binding_Modes_Markov_States**: A folder that contains the figures of the Markov state figures and a 2D depiction figure of the top 10 occurring binding modes in addition to PDB and PML files of the representative frames of the binding modes.


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
     - Markov chain Figure displaying binding modes with a minimum transition of 1% of the frames between the binding modes.
   * - markov_chain_plot_2.png
     - Markov chain Figure displaying binding modes with a minimum transition of 2% of the frames between the binding modes.
   * - markov_chain_plot_5.png
     - Markov chain Figure displaying binding modes with a minimum transition of 5% of the frames between the binding modes.
   * - markov_chain_plot_10.png
     - Markov chain Figure displaying binding modes with a minimum transition of 10% of the frames between the binding modes.

RMSD
------------------------------
**RMSD**: A folder that contains the files for the RMSD calculation during the simulation. The RMSD over time is always calculated, while the RMSD between the frames is optional.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - RMSD_over_time.csv
     - CSV file with the values of the RMSD  of the protein/backbone/ligand during the simulation time.
   * - RMSD_over_time.png
     - Figure of the RMSD of the protein/backbone/ligand during the simulation time.
   * - RMSD_between_the_frames.png*
     - Matrix figure displaying the RMSD between each consecutive frame.

Visualization
------------------------------
**Visualization**: A folder that contains the figures for the binding modes, the pdb files and pml files of the representative frames of the binding modes and the markov chain figure with the transitions between the binding modes.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Name
     - Description
   * - clouds.json
     - Topology file of the centered protein in PDB format.
   * - interacting_waters.dcd
     - Topology file of the centered protein in Gromacs GRO format.
   * - interacting_waters.pdb
     - Trajectory file of the centered protein in DCD format.
   * - interacting_waters.pkl
     - Trajectory file of the centered protein in XTC format.

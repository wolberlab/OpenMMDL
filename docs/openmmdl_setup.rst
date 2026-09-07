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

Interactive tutorial mode
-------------------------

OpenMMDL Setup includes guided tutorial mode. Click the ``Tutorials`` button in the Setup header, or open the tutorial section on the first page, to reach the tutorial selection page. Currently the available workflow is ``PDB small molecule``. Starting it opens an in-app tutorial panel that highlights the relevant controls, explains the current step, links back to this documentation, and can apply the tutorial choices for you: ``AMBER19`` with its ``OPC3`` water model, ``Single complex`` with a second ligand copy as an additional molecule, ``SMIRNOFF`` for the small molecules, pH ``7.4``, a water box with ionic strength ``0.15 M``, and a ``10 ns`` simulation length. The :doc:`PDB Path </pdb_path>` page describes the same pages with their default settings and all available options, and points out where the interactive tutorial deviates from the defaults. How the tutorial panel works and which tutorials exist is described on the :doc:`Interactive Tutorials </tutorials>` page.

The browser cannot select local files automatically. Use the ``Download files`` button in Setup to retrieve the bundled example files: the MOE-processed 5WYZ receptor (``5wyz_tutorial_receptor.pdb``, with the crystal N-glycans kept so the chain step has sugar chains to deselect) and the two copies of the CU-CPT9b ligand (``7VF``) taken from the crystal structure as separate files, ``7VF_A.sdf`` and ``7VF_B.sdf``. When the tutorial asks for them, upload the receptor, then ``7VF_A.sdf`` as the ligand and ``7VF_B.sdf`` as the additional molecule.


There are two possible options to create the input files for **OpenMMDL Simulation**:

1. The PDBFixer path, where a `pdb` file of the protein is used as an input for the preparation and simulation.
The PDBFixer path is described page by page :doc:`here </pdb_path>`.

Here is the table of the currently available forcefields and water models for the PDBFixer path: 

.. list-table:: PDBFixer path AMBER forcefields
   :header-rows: 1
   :widths: 26 34 34

   * - Water model
     - AMBER14 / AMBER19
     - AMBER99SB / AMBER99SB-ILDN / AMBER03 / AMBER10
   * - TIP3P
     - ✓
     - ✓
   * - TIP3P-FB
     - ✓
     - ✓
   * - SPC/E
     - ✓
     - ✓
   * - TIP4P-Ew
     - ✓
     - ✓
   * - TIP4P-FB
     - ✓
     - —
   * - TIP5P
     - —
     - ✓
   * - OPC3
     - ✓
     - ✓
   * - OPC
     - ✓
     - ✓


.. list-table:: PDBFixer path CHARMM forcefields
   :header-rows: 1
   :widths: 34 34 34

   * - Water model
     - CHARMM36
     - CHARMM36 2024
   * - SPC/E
     - ✓
     - ✓
   * - TIP4P-Ew
     - ✓
     - ✓
   * - TIP5P
     - ✓
     - ✓
   * - CHARMM default
     - ✓
     - ✓
   * - TIP3P-PME-B
     - ✓
     - ✓
   * - TIP3P-PME-F
     - ✓
     - ✓
   * - TIP4P-2005
     - ✓
     - ✓
   * - TIP5P-Ew
     - ✓
     - ✓

2. The Amber path, where `prmtop` and `inpcrd` files are used the preparation and simulation. This path allows us to either use already prepared `prmtop` and `inpcrd` as an input or create the `prmtop` and `inpcrd` from PDB files of the receptor and ligand.
The Amber path is described page by page :doc:`here </amber_path>`.

.. figure:: /_static/images/amber_ff.png
   :figwidth: 725px
   :align: center

In the table, the first row is the default setting, and the term `other` allows users to type their desired forcefields from those accessible in AmberTools 22.0 into the designated textbox.

**Glycoprotein receptors (Amber path)**

The Amber path additionally offers a *Glycoprotein* receptor type for proteins with covalently attached N- or O-linked glycans. When selected, OpenMMDL Setup:

1. Detects sugar residues by PDB residue name (NAG, NDG, BMA, MAN, FUC, FUL) and finds glycosidic linkages from CONECT records (falling back to a distance scan if absent).
2. Renames residues and atoms to GLYCAM-06j-1 conventions (e.g. NAG → 0YB / 4YB, MAN → 0MA / 3MA / 6MA, with the trimannose branch point as VMB).
3. Emits explicit ``bond`` statements for every protein-glycan and glycan-glycan covalent link, injected into the generated ``tleap.in`` after ``loadpdb``.
4. Sources both the chosen protein force field (e.g. ``leaprc.protein.ff14SB``) and the glycan force field (``leaprc.GLYCAM_06j-1`` by default).

Only the common N-glycan core sugars listed above are currently supported. Unrecognised PDB sugar codes or non-canonical linkage patterns raise a descriptive error; see ``openmmdl/openmmdl_setup/glycoprotein.py`` for the verified GLYCAM tables.

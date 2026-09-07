**PDB Path**
==============

This page walks through the PDB path of **OpenMMDL Setup** page by page and explains what each page shows and
which options can be changed. The example system is the human Toll-like receptor 8 (TLR8) in complex with the
antagonist CU-CPT9b, PDB entry `5WYZ <https://www.rcsb.org/structure/5WYZ>`_. TLR8 is a dimer and binds one ligand
per protomer, which is used here to show how a second molecule is added to a system.

The screenshots show the pages with their default settings. The interactive tutorial inside **OpenMMDL Setup**
(``Tutorials`` in the Setup header, then ``PDB small molecule``) uses the same files and pages, but deliberately
changes a few of the options so that their effect becomes visible: it switches the protein force field to
``AMBER19``, the small-molecule force field to ``SMIRNOFF``, the pH to ``7.4`` and the simulation length to ``10 ns``.
Whenever such a change is made in the interactive tutorial, the corresponding option is pointed out below.

Example files
------------------------------

The example files ship with OpenMMDL. Download them with the ``Download files`` button on the Setup
tutorial page, or copy them from ``openmmdl/openmmdl_simulation/tutorial_systems/pdb_path/5wyz_solvent`` in the repository:

* ``5wyz_tutorial_receptor.pdb``: the TLR8 ectodomain dimer (chains A and B) from 5WYZ, prepared in MOE.
  Hydrogens were added and the unresolved loops of the crystal structure were capped with ACE and NME groups
  instead of being rebuilt. The N-linked glycans of the crystal structure are still present, so the chain
  selection page has something to remove.
* ``7VF_A.sdf`` and ``7VF_B.sdf``: the two copies of CU-CPT9b (PDB residue ``7VF``), one bound to chain A and one to
  chain B, with the crystal coordinates. They contain only heavy atoms; **OpenMMDL Simulation** adds the hydrogens.

Starting OpenMMDL-Setup
------------------------------

We start by creating a folder where we copy the example files and later run the MD simulation.

To create a new folder we enter the following command lines:


.. code-block:: text

    mkdir openmmdl_pdb_tutorial


To start the OpenMMDL-Setup we need to activate the openmmdl environment. To do this we have to enter the following command lines:

.. code-block:: text

    conda activate openmmdl

Now that we have activated the openmmdl environment we can start OpenMMDL Setup. To do this you need to type the following:

.. code-block:: text

    openmmdl setup


Selecting Input Files
------------------------------

This will open the OpenMMDL-Setup, which you can use for the creation of the input files for OpenMMDL Simulation.



.. figure:: /_static/images/tutorials/PDB_Path/Selectfiletype.png
   :figwidth: 725px
   :align: center


The first page offers the two preparation paths. Here we select the ``PDB / mmCIF`` option and press ``Continue``.
The Amber path is described on its own page, :doc:`Amber Path </amber_path>`.


This leads us to the input page, where the receptor and the ligand files are selected and the way the system is modelled is chosen.


.. figure:: /_static/images/tutorials/PDB_Path/Selectfiles1.png
   :figwidth: 725px
   :align: center


**PDB File**: press ``Browse...`` (or drag the file onto the field) and select ``5wyz_tutorial_receptor.pdb``.

**Force Field**: the protein force field. The AMBER family (``AMBER14``, ``AMBER19``, ``AMBER99SB``, ``AMBER99SB-ILDN``,
``AMBER03``, ``AMBER10``) supports ligands; the CHARMM options (``CHARMM36``, ``CHARMM36 2024``) do not. ``AMBER19``
corresponds to the ff19SB protein force field, which was developed together with the OPC family of water models;
selecting it switches the **Water Model** to ``OPC3`` automatically. The interactive tutorial makes this switch to show it.

**Water Model**: the water models available for the chosen force field, see the tables in :doc:`Running OpenMMDL Setup </openmmdl_setup>`.

**Small Molecule Mode**: ``No small molecules`` prepares the receptor alone. ``Single complex`` reveals the ``Ligand File``
and ``Topology Code`` fields for one ligand. ``High-throughput / Multiple poses`` creates one simulation setup per entry
of a multi-molecule SDF file.

**Ligand File**: press ``Browse...`` and select ``7VF_A.sdf``. The **Topology Code** is the three-character residue name the
ligand receives in the generated topology and is used later during postprocessing and analysis. The default is ``UNK``.

**Additional Molecules**: ``+ Add additional molecule`` adds a row for a further SDF, MOL or MOL2 file. Select ``7VF_B.sdf``
here. Every additional molecule is parametrised on its own and included in the generated system, which is how cofactors,
a second ligand copy or other small molecules enter the simulation. Each molecule needs a unique topology code; the
defaults are ``L01``, ``L02`` and so on.

**Small Molecule Force Field**: ``GAFF`` is the General Amber Force Field, which assigns parameters through atom types, and
offers a choice of GAFF versions. ``SMIRNOFF`` is the format of the Open Force Field Initiative, which assigns parameters
by direct chemical perception with SMARTS patterns; it offers a ``Constrained`` or ``Unconstrained`` mode (the constrained
variants are intended for simulations that constrain bonds involving hydrogen) and a choice of OpenFF versions. Both are
provided through openmmforcefields. The interactive tutorial switches to ``SMIRNOFF`` to show these fields.

**Ligand Minimization** minimises the ligand geometry with MMFF94s before the simulation, and **Ligand Sanitization**
lets RDKit sanitise the ligand, which some ligands require and others do not tolerate. Both are off by default.

With the example files and the default options the page looks like this:


.. figure:: /_static/images/tutorials/PDB_Path/Selectfiles2.png
   :figwidth: 725px
   :align: center

``Continue`` becomes available once a receptor file is selected and all topology codes are valid and unique.

Selecting Chains
------------------------------

The next page lists the chains found in the PDB file next to a viewer showing the structure.


.. figure:: /_static/images/tutorials/PDB_Path/Modifychains1.png
   :figwidth: 725px
   :align: center


The first two rows are the protein chains A and B. The remaining rows are the N-linked glycans of the crystal structure:
the branched glycans as chains C to H and, at the end of the table, the single NAG residues that are listed under the chain
letters A and B. Their content column shows the sugar residue names (``NAG``, ``BMA``, ``MAN``) instead of ``Protein``.
Unticking a row removes that chain from the system.

For this system only the chains A and B are needed, so every glycan row is deselected. The AMBER protein force fields of
the PDB path carry no parameters for the sugars, and they are not needed for the ligand binding site. Simulations with
glycans are possible in the Amber path with the ``Glycoprotein`` receptor type described in :doc:`Running OpenMMDL Setup </openmmdl_setup>`.

The heterogen option below the table decides what happens to residues that are neither amino acids nor nucleotides.
Keep it at ``Keep all heterogens``: the ACE and NME groups that cap the chain breaks of the prepared receptor count as
heterogens, and the other two options would delete them.


.. figure:: /_static/images/tutorials/PDB_Path/Modifychains2.png
   :figwidth: 725px
   :align: center

With only the chains A and B selected we press ``Continue``.

Optional: Adding Residues
------------------------------

On this page PDBFixer lists sequence stretches that are present in the sequence records of the PDB file but missing
from the atom coordinates, and offers to rebuild them.

The example receptor was prepared in MOE and its chain breaks are capped, so PDBFixer finds nothing to add and
OpenMMDL Setup skips this page. If you start from the unprocessed entry ``5wyz.pdb`` from the RCSB PDB instead, the page
looks like this:

.. figure:: /_static/images/tutorials/PDB_Path/Addresidues1.png
   :figwidth: 725px
   :align: center

In that case the missing residues 101 to 111 and 435/438 to 460 of the chains A and B would be selected and rebuilt.

.. figure:: /_static/images/tutorials/PDB_Path/Addresidues2.png
   :figwidth: 725px
   :align: center

Optional: Adding Heavy Atoms
------------------------------

This page lists residues with missing heavy atoms, which PDBFixer adds automatically. It is skipped for the example
receptor, because its residues are complete. With an unprocessed PDB file it looks like this:

.. figure:: /_static/images/tutorials/PDB_Path/Addheavyatoms.png
   :figwidth: 725px
   :align: center

Adding Hydrogens and Water Box
------------------------------

The next page adds the missing hydrogens and, optionally, a water box or a membrane.

.. figure:: /_static/images/tutorials/PDB_Path/Addwater1.png
   :figwidth: 725px
   :align: center

**Add hydrogens**: PDBFixer assigns the protonation states of the titratable residues according to the pH entered in the
field. The default is ``7.0``; the interactive tutorial changes it to ``7.4`` for physiological conditions.

**Add water box**: the box can be given as explicit dimensions or as a padding distance. With ``Specify a padding distance``
the box is built around a sphere that contains the whole model plus the padding (default ``1.0`` nm), in the chosen
**Box Shape** (cube, truncated octahedron or rhombic dodecahedron).

**Add membrane and water** instead embeds the protein in a lipid bilayer of the chosen lipid type; the protein has to be
oriented for that, for example with a structure from the OPM database.

**Ionic strength**: ions are first added to neutralise the system and then up to the requested concentration, ``0.15``
molar with ``Na+`` and ``Cl-`` by default, which mimics physiological conditions.

.. figure:: /_static/images/tutorials/PDB_Path/Addwater2.png
   :figwidth: 725px
   :align: center

Pressing ``Continue`` takes a moment, because OpenMMDL adds the hydrogens, the water box and the ions before the next
page opens.

Simulation Options
------------------------------

The last page collects the simulation settings in tabs on the left and shows the generated simulation script on the
right. Every change updates the script immediately.

.. figure:: /_static/images/tutorials/PDB_Path/Simulationoptions.png
   :figwidth: 725px
   :align: center

* **System**: the **Simulation Length** (default ``50`` ns; the interactive tutorial sets ``10`` ns to keep the run short),
  the minimisation and equilibration protocol, the OpenMM platform (``CUDA`` for NVIDIA GPUs, ``OpenCL``, ``CPU`` or
  ``Reference``) and the numerical precision.
* **Output**: trajectory, log, checkpoint, XML and final-state files and their intervals.
* **Postprocessing**: centring of the complex and the topology and trajectory formats written after the simulation.
* **Analysis**: whether **OpenMMDL Analysis** runs after the simulation and with which settings.
* **Simulation** and **Integrator**: nonbonded method and cutoff, constraints, hydrogen mass repartitioning, time step,
  ensemble, temperature, friction and pressure.

The buttons above the script save the result. ``Save Script`` downloads only ``OpenMMDL_Simulation.py``, and
``Save All Files`` downloads ``openmmdl_simulation.zip``
with the script, the prepared receptor ``5wyz_tutorial_receptor-processed_openMMDL.pdb`` (with hydrogens, water box and
ions added) and the two ligand files ``7VF_A.sdf`` and ``7VF_B.sdf``.

Running OpenMMDL-Simulation
------------------------------

Unpack ``openmmdl_simulation.zip``. The folder ``openmmdl_simulation`` contains the four files listed above; the ligand
SDF files are the same files that were used as input for OpenMMDL Setup.

Now that we have the files in one folder we can start the MD simulation.

For this we start by activating the environment

.. code-block:: text

    conda activate openmmdl

Now that we activated the environment we start the simulation from inside the ``openmmdl_simulation`` folder.

For this enter the following command

.. code-block:: text

    openmmdl simulation -f tutorial_simulation -s OpenMMDL_Simulation.py -t 5wyz_tutorial_receptor-processed_openMMDL.pdb -l 7VF_A.sdf 7VF_B.sdf

By entering the command we create a folder called ``tutorial_simulation``, where the output of the MD simulation will appear.

As the input for the MD simulation we used ``-t`` to select the prepared receptor as the topology file, ``-l`` to pass both
ligand files in the same order as in OpenMMDL Setup (``7VF_A.sdf`` with the topology code ``UNK`` first, then ``7VF_B.sdf``
with ``L01``), and ``-s`` to specify the ``OpenMMDL_Simulation.py`` script that will run the MD simulation.

.. figure:: /_static/images/tutorials/PDB_Path/Outputfiles1.png
   :figwidth: 725px
   :align: center

During and after the simulation you can open the folder to see the progress.

After the simulation is finished the ``tutorial_simulation`` folder should look like the picture below.

.. figure:: /_static/images/tutorials/PDB_Path/Outputfiles2.png
   :figwidth: 725px
   :align: center

If there are files or folders missing, repeat the MD simulation.

To see what the separate files in the output represent follow this page:

* :doc:`MD Simulation Output </simulation_output>`

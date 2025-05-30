Frequently asked questions
==========================

.. contents::
   :local:

..


OpenMMDL
------------------------------------

What OS does OpenMMDL support?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**OpenMMDL** is available for Linux and Mac systems. **OpenMMDL Setup** and **OpenMMDL Simulation** are currently incompatible with Windows systems, thus only **OpenMMDL Analysis** of the **OpenMMDL** package can be used on Windows systems.

How can I install OpenMMDL?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**OpenMMDL** is available as a conda package. To install it, you can use the following command:
.. code-block:: text
   
   conda install -c conda-forge openmmdl

How do I cite OpenMMDL?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Please refer to the :doc:`citation` page.

OpenMMDL Setup
------------------------------------

Do I need to prepare the protein before using OpenMMDL Setup?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**OpenMMDL Setup** can prepare the protein for you in the PDBFixer path and add the missing residues in addition to protonating the protein, but it is not able to cap the residues and the preferred input for OpenMMDL Setup is a protein, which is already prepared by a user with having the missing residues modelled and the termini capped.

I have a special ligand in my system, what Path should I use in OpenMMDL Setup?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the case of a special ligand in the systems use the Amber Path and prepare the files there.

Can I use the CHARMM forcefield in the PDBFixer path for the simulations with my ligand?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
No, currently only the Amber force fields are compatible with protein-ligand simulations.

OpenMMDL Analysis
------------------------------------

How can I ensure that the residue numbering in the output files matches the original structure file for the convenience of post-MD analysis?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The '-ref' flag can be employed to align the residue numbers in the output topology file of the MD simulation with a reference protein.

Why can't I get the PDB or PML file for the most occurring poses in the 'Binding_Modes_Markov_States' folder ?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The cause could be many. One possibility, as revealed by our internal testing, could be a problematic SDF file of the ligand. By default, **OpenMMDL Analysis** uses `obabel` to directly extract the ligand from the topology file of the MD trajectory and subsequently generate the ligand's SDF file. However, `obabel` may struggle to correctly process certain small molecules, such as charged molecules, as observed in our internal testing.

A viable solution to this problem is to prepare the ligand SDF file using third party molecular tools such as MOE. Once prepared, the '-l' flag can be added to explicitly specify the correct SDF file for **OpenMMDL Analysis**.
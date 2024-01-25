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

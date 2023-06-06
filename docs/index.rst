.. openmmdl documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

OpenMMDL  @ Molecular Design Lab
=========================================================

.. image::
   https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT

.. image:: 
   https://github.com/talagayev/OpenMMDL/blob/main/openmmdl_setup/static/OpenMMDL_logo.svg
   :width: 350

**OpenMMDL** is an workflow, that allows an easy setup of **OpenMM** molecular dynamic simulations of protein-ligand complexes. **OpenMMDL** consists of:

*OpenMMDL_Setup* - A webserver tool that allows to modify the Protein and generate Input Files for MD Simulations.

*OpenMMDL_Simulation* - A script that performs a OpenMM MD Simulation with your Input files and further postprocesses the MD Simulation with MDTraj and MDAnalysis.



.. toctree::
   :maxdepth: 1
   :caption: User guide:

   installation

.. toctree::
   :maxdepth: 1
   :caption: OpenMMDL-Setup
   
   openmmdl_setup

.. toctree::
   :maxdepth: 1
   :caption: OpenMMDL-Simulation
   
   openmmdl_simulation

.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   
   tutorial_pdb_path
   tutorial_amber_path


.. toctree::
   :maxdepth: 1
   :caption: Developers:

   api




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

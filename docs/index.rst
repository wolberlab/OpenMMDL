.. openmmdl documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. figure:: /_static/images/OpenMMDL_logo_2.png
   :figwidth: 600px
   :align: center


**OpenMMDL** is a workflow that enables the setup, simulation, and analysis of
**OpenMM** molecular dynamics simulations of protein-ligand complexes.
**OpenMMDL** consists of:

*OpenMMDL Setup* - A webserver tool that allows users to modify the protein and
generate input files for MD simulations.

*OpenMMDL Simulation* - Runs an OpenMM MD simulation using the prepared input
files and postprocesses the trajectory with MDTraj and MDAnalysis.

*OpenMMDL Analysis* - Analyzes protein-ligand interactions in OpenMM MD
simulation trajectories.



.. toctree::
   :maxdepth: 1
   :caption: User guide:

   installation
   releases
   citation
   faq

.. toctree::
   :maxdepth: 1
   :caption: OpenMMDL-Setup
   
   openmmdl_setup

.. toctree::
   :maxdepth: 4
   :caption: OpenMMDL-Simulation
   
   openmmdl_simulation
   simulation_output

.. toctree::
   :maxdepth: 1
   :caption: OpenMMDL-Analysis
   
   openmmdl_analysis
   analysis_output

.. toctree::
   :maxdepth: 1
   :caption: Tutorials
   
   tutorial_pdb_path
   tutorial_amber_path


.. toctree::
   :maxdepth: 1
   :caption: API Docs
   
   api


Introduction
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

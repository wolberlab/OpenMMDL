.. openmmdl documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. figure:: /_static/images/OpenMMDL_logo_2.png
   :figwidth: 600px
   :align: center


**OpenMMDL** is an workflow, that allows an easy setup of **OpenMM** molecular dynamic simulations of protein-ligand complexes. **OpenMMDL** consists of:

*OpenMMDL Setup* - A webserver tool that allows to modify the Protein and generate Input Files for MD Simulations.

*OpenMMDL Simulation* - Performs a OpenMM MD Simulation with your Input files and further postprocesses the MD Simulation with MDTraj and MDAnalysis.

*OpenMMDL Analysis* - Analyzes the protein-ligand interactions in the OpenMM MD Simulation trajectories.



.. toctree::
   :maxdepth: 1
   :caption: User guide:

   installation
   changelog
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
   
   modules/openmmdl_analysis/analysis/markovchains
   modules/openmmdl_analysis/core/trajectories


Introduction
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Release Notes
====================

This page summarizes the published OpenMMDL releases and
release notes. For the full change history, see the
`CHANGELOG <https://github.com/wolberlab/OpenMMDL/blob/main/CHANGELOG.md>`_.

Release 1.3.0
------------------------------

This release focuses on simulation robustness, workflow improvements, and
postprocessing updates across the OpenMMDL pipeline.

**Major features & changes**

* Added an equilibration protocol to ``OpenMMDL Simulation``.
* Added support for multiple ligands in ``OpenMMDL Simulation``.
* Added restart handling for NaN-coordinate failures in simulation runs.
* Added a high-throughput simulation option in ``OpenMMDL Setup``.
* Added ligand ``mol2`` export with partial charges.
* Improved cleanup and postprocessing behavior.
* Added additional accepted boolean flag values in ``OpenMMDL Analysis``.

Release 1.2.0
------------------------------

This release focused on improvements across the workflow, including the CLI,
interaction analysis, simulation postprocessing, and broader force-field
support.

**Major features & changes**

* Added ``ProLIF`` as an option for protein-ligand interaction calculation.
* Added a command-line interface entry point for ``openmmdl``.
* Updated the ``OpenMMDL Setup`` layout.
* Added version selection for both ``GAFF`` and ``SMIRNOFF`` force fields.
* Added new force fields and water models, including ``CHARMM2024``.
* Fixed several simulation postprocessing issues.

Release 1.1.1
------------------------------

This was a short release focused on extending force-field and water-model
support.

**Major features & changes**

* Added the ``Amber19`` force field.
* Added ``OPC`` and ``OPC3`` water models.
* Updated the ``OpenFF`` force-field version and fixed related ``GAFF`` type
  issues.

Release 1.1.0
------------------------------

This release added major analysis and visualization functionality and included
substantial code refactoring.

**Major features & changes**

* Added support for selecting the final frame in ``OpenMMDL Analysis``.
* Added ``MDonatello`` visualization support.
* Added ``SMIRNOFF`` support for small-molecule force fields.
* Added ``PyMOL`` support for visualization.
* Added a ``Dockerfile`` and Docker installation path.
* Added the citation page in the documentation.
* Refactored ``OpenMMDL Analysis`` into clearer analysis, core, and
+  visualization sections.

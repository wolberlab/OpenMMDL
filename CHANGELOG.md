# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

<!--
The rules for this file:
  * entries are sorted newest-first.
  * summarize sets of changes - don't reproduce every git log comment here.
  * don't ever delete anything.
  * keep the format consistent:
    * do not use tabs but use spaces for formatting
    * 79 char width
    * YYYY-MM-DD date format (following ISO 8601)
  * accompany each entry with github issue/PR number (Issue #xyz)
-->

## [Unreleased]

### Authors
talagayev, NDoering99

### Added
- Changed to the use of RDKitConverter for Ligand recognition without -l flag in openmmdl_analysis (2025-04-10)
- Addition of MDontallo Visualization (2025-04-11)
- Addition of SMIRNOFF small molecule force field (2025-04-10, Issue #76)
- Addition of PyMOL support for visualization (2025-04-10)
- Added Dockerfile for image creation  (2025-04-10)
- Addition of citation page in documentation (2025-04-10)

### Fixed
- Fixed wrong color of the positive and negative ionizable interactions (2025-04-21, PR #137)
- Fixed the XMLSerializer error appearing during simulation (2025-04-10, Issue #122)
- Fixed the sanitization implementation in OpenMMDL Setup (2025-04-09)

### Changed
<!-- Changes in existing functionality -->

### Deprecated
<!-- Soon-to-be removed features -->

### Removed
<!-- Removed features -->

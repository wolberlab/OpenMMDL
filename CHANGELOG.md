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

## Version 1.1.0

### Authors
talagayev, NDoering99

### Added
- Addition of option to select final frame in OpenMMDL Analysis (2025-04-23, Issue #136, PR #140)
- Addition of `MDontallo` Visualization (2025-04-11)
- Addition of `SMIRNOFF` small molecule force field (2025-04-10, Issue #76)
- Addition of `PyMOL` support for visualization (2025-04-10)
- Added `Dockerfile` for image creation  (2025-04-10)
- Addition of citation page in documentation (2025-04-10)

### Fixed
- Fixed error cases, where `ligand name` was not specified (2025-04-22, Issue #138, PR #139)
- Fixed wrong color of the positive and negative ionizable interactions (2025-04-21, PR #137)
- Fixed the `XMLSerializer` error appearing during simulation (2025-04-10, Issue #122)
- Fixed the `sanitization` implementation in OpenMMDL Setup (2025-04-09)

### Changed
- Classes `ImageMerger` and `FigureArranger` were moved to `image_handler.py` (2025-04-30, Issue #141, PR #142)
- The class `TrajectorySaver` was moved to `trajectory_saving.py` (2025-04-24, Issue #141, PR #142)
- Changed to the use of `RDKitConverter` for Ligand recognition
  without `-l` flag in `OpenMMDL Analysis` (2025-04-10)

### Deprecated
<!-- Soon-to-be removed features -->

### Removed
<!-- Removed features -->

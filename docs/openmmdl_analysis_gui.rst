**OpenMMDL Analysis GUI**
=========================

The **OpenMMDL Analysis GUI** is a browser interface for **OpenMMDL Analysis**. It covers both
sides of the analysis:

* **Run a new analysis**: pick the topology, trajectory and ligand, set the options and start
  ``openmmdl analysis`` with a live log.
* **Open analysis results**: explore the folder of a finished analysis (one started from the GUI
  or from the command line) interactively.

Start it with:

.. code-block:: text

    openmmdl analysis-gui

The GUI runs on your own computer (``http://127.0.0.1:5050``) and opens in your browser. Nothing
is uploaded: files are read where they are, so large trajectories are no problem.

.. code-block:: text

    openmmdl analysis-gui path/to/results   # open a results folder right away
    openmmdl analysis-gui --port 5051        # another port
    openmmdl analysis-gui --no-browser       # do not open a browser tab

Run a new analysis
------------------

The form asks only for what every analysis needs: topology, trajectory, ligand SDF, ligand residue
name and the output folder. All other options of ``openmmdl analysis`` are under *Advanced
options* (their defaults are those of the command line; hover the ⓘ next to a field for its help
text). The command that will run can be shown and copied next to the form.

* **Drag and drop**: drop the topology, trajectory and SDF onto the drop area (or onto a single
  field); they go to the right fields by file extension. Dropped files are copied into
  ``<output folder>/inputs`` (a folder under ``~/openmmdl_analysis`` is chosen when no output folder
  is set yet). *Browse* picks files where they are, without copying.
* **Simulation folder**: *Use a simulation folder…* fills in the topology, trajectory and ligand SDF
  from the ``Final_Output`` (or ``MD_Postprocessing``) folder of ``openmmdl simulation``.
* **Ligand**: when a PDB topology is selected, its residues are inspected and the ligand residue
  name is suggested.
* *Work on a copy of the topology* (under *Advanced options*) is on by default, because the
  analysis rewrites its topology file in place.
* On Windows, keep the output folder path short (or enable long paths): the analysis writes files
  up to about 150 characters deep inside it, and Windows limits paths to 260 characters.

While the analysis runs, the GUI shows the overall progress with an estimate of the time left, what
is running right now (for example frames done and time per frame), and every step of the run with
the time it took: the steps follow the chosen options, and each binding-mode pass (by residue, by
ligand atom) lists its own steps. The full log is shown next to it; Python warnings of the
libraries are hidden unless *Show Python warnings* is ticked. A run can be cancelled. When it
finishes, *Open results* opens the explorer.

Explore results
---------------

Open a results folder with *Open results*. The bar on top switches between the two partner groupings
written by the analysis (*Residue*: ``BindingModes_residue``; *Ligand atom*:
``BindingModes_ligand``) and between three layouts:

* **Docked 3D**: the analysis tabs with the 3D view docked on the right (the divider can be dragged).
* **Full width**: each tab uses the full width; the 3D view is its own tab and appears next to the
  2D view in *Binding modes*.
* **Report**: all sections on one page, for reading or printing.

Tabs:

* **Overview**: key numbers, the most frequent interactions, the settings of the run and warnings
  from ``openmmdl.log``.
* **Interactions**: the barcode of every interaction over the frames. Filter by type, search, sort
  and set a minimal occurrence. Click an interaction to show it in 3D, or click a barcode to jump
  to that frame.
* **Binding modes**: the binding modes are recomputed live from ``df_all.csv``; the threshold
  slider has the same meaning as ``openmmdl analysis -b``. For the selected mode:

  * the ligand in 2D, either as a depiction with coloured halos on the interacting atoms or as the
    ProLIF ligand network;
  * chips to switch single interactions, or all interactions of one type, on and off. Clicking an
    atom in the depiction, or a residue or edge in the network, does the same. The 2D views and
    the 3D view always show the same interactions;
  * the representative frame (the frame whose contacts are closest to the mode average) and where
    the mode occurs along the trajectory.

* **Markov chain**: the transitions between binding modes as an interactive graph, recomputed with
  the rules of the analysis (minimal transition share as ``-m``, node colours by the third of the
  trajectory a mode mostly occurs in). Hover a node or edge for numbers, click a node to select
  the mode, drag nodes, scroll to zoom and drag the background to pan. The strip below shows the
  binding mode of every frame.
* **RMSD**: RMSD over time as an interactive chart, plus the frame-to-frame RMSD figure if it was
  computed.
* **Waters**: water bridges over time and, if ``-w`` was used, the stable water clusters, which can
  be shown in 3D.
* **Files**: every file of the folder; figures are shown as a gallery. *All figures (.zip)*,
  *Everything (.zip)* and *Folder (.zip)* download them as one archive.

Exporting
---------

Every figure of the GUI can be downloaded as it is shown, with the current threshold, selection
and zoom:

* **Markov chain** and **RMSD**: SVG or PNG (RMSD also as CSV).
* **Binding modes**: the 2D depiction as SVG or PNG, the ProLIF network as PNG, and all binding
  modes of the current threshold as CSV (frames, representative frame, interactions).
* **Interactions**: the barcode chart as PNG and the interactions with their frames as CSV.
* **3D view**: *PNG* saves the current view at twice the screen resolution.

In the Markov chain, *Structures* draws every binding mode as the ligand with that mode's
interactions highlighted (the same colours as in *Binding modes*); hovering a node shows a larger
version.

The 3D view plays the trajectory from ``Visualization/interacting_waters.dcd`` frame by frame,
draws the interactions of the current frame as dashed lines, and can show the interaction point
clouds of ``clouds.json`` and the bridging waters.

.. note::

    The binding modes in the GUI use every interaction above the threshold. The figures and the
    ``top_10_binding_modes.csv`` written by ``openmmdl analysis`` leave out one interaction column
    when grouping frames, so for some systems their mode numbering differs; the GUI then says so
    above the list of modes.

For developers
--------------

The page is written in TypeScript in ``frontend/analysis_gui`` and built into
``openmmdl/openmmdl_analysis_gui/static``. The built files are committed, so users do not need
Node.js. After changing the frontend:

.. code-block:: text

    cd frontend/analysis_gui
    npm ci
    npm run build       # bundles the page and copies vis-network next to it
    npm test            # binding-mode and Markov derivation tests
    npm run typecheck

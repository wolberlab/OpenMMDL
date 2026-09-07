**Interactive Tutorials**
=========================

**OpenMMDL Setup** contains guided tutorials that run directly inside the Setup interface. A tutorial adds a panel to
every page of the workflow that explains the current step, highlights the control to use, and can apply the expected
choice for you. Nothing is simulated during a tutorial; it ends with the same simulation package that a normal Setup
run produces, ready for **OpenMMDL Simulation**.

Starting a tutorial
------------------------------

Start **OpenMMDL Setup** as usual:

.. code-block:: text

    conda activate openmmdl
    openmmdl setup

Click ``Tutorials`` in the Setup header, or use the ``Open tutorials`` button on the first page. The tutorial page lists the
available workflows.

.. figure:: /_static/images/tutorials/Tutorials_page.png
   :figwidth: 725px
   :align: center

Each card offers two buttons:

* ``Download files`` downloads a zip archive with the example files of that workflow and a ``README.txt`` that says which
  file goes where. Browsers cannot select local files on their own, so the files have to be unpacked once and uploaded
  by hand when the panel asks for them.
* ``Start tutorial`` opens the first page of the workflow with the tutorial panel attached.

Working with the panel
------------------------------

The panel sits in the lower right corner of every page and shows the current page and step, the progress on that page,
and a short explanation of what the step does.

.. figure:: /_static/images/tutorials/Tutorial_panel.png
   :figwidth: 725px
   :align: center

* The control that belongs to the step is highlighted on the page and carries a label with the expected input.
* Steps with an action button, for example ``Set AMBER19`` or ``Add water box``, apply the expected choice when
  clicked. You can also make the change yourself; the panel shows a green ``Done`` mark once the page matches the
  step, whichever way the change was made.
* ``Back`` and ``Next`` move between the steps of a page, and ``Restart page`` returns to the first step of the
  current page. Moving to the next page is done with the page's own ``Continue`` button, as in a normal Setup run.
* ``Open docs`` opens the matching section of this documentation, and ``Download files`` is offered again on the
  pages where files are uploaded.
* ``Exit`` in the panel, or ``Exit Tutorial`` in the header, ends the tutorial and returns to the first page. The
  tutorial stays active across ``Back``, ``Start Over`` and page reloads until it is ended this way. On the last step
  of the last page the button reads ``Finish tutorial``.

Pages that have nothing to do for the example system, such as the missing-residue page for an already prepared
receptor, are skipped by Setup exactly as in a normal run, so the page counter of the panel may jump.

The tutorial state is kept in the browser session of the running Setup. Opening Setup in a second browser, or
restarting ``openmmdl setup``, starts without an active tutorial; start it again from the ``Tutorials`` page in that case.

Available tutorials
------------------------------

**PDB small molecule**

Prepares the human Toll-like receptor 8 dimer from PDB entry 5WYZ with both copies of its antagonist CU-CPT9b through
the PDB path. The example files are the MOE-processed receptor ``5wyz_tutorial_receptor.pdb``, which still contains the
N-glycans of the crystal structure, and the two ligand files ``7VF_A.sdf`` and ``7VF_B.sdf``.

The tutorial follows the pages described in :doc:`PDB Path </pdb_path>`, but changes a few options on purpose so that
their effect becomes visible:

* the protein force field is switched to ``AMBER19``, which makes Setup select the ``OPC3`` water model,
* the second ligand copy is added as an additional molecule with its own topology code,
* the small-molecule force field is switched from ``GAFF`` to ``SMIRNOFF``,
* the glycan chains are deselected on the chain page while all heterogens are kept,
* hydrogens are added at pH ``7.4`` with a water box and ``0.15 M`` ions,
* the simulation length is set to ``10 ns`` instead of the default ``50 ns``.

The resulting package is run with:

.. code-block:: text

    openmmdl simulation -f tutorial_simulation -s OpenMMDL_Simulation.py -t 5wyz_tutorial_receptor-processed_openMMDL.pdb -l 7VF_A.sdf 7VF_B.sdf

Further tutorials, for example for the Amber path, can be added to the same tutorial page.

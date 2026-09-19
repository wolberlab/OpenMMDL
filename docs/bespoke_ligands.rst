**Bespoke Ligand Parameters**
=============================

Introduction
------------------

General small molecule force fields (GAFF, OpenFF Sage) assign torsion parameters
by chemical environment. For ligands whose rotatable bonds are poorly covered,
`OpenFF BespokeFit <https://github.com/openforcefield/openff-bespokefit>`_ fits
molecule-specific torsion parameters against quantum chemistry torsion scans and
writes them as a SMIRNOFF force field (``.offxml``) that contains the base force
field plus the bespoke torsions.

**OpenMMDL Setup** wraps this step for both paths: the ligand of a system can either
be fitted with BespokeFit when the downloaded script runs, or an existing ``.offxml``
(for example the ``final.offxml`` of an earlier fit) can be imported. Everything
else in the setup, the simulation and the analysis stays unchanged.

Requirements
------------------

The OpenMMDL environment needs ``openff-toolkit``, ``openff-interchange`` and
``parmed``, which are part of ``environment.yml``. The Amber path also needs
AmberTools (``tleap``) as usual.

Fitting needs BespokeFit and a quantum chemistry program. Both are best installed
in a separate conda environment on Linux or WSL::

    conda create -n bespokefit -c conda-forge openff-bespokefit psi4

Give the path to ``openff-bespoke`` of that environment (e.g.
``~/miniconda3/envs/bespokefit/bin/openff-bespoke``) as **BespokeFit Executable**
in the setup, or leave the field empty when it is on the ``PATH``. Importing an
``.offxml`` needs no BespokeFit installation.

Ligand file
------------------

Provide the ligand as a single-molecule SDF (or MOL) file with explicit hydrogens
in the intended protonation state and the bound pose. The fit is done on this
molecule; bespoke torsions match its exact chemistry, so an imported ``.offxml``
must have been fitted for the same molecule and protonation state.

Setup interface
------------------

The **Ligand Parameters** selection is shown below the ligand file in the PDB path
(single complex mode) and under **Normal Ligand** in the Amber path:

* **Standard small molecule force field**: the GAFF or SMIRNOFF selection as before.
* **Fit bespoke torsions with OpenFF BespokeFit**: choose the executable, the QC level
  (B3LYP-D3BJ/DZVP with Psi4 is BespokeFit's default; GFN2-xTB with xtb is a fast
  choice for testing) and the QC resources (workers, cores and memory per worker).
* **Import a SMIRNOFF force field (.offxml)**: upload the file.

Invalid input (e.g. a ligand file without explicit hydrogens, or a file that is not
a SMIRNOFF force field) is reported on the page. In the Amber path, the generated
bash script shows the message instead of the script until the input is corrected.

PDB path
------------------

The downloaded ``OpenMMDL_Simulation.py`` fits or imports the parameters into the
directory ``ligand_parameters`` before the force field is built:

* ``ligand_parameters/final.offxml`` is the force field used for the ligand,
* ``ligand_parameters/fit.log`` shows the progress of BespokeFit (a fit takes hours),
* ``ligand_parameters/manifest.json`` records the inputs. A directory with the same
  inputs is reused, so re-running the script does not fit again.

Additional molecules keep the selected standard force field. Run the script as
usual with **OpenMMDL Simulation**; the imported ``.offxml`` or a finished
``ligand_parameters`` directory is passed with ``-p``::

    openmmdl simulation -f run -t protein.pdb -s OpenMMDL_Simulation.py -l ligand.sdf -p final.offxml

Amber path
------------------

The generated ``run_ambertools.sh`` starts with ``openmmdl ligand fit`` or
``openmmdl ligand import`` with ``--target amber``. This writes the fitted force field
and a tleap library into ``bespoke_ligand/`` (``ligand.lib``, ``ligand.frcmod``,
``ligand.pdb``, plus the exact ``ligand.prmtop``). The ligand gets the residue code
``BSP``, and antechamber/parmchk2 are not run for it. Ligand charges come from the
force field's charge model (AM1-BCC for OpenFF Sage).

tleap keeps only one of the three improper torsion terms that SMIRNOFF force fields
use per planar center. After tleap has built the solvated system, the script
therefore runs ``openmmdl ligand restore``, which replaces the ligand residue in the
final ``prmtop``/``inpcrd`` with the exact parameters. Run the bash script, then the
simulation script, as in the normal Amber path.

Command line
------------------

The same steps are available without the setup interface::

    openmmdl ligand fit --ligand ligand.sdf --output ligand_parameters \
        --executable ~/miniconda3/envs/bespokefit/bin/openff-bespoke \
        --qc-spec psi4 b3lyp-d3bj dzvp --workers 1 --cores 4 --memory-gb 8

    openmmdl ligand import --ligand ligand.sdf --offxml final.offxml --output ligand_parameters

    openmmdl ligand fit --ligand ligand.sdf --output bespoke_ligand --target amber
    openmmdl ligand restore --prmtop system.prmtop --inpcrd system.inpcrd --ligand-dir bespoke_ligand

Limitations
------------------

* One bespoke ligand per system; additional molecules use the standard force field.
* The pose library (high-throughput) mode of the PDB path is not supported, since a
  fit belongs to one molecule.
* Covalently bound ligands are not supported in the Amber path.
* BespokeFit fits torsions only; charges and all other parameters come from its base
  force field (OpenFF Sage 2.2.0). Virtual sites are not supported.

"""Put the exact bespoke ligand parameters back into a topology built by tleap.

tleap keeps only one of the three SMIRNOFF improper torsions per center when it
reads the exported library. After the solvated system has been built, the ligand
residue is therefore replaced by the ligand.prmtop written by `export_amber`.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

from openmmdl.ligand_parameters.core import RESNAME


def restore_ligand_parameters(prmtop, inpcrd, ligand_dir, resname=RESNAME):
    """Replace the ligand residue in prmtop/inpcrd (in place) by ligand_dir/ligand.prmtop."""
    import numpy as np
    import parmed as pmd

    prmtop, inpcrd, ligand_dir = Path(prmtop), Path(inpcrd), Path(ligand_dir)
    system = pmd.load_file(str(prmtop), xyz=str(inpcrd))
    ligand = pmd.load_file(str(ligand_dir / "ligand.prmtop"), xyz=str(ligand_dir / "ligand.inpcrd"))

    residues = [residue for residue in system.residues if residue.name == resname]
    if len(residues) != 1:
        raise ValueError(f"Expected exactly one {resname} residue in {prmtop.name}, found {len(residues)}")
    indices = [atom.idx for atom in residues[0].atoms]
    start, stop = min(indices), max(indices) + 1
    if sorted(indices) != list(range(start, stop)):
        raise ValueError("The ligand residue must be a contiguous block of atoms")
    for bond in system.bonds:
        if (bond.atom1.idx < start or bond.atom1.idx >= stop) != (bond.atom2.idx < start or bond.atom2.idx >= stop):
            raise ValueError("Covalently bound ligands are not supported")
    by_name = {atom.name: atom for atom in residues[0].atoms}
    if by_name.keys() != {atom.name for atom in ligand.atoms} or len(by_name) != len(ligand.atoms):
        raise ValueError("The ligand atoms in the topology do not match the fitted ligand")

    ligand.coordinates = np.array([system.coordinates[by_name[atom.name].idx] for atom in ligand.atoms])
    ligand.residues[0].name = resname
    # Slicing keeps every parameter of the receptor, membrane, solvent and ions.
    environment = system.copy(pmd.Structure)
    combined = ligand.copy(pmd.Structure)
    if start:
        combined = environment[:start] + combined
    if stop < len(system.atoms):
        combined += environment[stop:]
    combined.box = system.box

    # Structure.copy can leave water extra points (OPC, TIP4P) with stale child references
    # that break the molecule ordering of AmberParm; drop those.
    def prune_children(structure):
        for atom in structure.atoms:
            atom.children[:] = [
                child for child in atom.children
                if 0 <= child.idx < len(structure.atoms) and structure.atoms[child.idx] is child
            ]

    prune_children(combined)
    restored = pmd.amber.AmberParm.from_structure(combined)
    prune_children(restored)

    from openmm import app

    with tempfile.TemporaryDirectory(dir=prmtop.parent) as directory:
        new_prmtop, new_inpcrd = Path(directory) / prmtop.name, Path(directory) / inpcrd.name
        restored.write_parm(str(new_prmtop))
        restored.save(str(new_inpcrd), overwrite=True)
        # Make sure OpenMM accepts the result before replacing the tleap files.
        check = app.AmberPrmtopFile(str(new_prmtop)).createSystem(nonbondedMethod=app.NoCutoff)
        if check.getNumParticles() != len(system.atoms) or len(app.AmberInpcrdFile(str(new_inpcrd)).positions) != len(system.atoms):
            raise ValueError("Restoring the ligand parameters changed the number of particles")
        os.replace(new_prmtop, prmtop)
        os.replace(new_inpcrd, inpcrd)

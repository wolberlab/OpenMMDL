import os
import numpy as np
import mdtraj as md

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis import transformations as trans
from MDAnalysis.lib import distances


def mdtraj_conversion(pdb_file, mdtraj_output):
    """
    Center on the protein COM (per frame) then image molecules.
    For multichain proteins, anchor ALL protein chains so they are imaged into the same periodic image.
    """
    traj = md.load_dcd("trajectory.dcd", top=pdb_file)
    top = traj.topology

    # Build anchor_molecules for ALL protein molecules (chains are separate molecules)
    protein_atoms = set(top.select("protein"))
    anchor_molecules = []
    for mol in top.find_molecules():
        # mol is iterable of mdtraj Atom objects
        if any(a.index in protein_atoms for a in mol):
            anchor_molecules.append(mol)

    def image_with_anchors():
        if anchor_molecules:
            traj.image_molecules(inplace=True, anchor_molecules=anchor_molecules, make_whole=True)
        else:
            traj.image_molecules(inplace=True)

    # If unit cell info is missing, fall back to imaging only.
    if traj.unitcell_lengths is None or np.any(np.isnan(traj.unitcell_lengths)):
        image_with_anchors()
    else:
        # Pre-image once so protein isn't split across PBC before COM calculation
        image_with_anchors()

        # Center on protein COM (per frame)
        prot_idx = top.select("protein")
        if prot_idx.size > 0:
            prot = traj.atom_slice(prot_idx)
            com = md.compute_center_of_mass(prot)   # (n_frames, 3)
            box = traj.unitcell_lengths             # (n_frames, 3)
            shift = (box * 0.5) - com               # move protein COM to box center
            traj.xyz = traj.xyz + shift[:, None, :]

        # Image molecules into the primary unit cell (anchored for multichain)
        image_with_anchors()

    # Write trajectory
    if "dcd" in mdtraj_output:
        traj.save_dcd("centered_old_coordinates.dcd")
    if "xtc" in mdtraj_output:
        traj.save_xtc("centered_old_coordinates.xtc")

    # Write topology snapshot (frame 0)
    first = traj[0:1]
    if "pdb" in mdtraj_output:
        first.save_pdb("centered_old_coordinates_top.pdb")
    if "gro" in mdtraj_output:
        first.save_gro("centered_old_coordinates_top.gro")


def MDanalysis_conversion(
    post_mdtraj_pdb_file,
    post_mdtraj_dcd_file,
    mda_output,
    output_selection,
    ligand_name=None,
    special_ligname=None,
):
    """
    PBC-safe: unwrap protein, reimage ligand (min-image) to protein, center complex, wrap.
    Multichain-safe: reimage protein segments (chains) together + wrap protein as segments.
    After alignment: rewrap + reimage again to prevent PBC artifacts.
    """

    def _safe_select(u, sel):
        try:
            return u.select_atoms(sel)
        except Exception:
            return u.atoms[0:0]

    def _ligand_group(u):
        parts = []
        if ligand_name:
            parts.append(f"resname {ligand_name}")
        if special_ligname:
            parts.append(f"resname {special_ligname}")
        if not parts:
            return u.atoms[0:0]
        return _safe_select(u, " or ".join(parts))

    def _lipid_group(u):
        # Optional keyword; if unsupported, returns empty.
        return _safe_select(u, "lipid")

    def _has_pbc(ts):
        dims = getattr(ts, "dimensions", None)
        if dims is None or len(dims) < 3:
            return False
        dims = np.asarray(dims[:3], dtype=float)
        return np.all(np.isfinite(dims)) and np.all(dims > 0.0)

    def _can_unwrap(ag):
        """unwrap() needs fragments -> fragments needs bonds."""
        if ag.n_atoms == 0:
            return False
        try:
            _ = ag.fragments  # triggers NoDataError if no bonds
            return True
        except Exception:
            return False

    def reimage_to_reference(mobile_ag, ref_ag):
        """
        Translate mobile_ag so its COM is minimum-image to ref_ag COM.
        """
        def _transform(ts):
            if mobile_ag.n_atoms == 0 or ref_ag.n_atoms == 0 or not _has_pbc(ts):
                return ts

            ref_com = ref_ag.center_of_mass()
            mob_com = mobile_ag.center_of_mass()
            delta = mob_com - ref_com
            delta_min = distances.minimize_vectors(delta.reshape(1, 3), ts.dimensions)[0]
            mobile_ag.translate(delta_min - delta)
            return ts

        return _transform

    def reimage_residues_to_reference(residue_ag, ref_ag):
        """
        Translate each residue so its COM is minimum-image to ref_ag COM.
        Useful for lipids to keep membrane around the protein image.
        """
        def _transform(ts):
            if residue_ag.n_atoms == 0 or ref_ag.n_atoms == 0 or not _has_pbc(ts):
                return ts

            ref_com = ref_ag.center_of_mass()
            for res in residue_ag.residues:
                mob_com = res.atoms.center_of_mass()
                delta = mob_com - ref_com
                delta_min = distances.minimize_vectors(delta.reshape(1, 3), ts.dimensions)[0]
                res.atoms.translate(delta_min - delta)
            return ts

        return _transform

    def reimage_segments_to_reference(protein_ag):
        """
        Keep all protein segments (chains) in the same periodic image.
        This is what prevents dimers/oligomers from splitting across PBC.
        """
        segs = list(protein_ag.segments)
        if len(segs) <= 1:
            def _noop(ts):
                return ts
            return _noop

        ref = segs[0].atoms

        def _transform(ts):
            if not _has_pbc(ts):
                return ts
            ref_com = ref.center_of_mass()
            for seg in segs[1:]:
                mob = seg.atoms
                delta = mob.center_of_mass() - ref_com
                delta_min = distances.minimize_vectors(delta.reshape(1, 3), ts.dimensions)[0]
                mob.translate(delta_min - delta)
            return ts

        return _transform

    def _apply_pre_alignment_transforms(u):
        prot = _safe_select(u, "protein")
        lig = _ligand_group(u)
        lip = _lipid_group(u)

        # Anchor for centering: protein + ligand + lipids (if present)
        anchor = prot if prot.n_atoms else u.atoms
        if lig.n_atoms:
            anchor = anchor + lig
        if lip.n_atoms:
            anchor = anchor + lip

        protein = prot
        non_protein = _safe_select(u, "not protein")

        transforms = []

        # Make protein whole, keep chains together
        if protein.n_atoms:
            if _can_unwrap(protein):
                transforms.append(trans.unwrap(protein))
            transforms.append(reimage_segments_to_reference(protein))

        # Keep ligand in same protein image
        if lig.n_atoms and protein.n_atoms:
            if _can_unwrap(lig):
                transforms.append(trans.unwrap(lig))
            transforms.append(reimage_to_reference(lig, protein))

        # Keep lipids around protein image (optional)
        if lip.n_atoms and protein.n_atoms:
            transforms.append(reimage_residues_to_reference(lip, protein))

        # Center then wrap: protein by segments; everything else by residues
        transforms.append(trans.center_in_box(anchor, center="mass", wrap=False))

        if protein.n_atoms:
            transforms.append(trans.wrap(protein, compound="segments"))
        transforms.append(trans.wrap(non_protein, compound="residues"))

        u.trajectory.add_transformations(*transforms)

    def _rewrap_after_alignment(top_file, aligned_traj_file, out_traj_file, write_group_sel):
        """
        Re-open aligned trajectory, then:
        - keep chains together (segments)
        - reimage ligand to protein
        - keep lipids around protein
        - wrap protein by segments, others by residues
        """
        u2 = mda.Universe(top_file, aligned_traj_file)

        prot2 = _safe_select(u2, "protein")
        lig2 = _ligand_group(u2)
        lip2 = _lipid_group(u2)

        protein2 = prot2
        non_protein2 = _safe_select(u2, "not protein")

        transforms = []

        if protein2.n_atoms:
            transforms.append(reimage_segments_to_reference(protein2))

        if lig2.n_atoms and protein2.n_atoms:
            transforms.append(reimage_to_reference(lig2, protein2))

        if lip2.n_atoms and protein2.n_atoms:
            transforms.append(reimage_residues_to_reference(lip2, protein2))

        if protein2.n_atoms:
            transforms.append(trans.wrap(protein2, compound="segments"))
        transforms.append(trans.wrap(non_protein2, compound="residues"))

        u2.trajectory.add_transformations(*transforms)

        write_ag = _safe_select(u2, write_group_sel)
        with mda.Writer(out_traj_file, write_ag.n_atoms) as w:
            for _ts in u2.trajectory:
                w.write(write_ag)

    # ---- Main logic ----
    u = mda.Universe(post_mdtraj_pdb_file, post_mdtraj_dcd_file)

    # Apply per-frame transforms BEFORE writing unaligned files
    _apply_pre_alignment_transforms(u)

    # Always write topology from frame 0
    u.trajectory[0]

    # -------------------- PDB/DCD outputs --------------------
    if "pdb" in mda_output:
        # All-atoms
        if output_selection != "mda_prot_lig":
            all_atoms = u.atoms
            all_atoms.write("centered_top.pdb")

            with mda.Writer("centered_traj_unaligned.dcd", all_atoms.n_atoms) as w:
                for _ts in u.trajectory:
                    w.write(all_atoms)

            ref = mda.Universe("centered_top.pdb")
            mob = mda.Universe("centered_top.pdb", "centered_traj_unaligned.dcd")
            tmp = "centered_traj_aligned_tmp.dcd"

            align.AlignTraj(
                mob, ref,
                select="protein and name CA",
                weights="mass",
                filename=tmp,
            ).run()

            _rewrap_after_alignment("centered_top.pdb", tmp, "centered_traj.dcd", "all")
            os.remove(tmp)

        # Protein+ligand
        if output_selection != "mda_all":
            prot_lig_sel = "protein"
            if ligand_name:
                prot_lig_sel += f" or resname {ligand_name}"
            if special_ligname:
                prot_lig_sel += f" or resname {special_ligname}"

            prot_lig = _safe_select(u, prot_lig_sel)
            prot_lig.write("prot_lig_top.pdb")

            with mda.Writer("prot_lig_traj_unaligned.dcd", prot_lig.n_atoms) as w:
                for _ts in u.trajectory:
                    w.write(prot_lig)

            ref = mda.Universe("prot_lig_top.pdb")
            mob = mda.Universe("prot_lig_top.pdb", "prot_lig_traj_unaligned.dcd")
            tmp = "prot_lig_traj_aligned_tmp.dcd"

            align.AlignTraj(
                mob, ref,
                select="protein and name CA",
                weights="mass",
                filename=tmp,
            ).run()

            _rewrap_after_alignment("prot_lig_top.pdb", tmp, "prot_lig_traj.dcd", prot_lig_sel)
            os.remove(tmp)

    # -------------------- GRO/XTC outputs --------------------
    if "gro" in mda_output:
        # All-atoms
        if output_selection != "mda_prot_lig":
            all_atoms = u.atoms
            all_atoms.write("centered_top.gro")

            with mda.Writer("centered_traj_unaligned.xtc", all_atoms.n_atoms) as w:
                for _ts in u.trajectory:
                    w.write(all_atoms)

            ref = mda.Universe("centered_top.gro")
            mob = mda.Universe("centered_top.gro", "centered_traj_unaligned.xtc")
            tmp = "centered_traj_aligned_tmp.xtc"

            align.AlignTraj(
                mob, ref,
                select="protein and name CA",
                weights="mass",
                filename=tmp,
            ).run()

            _rewrap_after_alignment("centered_top.gro", tmp, "centered_traj.xtc", "all")
            os.remove(tmp)

        # Protein+ligand
        if output_selection != "mda_all":
            prot_lig_sel = "protein"
            if ligand_name:
                prot_lig_sel += f" or resname {ligand_name}"
            if special_ligname:
                prot_lig_sel += f" or resname {special_ligname}"

            prot_lig = _safe_select(u, prot_lig_sel)
            prot_lig.write("prot_lig_top.gro")

            with mda.Writer("prot_lig_traj_unaligned.xtc", prot_lig.n_atoms) as w:
                for _ts in u.trajectory:
                    w.write(prot_lig)

            ref = mda.Universe("prot_lig_top.gro")
            mob = mda.Universe("prot_lig_top.gro", "prot_lig_traj_unaligned.xtc")
            tmp = "prot_lig_traj_aligned_tmp.xtc"

            align.AlignTraj(
                mob, ref,
                select="protein and name CA",
                weights="mass",
                filename=tmp,
            ).run()

            _rewrap_after_alignment("prot_lig_top.gro", tmp, "prot_lig_traj.xtc", prot_lig_sel)
            os.remove(tmp)

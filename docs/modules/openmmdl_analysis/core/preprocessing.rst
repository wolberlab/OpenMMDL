API Documentation for Preprocessing
============================


.. py:class:: Preprocessing()

    A class providing utilities for preparing and modifying PDB and trajectory files for molecular dynamics analysis.

    .. py:method:: renumber_protein_residues(input_pdb, reference_pdb, output_pdb)

        Renumber protein residues in a trajectory PDB file based on a reference PDB structure.

        :param str input_pdb: Path to the input PDB file to be renumbered.
        :param str reference_pdb: Path to the reference PDB file for correct residue numbering.
        :param str output_pdb: Path to the output PDB file with updated residue numbering.

    .. py:method:: increase_ring_indices(ring, lig_index)

        Adjusts ligand ring atom indices to align with full protein-ligand complex indices.

        :param list ring: List of ligand atom indices in the ring.
        :param int lig_index: Index offset from the full structure.
        :returns: List of adjusted atom indices.
        :rtype: list

    .. py:method:: process_pdb_file(input_pdb_filename)

        Modifies residue names in a PDB file to standard forms (e.g., water residues renamed to "HOH").

        :param str input_pdb_filename: Path to the PDB file to be processed.

    .. py:method:: extract_and_save_ligand_as_sdf(input_pdb_filename, output_filename, target_resname)

        Extracts a ligand from a receptor-ligand complex PDB and saves it in SDF format.

        :param str input_pdb_filename: Path to the complex PDB file.
        :param str output_filename: Path for the output SDF file.
        :param str target_resname: Residue name of the ligand.

    .. py:method:: renumber_atoms_in_residues(input_pdb_file, output_pdb_file, lig_name)

        Renames ligand atom names based on element and count within the residue for clarity.

        :param str input_pdb_file: Path to the original PDB file.
        :param str output_pdb_file: Path to save the updated PDB file.
        :param str lig_name: Name of the ligand residue.

    .. py:method:: replace_atom_type(data)

        Corrects atom type annotations in ligand lines marked with 'X'.

        :param str data: Contents of the PDB file as a string.
        :returns: Modified PDB file contents.
        :rtype: str

    .. py:method:: process_pdb(input_file, output_file)

        Wrapper function that processes and writes corrected PDB content to a file.

        :param str input_file: Path to the input PDB file.
        :param str output_file: Path to the output PDB file.

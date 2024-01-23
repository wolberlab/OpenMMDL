OpenMMDL simulation functions
=============================

This page displays all the functions of **OpenMMDL Simulation**.

openmmdl_simulation.scripts.cleaning_procedures
------------------------------

.. py:function:: cleanup(protein_name)
    
    Cleans up the PDB Reporter Output File and MDTraj Files of the performed simulation.
    
    :param str protein_name: Name of the protein PDB.

    :returns: None.
    :rtype: None
   

.. py:function:: create_directory_if_not_exists(directory_path)
    
    Create a directory if it doesn't exist, or overwrite it if already does.
    
    :param str directory_path: Path of the directory that you want to create.

    :returns: None.
    :rtype: None


.. py:function:: copy_file(src, dest)
    
    Copy a file to the destination path.
    
    :param str src: Path of the file that needs to be copied.
    :param str dest: Path of destination where the file needs to be copied to.

    :returns: None.
    :rtype: None


.. py:function:: organize_files(source, destination)
    
    Organizes the files and moves them from the source to the destination directory.
    
    :param str source: Path of the file that needs to be moved.
    :param str destination: Path of destination where the file needs to be moved to.

    :returns: None.
    :rtype: None




.. py:function:: post_md_file_movement(protein_name: str, prmtop: str = None, inpcrd: str = None, ligands: List[str] = None)
    
    Organizes and moves the files after the MD simulation to their respective directories.
    
    :param str protein_name: Name of the protein PDB.
    :param prmtop: Path to the AMBER topology file.
    :param inpcrd: Path to the AMBER coordinate file.
    :param ligands: List of paths to the ligand files.
    :type prmtop: Optional [str]
    :type inpcrd: Optional [str]
    :type ligands: Optional [List[str]]


    :returns: None.
    :rtype: None


openmmdl_simulation.scripts.forcefield_water
------------------------------

.. py:function:: ff_selection(ff)
    
    Selects the required XML forcefield file.
    
    :param str ff: Input forcefield.

    :returns: Selected XML forcefield file.
    :rtype: str


.. py:function:: water_forcefield_selection(water, forcefield_selection)
    
    Selects the required XML forcefield file.
    
    :param str water: The chosen water model.
    :param str forcefield_selection: The selected force field.

    :returns: The XML filename of the water forcefield.
    :rtype: str


.. py:function:: water_model_selection(water, forcefield_selection)
    
    Selects the required water model forcefield XML file according to water selection and previous force field selection.
    
    :param str water: Water model input.
    :param str forcefield_selection: Input of selected forcefield XML file.

    :returns: Water model forcefield XML file.
    :rtype: str


.. py:function:: generate_forcefield(protein_ff, solvent_ff, add_membrane, rdkit_mol=None)
    
    Generate an OpenMM Forcefield object and register a small molecule.
    
    :param str protein_ff: Input of selected forcefield XML File.
    :param str solvent_ff: Input of selected water model forcefield XML File.
    :param bool add_membrane: Selection if the system should be built with a membrane.
    :param rdkit.Chem.rdchem.Mol rdkit_mol: Small molecule to register in the force field.

    :returns: Forcefield with a registered small molecule.
    :rtype: simtk.openmm.app.Forcefield


.. py:function:: generate_transitional_forcefield(protein_ff, solvent_ff, add_membrane, rdkit_mol=None)
    
    Generate an OpenMM transitional forcefield object with TIP3P water model for membrane building and register a small molecule.
    
    :param str protein_ff: Name of the force field in XML format.
    :param str solvent_ff: Name of the water model force field in XML format.
    :param bool add_membrane: Selection if the system should be built with a membrane.
    :param rdkit.Chem.rdchem.Mol rdkit_mol: Small molecule to register in the force field.

    :returns: A transitional forcefield with TIP3P water and a registered small molecule.
    :rtype: simtk.openmm.app.Forcefield


openmmdl_simulation.scripts.post_md_conversions
------------------------------

.. py:function:: mdtraj_conversion(pdb_file, mdtraj_output)
    
    Recenter and apply periodic boundary conditions to the molecules in each frame of the trajectory, and save the centered trajectory and its first frame.
    
    :param str pdb_file: Name of the PDB file. This PDB file stores the extracted frames from the MD trajectory.
    :param str mdtraj_output: The selected format that will be used as an output of the topology and trajectory.

    :returns: None.
    :rtype: None


.. py:function:: MDanalysis_conversion(post_mdtraj_pdb_file, post_mdtraj_dcd_file, mda_output, output_selection, ligand_name=None, special_ligname=None)
    
    Translate the trajectory so that all frames coincide with its center of geometry.
    
    :param str post_mdtraj_pdb_file: Name of the post-MDtraj PDB file.
    :param str post_mdtraj_dcd_file: Name of the post-MDtraj DCD File.
    :param str ligand_name: Ligand name saved in the PDB file.
    :param str special_ligname: Special residue name saved in the PDB file.
    :param str mda_output: Selection of output formats.
    :param str output_selection: Selection of topologies with specific atom selections that will be created.

    :returns: None.
    :rtype: None

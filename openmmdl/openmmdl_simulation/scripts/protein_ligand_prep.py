import mdtraj as md
import numpy as np
import simtk.openmm.app as app
from rdkit import Chem
from openff.toolkit.topology import Molecule
from simtk.openmm.app import PDBFile
from simtk.openmm import unit
from simtk.openmm import Vec3


class MoleculePreparation:
    """
    A class for preparing and converting ligand molecules using RDKit and OpenMM.

    This class handles reading ligand files in different formats, adding hydrogens to the ligand structure,
    and optionally minimizing the molecule's energy. It also provides functionality to convert the prepared
    RDKit molecule into an OpenMM-compatible format and to merge the ligand with a protein structure.

    Attributes:
        ligand_file (str): Path to the ligand file (supports SDF, MOL, and MOL2 formats).
        minimize_molecule (bool): Flag to determine whether to minimize the ligand's energy using MMFF.
        ligand (rdkit.Chem.rdchem.Mol): The prepared RDKit molecule object.
    """

    def __init__(self, ligand_file, minimize_molecule=True):
        """
        Initialize the MoleculePreparation instance with the ligand file and minimization flag.
        """
        self.ligand_file = ligand_file
        self.minimize_molecule = minimize_molecule
        self.ligand = None

    def prepare_ligand(self):
        """Reads an SDF File into RDKit, adds hydrogens to the structure, minimizes it if selected, and creates an openforcefield Molecule object.
        Inspired by @teachopencadd T019.

        Returns:
            rdkitmolh (rdkit.Chem.rdchem.Mol): The prepared and converted ligand.
        """
        file_name = self.ligand_file.lower()
        if file_name.endswith(".sdf"):
            rdkit_mol = Chem.SDMolSupplier(self.ligand_file, sanitize=False)
            for mol in rdkit_mol:
                rdkit_mol = mol
        elif file_name.endswith(".mol") and not file_name.endswith(".mol2"):
            print(self.ligand_file)
            rdkit_mol = Chem.rdmolfiles.MolFromMolFile(self.ligand_file, sanitize=False)
        elif file_name.endswith(".mol2"):
            rdkit_mol = Chem.rdmolfiles.MolFromMol2File(
                self.ligand_file, sanitize=False
            )

        print("Adding hydrogens")
        rdkitmolh = Chem.AddHs(rdkit_mol, addCoords=True)
        Chem.AssignAtomChiralTagsFromStructure(rdkitmolh)

        if self.minimize_molecule:
            Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(
                mol=rdkitmolh, mmffVariant="MMFF94s", maxIters=2000
            )

        Molecule(rdkitmolh)
        self.ligand = rdkitmolh
        return rdkitmolh

    @staticmethod
    def rdkit_to_openmm(rdkit_mol, name):
        """Convert an RDKit molecule to an OpenMM molecule.
        Inspired by @teachopencadd T019, @hannahbrucemcdonald and @glass-w.

        Args:
            rdkit_mol (rdkit.Chem.rdchem.Mol): RDKit molecule to convert.
            name (str): Molecule name.

        Returns:
            omm_molecule (simtk.openmm.app.Modeller): OpenMM modeller object holding the molecule of interest.
        """
        off_mol = Molecule.from_rdkit(rdkit_mol)
        off_mol.name = name

        element_counter_dict = {}
        for off_atom, rdkit_atom in zip(off_mol.atoms, rdkit_mol.GetAtoms()):
            element = rdkit_atom.GetSymbol()
            if element in element_counter_dict.keys():
                element_counter_dict[element] += 1
            else:
                element_counter_dict[element] = 1
            off_atom.name = element + str(element_counter_dict[element])

        off_mol_topology = off_mol.to_topology()
        mol_topology = off_mol_topology.to_openmm()
        new_mol_positions = []

        for mol_position in off_mol.conformers[0]:
            new_mol_positions.append(mol_position.magnitude / 10.0)

        omm_mol = app.Modeller(mol_topology, new_mol_positions * unit.nanometers)
        return omm_mol

    @staticmethod
    def merge_protein_and_ligand(protein, ligand):
        """Merge two OpenMM objects. Inspired by @teachopencadd T019.

        Args:
            protein (pdbfixer.pdbfixer.PDBFixer): Protein to merge.
            ligand (simtk.openmm.app.Modeller): Ligand to merge.

        Returns:
            complex_topology (simtk.openmm.app.topology.Topology): The combined topology of the protein and ligand.
            complex_positions (simtk.unit.quantity.Quantity): The combined positions of the protein and ligand.
        """
        md_protein_topology = md.Topology.from_openmm(protein.topology)
        md_ligand_topology = md.Topology.from_openmm(ligand.topology)
        md_complex_topology = md_protein_topology.join(md_ligand_topology)

        complex_topology = md_complex_topology.to_openmm()

        total_atoms = len(protein.positions) + len(ligand.positions)
        complex_positions = unit.Quantity(
            np.zeros([total_atoms, 3]), unit=unit.nanometers
        )
        complex_positions[: len(protein.positions)] = protein.positions
        complex_positions[len(protein.positions) :] = ligand.positions

        return complex_topology, complex_positions


class WaterComplexBuilder:
    """
    A class for building solvent boxes around a protein complex using different water models.

    This class handles the addition of solvent to a protein structure, either by specifying padding distance
    or absolute dimensions of the solvent box. It also provides functionality for converting between different
    water models.

    Attributes:
        model_water (str): Selected water model (e.g., 'charmm', 'tip3pfb', 'tip4pew', etc.).
        forcefield (openmm.app.forcefield.ForceField): Selected Forcefield for the system.
        protein_pdb (pdbfixer.pdbfixer.PDBFixer): Protein structure as a pdbfixer object.
        protein_name (str): Name of the protein file (used for output file naming).
        water_positive_ion (str): Positive ion to add to the solvent box.
        water_negative_ion (str): Negative ion to add to the solvent box.
        water_ionicstrength (float): Ionic strength of the solvent in Molar.
        modeller (openmm.app.modeller.Modeller): OpenMM Modeller object for building and modifying the complex.
    """

    def __init__(
        self,
        model_water,
        forcefield,
        protein_pdb,
        protein_name,
        water_positive_ion,
        water_negative_ion,
        water_ionicstrength,
        modeller,
    ):
        """
        Initialize the ComplexBuilder with required parameters.
        """
        self.model_water = model_water
        self.forcefield = forcefield
        self.protein_pdb = protein_pdb
        self.protein_name = protein_name
        self.water_positive_ion = water_positive_ion
        self.water_negative_ion = water_negative_ion
        self.water_ionicstrength = water_ionicstrength
        self.modeller = modeller

    def water_padding_solvent_builder(self, water_padding_distance):
        """
        Build a solvent box with padding distance.

        Args:
            water_padding_distance (float): Solvent padding distance in nanometers.

        Returns:
            modeller (openmm.app.modeller.Modeller): The complex with solvent.
        """
        with open(f"prepared_no_solvent_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(
                self.protein_pdb.topology, self.protein_pdb.positions, outfile
            )

        if self.model_water in ["charmm", "tip3pfb", "tip3"]:
            self.modeller.addSolvent(
                self.forcefield,
                padding=water_padding_distance * unit.nanometers,
                positiveIon=self.water_positive_ion,
                negativeIon=self.water_negative_ion,
                ionicStrength=self.water_ionicstrength * unit.molar,
            )
        elif self.model_water == "charmm_tip4pew":
            self.protein_pdb.addSolvent(
                padding=water_padding_distance * unit.nanometers,
                positiveIon=self.water_positive_ion,
                negativeIon=self.water_negative_ion,
                ionicStrength=self.water_ionicstrength * unit.molar,
            )
        else:
            if self.model_water == "tip4pfb":
                self.model_water = "tip4pew"
            self.modeller.addSolvent(
                self.forcefield,
                model=self.model_water,
                padding=water_padding_distance * unit.nanometers,
                positiveIon=self.water_positive_ion,
                negativeIon=self.water_negative_ion,
                ionicStrength=self.water_ionicstrength * unit.molar,
            )

        with open(f"solvent_padding_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, outfile)
        print(f"Protein with buffer solvent prepared")

        return self.modeller

    def water_absolute_solvent_builder(self, water_box_x, water_box_y, water_box_z):
        """
        Build a solvent box with absolute dimensions.

        Args:
            water_box_x (float): Vector x of the solvent box in nanometers.
            water_box_y (float): Vector y of the solvent box in nanometers.
            water_box_z (float): Vector z of the solvent box in nanometers.

        Returns:
            modeller (openmm.app.modeller.Modeller): The complex with solvent.
        """
        with open(f"prepared_no_solvent_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(
                self.protein_pdb.topology, self.protein_pdb.positions, outfile
            )

        if self.model_water in ["charmm", "tip3pfb", "tip3"]:
            self.modeller.addSolvent(
                self.forcefield,
                boxSize=Vec3(water_box_x, water_box_y, water_box_z) * unit.nanometers,
                positiveIon=self.water_positive_ion,
                negativeIon=self.water_negative_ion,
                ionicStrength=self.water_ionicstrength * unit.molar,
            )
        elif self.model_water == "charmm_tip4pew":
            self.protein_pdb.addSolvent(
                boxSize=Vec3(water_box_x, water_box_y, water_box_z) * unit.nanometers,
                positiveIon=self.water_positive_ion,
                negativeIon=self.water_negative_ion,
                ionicStrength=self.water_ionicstrength * unit.molar,
            )
        else:
            if self.model_water == "tip4pfb":
                self.model_water = "tip4pew"
            self.modeller.addSolvent(
                self.forcefield,
                model=self.model_water,
                boxSize=Vec3(water_box_x, water_box_y, water_box_z) * unit.nanometers,
                positiveIon=self.water_positive_ion,
                negativeIon=self.water_negative_ion,
                ionicStrength=self.water_ionicstrength * unit.molar,
            )

        with open(f"solvent_absolute_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, outfile)
        print(f"Protein with absolute solvent prepared")

        return self.modeller

    def water_conversion(self):
        """
        Convert the water model of an OpenMM object.

        Args:
            model_water (str): The name of the preferred converted water model.

        Returns:
            modeller (openmm.app.modeller.Modeller): The converted object.
        """
        with open(f"pre_converted_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, outfile)

        self.modeller.convertWater(self.model_water)

        with open(f"converted_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, outfile)

        return self.modeller


class MembraneBuilder:
    """
    A class for building a membrane around a protein complex.

    This class facilitates the addition of a membrane to a protein structure, including specifying lipid types,
    membrane padding, and ionic conditions. It supports different forcefields and water models to accommodate
    various simulation setups.

    Attributes:
        ff (str): Selected forcefield as a string (e.g., 'CHARMM36').
        model_water (str): Selected water model (e.g., 'charmm', 'tip4pew', 'tip5p').
        forcefield (openmm.app.forcefield.ForceField): Selected forcefield for the simulation.
        transitional_forcefield (openmm.app.forcefield.ForceField): Transitional forcefield for specific water models.
        protein_pdb (pdbfixer.pdbfixer.PDBFixer): Protein structure as a pdbfixer object.
        modeller (openmm.app.modeller.Modeller): OpenMM Modeller object used to build the complex.
        membrane_lipid_type (str): Lipid type to be used in the membrane.
        membrane_padding (float): Minimum padding distance for the membrane in nanometers.
        membrane_positive_ion (str): Positive ion to be added to the membrane.
        membrane_negative_ion (str): Negative ion to be added to the membrane.
        membrane_ionicstrength (float): Ionic strength of the membrane in Molar.
        protein_name (str): Name of the protein file (used for naming output files).
    """

    def __init__(
        self,
        ff,
        model_water,
        forcefield,
        transitional_forcefield,
        protein_pdb,
        modeller,
        membrane_lipid_type,
        membrane_padding,
        membrane_positive_ion,
        membrane_negative_ion,
        membrane_ionicstrength,
        protein_name,
    ):
        """
        Initialize the MembraneBuilder with necessary parameters.
        """
        self.ff = ff
        self.model_water = model_water
        self.forcefield = forcefield
        self.transitional_forcefield = transitional_forcefield
        self.protein_pdb = protein_pdb
        self.modeller = modeller
        self.membrane_lipid_type = membrane_lipid_type
        self.membrane_padding = membrane_padding
        self.membrane_positive_ion = membrane_positive_ion
        self.membrane_negative_ion = membrane_negative_ion
        self.membrane_ionicstrength = membrane_ionicstrength
        self.protein_name = protein_name

    def build_membrane(self):
        """
        Build a membrane with minimum padding.

        Returns:
            modeller (openmm.app.modeller.Modeller): The complex with solvent.
        """
        # Writing out the protein without solvent
        with open(f"prepared_no_solvent_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(
                self.protein_pdb.topology, self.protein_pdb.positions, outfile
            )

        # Adds a membrane to the selected protein
        if self.ff == "CHARMM36":
            self.protein_pdb.addMembrane(
                lipidType=self.membrane_lipid_type,
                minimumPadding=self.membrane_padding * unit.nanometer,
                positiveIon=self.membrane_positive_ion,
                negativeIon=self.membrane_negative_ion,
                ionicStrength=self.membrane_ionicstrength * unit.molar,
            )
            self.modeller = app.Modeller(
                self.protein_pdb.topology, self.protein_pdb.positions
            )
        else:
            if self.model_water == "charmm":
                self.modeller.addMembrane(
                    self.forcefield,
                    lipidType=self.membrane_lipid_type,
                    minimumPadding=self.membrane_padding * unit.nanometer,
                    positiveIon=self.membrane_positive_ion,
                    negativeIon=self.membrane_negative_ion,
                    ionicStrength=self.membrane_ionicstrength * unit.molar,
                )
            else:
                if self.model_water in ["tip4pew", "tip5p"]:
                    self.modeller.addMembrane(
                        self.transitional_forcefield,
                        lipidType=self.membrane_lipid_type,
                        minimumPadding=self.membrane_padding * unit.nanometer,
                        positiveIon=self.membrane_positive_ion,
                        negativeIon=self.membrane_negative_ion,
                        ionicStrength=self.membrane_ionicstrength * unit.molar,
                    )
                else:
                    self.modeller.addMembrane(
                        self.forcefield,
                        lipidType=self.membrane_lipid_type,
                        minimumPadding=self.membrane_padding * unit.nanometer,
                        positiveIon=self.membrane_positive_ion,
                        negativeIon=self.membrane_negative_ion,
                        ionicStrength=self.membrane_ionicstrength * unit.molar,
                    )

        with open(f"membrane_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, outfile)

        print(f"Protein with Membrane {self.membrane_lipid_type} prepared")

        return self.modeller

import mdtraj as md
import numpy as np
import pdbfixer
import simtk.openmm.app as app
from rdkit import Chem
from rdkit.Chem import rdForceFieldHelpers
from openff.toolkit.topology import Molecule
from simtk.openmm.app import PDBFile
from simtk.openmm import unit
from simtk.openmm import Vec3


class PDBFormatter:
    def __init__(self, config_parser):
        self.config_parser = config_parser
        self.ligand = self.config_parser.get_variable('ligand', None)
        self.protein = self.config_parser.get_variable('protein')
        self.protein_pdb = self._get_protein_pdb()

    def _get_protein_pdb(self):
        if self.ligand is None:
            protein_pdb = PDBFile(self.protein)
        else:
            protein_pdb = pdbfixer.PDBFixer(self.protein)
        return protein_pdb

    def get_protein_pdb(self):
        return self.protein_pdb


class ComplexBuilder:
    def __init__(self, config_parser, protein_pdb, omm_ligand=None):
        """
        Initialize the ComplexBuilder with configuration and entities.

        Parameters:
        - config_parser (ConfigParser): The configuration parser instance.
        - protein_pdb (object): The protein PDB object.
        - molecule_preparation (MoleculePreparation): The molecule preparation instance.
        - omm_ligand (object, optional): The ligand in OpenMM format. Default is None.
        """
        self.config_parser = config_parser
        self.protein_pdb = protein_pdb
        self.omm_ligand = omm_ligand

    def build_complex(self):
        """
        Build the complex based on the presence of a ligand in the configuration.

        Returns:
        - tuple: (complex_topology, complex_positions) or (None, None)
        """
        ligand = self.config_parser.ligand

        if ligand is not None and self.omm_ligand is not None:
            print("Merging protein and ligand...")
            self.molecule_preparation = MoleculePreparation(self.config_parser.ligand)
            complex_topology, complex_positions = self.molecule_preparation.merge_protein_and_ligand(
                self.protein_pdb, self.omm_ligand
            )
            print("Complex construction complete.")
            return complex_topology, complex_positions
        else:
            print("No ligand specified or ligand not available; skipping complex construction.")
            return None, None


class ModellerCreator:
    def __init__(self, config_parser, protein_pdb, complex_topology=None, complex_positions=None):
        """
        Initialize the ModellerCreator with the necessary parameters.

        Parameters:
        - config_parser (ConfigParser): The configuration parser instance.
        - protein_pdb (app.PDBFile): The protein PDB object.
        - complex_topology (Topology, optional): The topology of the protein-ligand complex. Default is None.
        - complex_positions (Quantity, optional): The positions of the protein-ligand complex. Default is None.
        """
        self.config_parser = config_parser
        self.protein_pdb = protein_pdb
        self.complex_topology = complex_topology
        self.complex_positions = complex_positions

    def create_modeller(self):
        """
        Create and return an OpenMM Modeller object based on the presence of a ligand.

        Returns:
        - modeller (app.Modeller): The OpenMM Modeller object.
        """
        ligand = self.config_parser.ligand

        if ligand is None:
            print("Creating modeller with protein only...")
            modeller = app.Modeller(self.protein_pdb.topology, self.protein_pdb.positions)
        else:
            print("Creating modeller with protein-ligand complex...")
            modeller = app.Modeller(self.complex_topology, self.complex_positions)

        return modeller

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

    def __init__(self, ligand_file, minimize_molecule=False):
        """
        Initialize the MoleculePreparation instance with the ligand file and minimization flag.
        """
        self.ligand_file = ligand_file
        self.minimize_molecule = minimize_molecule
        self.ligand = None

    def prepare_ligand(self):
        """Reads an SDF File into RDKit, adds hydrogens to the structure, minimizes it if selected, and creates an openforcefield Molecule object. Inspired by @teachopencadd T019.

        Returns:
            rdkitmolh (rdkit.Chem.rdchem.Mol): The prepared and converted ligand.
        """
        file_name = self.ligand_file.lower()
        if file_name.endswith(".sdf"):
            rdkit_mol = Chem.SDMolSupplier(self.ligand_file, sanitize=False)
            for mol in rdkit_mol:
                rdkit_mol = mol
        elif file_name.endswith(".mol") and not file_name.endswith(".mol2"):
            rdkit_mol = Chem.rdmolfiles.MolFromMolFile(self.ligand_file, sanitize=False)
        elif file_name.endswith(".mol2"):
            rdkit_mol = Chem.rdmolfiles.MolFromMol2File(
                self.ligand_file, sanitize=False
            )

        print("Adding hydrogens")
        rdkitmolh = Chem.AddHs(rdkit_mol, addCoords=True)
        Chem.AssignAtomChiralTagsFromStructure(rdkitmolh)

        if self.minimize_molecule == True:
            rdForceFieldHelpers.MMFFOptimizeMolecule(
                mol=rdkitmolh, mmffVariant="MMFF94s", maxIters=2000
            )

        Molecule(rdkitmolh)
        self.ligand = rdkitmolh
        return rdkitmolh

    @staticmethod
    def rdkit_to_openmm(rdkit_mol, name):
        """Convert an RDKit molecule to an OpenMM molecule. Inspired by @teachopencadd T019, @hannahbrucemcdonald and @glass-w.

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


class WaterMembraneComplexBuilder:
    """
    A class for building solvent boxes or membranes around a protein complex using different water models.

    This class handles the addition of solvent to a protein structure, either by specifying padding distance
    or absolute dimensions of the solvent box. It also provides functionality for converting between different
    water models and adding membranes.

    Attributes:
        model_water (str): Selected water model (e.g., 'charmm', 'tip3pfb', 'tip4pew', etc.).
        forcefield (openmm.app.ForceField): Selected forcefield for the system.
        protein_pdb (pdbfixer.pdbfixer.PDBFixer): Protein structure as a pdbfixer object.
        protein_name (str): Name of the protein file (used for output file naming).
        modeller (openmm.app.Modeller): OpenMM Modeller object for building and modifying the complex.
        membrane_lipid_type (str, optional): Lipid type to be used in the membrane.
        membrane_padding (float, optional): Minimum padding distance for the membrane in nanometers.
        membrane_positive_ion (str, optional): Positive ion to add to the membrane.
        membrane_negative_ion (str, optional): Negative ion to add to the membrane.
        membrane_ionicstrength (float, optional): Ionic strength of the membrane in Molar.
        water_positive_ion (str, optional): Positive ion to add to the solvent box.
        water_negative_ion (str, optional): Negative ion to add to the solvent box.
        water_ionicstrength (float, optional): Ionic strength of the solvent in Molar.
    """

    def __init__(
        self,
        model_water,
        forcefield,
        protein_pdb,
        protein_name,
        modeller,
        membrane_lipid_type=None,
        membrane_padding=None,
        membrane_positive_ion=None,
        membrane_negative_ion=None,
        membrane_ionicstrength=None,
        water_positive_ion=None,
        water_negative_ion=None,
        water_ionicstrength=None,
    ):
        """
        Initialize the WaterComplexBuilder with required parameters.
        """
        self.model_water = model_water
        self.forcefield = forcefield
        self.protein_pdb = protein_pdb
        self.protein_name = protein_name
        self.modeller = modeller
        self.membrane_lipid_type = membrane_lipid_type
        if membrane_padding != None:
            self.membrane_padding = float(membrane_padding)

        else:
            self.membrane_padding = membrane_padding
        self.membrane_positive_ion = membrane_positive_ion
        self.membrane_negative_ion = membrane_negative_ion
        if membrane_ionicstrength != None:
            self.membrane_ionicstrength = float(membrane_ionicstrength)        
        else:
            self.membrane_ionicstrength = membrane_ionicstrength
        self.water_positive_ion = water_positive_ion
        self.water_negative_ion = water_negative_ion
        if water_ionicstrength != None:
            self.water_ionicstrength = float(water_ionicstrength)        
        else:
            self.water_ionicstrength = water_ionicstrength

    def build_water_padding_solvent(self, water_padding_distance):
        """
        Build a solvent box with padding distance.

        Args:
            water_padding_distance (float): Solvent padding distance in nanometers.

        Returns:
            modeller (openmm.app.Modeller): The complex with solvent.
        """
        with open(f"prepared_no_solvent_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.protein_pdb.topology, self.protein_pdb.positions, outfile)

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

    def build_water_absolute_solvent(self, water_box_x, water_box_y, water_box_z):
        """
        Build a solvent box with absolute dimensions.

        Args:
            water_box_x (float): Vector x of the solvent box in nanometers.
            water_box_y (float): Vector y of the solvent box in nanometers.
            water_box_z (float): Vector z of the solvent box in nanometers.

        Returns:
            modeller (openmm.app.Modeller): The complex with solvent.
        """
        with open(f"prepared_no_solvent_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.protein_pdb.topology, self.protein_pdb.positions, outfile)

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

    def convert_water(self):
        """
        Convert the water model of an OpenMM object.

        Returns:
            modeller (openmm.app.Modeller): The converted object.
        """
        with open(f"pre_converted_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, outfile)

        self.modeller.convertWater(self.model_water)

        with open(f"converted_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.modeller.topology, self.modeller.positions, outfile)

        return self.modeller

    def build_membrane(self):
        """
        Build a membrane with minimum padding.

        Returns:
            modeller (openmm.app.Modeller): The complex with membrane.
        """
        with open(f"prepared_no_solvent_{self.protein_name}", "w") as outfile:
            PDBFile.writeFile(self.protein_pdb.topology, self.protein_pdb.positions, outfile)

        if self.forcefield == "CHARMM36":
            self.protein_pdb.addMembrane(
                lipidType=self.membrane_lipid_type,
                minimumPadding=self.membrane_padding * unit.nanometer,
                positiveIon=self.membrane_positive_ion,
                negativeIon=self.membrane_negative_ion,
                ionicStrength=self.membrane_ionicstrength * unit.molar,
            )
            self.modeller = app.Modeller(self.protein_pdb.topology, self.protein_pdb.positions)
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



class EnvironmentBuilder:
    def __init__(self, config_parser, forcefield_name, model_water, forcefield, transitional_forcefield, protein_pdb, modeller, protein, ligand=None):
        """
        Initialize the EnvironmentBuilder with the necessary parameters.

        Parameters:
        - config_parser (ConfigParser): The configuration parser instance.
        - forcefield_name: Forcefield object.
        - model_water: Water model.
        - forcefield: Selected forcefield.
        - transitional_forcefield: Transitional forcefield for membrane.
        - protein_pdb (app.PDBFile): The protein PDB object.
        - modeller: OpenMM Modeller object.
        - protein: Protein object or identifier.
        - ligand (object, optional): The ligand in OpenMM format. Default is None.
        """
        self.config_parser = config_parser
        self.forcefield_name = forcefield_name
        self.model_water = model_water
        self.forcefield = forcefield
        self.transitional_forcefield = transitional_forcefield
        self.protein_pdb = protein_pdb
        self.modeller = modeller
        self.protein = protein
        self.ligand = ligand

        self.water_complex_builder = WaterMembraneComplexBuilder(
            model_water=self.model_water,
            forcefield=self.forcefield,
            protein_pdb=self.protein_pdb,
            protein_name=protein,
            modeller=self.modeller,
            membrane_lipid_type=self.config_parser.membrane_lipid_type,
            membrane_padding=self.config_parser.membrane_padding,
            membrane_positive_ion=self.config_parser.membrane_positive_ion,
            membrane_negative_ion=self.config_parser.membrane_negative_ion,
            membrane_ionicstrength=self.config_parser.membrane_ionicstrength,
            water_positive_ion=self.config_parser.water_positive_ion,
            water_negative_ion=self.config_parser.water_negative_ion,
            water_ionicstrength=self.config_parser.water_ionicstrength
        )

    def build_environment(self):
        """
        Build the environment based on the presence of a ligand and the configuration settings.
        """
        add_membrane = self.config_parser.add_membrane
        water_box = self.config_parser.Water_Box

        if self.ligand is not None:
            self._handle_with_ligand(add_membrane, water_box)
        else:
            self._handle_without_ligand(add_membrane, water_box)

    def _handle_with_ligand(self, add_membrane, water_box):
        """
        Handle environment setup when a ligand is present.
        """
        if add_membrane == "True":
            self._membrane_builder()
        else:
            if water_box == "Buffer":
                self._water_padding_solvent_builder()
            elif water_box == "Absolute":
                self._water_absolute_solvent_builder()

        if add_membrane and (self.model_water == 'tip4pew' or self.model_water == 'tip5p'):
            self._water_conversion()

    def _handle_without_ligand(self, add_membrane, water_box):
        """
        Handle environment setup when no ligand is present.
        """
        if add_membrane == "True":
            self._membrane_builder()
        else:
            if water_box == "Buffer":
                self._water_padding_solvent_builder()
            elif water_box == "Absolute":
                self._water_absolute_solvent_builder()

        if add_membrane and (self.model_water == 'tip4pew' or self.model_water == 'tip5p'):
            self._water_conversion()

    def _membrane_builder(self):
        """
        Build the membrane environment.
        """
        self.water_complex_builder.build_membrane()

    def _water_padding_solvent_builder(self):
        """
        Build the water padding solvent environment.
        """
        water_padding_distance = float(self.config_parser.water_padding_distance)
        self.water_complex_builder.build_water_padding_solvent(water_padding_distance)

    def _water_absolute_solvent_builder(self):
        """
        Build the water absolute solvent environment.
        """
        water_box_x = float(self.config_parser.water_box_x)
        water_box_y = float(self.config_parser.water_box_y)
        water_box_z = float(self.config_parser.water_box_t)
        self.water_complex_builder.build_water_absolute_solvent(water_box_x, water_box_y, water_box_z)

    def _water_conversion(self):
        """
        Convert water model if necessary.
        """
        self.water_complex_builder.convert_water()


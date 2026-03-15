import os
import parmed as pmd
import simtk.openmm.app as app
from simtk.openmm.app import PDBFile
from simtk.openmm import unit
from openmm.openmm import XmlSerializer
import pdbfixer

from openmmdl.openmmdl_simulation.scripts.forcefield_water import (
    ff_selection,
    water_forcefield_selection,
    water_model_selection,
    generate_forcefield,
    generate_transitional_forcefield,
)
from openmmdl.openmmdl_simulation.scripts.protein_ligand_prep import (
    prepare_ligand,
    rdkit_to_openmm,
    merge_protein_and_ligand,
    water_padding_solvent_builder,
    water_absolute_solvent_builder,
    membrane_builder,
    water_conversion,
)


class PreparedSystemBuilder:
    def __init__(self, config):
        self.config = config

    def build(self):
        protein = self.config["protein"]
        protein_name = os.path.basename(protein)
        ligand = self.config.get("ligand")
        ligand_name = self.config.get("ligand_name", "UNK")
        minimization = self.config.get("ligand_minimization", False)
        sanitization = self.config.get("ligand_sanitization", False)

        ff = self.config["forcefield"]
        water = self.config.get("water_model")
        add_membrane = self.config.get("add_membrane", False)
        use_solvent = self.config.get("use_solvent", False)

        small_molecule_ff = self.config.get("small_molecule_forcefield")
        small_molecule_ff_version = self.config.get("small_molecule_forcefield_version")

        forcefield_selected = ff_selection(ff)
        water_selected = water_forcefield_selection(
            water=water,
            forcefield_selection=ff_selection(ff),
        )
        model_water = water_model_selection(
            water=water,
            forcefield_selection=ff_selection(ff),
        )

        if ligand:
            ligand_prepared = prepare_ligand(
                ligand,
                minimize_molecule=minimization,
            )
            omm_ligand = rdkit_to_openmm(ligand_prepared, ligand_name)
            protein_pdb = pdbfixer.PDBFixer(str(protein))

            if add_membrane:
                transitional_forcefield = generate_transitional_forcefield(
                    protein_ff=forcefield_selected,
                    solvent_ff=water_selected,
                    add_membrane=add_membrane,
                    smallMoleculeForceField=small_molecule_ff,
                    smallMoleculeForceFieldVersion=small_molecule_ff_version,
                    rdkit_mol=ligand_prepared,
                )

            forcefield = generate_forcefield(
                protein_ff=forcefield_selected,
                solvent_ff=water_selected,
                add_membrane=add_membrane,
                smallMoleculeForceField=small_molecule_ff,
                smallMoleculeForceFieldVersion=small_molecule_ff_version,
                rdkit_mol=ligand_prepared,
            )

            complex_topology, complex_positions = merge_protein_and_ligand(
                protein_pdb,
                omm_ligand,
            )
            modeller = app.Modeller(complex_topology, complex_positions)

        else:
            protein_pdb = PDBFile(protein)

            if water_selected is not None:
                forcefield = generate_forcefield(
                    protein_ff=forcefield_selected,
                    solvent_ff=water_selected,
                    add_membrane=add_membrane,
                    smallMoleculeForceField=small_molecule_ff,
                    smallMoleculeForceFieldVersion=small_molecule_ff_version,
                    rdkit_mol=None,
                )
            else:
                forcefield = app.ForceField(forcefield_selected)

            if add_membrane:
                transitional_forcefield = generate_transitional_forcefield(
                    protein_ff=forcefield_selected,
                    solvent_ff=water_selected,
                    add_membrane=add_membrane,
                    smallMoleculeForceField=small_molecule_ff,
                    smallMoleculeForceFieldVersion=small_molecule_ff_version,
                    rdkit_mol=None,
                )

            modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)

        if use_solvent:
            if add_membrane:
                membrane_builder(
                    ff,
                    model_water,
                    forcefield,
                    transitional_forcefield,
                    protein_pdb,
                    modeller,
                    self.config["membrane_lipid_type"],
                    self.config["membrane_padding"],
                    self.config["membrane_positive_ion"],
                    self.config["membrane_negative_ion"],
                    self.config["membrane_ionicstrength"],
                    protein_name,
                )
            else:
                water_box_mode = self.config["water_box_mode"]
                if water_box_mode == "Buffer":
                    water_padding_solvent_builder(
                        model_water,
                        forcefield,
                        self.config["water_padding_distance"],
                        protein_pdb,
                        modeller,
                        self.config["water_positive_ion"],
                        self.config["water_negative_ion"],
                        self.config["water_ionicstrength"],
                        protein_name,
                    )
                elif water_box_mode == "Absolute":
                    water_absolute_solvent_builder(
                        model_water,
                        forcefield,
                        self.config["water_box_x"],
                        self.config["water_box_y"],
                        self.config["water_box_z"],
                        protein_pdb,
                        modeller,
                        self.config["water_positive_ion"],
                        self.config["water_negative_ion"],
                        self.config["water_ionicstrength"],
                        protein_name,
                    )

        if add_membrane and model_water in ("tip4pew", "tip5p"):
            water_conversion(model_water, modeller, protein_name)

        topology = modeller.topology
        positions = modeller.positions

        nonbonded_method = getattr(app, self.config["nonbonded_method"])
        kwargs = {
            "nonbondedMethod": nonbonded_method,
            "constraints": self.config["constraints"],
            "rigidWater": self.config["rigid_water"],
        }

        if self.config.get("nonbonded_cutoff") is not None:
            kwargs["nonbondedCutoff"] = self.config["nonbonded_cutoff"] * unit.nanometers

        if self.config.get("ewald_error_tolerance") is not None:
            kwargs["ewaldErrorTolerance"] = self.config["ewald_error_tolerance"]

        if self.config.get("hydrogen_mass") is not None:
            kwargs["hydrogenMass"] = self.config["hydrogen_mass"] * unit.amu

        system = forcefield.createSystem(topology, **kwargs)

        return topology, system, positions


class PreparedSystemExporter:
    def __init__(self, topology, system, positions, output_prefix="prepared_system", output_dir="."):
        self.topology = topology
        self.system = system
        self.positions = positions
        self.output_prefix = output_prefix
        self.output_dir = output_dir
        self._structure = None

    def _path(self, suffix):
        return os.path.join(self.output_dir, f"{self.output_prefix}{suffix}")

    def build_parmed_structure(self):
        if self._structure is None:
            self._structure = pmd.openmm.load_topology(
                self.topology,
                self.system,
                xyz=self.positions,
            )
        return self._structure

    def write_processed_pdb(self):
        path = self._path(".pdb")
        with open(path, "w") as handle:
            PDBFile.writeFile(self.topology, self.positions, handle)
        return path

    def write_amber(self):
        structure = self.build_parmed_structure()
        prmtop = self._path(".prmtop")
        inpcrd = self._path(".inpcrd")
        structure.save(prmtop, overwrite=True)
        structure.save(inpcrd, overwrite=True)
        return [prmtop, inpcrd]

    def write_gromacs(self):
        structure = self.build_parmed_structure()
        top = self._path(".top")
        gro = self._path(".gro")
        structure.save(top, overwrite=True)
        structure.save(gro, overwrite=True)
        return [top, gro]

    def write_psf(self):
        structure = self.build_parmed_structure()
        psf = self._path(".psf")
        pdb = self._path(".pdb")
        structure.save(psf, overwrite=True)
        structure.save(pdb, overwrite=True)
        return [psf, pdb]

    def write_openmm_xml(self):
        xml_path = self._path(".xml")
        with open(xml_path, "w") as handle:
            handle.write(XmlSerializer.serialize(self.system))
        return xml_path

    def export_selected(self, formats):
        written_files = []

        if "processed_pdb" in formats:
            written_files.append(self.write_processed_pdb())

        if "amber" in formats:
            written_files.extend(self.write_amber())

        if "gromacs" in formats:
            written_files.extend(self.write_gromacs())

        if "psf" in formats:
            written_files.extend(self.write_psf())

        if "openmm_xml" in formats:
            written_files.append(self.write_openmm_xml())

        return written_files
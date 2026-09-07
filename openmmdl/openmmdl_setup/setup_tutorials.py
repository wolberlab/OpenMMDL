"""Step definitions and bundled files for the guided tutorials shown inside OpenMMDL Setup.

Each page entry is keyed by the ``tutorial_page`` name a template is rendered with. A step
highlights ``target`` (a CSS selector), may offer an ``action`` button handled in
``static/openmmdl_tutorial.js``, and may declare a ``done`` condition that shows a check mark
once the user (or the action) has made the expected change.
"""

from pathlib import Path

PDB_TUTORIAL_ID = "pdb-small-molecule"

# Files served by the "Download files" button, taken from the bundled tutorial system.
PDB_TUTORIAL_DIR = (
    Path(__file__).resolve().parents[1]
    / "openmmdl_simulation"
    / "tutorial_systems"
    / "pdb_path"
    / "5wyz_solvent"
)

PDB_TUTORIAL_FILES = (
    "5wyz_tutorial_receptor.pdb",
    "7VF_A.sdf",
    "7VF_B.sdf",
)

PDB_TUTORIAL_README = """OpenMMDL PDB small-molecule tutorial files (bundled with the openmmdl package).

5wyz_tutorial_receptor.pdb
    Human TLR8 ectodomain dimer from PDB entry 5WYZ, processed in MOE: chains A and B,
    hydrogens added, chain breaks capped with ACE/NME. The N-linked glycans of the
    crystal structure (NAG/BMA/MAN, chains C-H and single NAGs on A/B) are kept so the
    chain-selection step has sugar chains to remove. The coordinates are in the
    original crystal frame.

7VF_A.sdf, 7VF_B.sdf
    The two copies of the CU-CPT9b ligand (PDB residue 7VF) from 5WYZ, one per binding
    site, exported from MOE with the crystal heavy-atom coordinates. Hydrogens are added
    by OpenMMDL Simulation.

In the tutorial, upload the PDB as the receptor, 7VF_A.sdf as the ligand
(topology code UNK) and 7VF_B.sdf as the additional molecule (topology code L01).
"""


PDB_TUTORIAL_PAGES = {
    "select_file_type": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Choose the PDB path",
        "pageNumber": 1,
        "totalPages": 7,
        "pageKey": "select_file_type",
        "docsUrl": "pdb_path.html",
        "showFiles": True,
        "steps": [
            {
                "target": "input[name='type'][value='pdb']",
                "expectedInputLabel": "Select the PDB / mmCIF input type",
                "title": "Start with a PDB/mmCIF system",
                "body": "The PDB path prepares a protein structure with PDBFixer before OpenMMDL generates the simulation files. Select this option for the small-molecule PDB tutorial.",
                "action": "selectPdbPath",
                "actionLabel": "Select PDB path",
                "done": {
                    "selector": "input[name='type'][value='pdb']",
                    "condition": "checked",
                },
            },
            {
                "target": "input[type='submit']",
                "title": "Continue to file selection",
                "body": "Continue opens the page where you provide the receptor structure and optional ligand files.",
            },
        ],
    },
    "configure_pdb_file": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Select the PDB and ligand files",
        "pageNumber": 2,
        "totalPages": 7,
        "pageKey": "configure_pdb_file",
        "docsUrl": "pdb_path.html#selecting-input-files",
        "showFiles": True,
        "steps": [
            {
                "target": ".om-filedrop[data-file-input='file']",
                "expectedInputLabel": "Upload tutorial receptor PDB file",
                "done": {
                    "selector": "#file",
                    "condition": "hasValue",
                },
                "title": "Upload the receptor PDB",
                "body": "Use Browse to select the tutorial receptor PDB. Browsers cannot select local files automatically, so this step remains manual.",
            },
            {
                "target": "#forcefield",
                "title": "Switch the force field to AMBER19",
                "body": "Choose AMBER19. OpenMMDL will automatically use OPC3 as the default water model for AMBER19.",
                "action": "setForceFieldAmber19",
                "actionLabel": "Set AMBER19",
                "done": {
                    "selector": "#forcefield",
                    "condition": "valueEquals",
                    "value": "AMBER19",
                },
                "expectedOptionLabel": "AMBER19",
            },
            {
                "target": "#smallMoleculeMode",
                "title": "Enable a single ligand complex",
                "body": "The small-molecule PDB tutorial uses one ligand, so choose Single complex. This reveals the ligand upload and topology-code fields.",
                "action": "configurePdbTutorialDefaults",
                "actionLabel": "Choose single complex",
                "done": {
                    "selector": "#smallMoleculeMode",
                    "condition": "valueEquals",
                    "value": "single",
                },
                "expectedOptionLabel": "Single complex",
            },
            {
                "target": ".om-filedrop[data-file-input='sdfFile']",
                "expectedInputLabel": "Upload 7VF_A.sdf",
                "done": {
                    "selector": "#sdfFile",
                    "condition": "hasValue",
                },
                "title": "Upload the ligand file",
                "body": "Select 7VF_A.sdf from the downloaded tutorial files: the CU-CPT9b ligand bound in the chain A site. Keep the topology code as UNK unless you need a different three-character residue code.",
            },
            {
                "target": "#addCompanionMolecule",
                "expectedInputLabel": "Add one additional molecule",
                "title": "Add an additional molecule",
                "body": "TLR8 is a dimer with one ligand in each protomer. Click + Add additional molecule to include the second CU-CPT9b copy, with its own topology code, in the generated system.",
                "action": "addTutorialCompanionMolecule",
                "actionLabel": "Add additional molecule",
                "done": {
                    "selector": "[data-companion-row]",
                    "condition": "exists",
                },
            },
            {
                "target": ".om-companion-upload",
                "expectedInputLabel": "Upload 7VF_B.sdf",
                "done": {
                    "selector": "input[name='companionFile']",
                    "condition": "hasValue",
                },
                "title": "Upload the additional molecule",
                "body": "Use Browse in the new row to upload 7VF_B.sdf, the ligand copy bound in the chain B site. Keep or edit its topology code (default L01).",
            },
            {
                "target": "#smallMoleculeForceField",
                "title": "Switch the ligand force field to SMIRNOFF",
                "body": "Change the small-molecule force field from GAFF to SMIRNOFF.",
                "action": "setLigandForceFieldSmirnoff",
                "actionLabel": "Set SMIRNOFF",
                "done": {
                    "selector": "#smallMoleculeForceField",
                    "condition": "valueEquals",
                    "value": "smirnoff",
                },
                "expectedOptionLabel": "SMIRNOFF",
            },
            {
                "target": "#continue",
                "title": "Proceed after both files are selected",
                "body": "Continue remains disabled until a receptor file is selected and all visible topology codes are valid and unique.",
            },
        ],
    },
    "select_chains": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Keep the receptor chains",
        "pageNumber": 3,
        "totalPages": 7,
        "pageKey": "select_chains",
        "docsUrl": "pdb_path.html#selecting-chains",
        "steps": [
            {
                "target": "#table",
                "expectedInputLabel": "Keep chains A and B, deselect the glycan chains",
                "title": "Review detected chains",
                "body": "The receptor contains the two TLR8 protomers, chains A and B, plus the N-linked glycans of the crystal structure (NAG, BMA, MAN) as separate chains. Deselect every glycan chain so PDBFixer removes the sugars, and keep all heterogens: the ACE/NME caps at the chain breaks count as heterogens.",
                "action": "selectFirstTwoChains",
                "actionLabel": "Keep chains A and B",
            },
            {
                "target": "#continue",
                "title": "Continue to structure checks",
                "body": "PDBFixer next checks for missing residues and heavy atoms. Pages with nothing to add are skipped automatically, so the tutorial may jump straight to the hydrogens and water step.",
            },
        ],
    },
    "add_residues": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Add selected missing residues",
        "pageNumber": 4,
        "totalPages": 7,
        "pageKey": "add_residues",
        "docsUrl": "pdb_path.html#optional-adding-residues",
        "steps": [
            {
                "target": "#table",
                "expectedInputLabel": "Select all missing residue spans",
                "title": "Inspect missing sequence spans",
                "body": "For the PDB tutorial, add all missing residue spans so PDBFixer rebuilds them. Use Select all missing residues, then continue.",
                "action": "selectAllMissingResidues",
                "actionLabel": "Select all missing residues",
            },
            {
                "target": "input[type='submit']",
                "title": "Continue after selecting spans",
                "body": "The next pages convert non-standard residues if needed and add missing heavy atoms.",
            },
        ],
    },
    "convert_residues": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Convert non-standard residues if present",
        "pageNumber": 4,
        "totalPages": 7,
        "pageKey": "convert_residues",
        "docsUrl": "pdb_path.html",
        "steps": [
            {
                "target": "#table",
                "expectedInputLabel": "Keep the proposed standard replacements",
                "title": "Review proposed residue conversions",
                "body": "OpenMMDL shows non-standard residues only when they are detected. Keep the proposed standard replacements unless you have a specific modeling reason to change them.",
            },
        ],
    },
    "add_heavy_atoms": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Add missing heavy atoms",
        "pageNumber": 5,
        "totalPages": 7,
        "pageKey": "add_heavy_atoms",
        "docsUrl": "pdb_path.html#optional-adding-heavy-atoms",
        "steps": [
            {
                "target": "#table",
                "expectedInputLabel": "Review the heavy atoms to be added",
                "title": "Review atoms that will be rebuilt",
                "body": "This table is informational. PDBFixer will add the listed heavy atoms when you continue.",
            },
            {
                "target": "input[value='Continue']",
                "title": "Continue to solvent setup",
                "body": "The next page sets protonation, water, membrane, and ion options.",
            },
        ],
    },
    "add_hydrogens": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Add hydrogens and a water box",
        "pageNumber": 6,
        "totalPages": 7,
        "pageKey": "add_hydrogens",
        "docsUrl": "pdb_path.html#adding-hydrogens-and-water-box",
        "steps": [
            {
                "target": "#phfield",
                "expectedInputLabel": "Set the protonation pH to 7.4",
                "title": "Set the protonation pH",
                "body": "The PDB tutorial uses pH 7.4 for adding hydrogens.",
                "action": "setTutorialPh",
                "actionLabel": "Set pH to 7.4",
            },
            {
                "target": "#addWater",
                "done": {
                    "selector": "#addWater",
                    "condition": "checked",
                },
                "expectedInputLabel": "Enable the explicit water box",
                "title": "Add an explicit water box",
                "body": "Select Add water box for the solvent tutorial. The ionic strength should remain 0.15 M.",
                "action": "enableTutorialWaterBox",
                "actionLabel": "Add water box",
            },
            {
                "target": "input[value='Continue']",
                "title": "Process the solvated structure",
                "body": "Continue can take some time because OpenMMDL adds atoms, hydrogens, solvent, and ions before moving to simulation options.",
            },
        ],
    },
    "simulation_options": {
        "id": "pdb-small-molecule",
        "title": "PDB small molecule tutorial",
        "pageTitle": "Generate the simulation package",
        "pageNumber": 7,
        "totalPages": 7,
        "pageKey": "simulation_options",
        "docsUrl": "pdb_path.html#simulation-options",
        "steps": [
            {
                "target": "#sim_length",
                "expectedInputLabel": "Set the simulation length to 10 ns",
                "title": "Set the tutorial simulation length",
                "body": "The PDB tutorial uses 10 ns. You can shorten this further for a test run if needed.",
                "action": "setSimulationLength",
                "actionLabel": "Set to 10 ns",
            },
            {
                "target": "button[onclick*=downloadPackage]",
                "expectedInputLabel": "Save all generated files",
                "title": "Save all generated files",
                "body": "Save All Files downloads openmmdl_simulation.zip with everything OpenMMDL Simulation needs: the script OpenMMDL_Simulation.py, the prepared receptor 5wyz_tutorial_receptor-processed_openMMDL.pdb (chains A and B, hydrogens, water box and ions), and the two ligands 7VF_A.sdf (topology code UNK) and 7VF_B.sdf (topology code L01), parametrised with SMIRNOFF in an AMBER19/OPC3 system for a 10 ns run. Unpack it and start the run with: openmmdl simulation -f tutorial_simulation -s OpenMMDL_Simulation.py -t 5wyz_tutorial_receptor-processed_openMMDL.pdb -l 7VF_A.sdf 7VF_B.sdf",
            },
        ],
    },
}

from __future__ import annotations

import argparse
import os


TRUE_BOOL_VALUES = {"1", "true", "t", "yes", "y", "on"}
FALSE_BOOL_VALUES = {"0", "false", "f", "no", "n", "off"}


def parse_bool_flag(value):
    """Parse flexible CLI boolean values like true/false, yes/no, y/n, 1/0."""
    if isinstance(value, bool):
        return value

    normalized_value = str(value).strip().lower()
    if normalized_value in TRUE_BOOL_VALUES:
        return True
    if normalized_value in FALSE_BOOL_VALUES:
        return False

    raise argparse.ArgumentTypeError(
        "Boolean value expected. Use one of: true/false, yes/no, y/n, 1/0, on/off."
    )


LOGO = "\n".join(
    [
        r"     ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      ",
        r"   .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      ",
        r"  / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      ",
        r" ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    ",
        r" |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    ",
        r" : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    ",
        r"  \ `_/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  ",
        r"   '. \_/``'.'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ ",
        r"     '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` ",
        r"              Prepare and Perform OpenMM Protein-Ligand MD Simulations                                 ",
        r"                                     Version 1.3.0                                                     ",
    ]
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openmmdl analysis",
        description=LOGO,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("-t", dest="topology", help="Topology File after MD Simulation", required=True)
    parser.add_argument("-d", dest="trajectory", help="Trajectory File in DCD Format", required=True)
    parser.add_argument("-n", dest="ligand_name", help="Ligand Name (3 Letter Code in PDB)", default=None)
    parser.add_argument("-l", dest="ligand_sdf", help="Ligand in SDF Format", default=None)
    parser.add_argument("-b", dest="binding", help="Binding Mode Treshold for Binding Mode in %%", default=40)
    parser.add_argument(
        "-df",
        dest="dataframe",
        help='Dataframe (use if the interactions were already calculated, default name would be "df_all.csv")',
        default=None,
    )
    parser.add_argument(
        "-f",
        dest="frames",
        help="final frame of the analysis (if you want to analyze only a certain part of the trajectory)",
        default=None,
    )
    parser.add_argument(
        "-m",
        dest="min_transition",
        help="Minimal Transition percentage for Markov State Model",
        default=1,
    )
    parser.add_argument(
        "-c",
        dest="cpu_count",
        help="cores, specify how many PC cores should be used, default is half of the PC cores",
        default=os.cpu_count() // 2,
    )
    parser.add_argument(
        "-p",
        dest="generate_pml",
        help="Generate .pml files for pharmacophore visualization (accepts true/false, yes/no, y/n)",
        default=False,
        type=parse_bool_flag,
    )
    parser.add_argument(
        "-r",
        dest="frame_rmsd",
        help="Calculate RMSD differences between frames (accepts true/false, yes/no, y/n)",
        default=False,
        type=parse_bool_flag,
    )
    parser.add_argument(
        "-nuc",
        dest="receptor_nucleic",
        help="Treat nucleic acids as receptor (accepts true/false, yes/no, y/n)",
        default=False,
        type=parse_bool_flag,
    )
    parser.add_argument(
        "-s",
        dest="special_ligand",
        help="Calculate interactions with special ligands",
        default=None,
    )
    parser.add_argument(
        "-pep",
        dest="peptide",
        help="Calculate interactions with peptides. Give the peptides chain name as input. Defaults to None",
        default=None,
    )
    parser.add_argument(
        "-ref",
        dest="reference",
        help="Add a reference PDB to renumber the residue numbers",
        default=None,
    )
    parser.add_argument(
        "-w",
        dest="stable_water_analysis",
        help="Should stable water analysis be performed? (accepts true/false, yes/no, y/n)",
        default=False,
        type=parse_bool_flag,
    )
    parser.add_argument(
        "-rep",
        dest="representative_frame",
        help="Calculate the representative frame for each binding mode (accepts true/false, yes/no, y/n)",
        default=False,
        type=parse_bool_flag,
    )
    parser.add_argument(
        "--watereps",
        dest="water_eps",
        help="Set the Eps for clustering, this defines how big clusters can be spatially in Angstrom",
        default=1.0,
    )
    parser.add_argument(
        "--figure",
        dest="figure_type",
        help="File type for the figures, default is png. Can be changed to all file types supported by matplotlib.",
        default="png",
    )
    parser.add_argument(
        "--interaction_package",
        dest="interaction_package",
        choices=["plip", "prolif"],
        default="plip",
        help=(
            "Protein-ligand interaction engine. "
            "'plip' uses PLIP for interaction calculation. "
            "'prolif' uses ProLIF for interaction calculation."
        ),
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging with module names and debug messages.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    from openmmdl.openmmdl_analysis.openmmdlanalysis import run_analysis

    return run_analysis(args)


if __name__ == "__main__":
    raise SystemExit(main())
from __future__ import annotations

import argparse


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
        prog="openmmdl simulation",
        description=LOGO,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-f",
        dest="folder",
        type=str,
        help="Folder Name for MD Simulation",
        required=True,
    )
    parser.add_argument(
        "-s",
        dest="script",
        type=str,
        help="MD Simulation script",
        required=True,
    )
    parser.add_argument(
        "-t",
        dest="topology",
        help="Protein Topology PDB/Amber File",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--ligand",
        dest="ligands",
        action="append",
        nargs="+",
        help=(
            "Ligand/cofactor/additional-molecule file(s) in SDF, MOL or MOL2 format. "
            "Repeat -l or pass multiple files after one -l."
        ),
        default=None,
    )
    parser.add_argument(
        "-c",
        dest="coordinate",
        help="Amber coordinates file",
        default=None,
    )
    parser.add_argument(
        "--failure-retries",
        dest="failure_retries",
        type=int,
        default=10,
        help=(
            "Number of reruns if OpenMM fails with 'Particle coordinate is NaN' "
            "(default: 10)"
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    from openmmdl.openmmdl_simulation.openmmdlsimulation import run_simulation

    return run_simulation(args)


if __name__ == "__main__":
    raise SystemExit(main())
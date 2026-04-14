"""
Perform Simulations of Protein-ligand complexes with OpenMM
"""

import argparse
import os
import shutil
import subprocess
import sys

parser = argparse.ArgumentParser()


logo = "\n".join(
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


def main():
    parser = argparse.ArgumentParser(
        prog="openmmdl_simulation",
        description=logo,
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
    parser.add_argument("-t", dest="topology", help="Protein Topology PDB/Amber File", required=True)
    parser.add_argument(
        "-l",
        "--ligand",
        dest="ligands",
        action="append",
        nargs="+",
        help="Ligand/cofactor/additional-molecule file(s) in SDF, MOL or MOL2 format. Repeat -l or pass multiple files after one -l.",
        default=None,
    )
    parser.add_argument("-c", dest="coordinate", help="Amber coordinates file", default=None)
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
    input_formats = [".py", ".pdb", ".sdf", ".mol", ".prmtop", ".inpcrd", ".mol2"]
    args = parser.parse_args()

    ligand_files = []
    if args.ligands is not None:
        for group in args.ligands:
            ligand_files.extend(group)

    if not os.path.exists(args.folder):
        os.mkdir(args.folder)
    else:
        shutil.rmtree(args.folder)
        os.mkdir(args.folder)
    if os.path.exists(args.folder):
        if input_formats[0] in args.script:
            if os.path.exists(args.script):
                shutil.copy(args.script, args.folder)
            else:
                print("Wrong python script path, try the absolute path")
        if input_formats[1] in args.topology:
            if os.path.exists(args.topology):
                shutil.copy(args.topology, args.folder)
            else:
                print("Wrong topology file path, try the absolute path")
        elif input_formats[4] in args.topology:
            if os.path.exists(args.topology):
                shutil.copy(args.topology, args.folder)
            else:
                print("Wrong topology file path, try the absolute path")
        else:
            print("Wrong Format, don't forget the .pdb/.prmtop of the file")
        for ligand_file in ligand_files:
            if any(ligand_file.endswith(ext) for ext in (input_formats[2], input_formats[3], input_formats[6])):
                if os.path.exists(ligand_file):
                    shutil.copy(ligand_file, args.folder)
                else:
                    print(f"Wrong ligand file path, try the absolute path: {ligand_file}")
            else:
                print(f"Wrong ligand format for {ligand_file}, use .sdf, .mol or .mol2")
        if args.coordinate is not None:
            if input_formats[5] in args.coordinate:
                if os.path.exists(args.coordinate):
                    shutil.copy(args.coordinate, args.folder)
                else:
                    print("Wrong coordinates file path, try the absolute path")
            else:
                print("Wrong Format, don't forget the .inpcrd of the coordinate file")
        os.chdir(args.folder)

        script_name = os.path.basename(args.script)
        keep_files = {script_name, os.path.basename(args.topology)}

        for ligand_file in ligand_files:
            keep_files.add(os.path.basename(ligand_file))
        if args.coordinate is not None:
            keep_files.add(os.path.basename(args.coordinate))

        nan_message = "Particle coordinate is NaN"

        for attempt in range(args.failure_retries + 1):
            process = subprocess.Popen(
                [sys.executable, "-u", script_name],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
            )

            output_lines = []
            for line in process.stdout:
                print(line, end="")          # shows output live
                output_lines.append(line)    # also stores it for error check

            process.wait()
            combined_output = "".join(output_lines)

            if process.returncode == 0:
                break

            if nan_message not in combined_output or attempt == args.failure_retries:
                raise SystemExit(process.returncode)

            print(f"Detected OpenMM NaN error. Retrying ({attempt + 1}/{args.failure_retries})...")

            for entry in os.listdir("."):
                if entry not in keep_files:
                    path = os.path.join(".", entry)
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)


if __name__ == "__main__":
    main()


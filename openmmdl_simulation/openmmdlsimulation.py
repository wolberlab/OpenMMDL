"""
mmdl_simulation.py
Perform Simulations of Protein-ligand complexes with OpenMM
"""
import argparse
import sys
import os
import shutil
import argparse
parser = argparse.ArgumentParser()


logo = '\n'.join(["     ,-----.    .-------.     .-''-.  ,---.   .--.,---.    ,---.,---.    ,---. ______       .---.      ",
                  "   .'  .-,  '.  \  _(`)_ \  .'_ _   \ |    \  |  ||    \  /    ||    \  /    ||    _ `''.   | ,_|      ",
                  "  / ,-.|  \ _ \ | (_ o._)| / ( ` )   '|  ,  \ |  ||  ,  \/  ,  ||  ,  \/  ,  || _ | ) _  \,-./  )      ",
                  " ;  \  '_ /  | :|  (_,_) /. (_ o _)  ||  |\_ \|  ||  |\_   /|  ||  |\_   /|  ||( ''_'  ) |\  '_ '`)    ",
                  " |  _`,/ \ _/  ||   '-.-' |  (_,_)___||  _( )_\  ||  _( )_/ |  ||  _( )_/ |  || . (_) `. | > (_)  )    ",
                  " : (  '\_/ \   ;|   |     '  \   .---.| (_ o _)  || (_ o _) |  || (_ o _) |  ||(_    ._) '(  .  .-'    ",
                  "  \ `_/  \  ) / |   |      \  `-'    /|  (_,_)\  ||  (_,_)  |  ||  (_,_)  |  ||  (_.\.' /  `-'`-'|___  ",
                  "   '. \_/``'.'  /   )       \       / |  |    |  ||  |      |  ||  |      |  ||       .'    |        \ ",
                  "     '-----'    `---'        `'-..-'  '--'    '--''--'      '--''--'      '--''-----'`      `--------` ",
                  "              Prepare and Perform OpenMM Protein-Ligand MD Simulations                                 ",
                  "                                     Alpha Version                                                     "])




if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='openmmdl', description=logo, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', dest='folder', type=str, help='Folder Name for MD Simulation', required=True,)
    parser.add_argument('-s', dest='script', type=str, help='MD Simulation script', required=True,)
    parser.add_argument('-p', dest='protein', help='Protein PDB File', required=True)
    parser.add_argument('-l', dest='ligand', help='SDF File of Ligand', default=None)
    input_formats = ['.py', '.pdb', '.sdf'] 
    args = parser.parse_args()
    if not os.path.exists(args.folder):
        os.mkdir(args.folder)
    else:
        shutil.rmtree(args.folder)
        os.mkdir(args.folder)
    script_dir = os.path.abspath( os.path.dirname( __file__ ))
    if os.path.exists(args.folder):
        if input_formats[0] in args.script:
            if os.path.exists(args.script):
                shutil.copy(args.script, args.folder)
            else:
                print("Wrong python script path, try the absolute path")
        if input_formats[1] in args.protein:
            if os.path.exists(args.protein):
                shutil.copy(args.protein, args.folder)
            else:
                print("Wrong pdb file path, try the absolute path")
        else:
            print("Wrong Format, don't forget the .pdb of the pdb file")
        if args.ligand != None:
            if input_formats[2] in args.ligand:
                if os.path.exists(args.ligand):
                    shutil.copy(args.ligand, args.folder)
                else:
                    print("Wrong pdb file path, try the absolute path")
            else:
                print("Wrong Format, don't forget the .pdb of the pdb file")
        shutil.copytree(f"{script_dir}/scripts", f"{args.folder}/scripts")
        os.system(f"python3 {args.folder}/{args.script}")


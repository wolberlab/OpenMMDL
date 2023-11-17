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




def main():
    parser = argparse.ArgumentParser(prog='openmmdl_simulation', description=logo, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', dest='folder', type=str, help='Folder Name for MD Simulation', required=True,)
    parser.add_argument('-s', dest='script', type=str, help='MD Simulation script', required=True,)
    parser.add_argument('-t', dest='topology', help='Protein Topology PDB/Amber File', required=True)
    parser.add_argument('-l', dest='ligand', help='SDF File of Ligand', default=None)
    parser.add_argument('-c', dest='coordinate', help='Amber coordinates file', default=None)
    input_formats = ['.py', '.pdb', '.sdf', '.mol', '.prmtop', '.inpcrd'] 
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
        if args.ligand != None:
            if input_formats[2] in args.ligand or input_formats[3] in args.ligand:
                if os.path.exists(args.ligand):
                    shutil.copy(args.ligand, args.folder)
                else:
                    print("Wrong ligand file path, try the absolute path")
            else:
                print("Wrong Format, don't forget the .sdf of the ligand file")
        if args.coordinate != None:
            if input_formats[5] in args.coordinate:
                if os.path.exists(args.coordinate):
                    shutil.copy(args.coordinate, args.folder)
                else:
                    print("Wrong coordinates file path, try the absolute path")
            else:
                print("Wrong Format, don't forget the .inpcrd of the coordinate file")
        os.chdir(args.folder)
        os.system(f"python3 *.py")
        
if __name__ == "__main__":
    main()


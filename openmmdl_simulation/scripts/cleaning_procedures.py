import os
import shutil

def cleanup(protein_name):
    """
    Cleans up the PDB Reporter Output File and MDTraj Files of the performed simulation

    Parameters
    ----------
    protein_name: str
        Name of the protein pdb.    
    """ 
    print("Cleaning Up :)")
    os.remove(f'output_{protein_name}')
    os.remove(f'centered_old_coordinates.pdb')
    os.remove(f'centered_old_coordinates.dcd')
    print("It's done")
    
def post_md_file_movement(protein_name, ligand=None):
    """
    Moves the files to the required folders

    Parameters
    ----------
    protein_name: str
        Name of the protein pdb.    
    ligand: str
        Name of the ligand sdf.
    """ 
    # Creation of necessary Folders + Movement of Files to Folders
    # Input Files
    
    if not os.path.exists("Input_Files"):
        os.mkdir("Input_Files")    
    else:
        shutil.rmtree("Input_Files")
        os.mkdir("Input_Files")

    if   os.path.exists(f"{ligand}"):
            shutil.copy(f"{ligand}", f"Input_Files/{ligand}")
    shutil.copy(f"{protein_name}", f"Input_Files/{protein_name}")
    
    # MD Simulation Files
    # Pre Simulation Files
    
    if not os.path.exists("MD_Files"):
        os.mkdir("MD_Files")    
    else:
        shutil.rmtree("MD_Files")
        os.mkdir("MD_Files")
        
    if not os.path.exists("MD_Files/Pre_MD"):
        os.mkdir("MD_Files/Pre_MD")    
    else:
        shutil.rmtree("MD_Files/Pre_MD")
        os.mkdir("MD_Files/Pre_MD")
    
    os.rename(f'prepared_no_solvent_{protein_name}', f"MD_Files/Pre_MD/prepared_no_solvent_{protein_name}")
    
    if   os.path.exists(f'solvent_padding_{protein_name}'):
         os.rename(f'solvent_padding_{protein_name}', f'MD_Files/Pre_MD/solvent_padding_{protein_name}')
    elif os.path.exists(f'solvent_absolute_{protein_name}'):
         os.rename(f'solvent_absolute_{protein_name}', f'MD_Files/Pre_MD/solvent_absolute_{protein_name}')
    elif os.path.exists(f'membrane_{protein_name}'):
         os.rename(f'membrane_{protein_name}', f'MD_Files/Pre_MD/membrane_{protein_name}')
    
    
    # Topology files after Minimization and Equilibration
    if not os.path.exists("MD_Files/Minimization_Equilibration"):
        os.mkdir("MD_Files/Minimization_Equilibration")    
    else:
        shutil.rmtree("MD_Files/Minimization_Equilibration")
        os.mkdir("MD_Files/Minimization_Equilibration")
    
    os.rename(f'Energyminimization_{protein_name}', f"MD_Files/Minimization_Equilibration/Energyminimization_{protein_name}")
    os.rename(f'Equilibration_{protein_name}', f"MD_Files/Minimization_Equilibration/Equilibration_{protein_name}")
    
    
    # Simulation Output Files
    if not os.path.exists("MD_Files/MD_Output"):
        os.mkdir("MD_Files/MD_Output")    
    else:
        shutil.rmtree("MD_Files/MD_Output")
        os.mkdir("MD_Files/MD_Output")    
    
    if  os.path.exists(f'output_{protein_name}'):
        os.rename(f'output_{protein_name}', f'MD_Files/MD_Output/output_{protein_name}')

    if  os.path.exists(f'trajectory.dcd'):
        os.rename(f'trajectory.dcd', f'MD_Files/MD_Output/trajectory.dcd')
    
    # Post Processing Files
        
    if not os.path.exists("MD_Postprocessing"):
        os.mkdir("MD_Postprocessing")    
    else:
        shutil.rmtree("MD_Postprocessing")
        os.mkdir("MD_Postprocessing")
    
    # MDtraj Files
    if  os.path.exists(f'centered_old_coordinates.pdb'):
        if not os.path.exists("MD_Postprocessing/1_MDTraj"):
            os.mkdir("MD_Postprocessing/1_MDTraj")    
        else:
            shutil.rmtree("MD_Postprocessing/1_MDTraj")
            os.mkdir("MD_Postprocessing/1_MDTraj")
        os.rename(f'centered_old_coordinates.pdb', f'MD_Postprocessing/1_MDTraj/centered_old_coordinates.pdb')
    
    if  os.path.exists(f'centered_old_coordinates.dcd'):
        os.rename(f'centered_old_coordinates.dcd', f'MD_Postprocessing/1_MDTraj/centered_old_coordinates.dcd')
    
    # MDAnalysis Files
    if not os.path.exists("MD_Postprocessing/2_MDAnalysis"):
        os.mkdir("MD_Postprocessing/2_MDAnalysis")    
    else:
        shutil.rmtree("MD_Postprocessing/2_MDAnalysis")
        os.mkdir("MD_Postprocessing/2_MDAnalysis")

    os.rename(f'centered_top.pdb', f'MD_Postprocessing/2_MDAnalysis/centered_top.pdb')
    os.rename(f'centered_traj.dcd', f'MD_Postprocessing/2_MDAnalysis/centered_traj.dcd')
    os.rename(f'prot_lig_top.pdb', f'MD_Postprocessing/2_MDAnalysis/prot_lig_top.pdb')
    os.rename(f'prot_lig_traj.dcd', f'MD_Postprocessing/2_MDAnalysis/prot_lig_traj.dcd')
    
    # Checkpoints
    if not os.path.exists("Checkpoints"):
        os.mkdir("Checkpoints")    
    else:
        shutil.rmtree("Checkpoints")
        os.mkdir("Checkpoints") 
    
    if os.path.exists(f'checkpoint.chk'):
        os.rename(f'checkpoint.chk', f'Checkpoints/checkpoint.chk')
    if os.path.exists(f'10x_checkpoint.chk'):
        os.rename(f'10x_checkpoint.chk', f'Checkpoints/10x_checkpoint.chk')
    if os.path.exists(f'100x_checkpoint.chk'):
        os.rename(f'100x_checkpoint.chk', f'Checkpoints/100x_checkpoint.chk')  
    
    # Analysis

    if os.path.exists(f'RMSD_over_time.csv'):
        if not os.path.exists("Analysis"):
            os.mkdir("Analysis")    
        else:
            shutil.rmtree("Analysis")
            os.mkdir("Analysis") 

    if os.path.exists(f'RMSD_over_time.csv'):
        os.rename(f'RMSD_over_time.csv', f'Analysis/RMSD_over_time.csv')
        os.rename(f'RMSD_between_the_frames.png', f'Analysis/RMSD_between_the_frames.png')
        os.rename(f'RMSD_over_time.png', f'Analysis/RMSD_over_time.png')
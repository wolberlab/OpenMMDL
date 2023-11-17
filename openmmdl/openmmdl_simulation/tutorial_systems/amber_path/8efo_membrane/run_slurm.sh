#!/bin/bash
#SBATCH --job-name=8efo_8qy_job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=all-gpu	
#SBATCH --gres=gpu:1
#SBATCH --exclude=n10,n9,n8,n17,n19
#SBATCH --array=1-5		# adjust acording to the amount of repilicas needed


# simulation variables
workdir=/home/yuchen/or/OpenMMDL/8efo

# environment variables
eval "$(/home/yuchen/miniconda3/bin/conda shell.bash hook)" 		# set this to where your conda is
conda activate openmmdl   						# change if openmmdl is installed in different env

# execute simulation
cd ${workdir}
mkdir ./${SLURM_ARRAY_TASK_ID} 
cp ./OpenMMDL_Simulation.py ./${SLURM_ARRAY_TASK_ID}/
cp ./system.opc.inpcrd ./${SLURM_ARRAY_TASK_ID}/
cp ./system.opc.prmtop ./${SLURM_ARRAY_TASK_ID}/
cp ./8QY.sdf ./${SLURM_ARRAY_TASK_ID}/
cd ./${SLURM_ARRAY_TASK_ID}/

python3 OpenMMDL_Simulation.py

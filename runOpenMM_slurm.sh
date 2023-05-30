#!/bin/bash
# Leons wrap for openmm-run
# currently quite hacky
# the paths must be absolute paths.


usage() {
        echo ""
        echo "Please make sure all required parameters are given"
        echo "Usage: $0 <OPTIONS>"
        echo "Required Parameters are EITHER:"
        echo "-t <topology>         Path to the input topology file (meaning the protein, e.g. PDB)"
        echo "-f <output_dir>       Path to a directory that will store the results!"
        echo "-s <simulationscript> Path to the simulationscript. This is the python script which is the output of the preparational step"
        echo "OR instead of -f -s -t:"
        echo "-i <inputdir> 	    A directory containing the topology file as a pdb, the script file (python), and possibly a ligand file!"
        echo "OPTIONAL STUFF:"
        echo "-l <ligand> 	    Path to the ligand SDF file in case you are simulating a ligand in the complex!"
        echo "-n <REPLICA> 	    Number of replicas to be run (default = 1)"
        echo "-g <GPUs> 	    specify the GPUs that should be used (default is all-gpu)"
        echo "-c 	 	    If -c is added to the command line, the cn-gpus are used by default"
        echo ""
        echo "an example command using the -i option instead of the -f -l -s -t option could look like: "
        echo "bash runOpenMM.sh -i /mdspace/leon-moveWIP/Ligand-search-stability-check/ALk2R206H-D207Q-backmut/simulation -n 5 -c"
        echo ""
        exit 1
}



while getopts ":n:t:f:i:s:l:g:c" i; do
        case "${i}" in
        t)
                topology=$OPTARG
        ;;
        f)
                output_dir=$OPTARG
        ;;
        s)
                simulationscript=$OPTARG
        ;;
        l)
                ligand=$OPTARG
        ;;
        n)
                REPLICA=$OPTARG
        ;;
        g)
                GPUs=$OPTARG
        ;;
        c)
                GPUs=cn-gpu
        ;;
        i)
                inputfolder=$OPTARG
                topology=$inputfolder/*.pdb
                simulationscript=$inputfolder/OpenMMDL_Simulation.py
                ligand=$inputfolder/*.sdf
                mkdir $inputfolder/outputs
                output_dir=$inputfolder/outputs
        ;;
        esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
echo $GPUs
if [[ "$topology" == "" || "$output_dir" == "" || "$simulationscript" == "" ]] ; then
    usage
fi

if [[ "$GPUs" == "" ]] ; then
    GPUs=all-gpu
fi
echo $GPUs
if [[ "$REPLICA" == "" ]] ; then
    REPLICA=1
fi
# Run AlphaFold with required parameters
cd $output_dir
for counter in $(seq 0 $((REPLICA - 1)))
do
  if [ -d "$counter" ]; then
    rm -r $counter
  fi
  mkdir $counter
  echo "Submitting replica ${counter} ..."
  sbatch -p $GPUs --nodes=1 --gres=gpu:1 $SCRIPT_DIR/SlurmWrap.sh -t $topology -f $output_dir/$counter -s $simulationscript -l $ligand
done


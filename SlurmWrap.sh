#!/bin/bash
# Leons simple wrap for openmm-run

while getopts ":t:f:s:l:" i; do
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
        esac
done
$(python3 ~/OpenMMDL/openmmdl_simulation/openmmdlsimulation.py -t $topology -f $output_dir -s $simulationscript -l $ligand)

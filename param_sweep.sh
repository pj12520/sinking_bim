#!/bin/bash 
#
#PBS -l walltime=120:00:00

# Add module--------------------------------------

module add gcc/4.4.6

# Set the working directory for this job----------

# export RUNDIR="${HOME}/param_sweep2"
cd $RUNDIR
# echo $RUNDIR
# Name of application----------------------------

export APPLICATION="${HOME}/sweep1/impact"

# Execute the code
# echo $INPUT
 $APPLICATION $INPUT




#!/bin/bash 
#
#PBS -l walltime=120:00:00

# Add module--------------------------------------

module add languages/gcc-4.7
module add gnu_builds/gsl-1.16

# Set the working directory for this job----------

# export RUNDIR="${HOME}/param_sweep2"
cd $RUNDIR
# echo $RUNDIR
# Name of application----------------------------

export APPLICATION="${HOME}/sinking_bim/trunk/impact"

# Execute the code
# echo $INPUT
 $APPLICATION $INPUT




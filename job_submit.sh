#!/bin/bash

for i in $(ls -d */); do

    cd ${i}

    for j in $(ls -d */); do

	cd ${j}

	for k in $(ls -d */); do

	    cd ${k}

	    export INPUT=${HOME}/sweep1/${i}/${j}/${k}/dimensionless_input.dat

	    export RUNDIR=${HOME}/sweep1/${i}/${j}/${k}

	    qsub -V ../../../param_sweep.sh
# echo $INPUT
# echo $RUNDIR
	    sleep 10

	    cd ..
	done

	cd ..

    done

    cd ..
   
done
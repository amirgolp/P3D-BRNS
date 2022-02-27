#!/bin/bash
#$ -cwd
#$ -j y
#$ -N bioGrowth
#$ -o $JOB_NAME.$JOB_ID.output
#$ -M amir.golparvar@ufz.de
#$ -m beas
#$ -l h_rt=24:0:00
#$ -l h_vmem=8G
#$ -pe openmpi-orte 4
##$ -l cpu_model=E5-2690v4
module load foss/2018a
#module load OpenFOAM/6
module load OpenFOAM-Extend/4.0
echo "Executing job commands, current working directory is `pwd`"
source $FOAM_BASH
mpiexec -np 4 bioGrowth -parallel

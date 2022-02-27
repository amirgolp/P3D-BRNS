#!/bin/bash
#$ -cwd
#$ -j y
#$ -N corner
#$ -o $JOB_NAME.$JOB_ID.output
#$ -M amir.golparvar@ufz.de
#$ -m beas
#$ -l h_rt=48:0:00
#$ -l h_vmem=4G
#$ -pe openmpi-orte 16
##$ -l cpu_model=E5-2690v4
module load foss/2018a
#ml foss/2019b
#module load OpenFOAM/7
module load OpenFOAM-Extend/4.0
echo "Executing job commands, current working directory is `pwd`"
source $FOAM_BASH
mpiexec -np 16 interAMFoam -parallel


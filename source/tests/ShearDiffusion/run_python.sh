# This will run matlab scripts using this bash script

#!/bin/bash

#PBS -V
#PBS -l nodes=1:ppn=8
#PBS -l walltime=2000:00:00
#PBS -N TGVortexTestError


cd $PBS_O_WORKDIR

python ShearDiffusion.py

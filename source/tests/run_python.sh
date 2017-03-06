# This will run matlab scripts using this bash script

#!/bin/bash

#PBS -V
#PBS -l nodes=1:ppn=24
#PBS -l walltime=20:00:00
#PBS -N TGVortexTestError


cd $PBS_O_WORKDIR

python TGVortexTest.py >> TGVortexTest_out.txt

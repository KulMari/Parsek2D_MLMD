#!/bin/bash

#MSUB -r TestNewInputCentered               # Request name
#MSUB -n 128  #512  #128 #128 #2048 #512                        # Number of tasks to use
#MSUB -T 1800# 72000 #1800 #86300 #72000 #1800                       # Elapsed time limit in seconds
#MSUB -o $SCRATCHDIR/MLMD/TestNewInput/prova/TestNewInputCentered.o            # Standard output. %I is the job id
#MSUB -e $SCRATCHDIR/MLMD/TestNewInput/prova/TestNewInputCentered.e        # Error output. %I is the job id
#MSUB -q standard #xlarge
#MSUB -Q test
#MSUB -@ mariaelena.innocenti@wis.kuleuven.be: begin,end
#MSUB -A ra1254 #new allocation

BRIDGE_MSUB_PWD='/ccc/cont005/home/ra0747/innocenm/Parsek2D-MLMD_DEV/Parsek2D_MLMD'
export BRIDGE_MSUB_PWD
#module load valgrind
set -x
cd ${BRIDGE_MSUB_PWD}
##ccc_mprun ./iPIC3D inputfiles/GEM.inp   #Test1.inp
ccc_mprun ./ParsekEM  INPUTFILE_Mari_NewFormat/inputParsek2D_MLMD_NewFormat 
#ccc_mprun -d ddt ./ParsekEM  INPUTFILE_Mari/inputfile.SmoothTest_1   

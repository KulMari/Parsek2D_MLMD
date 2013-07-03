#!/bin/bash

#MSUB -r TwoLevels_16x16_ddt               # Request name
#MSUB -n 512  #512  #128 #128 #2048 #512                        # Number of tasks to use
#MSUB -T 6000# 72000 #1800 #86300 #72000 #1800                       # Elapsed time limit in seconds
#MSUB -o $SCRATCHDIR/MLMD/Smoothing/TwoLevels_8x8_ddt/TwoLevels_16x16_ddt.o            # Standard output. %I is the job id
#MSUB -e $SCRATCHDIR/MLMD/Smoothing/TwoLevels_8x8_ddt/TwoLevels_16x16_ddt.e        # Error output. %I is the job id
#MSUB -q standard #xlarge
#####MSUB -Q test
#MSUB -@ mariaelena.innocenti@gmail.com: begin,end
#MSUB -A ra1254 #new allocation

BRIDGE_MSUB_PWD='/ccc/cont005/home/ra0747/innocenm/Parsek2D-multilevel_TimeStamping'
export BRIDGE_MSUB_PWD
#module load valgrind
set -x
cd ${BRIDGE_MSUB_PWD}
##ccc_mprun ./iPIC3D inputfiles/GEM.inp   #Test1.inp
###ccc_mprun ./ParsekEM  INPUTFILE_Mari/inputfile.SmoothTest_1
ccc_mprun -d ddt ./ParsekEM  INPUTFILE_Mari/inputfile.SmoothTest_1   

#!/bin/bash

#MSUB -r TestGit               # Request name
#MSUB -n 256  #512  #128 #128 #2048 #512                        # Number of tasks to use
#MSUB -T 1800                       # Elapsed time limit in seconds
#MSUB -o $SCRATCHDIR/MLMD/PRACE/TestSuite/TestGit/TestGit.o            # Standard output. %I is the job id
#MSUB -e $SCRATCHDIR/MLMD/PRACE/TestSuite/TestGit/TestGit.e        # Error output. %I is the job id
#MSUB -q standard #xlarge
####MSUB -Q test
####MSUB -@ mariaelena.innocenti@wis.kuleuven.be: begin,end
####MSUB -A ra1254 #new allocation
#MSUB -A pa1802
#####MSUB -X 


#BRIDGE_MSUB_PWD='/ccc/cont005/home/ra0747/innocenm/Parsek2D-MLMD_DEV/Parsek2D_MLMD'
BRIDGE_MSUB_PWD='.'
export BRIDGE_MSUB_PWD
#module load valgrind
set -x
cd ${BRIDGE_MSUB_PWD}

ccc_mprun ./ParsekEM  INPUTFILE_Mari_NewFormat/PRACE_Input_TestSuite 


#!/bin/bash

#used for PRACE scalasca tests

#MSUB -r ScalascaTest1               # Request name
#MSUB -n 256  #512  #128 #128 #2048 #512                        # Number of tasks to use
#MSUB -T 1800                       # Elapsed time limit in seconds
#MSUB -o $SCRATCHDIR/MLMD/PRACE/ScalascaTest1/ScalascaTest1.o            # Standard output. %I is the job id
#MSUB -e $SCRATCHDIR/MLMD/PRACE/ScalascaTest1/ScalascaTest1.e        # Error output. %I is the job id
#MSUB -q standard #xlarge
####MSUB -Q test
####MSUB -@ mariaelena.innocenti@wis.kuleuven.be: begin,end
####MSUB -A ra1254 #new allocation
#MSUB -A pa1802
#####MSUB -X 


BRIDGE_MSUB_PWD='/ccc/cont005/home/ra0747/innocenm/Parsek2D-MLMD_DEV/Parsek2D_MLMD'
export BRIDGE_MSUB_PWD
#module load valgrind
set -x
cd ${BRIDGE_MSUB_PWD}

export ELG_BUFFER_SIZE=10000000000000   #with 1000000000 intermediate buffer, then crashes
export EPK_FILTER=$HOME/Parsek2D-MLMD_DEV/Parsek2D_MLMD/scala_filter.txt
#export EPK_GDIR= $SCRATCHDIR/MLMD/PRACE/ScalascaTest1
#export EPK_LDIR= $SCRATCHDIR/MLMD/PRACE/ScalascaTest1

module load scalasca/1.4.3
#####module load ddt
##ccc_mprun ./iPIC3D inputfiles/GEM.inp   #Test1.inp
scalasca -analyze -s ccc_mprun ./ParsekEM  INPUTFILE_Mari_NewFormat/PRACE_Input 
#ccc_mprun ./ParsekEM  INPUTFILE_Mari_NewFormat/PRACE_Input
#ccc_mprun -d ddt ./ParsekEM  INPUTFILE_Mari/inputfile.SmoothTest_1   

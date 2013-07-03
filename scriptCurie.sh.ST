#!/bin/bash

#MSUB -r DoubleHarris_ratio12_correct.ST               # Request name
#MSUB -n 512  #128 #128 #2048 #512                        # Number of tasks to use
##MSUB -n 32
#####MSUB -N 2 
#MSUB -T 320 #1800 #1800                      # Elapsed time limit in seconds
#MSUB -o $SCRATCHDIR/MLMD/TESI/DoubleHarris_ratio12_correct/ratio12_corrST.o            # Standard output. %I is the job id
#MSUB -e $SCRATCHDIR/MLMD/TESI/DoubleHarris_ratio12_correct/ratio12_corrST.e       # Error output. %I is the job id
#MSUB -q xlarge
####MSUB -Q test
#MSUB -@ mariaelena.innocenti@gmail.com: begin,end
#MSUB 

BRIDGE_MSUB_PWD='/ccc/cont005/home/ra0747/innocenm/Parsek2D-multilevel_TimeStamping'
export BRIDGE_MSUB_PWD
#module load valgrind
set -x
cd ${BRIDGE_MSUB_PWD}
##ccc_mprun ./iPIC3D inputfiles/GEM.inp   #Test1.inp
ccc_mprun ./ParsekEM INPUTFILE_Mari/inputfile.DoubleHarrisNew_ratio12highmassratio_CORRECT_ST
#ccc_mprun -d ddt ./ParsekEM INPUTFILE_Mari/inputfile.maxwell2grids
#°ccc_mprun ./ParsekEM INPUTFILE_Mari/inputfile.maxwellTimestamp
#ccc_mprun ./ParsekEM INPUTFILE_Mari/inputfile.vacuum-2grids 
#valgrind --tool=memcheck --track-origins=yes  --leak-check=full ccc_mprun ./ParsekEM INPUTFILE_Mari/inputfile.maxwell2grids #put a comment or a new line, otherwise a string is added here and it makes you go into restart
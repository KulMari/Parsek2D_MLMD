#!/bin/bash                                                                     
#MSUB -r MLMD_visu                # Request name                                
#MSUB -n 8 #4 #6 #7 #50                      # Number of tasks to use             
### use dimension of the field * number of species                              
#MSUB -T 6000                      # Elapsed time limit in seconds      
#MSUB -o $HOME/MLMD_visu_python.o              # Standard output. %I is the job id   
#MSUB -e $HOME/MLMD_visu_python.e              # Error output. %I is the job id                 
#MSUB -q standard                                                               
#####MSUB -@ mariaelena.innocenti@gmail.com["begin,end"]                        
######MSUB -Q test                                                                   
#MSUB -A ra1254 #new allocation                                                 
module load python

BRIDGE_MSUB_PWD='/ccc/cont005/home/ra0747/innocenm/MLMD_visu_python'
export BRIDGE_MSUB_PWD

set -x
cd ${BRIDGE_MSUB_PWD}
#ccc_mprun python  read_proc_para_total2D_WRONGRECRATE.py                       
ccc_mprun python proc2vtrMLMD.py
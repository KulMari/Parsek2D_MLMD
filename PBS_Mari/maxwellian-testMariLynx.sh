/Users/Mari/Library/Mail\ Downloads/Lynx\ -\ ExaWiki.pdf 
#PBS -o out.txt

#module load icc-12.1.0
#module load mvapich2
module load hdf5-1.6.8
#module load mpich2-1.5
#module load openmpi
module load openmpi-1.5.5

echo Start: host `hostname`, date `date`
echo Assigned nodes: `cat $PBS_NODEFILE`

cd $PBS_O_WORKDIR

#Mpich2 asking for a total number of process#

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/openmpi-1.5.4/lib:/opt/mpich2-trunk-r9185/lib:/opt/hdf5-1.6.8/lib

echo $LD_LIBRARY_PATH
#export NPROCS=32  #added for the nproc mess
echo $NPROCS

#echo $PBS_NODEFILE
#echo "contents of nodefile is"
#cat $PBS_NODEFILE

#with mpich2-1.5
#mpirun -f $PBS_NODEFILE -n 32 /home/mariai/Parsek_AMR/Parsek_AMR/AMR_In_Progress/Parsek2D-amr/ParsekEM /home/mariai/Parsek_AMR/Parsek_AMR/AMR_In_Progress/Parsek2D-amr/INPUTFILE_Mari/inputfile.maxwell2grids8ppg
#with openmpi
#ulimit -c unlimited
#mpirun -hostfile $PBS_NODEFILE -n 16 -x LD_LIBRARY_PATH ../ParsekEM ../INPUTFILE_Mari/inputfile.maxwell2grids8ppg
mpirun -hostfile $PBS_NODEFILE -n 8 -x LD_LIBRARY_PATH ../ParsekEM ../INPUTFILE_Mari/inputfile.maxwell1gridArnaud

#mpirun -hostfile $PBS_NODEFILE -n 32 -x LD_LIBRARY_PATH ../ParsekEM ../INPUTFILE_Mari/inputfile.maxwell1gridArnaudFGI
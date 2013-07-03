
#module load icc-12.1.0
##module load mvapich2
module load hdf5-1.8.9
##module load mpich2-1.5
module load openmpi

echo Start: host `hostname`, date `date`
echo Assigned nodes: `cat $PBS_NODEFILE`

cd $PBS_O_WORKDIR

#Mpich2 asking for a total number of process#
#mpirun -f $PBS_NODEFILE -n 16 /home/beck/Parsek2D-multilevel/ParsekEM /home/beck/Parsek2D-multilevel/inputfiles/inputfile.maxwell1grid
#Mpich2 asking for a number of process per node
##mpirun -f $PBS_NODEFILE -ppn 8 /home/beck/Parsek2D-multilevel/ParsekEM /home/beck/Parsek2D-multilevel/inputfiles/inputfile.maxwellian-test-lynx
#Openmpi asking for a total number of process
mpirun -hostfile $PBS_NODEFILE -n 32 -x LD_LIBRARY_PATH /home/beck/Parsek2D-multilevel/ParsekEM /home/beck/Parsek2D-multilevel/inputfiles/inputfile.maxwell2grids

##echo $PATH

##NPROCS=`wc -l < $PBS_NODEFILE`
##n_node=$(cat $PBS_NODEFILE | uniq | wc -l)
##NODEFILE_UNIQ=/tmp/`basename ${PBS_NODEFILE}`.uniq
##cat $PBS_NODEFILE | uniq > $NODEFILE_UNIQ
### Boot the MPI2 engine.
##mpdboot --rsh=ssh -n $n_node --file=${NODEFILE_UNIQ} --verbose
####mpdboot --rsh=ssh -n $n_node --file=$PBS_NODEFILE --verbose

##echo $NPROCS

##mpiexec -l -n $NPROCS -env I_MPI_DEBUG 4 /user/leuven/304/vsc30483/Parsek2D-amr/ParsekEM /user/leuven/304/vsc30483/Parsek2D-amr/inputfiles/inputfile.maria-elena_firststeps
##mpiexec -l -n $NPROCS /user/leuven/304/vsc30483/Parsek2D-amr/ParsekEM /user/leuven/304/vsc30483/Parsek2D-amr/inputfiles/inputfile.maria-elena_firststeps
##mpdallexit


#module load icc-12.1.0
#module load mvapich2
module load mpich2-1.5
module load hdf5-1.6.8

cd $PBS_O_WORKDIR

mpirun -f $PBS_NODEFILE -ppn 8 ./ParsekEM inputfiles/inputfile.maxwellian-test-lynx


import os
import scipy
from pylab import *
import tables
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD
rankproc = comm.Get_rank() # Proc rank
numproc = comm.Get_size() # Total number of procs
figure(0,figsize=(8, 12))

######### INPUT THE FOLLOWING PARAMETERS ############################################
directory = "/scratch/leuven/304/vsc30483/results/maxwell-amr/" #Directory of the simulation data
####################################################################################


if (matplotlib.pyplot.isinteractive()):
    matplotlib.pyplot.ioff()

###Gathering in formation and initialization
setting_filename=directory + "settings.hdf" # Setting file name
h5setting_file = tables.openFile(setting_filename, mode = "r", title = "Setting_file")
collective_group = h5setting_file.root.collective
topology_group = h5setting_file.root.topology
Nprocs = topology_group._f_getChild("Nprocs").read()[0].astype('int32')
dx = collective_group._f_getChild("Dx").read()[0]
dy = collective_group._f_getChild("Dy").read()[0]
Nxc = collective_group._f_getChild("Nxc").read()[0]
Nyc = collective_group._f_getChild("Nyc").read()[0]
Ncycle = collective_group._f_getChild("Ncycles").read()[0]
Ngrids = collective_group._f_getChild("Ngrids").read()[0]
Ratio = collective_group._f_getChild("Ratio").read()[0]
XLEN = topology_group._f_getChild("XLEN").read()[0]
YLEN = topology_group._f_getChild("YLEN").read()[0]
h5setting_file.close()




proc_filename = directory + "proc"+repr(rankproc)+".hdf" # Proc file name
h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
list_group = h5proc_file.root.fields._f_listNodes("Group")
list_array=list_group[0]._f_listNodes("Array") 

n_time=len(list_array)
shape=scipy.shape(list_array[0].read()[:,:,0])
h5proc_file.close()

list_cycle=range(n_time)

Nbyproc=(XLEN*YLEN)/numproc #Minimum Nomber of Y you have to read
reste=XLEN*YLEN-Nbyproc*numproc
print Nbyproc, reste
if (rankproc<reste):
    Nrange=scipy.arange(rankproc*(Nbyproc+1),(rankproc+1)*(Nbyproc+1))
else:
    Nrange=scipy.arange((reste*(Nbyproc+1)+(rankproc-reste)*Nbyproc),(reste*(Nbyproc+1)+(rankproc-reste+1)*Nbyproc))

list_read_proc=scipy.arange(Nrange[0],Nrange[-1]+1)

for grid in range(Ngrids):
    Vcell = dx * dy / Ratio**(2*grid)
    Eelectric = scipy.zeros((Ncycle))
    Emagnetic = scipy.zeros((Ncycle))
    Ekinetic = scipy.zeros((Ncycle))
    list_read_proc=scipy.array(list_read_proc)+ grid*Nprocs
    print list_read_proc
    
    
            ### Start reading #######
    for read_proc in list_read_proc:
        proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name
        
        h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
        
        list_energy = h5proc_file.root.energy.electric._f_listNodes()
        for energy_cycle in list_energy:
            num = repr(energy_cycle).split()[0][23:]
            Eelectric[int(num)] += energy_cycle.read()[0]
        list_energy = h5proc_file.root.energy.magnetic._f_listNodes()
        for energy_cycle in list_energy:
            num = repr(energy_cycle).split()[0][23:]
            Emagnetic[int(num)] += energy_cycle.read()[0]
        list_energy = h5proc_file.root.energy.kinetic.species_0._f_listNodes()
        for energy_cycle in list_energy:
            num = repr(energy_cycle).split()[0][32:]
            Ekinetic[int(num)] += energy_cycle.read()[0]
        list_energy = h5proc_file.root.energy.kinetic.species_1._f_listNodes()
        for energy_cycle in list_energy:
            num = repr(energy_cycle).split()[0][32:]
            Ekinetic[int(num)] += energy_cycle.read()[0]

        h5proc_file.close()

    Eelectric = Vcell * Eelectric[Eelectric > 0]
    Emagnetic = Vcell * Emagnetic[Emagnetic > 0]
    Ekinetic  = Ekinetic [Ekinetic  > 0]
    count = scipy.size(Eelectric)
    print "grid = ", grid, " count =", count
    if rankproc  > 0 :
        comm.Reduce([Eelectric,count,MPI.DOUBLE],None, op = MPI.SUM,root=0)
        comm.Reduce([Emagnetic,count,MPI.DOUBLE],None, op = MPI.SUM,root=0)
        comm.Reduce([Ekinetic,count,MPI.DOUBLE],None, op = MPI.SUM,root=0)
    if rankproc == 0 :               
        Eelectricrecv = scipy.zeros((count))
        Emagneticrecv = scipy.zeros((count))
        Ekineticrecv = scipy.zeros((count))
        comm.Reduce([Eelectric,count,MPI.DOUBLE],[Eelectricrecv,count,MPI.DOUBLE], op = MPI.SUM,root=0)
        comm.Reduce([Emagnetic,count,MPI.DOUBLE],[Emagneticrecv,count,MPI.DOUBLE], op = MPI.SUM,root=0)
        comm.Reduce([Ekinetic,count,MPI.DOUBLE],[Ekineticrecv,count,MPI.DOUBLE], op = MPI.SUM,root=0)
        norm = Eelectricrecv[0]+Emagneticrecv[0]+Ekineticrecv[0]
        figure(0)
        subplot(Ngrids,1,grid+1)
        semilogy(Eelectricrecv/norm,label="Electric Energy",linewidth=2)
        semilogy(Emagneticrecv/norm,label="Magnetic Energy",linewidth=2)
        semilogy(Ekineticrecv/norm,label="Kinetic Energy",linewidth=2)
        title("grid "+ repr(grid))
        legend()
        ylim(10**(-8),10)
        
# End of grid loop
if rankproc == 0 :
    figure(0)

    savefig(directory+'Energy.png',format="png")
print "Done" 

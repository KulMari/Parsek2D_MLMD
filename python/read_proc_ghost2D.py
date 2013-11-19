import os
import scipy
from pylab import *
##import matplotlib
##from matplotlib import pyplot
import tables
import sys
from mpi4py import MPI

comm = MPI.COMM_WORLD
rankproc = comm.Get_rank() # Proc rank
numproc = comm.Get_size() # Total number of procs
figure(0,figsize=(8, 12))

######### INPUT THE FOLLOWING PARAMETERS ############################################
directory = "/home/beck/results/maxwell2grids/" #Directory of the simulation data
list_field=['B'] # List of the fields you want to extract
cycle_step = 1 # step between two displayed cycles
valuemax = 4e-3 #max threshold. Put 0 for non threshold.
level = 0 #Grid level you want to analyze
outputcycle = 1 # outputcycle as input in the inputfile
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
print "shape = ", shape
onebuffersize = 2*shape[0]+2*shape[1]+4
h5proc_file.close()


list_cycle=range(n_time)[::cycle_step]

list_read_proc=scipy.arange(rankproc,XLEN*YLEN,numproc)

grid=1
list_read_proc=scipy.array(list_read_proc)+ level*grid*Nprocs

Fields_local = scipy.zeros((3,shape[0]+2,shape[1]+2),dtype="float32")

field=list_field[0]
if(field=='E'):
    init=3#To keep E fields only
if(field=='B'):
    init=0#To keep B field only


for cycle in list_cycle:    
    MPI.COMM_WORLD.Barrier()
    ### Start reading #######
    for read_proc in list_read_proc:
        proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name
        
        h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
        list_group = h5proc_file.root.fields._f_listNodes("Group")
        list_array=list_group[0]._f_listNodes("Array") 
        
        list_group=list_group[init:init+3]# Keep only E or B depending on init=3 or 0.
        k=0 #Counts the 3 components
        for group in list_group:
            list_array=group._f_listNodes("Array") #List of all recorded arrays in "group"
            Fields_local[k,1:-1,1:-1]=list_array[cycle].read()[:,:,0]
            num = repr(list_array[cycle]).split()[0][17:]
            ## Filling the ghost nodes ##
            if k==0 and int(num)>0: #Because outputghost gives Ez
                ghostcycle=int(num)-1 #Correction needed when outputghost is done AFTER calculate field
                ghostbuffer = scipy.fromfile(directory+"ghost"+repr(read_proc),sep="\n")[ghostcycle/outputcycle*onebuffersize:(ghostcycle/outputcycle+1)*onebuffersize]
                Fields_local[k,:,0] = ghostbuffer[:shape[0]+2]
                Fields_local[k,:,-1] = ghostbuffer[shape[0]+2:2*shape[0]+4]
                Fields_local[k,0,1:-1] = ghostbuffer[2*shape[0]+4:2*shape[0]+shape[1]+4]
                Fields_local[k,-1,1:-1] = ghostbuffer[2*shape[0]+shape[1]+4:2*shape[0]+2*shape[1]+4]

            #  Plot single frames 
            figure(1)
            if valuemax > 0. :
                #imshow(Fields_local[k,:,:].T,origin=0,vmin = 0,vmax=valuemax) 
                pcolor(Fields_local[k,:,:].T,vmin = -valuemax,vmax=valuemax) 
            else:
                imshow(Fields_local[k,:,:].T,origin=0) 
            colorbar(format='%5.3e')
            xlabel('x', fontsize=10)
            ylabel('y', fontsize=10)
            if k==0:
                indice='x'
            if k==1:
                indice='y'
            if k==2:
                indice='z'
            Title = field+indice+" cycle "+num+ " proc "+repr(read_proc)
            title(Title, fontsize=10) # Titre normal
            savefig(directory+field+indice+repr(read_proc)+"_cycle"+num.rjust(8,'0')+'.png',format="png")
            close(1)
            k=k+1

        h5proc_file.close()
   #End of read_proc loop 


#End of cycle loop
if rankproc ==0:
    for read_proc in range(level*XLEN*YLEN,XLEN*YLEN*(level+1)):
        os.system("convert -delay 10  -loop 0 "+directory+field+"x"+repr(read_proc)+"_*.png "+ directory+field+"x"+repr(read_proc)+".gif")
        os.system("convert -delay 10  -loop 0 "+directory+field+"y"+repr(read_proc)+"_*.png "+ directory+field+"y"+repr(read_proc)+".gif")
        os.system("convert -delay 10  -loop 0 "+directory+field+"z"+repr(read_proc)+"_*.png "+ directory+field+"z"+repr(read_proc)+".gif")
    os.system("rm -f "+directory+"*"+"*_cycle*.png")
print "Done" 



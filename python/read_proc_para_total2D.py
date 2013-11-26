import os
import scipy
import matplotlib
matplotlib.use('Agg')
#from pylab import *
from matplotlib import pyplot
import tables
import sys
from mpi4py import MPI
import math
print "Start"
sys.stdout.flush()

comm = MPI.COMM_WORLD
rankproc = comm.Get_rank() # Proc rank
numproc = comm.Get_size() # Total number of procs

######### INPUT THE FOLLOWING PARAMETERS ############################################
directory = "/home/arnaud/data/results/mlmd/vacuum-2grids_withoutRTC/" #Directory of the simulation data
list_field=['E','B'] #List of the fields you want to extract (E,B,J,rhoN,PN where N is the number of the species)
first_cycle = 0
last_cycle = 495
cycle_step = 5 # step between two displayed cycles (minimum is given by outputcyle of the inputfile)
valuemax = 0 #max threshold. Put 0 for non threshold.
####################################################################################

print "Start reading collectives"
sys.stdout.flush()

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
numgrid = rankproc*Ngrids/numproc #Number of the grid this proc is working on
gridcomm = MPI.COMM_WORLD.Split(numgrid, 1)
rankproc = gridcomm.Get_rank() # Proc rank in the new comm
numproc = gridcomm.Get_size() # Total number of procs in the new comm

##Parameters for som theoretical comparisons ##
##Lx=1.+2*dx           
##Ly=1.+2*dy           
##dt = 0.01125
##E0 = 0.5
##kx = math.pi/2./Lx
##ky = math.pi/2./Ly
##k_tot= scipy.sqrt(kx**2+ky**2)
##startt=0.25/k_tot
##x=scipy.linspace(0,Lx,Nxc+3)
##y=scipy.linspace(0,Ly,Nyc+3)
##X,Y=scipy.meshgrid(x,y)
##############################################

proc_filename = directory + "proc"+repr(rankproc)+".hdf" # Proc file name
h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
list_group = h5proc_file.root.fields._f_listNodes("Group")
list_array=list_group[0]._f_listNodes("Array") 

n_time=len(list_array)
shape=scipy.shape(list_array[0].read()[:,:,0])
h5proc_file.close()


list_cycle=range(first_cycle,last_cycle+1,cycle_step) 
#singlepointevolution = scipy.zeros(Ncycle)
#singlepointtheory = scipy.zeros(Ncycle)

if rankproc < 6:
    print "Reading directory " + directory
    Fields_global=scipy.zeros((Nxc+1,Nyc+1),dtype="float32") #This assumes same size for each proc TEMPORARY !!
    theory = Fields_global
    print "Field_global size ", scipy.shape(Fields_global)
else:
    Fields_global=None
num_send = scipy.zeros((1)).astype(int)
sys.stdout.flush()

Xbyproc=XLEN/numproc #Minimum Nomber of X you have to read
reste=XLEN-Xbyproc*numproc
if (rankproc<reste):
    Xrange=scipy.arange(rankproc*(Xbyproc+1),(rankproc+1)*(Xbyproc+1))
else:
    Xrange=scipy.arange((reste*(Xbyproc+1)+(rankproc-reste)*Xbyproc),(reste*(Xbyproc+1)+(rankproc-reste+1)*Xbyproc))

xread_proc=len(Xrange) #Number of proc file you want to read in the y direction
list_read_proc=scipy.zeros((len(Xrange)*YLEN))
if len(Xrange) > 0:
    list_read_proc=scipy.arange(Xrange[0]*YLEN,(Xrange[-1]+1)*YLEN)

list_read_proc_init = list_read_proc

grid = numgrid
list_read_proc=scipy.array(list_read_proc_init)+ numgrid*Nprocs

Fields_local = scipy.zeros((6,xread_proc*(shape[0]-1)+1,Nyc+1),dtype="float32")

for field in list_field:
    
    if   (field=='B'):
        init=0
        K = 3 #Dimensionality
        quantitytype = "fields"
        quantity = "B"
    elif (field=='E'):
        init=3
        K = 3
        quantitytype = "fields"
        quantity = "E"
    elif (field[0]=="r"):
        init=9
        K=1
        quantitytype = "moments"
        quantity = "rho"
        name = 'species_'+field[3]
    elif (field[0]=="J"):
        init=0
        K=3
        quantitytype = "moments"
        quantity = "J"
        name = 'species_'+field[1]
    elif (field[0]=="P"):
        init=3
        K=6
        quantitytype = "moment"
        species = field[1]
    else:
        sys.exit("Field not supported")
    
    for cycle in list_cycle:    
        print field, cycle
        num = repr(cycle)
        namecycle = 'cycle_'+num
        sys.stdout.flush()

        gridcomm.Barrier()
        ### Start reading #######
        for read_proc in list_read_proc:
            proc_filename = directory + "proc"+repr(read_proc)+".hdf" # Proc file name
            
            h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
            if quantitytype == "moments":
                base = h5proc_file.root.moments._f_getChild(name)
            else: 
                base = h5proc_file.root.fields
            if quantity == 'rho':
                list_group=[base.rho]
            elif K==3:
                list_group=[base._f_getChild(quantity+'x'),base._f_getChild(quantity+'y'),base._f_getChild(quantity+'z')]
                
            p_coordinate=h5proc_file.root.topology.cartesian_coord.read() 
            p_coordinate[0]=p_coordinate[0]-(list_read_proc[0]-numgrid*Nprocs)/YLEN
            x_start=p_coordinate[0]*(shape[0]-1) #-1 to account for the 1 node overlapping
            x_end=x_start+shape[0]
            y_start=p_coordinate[1]*(shape[1]-1)
            y_end=y_start+shape[1]
            k=0 #Counts the components
            for group in list_group:
                Fields_local[k,x_start:x_end,y_start:y_end]=group._f_getChild(namecycle).read()[:,:,0]
                k=k+1
            h5proc_file.close()
        #End of read_proc loop 
        #Reassemble data from each processors on procs 0,1 and 2 ####
        for k in range(K):
            if (rankproc > 0 ):
                buffer=Fields_local[k,1:,:].astype("float32") #Account for overlapping
            else:
                buffer=Fields_local[k,:,:].astype("float32") 
            vec_count=scipy.arange(numproc,dtype='int32') #Initialize for all procs
            matsize=scipy.int32(buffer.size)
            gridcomm.Allgather([matsize,1,MPI.INT],[vec_count,1,MPI.INT])
            vec_stride = [vec_count[:j].sum() for j in range(numproc)]

            gridcomm.Gatherv([buffer,buffer.size,MPI.FLOAT],[Fields_global,vec_count,vec_stride,MPI.FLOAT],root=k)
        #End of the field gathering loop ############
        #################################################################
        #  Plot single frames by procs 0,1 and 2 #######################
        if rankproc < K  :
            if rankproc==0 and K==3:
                indice='x'
            elif rankproc==1 and K==3:
                indice='y'
            elif rankproc ==2 and K==3:
                indice='z'
                #if field == 'B' and grid == 0:
                #    singlepointevolution[int(num)]=Fields_global[25,25]#point [5,5] is taken arbitrarily
                #    time=int(num)*dt+startt
                #    singlepointtheory[int(num)]= E0*k_tot*scipy.cos(kx*26*dx)*scipy.sin(ky*26*dy)*scipy.cos(k_tot*time)/kx
                #    print "num = ", int(num)
                #    theory = E0*k_tot*scipy.cos(kx*X)*scipy.sin(ky*Y)*scipy.cos(k_tot*time)/kx 
            elif rankproc==0 and K==6:
                indice='xx'
            elif rankproc==1 and K==6:
                indice='xy'
            elif rankproc ==2 and K==6:
                indice='xz'
            elif rankproc==3 and K==6:
                indice='yy'
            elif rankproc==4 and K==6:
                indice='yz'
            elif rankproc ==5 and K==6:
                indice='zz'
            else:
                indice = ''
                    

            h5proc_file = tables.openFile(proc_filename, mode = "r", title = "Proc_file")
            Ox=h5proc_file.root.topology.Ox.read()[0] 
            Oy=h5proc_file.root.topology.Oy.read()[0] 
            h5proc_file.close()
            #difference = Fields_global.T - theory
            #pyplot.figure(1,figsize=(8, 12))
            pyplot.figure()
            if valuemax > 0. :
                pyplot.imshow(Fields_global.T,origin=0,extent=[Ox,dx*Nxc/Ratio**numgrid+Ox,Oy,dy*Nyc/Ratio**grid+Oy],vmin = -valuemax,vmax=valuemax) 
            else:
                pyplot.imshow(Fields_global[:,:].T,origin=0,extent=[Ox,dx*Nxc/Ratio**grid+Ox,Oy,dy*Nyc/Ratio**grid+Oy]) 
               #pyplot.imshow(Fields_global[:,:].T,origin=0,extent=[Ox,dx*Nxc/Ratio**grid+Ox,Oy,dy*Nyc/Ratio**grid+Oy],vmin = valuemin,vmax=valuemax) 
                 # Fields_global[15:-16,15:-16] = 0
                 # pcolor(Fields_global[:,:].T) 
            pyplot.colorbar(format='%5.3e')
            pyplot.xlabel('x', fontsize=10)
            pyplot.ylabel('y', fontsize=10)
            Title = field+indice+" cycle "+num
            pyplot.title(Title, fontsize=10) # Titre normal
            pyplot.savefig(directory+field+indice+"grid_"+repr(numgrid)+"_cycle"+num.rjust(8,'0')+'.png',format="png")
            pyplot.close(1)

    #End of cycle loop
    if rankproc < K:
        os.system("convert -delay 10  -loop 0 "+directory+field+indice+"grid_"+repr(numgrid)+"*.png "+ directory+field+indice+"_grid"+repr(numgrid)+".gif")
        os.system("rm -f "+directory+"*"+indice+"grid_"+repr(numgrid)+"_cycle*.png")
# End of field loop      
#if rankproc == 2:
#    singlepointevolution[singlepointevolution!=0].tofile(directory+"singlepointdata.txt",sep="\n")        
#    singlepointtheory[singlepointtheory!=0].tofile(directory+"singlepointtheory.txt",sep="\n")        
print "Done" 



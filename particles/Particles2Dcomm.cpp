/*******************************************************************************************
Particles2Dcomm.cpp  -  Class for particles of the same species, in a 2D space and 3component velocity
-------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess, Maria Elena Innocenti
********************************************************************************************/

#include "mpi.h"
#include <iostream>
#include <math.h>
#include <float.h>   // for DBL_EPSILON
#include "../processtopology/VirtualTopology.h"
#include "../processtopology/VCtopology.h"
#include "../inputoutput/CollectiveIO.h"
#include "../inputoutput/Collective.h"
#include "../communication/ComParticles.h"
#include "../utility/Alloc.h"
#include "../mathlib/Basic.h"
#include "../mathlib/Bessel.h"
#include "../bc/BcParticles.h"
#include "../grids/Grid.h"
#include "../grids/Grid2DCU.h"
#include "../fields/Field.h"

#include "Particles2Dcomm.h"

#include "hdf5.h"
#include <vector>
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-32
/**
 * 
 * Class for particles of the same species, in a 2D space and 3component velocity
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Enrico Camporeale, Enrico Camporeale, David Burgess
 * @version 2.0
 *
 */

/** constructor */
Particles2Dcomm::Particles2Dcomm(){
  // see allocate(int species, CollectiveIO* col, VirtualTopology* vct, Grid* grid)
	
}
/** deallocate particles */
Particles2Dcomm::~Particles2Dcomm(){
  delete[] x;
  delete[] y;
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] q;
  delete[] xptilde;
  delete[] yptilde;
  delete[] uptilde;
  delete[] vptilde;
  delete[] wptilde;
  // deallocate buffers
  delete[] b_XDX;
  delete[] b_XSN;
  delete[] b_YDX;
  delete[] b_YSN;
  delete[] MIN_VAL_VEC_COMM;

  // AMR, ME
  // deallocate buffer for particle repop
  delete[] MIN_VAL_VEC;

    
  delete[] REPOP_b_BOTTOM;
  delete[] REPOP_b_TOP;
  delete[] REPOP_b_LEFT;
  delete[] REPOP_b_RIGHT;

  delete []targetBC;
  delete []BCSide;
  delete []fromBC;
  delete []BCSidecu;
  // end these buffers do not exist for the finest grid 
    
  // this buffer exists only if level >0
  delete[] REPOP_receive_b;
  delete[] SplittedParticles_Comm_BOTTOM;
  delete[] SplittedParticles_Comm_TOP;
  delete[] SplittedParticles_Comm_LEFT;
  delete[] SplittedParticles_Comm_RIGHT;
  delete[] MIN_VAL_VEC_SP;

  delete [] OSParticles_Comm_BOTTOM;
  delete [] OSParticles_Comm_TOP;
  delete [] OSParticles_Comm_LEFT;
  delete [] OSParticles_Comm_RIGHT;
  delete [] MIN_VAL_VEC_OS;
  // end this buffer exists only if level >0
  // this buffer exists only if the level is a coarse one
  delete [] AlreadyAccumulated;

  delete[] RP_x;
  delete[] RP_y;
  delete[] RP_u;
  delete[] RP_v;
  delete[] RP_w;
  delete[] RP_q;
  delete[] RP_ParticleID;
 }
/** constructors fo a single species*/
void Particles2Dcomm::allocate(int species, CollectiveIO* col, VirtualTopology* vct, Grid* grid){

  bool PrintSize= false; // memory prints useful for estimating memory consumption

  // info from collectiveIO
  ns = species;
  npcel  = col->getNpcel(species);
  npcelx = col->getNpcelx(species);
  npcely = col->getNpcely(species);
  nop   = col->getNp(species)/(vct->getNprocs());  
  //now, the division by the number of cores is done in Collective already to prevent the int to go out of boundaries
  //npmax =  col->getNpMax(species)/(vct->getNprocs());  
  npmax =  (int) col->getNpMax(species); 
  //cout <<"col->getNpMax(species) " << col->getNpMax(species)  << " npmax " << npmax <<endl;
  //npmax =  2*col->getNpMax(species); // for maxwellian test
  qom   = col->getQOM(species);
  uth   = col->getUth(species);
  vth   = col->getVth(species);
  wth   = col->getWth(species);
  u0    = col->getU0(species);
  v0    = col->getV0(species);
  w0    = col->getW0(species);
  dt    = col->getDt();
  int level= grid->getLevel();
  double ratio= grid->getRatio();
  int TimeRatio= col->getTimeRatio();
  SubCycling= col->getSubCycling();
  if (level and SubCycling)
    dt= dt/(double)pow(TimeRatio,level);
  //cout << "level " << level << ", dt " <<dt <<endl;

  Lx     = col->getLx()/pow(col->getRatio(),grid->getLevel());
  Ly     = col->getLy()/pow(col->getRatio(),grid->getLevel());
  delta  = col->getDelta();
  TrackParticleID =col-> getTrackParticleID(species);
  c = col->getC();
  // info for mover
  NiterMover = col->getNiterMover();
  // velocity of the injection from the wall
  Vinj = col->getVinj();
  // info from Grid
  xstart = grid->getXstart();
  xend   = grid->getXend();
  ystart = grid->getYstart();
  yend   = grid->getYend();
  // AMR, ME
  Modified_xstart= grid->getmodifiedXstart(vct);
  Modified_xend= grid->getmodifiedXend(vct);
  Modified_ystart= grid->getmodifiedYstart(vct);
  Modified_yend= grid->getmodifiedYend(vct);
  // end AMR, ME
  nxn = grid->getNXN();
  nyn = grid->getNYN();
  dx  = grid->getDX();
  dy  = grid->getDY();
  invVOL = grid->getInvVOL();
  // info from VirtualTopology
  cVERBOSE = vct->getcVERBOSE();
    
  // boundary condition for particles
  bcPfaceXright = col->getBcPfaceXright();
  bcPfaceXleft = col->getBcPfaceXleft();
  bcPfaceYright = col->getBcPfaceYright();
  bcPfaceYleft = col->getBcPfaceYleft();
    
  ////////////////////////////////////////////////////////////////
  ////////////////     ALLOCATE ARRAYS   /////////////////////////
  ////////////////////////////////////////////////////////////////

  //if (vct->getCartesian_rank_COMMTOTAL() == 0 || vct->getCartesian_rank_COMMTOTAL()==3000)
  //  cout <<"Starting allocating particle vectors" <<endl;

  // for the mover                                                                        
  allocArr3(&XN, nxn, nyn, 1);
  allocArr3(&YN, nxn, nyn, 1);
  XN= grid->xn;
  YN=grid->yn;


  // positions
  x = new double[npmax];
  y = new double[npmax];
  xptilde = new double[npmax];
  yptilde = new double[npmax];
  // velocities
  u = new double[npmax];
  v = new double[npmax];
  w = new double[npmax];
  uptilde = new double[npmax];
  vptilde = new double[npmax];
  wptilde = new double[npmax];
  // charge
  q = new double[npmax];
  // MLMD
  // for coarse particles collection method, with PRACollectionMethod= 0  
  if (grid->getLevel()< col->getNgrids()-1  && PRACollectionMethod==0) //needed only for one (PRACollectionMethod==0 coarse particles collction method)                                                                                                                                 
    {
      AlreadyAccumulated = new bool[npmax];
    }
  // end MLMD
  //ID
  if (TrackParticleID){
    ParticleID= new unsigned long[npmax];
    BirthRank[0]=vct->getCartesian_rank();
    if (vct->getNprocs()>1) 
      BirthRank[1]= (int) ceil(log10((double) (vct->getNprocs()))); // Number of digits needed for # of process in ID
    else BirthRank[1]=1;
    if (BirthRank[1]+ (int) ceil(log10((double) (npmax)))>10 && BirthRank[0] == 0 ) {
      cerr<< "Error: can't Track particles in Particles2Dcomm::allocate"<<endl;
      cerr<< "Unsigned long 'ParticleID' cannot store all the particles"<<endl;
      return ;
    }
  }

  //Debug1-1
  if (PrintSize)
    {
      MPI_Barrier(vct->getCART_COMM());
      if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN()))   )
	cout <<"Each core of level " << grid->getLevel() <<" allocated 13 vectors of size npmax= " <<npmax << " for basic particle info"<<endl; 
      if (vct->getCartesian_rank_COMMTOTAL() == 0 || vct->getCartesian_rank_COMMTOTAL()==3000)
	cout <<"Each core of level " << grid->getLevel() <<" allocated 13 vectors of size npmax= " <<npmax << " for basic particle info"<<endl;
    }

  // BUFFERS
  // the buffer size should be decided depending on number of particles
  if (TrackParticleID)
    nVar=12;
  else 
    nVar=11;
  
  buffer_size=(int) (0.05*nop*nVar+1); 
  MAX_BUFFER_SIZE = (int)(1*nop*nVar+1); //new

  // if they try to resize beyond this, the simulation is terminated
  b_XDX = new double[buffer_size];
  b_XDX_ptr = b_XDX; // alias to make the resize
  b_XSN = new double[buffer_size];
  b_XSN_ptr = b_XSN; // alias to make the resize
  b_YDX = new double[buffer_size];
  b_YDX_ptr = b_YDX; // alias to make the resize
  b_YSN = new double[buffer_size];
  b_YSN_ptr = b_YSN; // alias to make the resize
  // to be used in setToMINVAL_comm
  MIN_VAL_VEC_COMM = new double[MAX_BUFFER_SIZE];
  for (int i=0; i< MAX_BUFFER_SIZE; i++)
    {MIN_VAL_VEC_COMM[i]=MIN_VAL;}
  
  //Debug1-2
  if (PrintSize)
    {
      MPI_Barrier(vct->getCART_COMM());
      if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN()))   )
        cout <<"Each core of level " << grid->getLevel() <<" allocated 4 vectors of size buffer_size= " <<buffer_size << " and 1 vector of size MAX_BUFFER_SIZE= " << MAX_BUFFER_SIZE << " for communication within the grid; alias ptr for resizing"<<endl;
    }

  // AMR, ME
  // allocate buffers for the AMR repopulation of particles
  // with size MAX_REPOP_SIZE

  //MAX_NP_REPOP_SIZE= (int) (0.1*nop);  // attempt9, working
  //MAX_NP_REPOP_SIZE= (int) (nop); 
  MAX_NP_REPOP_SIZE= (int) (5*nop);  
  //cout << "MAX_NP_REPOP_SIZE: " << MAX_NP_REPOP_SIZE <<endl;
  //cout << "npmax " <<npmax <<endl;

  // this for all levels
  MIN_VAL_VEC= new double[MAX_NP_REPOP_SIZE* nVar];
  for (int i=0; i<MAX_NP_REPOP_SIZE* nVar; i++ )
    {MIN_VAL_VEC[i]= MIN_VAL;}
	
  REPOP_b_BOTTOM = new double[MAX_NP_REPOP_SIZE* nVar];
  REPOP_b_BOTTOM_ptr = REPOP_b_BOTTOM; // alias for resize
  np_REPOP_b_BOTTOM = 0; // to do also at the beginning of each time step
	
  REPOP_b_TOP = new double[MAX_NP_REPOP_SIZE* nVar];
  REPOP_b_TOP_ptr = REPOP_b_TOP; // alias for resize                     
  np_REPOP_b_TOP = 0; // to do also at the beginning of each time step
	
  REPOP_b_LEFT = new double[MAX_NP_REPOP_SIZE* nVar];
  REPOP_b_LEFT_ptr = REPOP_b_LEFT; // alias for resize                         
  np_REPOP_b_LEFT = 0; // to do also at the beginning of each time step
	
  REPOP_b_RIGHT = new double[MAX_NP_REPOP_SIZE* nVar];
  REPOP_b_RIGHT_ptr = REPOP_b_RIGHT; // alias for resize                         
  np_REPOP_b_RIGHT = 0; // to do also at the beginning of each time step

  //Debug1-3
  if (PrintSize)
    {
      MPI_Barrier(vct->getCART_COMM());
      if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN()))   )
        cout <<"Each core of level " << grid->getLevel() <<" allocated 5 vectors of size MAX_NP_REPOP_SIZE* nVar= " <<MAX_NP_REPOP_SIZE* nVar << " for particle repopulation; alias for resizing"<<endl;
    }


  // if FinerLevel_PRAOps==1,  PRA particles are stored to be passed to finer grids
  if (grid->getLevel()< col->getNgrids()-1)
    {FinerLevels_PRAOps=1;}
  else
    {FinerLevels_PRAOps=0;}

  // to receive and split PRA particles, if level >0
  //max_np_SplitPartComm= (int) (nop*0.5); //attempt9, working
  max_np_SplitPartComm= (int) (nop);
  //cout << "max_np_SplitPartComm " << max_np_SplitPartComm <<", nop " << nop <<endl;
  MAX_NP_SPLIPARTCOMM =  max_np_SplitPartComm* 10; // *10; 
                                                     
  if (grid->getLevel()>0)
    {
      REPOP_receive_b = new double[MAX_NP_REPOP_SIZE* nVar];
      
      SplittedParticles_Comm_BOTTOM = new double [max_np_SplitPartComm* nVar];
      SplittedParticles_Comm_BOTTOM_ptr= SplittedParticles_Comm_BOTTOM; // alias for resize
      SplittedParticles_Comm_TOP = new double [max_np_SplitPartComm* nVar];
      SplittedParticles_Comm_TOP_ptr= SplittedParticles_Comm_TOP; // alias for resize
      SplittedParticles_Comm_LEFT = new double [max_np_SplitPartComm* nVar];
      SplittedParticles_Comm_LEFT_ptr= SplittedParticles_Comm_LEFT; // alias for resize 
      SplittedParticles_Comm_RIGHT = new double [max_np_SplitPartComm* nVar];
      SplittedParticles_Comm_RIGHT_ptr= SplittedParticles_Comm_RIGHT; // alias for resize 

      MIN_VAL_VEC_SP = new double [MAX_NP_SPLIPARTCOMM* nVar];
      for (int i=0; i<MAX_NP_SPLIPARTCOMM* nVar; i++ )
	{MIN_VAL_VEC_SP[i]= MIN_VAL;}

    }
  
  if (PrintSize)
    {
      MPI_Barrier(vct->getCART_COMM());
      if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN()))   )
	{
	  cout <<"Each core of the refined grid allocated 4 vectors of size max_np_SplitPartComm* nVar= " <<max_np_SplitPartComm* nVar <<" and 1 vector of size MAX_NP_SPLIPARTCOMM* nVar= " << MAX_NP_SPLIPARTCOMM* nVar << " for split particles, alias for resize "<<endl;
	  cout << "I am level " << grid->getLevel()<<endl;
	}
    }

  // to save particle with contribute to the first/last ghost node, if level>0 && RefLevelAdj ==1
  if  (grid->getLevel()>0 && vct->getRefLevelAdj()==1)
    {
      // buffers for hosting OS particles, dim is just an educated guess
      npmax_OS = (int) (2*(col->getNxc()+col->getNyc())*npcelx * npcely); 
      cout << "max OS particles: " << npmax_OS << endl;
      
      OS_x= new double [npmax_OS];
      OS_y= new double [npmax_OS];
      OS_u= new double [npmax_OS];
      OS_v= new double [npmax_OS];
      OS_w= new double [npmax_OS];
      OS_q= new double [npmax_OS];
      if (TrackParticleID){
	OS_ParticleID= new unsigned long[npmax_OS];
      }
            
      // buffers for communicating OS particle
      if (TrackParticleID)
	nVarOS=7;
      else
	nVarOS=6;

      max_np_OsPartComm= (int) (npmax_OS*0.01);
      MAX_NP_OSPARTCOMM= npmax_OS*  col->getRatio();   //probably overdim; do some tests
      OSParticles_Comm_BOTTOM = new double [max_np_OsPartComm* nVarOS];
      OSParticles_Comm_BOTTOM_ptr= OSParticles_Comm_BOTTOM; //alias for the resize
      OSParticles_Comm_TOP = new double [max_np_OsPartComm* nVarOS];
      OSParticles_Comm_TOP_ptr= OSParticles_Comm_TOP; //alias for the resize  
      OSParticles_Comm_LEFT = new double [max_np_OsPartComm* nVarOS];
      OSParticles_Comm_LEFT_ptr= OSParticles_Comm_LEFT; //alias for the resize 
      OSParticles_Comm_RIGHT = new double [max_np_OsPartComm* nVarOS];
      OSParticles_Comm_RIGHT_ptr= OSParticles_Comm_RIGHT; //alias for the resize   
      
      MIN_VAL_VEC_OS = new double [MAX_NP_OSPARTCOMM* nVarOS];
      for (int i=0; i<MAX_NP_OSPARTCOMM* nVarOS; i++ )
        {MIN_VAL_VEC_OS[i]= MIN_VAL;}
    }

  // to host Repopulated Particles when subcycling
  SizeRP_Sub= (int) (nop);
  if (grid->getLevel() and SubCycling)
    {
      RP_x= new double[SizeRP_Sub];
      RP_y= new double[SizeRP_Sub];
      RP_u= new double[SizeRP_Sub];
      RP_v= new double[SizeRP_Sub];
      RP_w= new double[SizeRP_Sub];
      RP_q= new double[SizeRP_Sub];
      if (TrackParticleID)
	RP_ParticleID= new unsigned long[SizeRP_Sub];
    }

  // if RESTART is true initialize the particle in allocate method
  restart = col->getRestart_status();
  if (restart!=0){
    if (vct->getCartesian_rank()==0 && ns==0)
      cout << "LOADING PARTICLES FROM RESTART FILE in " + col->getRestartDirName() + "/restart.hdf" << endl;
    stringstream ss;
    //ss << vct->getCartesian_rank(); // this before MLMD
    ss << vct->getCartesian_rank_COMMTOTAL(); 
    string name_file = col->getRestartDirName() + "/restart" + ss.str() + ".hdf";
    //cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "is opening " << name_file << endl;
    // hdf stuff 
    hid_t    file_id, dataspace;
    hid_t    datatype, dataset_id;
    herr_t   status;
    size_t   size;
    hsize_t     dims_out[1];           /* dataset dimensions */
    int status_n;
		
    // open the hdf file
    file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 0){
      cout << "particelle" << endl;
      //cout << "couldn't open file: " << name_file << endl;
      //cout << "RESTART NOT POSSIBLE" << endl;
    }
		
    stringstream species_name;
    species_name << ns;
    // the cycle of the last restart is set to 0
    string name_dataset = "/particles/species_" + species_name.str() + "/x/cycle_0";
    dataset_id = H5Dopen1(file_id,name_dataset.c_str());
    datatype  = H5Dget_type(dataset_id);  
    size  = H5Tget_size(datatype);
    dataspace = H5Dget_space(dataset_id);    /* dataspace handle */
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
		
    // get how many particles there are on this processor for this species
    status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    nop = dims_out[0]; // this the number of particles on the processor!
    // get x
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,x);
    // close the data set
    status = H5Dclose(dataset_id);
		
    // get y
    name_dataset = "/particles/species_" + species_name.str() + "/y/cycle_0";
    dataset_id = H5Dopen1(file_id, name_dataset.c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,y);
    status = H5Dclose(dataset_id);
		
    // get u
    name_dataset = "/particles/species_" + species_name.str() + "/u/cycle_0";
    dataset_id = H5Dopen1(file_id, name_dataset.c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,u);
    status = H5Dclose(dataset_id);
    // get v
    name_dataset = "/particles/species_" + species_name.str() + "/v/cycle_0";
    dataset_id = H5Dopen1(file_id, name_dataset.c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,v);
    status = H5Dclose(dataset_id);
    // get w
    name_dataset = "/particles/species_" + species_name.str() + "/w/cycle_0";
    dataset_id = H5Dopen1(file_id, name_dataset.c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,w);
    status = H5Dclose(dataset_id);
    // get q
    name_dataset = "/particles/species_" + species_name.str() + "/q/cycle_0";
    dataset_id = H5Dopen1(file_id, name_dataset.c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,q);
    status = H5Dclose(dataset_id);
    // ID	
    if (TrackParticleID){
      herr_t (*old_func)(void*);
      void *old_client_data;
      H5Eget_auto1(&old_func, &old_client_data);
      /* Turn off error handling */
      H5Eset_auto1(NULL, NULL);
      name_dataset = "/particles/species_" + species_name.str() + "/ID/cycle_0";
      dataset_id = H5Dopen1(file_id, name_dataset.c_str());
			
      H5Eset_auto1(old_func, old_client_data);
      if (dataset_id>0)
	status = H5Dread(dataset_id, H5T_NATIVE_ULONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,ParticleID);
      else{ 
	for (register int counter=0; counter<nop; counter++)
	  ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];}
    }
    // close the hdf file
    status = H5Fclose(file_id);
		
  }
    
}


/** Interpolation Particle --> Grid */
void Particles2Dcomm::interpP2G(Field* EMf, Grid *grid, VirtualTopology* vct){
  double*** weight;// = newArr3(double,2,2,1);
  allocArr3(&weight, 2, 2, 1);
  double*** temp;// = newArr3(double,2,2,1);
  allocArr3(&temp, 2, 2, 1);

  int ix,iy, temp2,temp1;
  double inv_dx, inv_dy;
  inv_dx = 1.0/dx;
  inv_dy = 1.0/dy;

  //cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << "nop before interpP2G: " << nop<<endl;

  for (register int i=0; i < nop; i++){

    if (i> npmax)
      {
	cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << " in interpP2G, i " << i << " npmax " <<npmax <<endl;
      }


    ix = 2 +  int(floor((x[i]-xstart)*inv_dx));
    iy = 2 +  int(floor((y[i]-ystart)*inv_dy));
 
    //cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << "in interpP2G x[i] " << x[i] << ", y[i] " << y[i] <<endl;

    if (ix-1<0 || iy-1<0 || ix> grid->getNXN()-1 || iy> grid->getNYN()-1)
      {
	cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()  << "Particle mess in P2G, exiting " << endl;
	cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << " ix " << ix << " iy " << iy  << " x: " << x[i] << ", y: " << y[i] << ", qom: " << qom << ", ID: " << ParticleID[i]<< ", PRA_oxStartLeft: " << PRA_oxStartLeft <<", PRA_oyStartLeft: " <<  PRA_oyStartLeft<< ", xstart: " << xstart << ", xend: " << xend << ", ystart: " << ystart << ", yend: " << yend<<", xstart-dx: " << xstart-dx <<",ystart-dy: " << ystart-dy<<", Lx: " <<Lx <<", Ly: " <<Ly  << ", i: " <<i << ", nop: "<<nop  <<  endl;
	fflush(stdout);
	//exit(-1);
	continue;
      }
        
    weight[1][1][0] = ((x[i] - grid->getXN(ix-1,iy-1,0))*inv_dx)*((y[i] - grid->getYN(ix-1,iy-1,0))*inv_dy);
    weight[1][0][0] = ((x[i] - grid->getXN(ix-1,iy,0))*inv_dx)*((grid->getYN(ix-1,iy,0) - y[i])*inv_dy);
    weight[0][1][0] = ((grid->getXN(ix,iy-1,0) - x[i])*inv_dx)*((y[i] - grid->getYN(ix,iy-1,0))*inv_dy);
    weight[0][0][0] = ((grid->getXN(ix,iy,0) - x[i])*inv_dx)*((grid->getYN(ix,iy,0) - y[i])*inv_dy);

    if (0 && vct->getCartesian_rank_COMMTOTAL()== 496)
      {
	cout << "x[i] " << x[i] << " y[i] " << y[i] << " ix " << ix << " iy " << iy  << " xstart "<< xstart << " xend " << xend << " grid->getXN(ix,iy,0) " << grid->getXN(ix,iy,0) << " grid->getYN(ix,iy,0) " << grid->getYN(ix,iy,0) << " weight[1][1][0] " << weight[1][1][0] <<" weight[0][0][0] " <<weight[0][0][0] <<" weight[1][0][0] " <<weight[1][0][0] <<" weight[0][1][0] " <<weight[0][1][0] << " dx " << grid->getDX()<< endl;
      }
    scale(weight,q[i],2,2);
    // add charge density
    EMf->addRho(weight,ix,iy,0,ns);
    // add current density - X
    eqValue(0.0,temp,2,2);
    addscale(u[i],temp,weight,2,2);
    EMf->addJx(temp,ix,iy,0,ns);
    // add current density - Y
    eqValue(0.0,temp,2,2);
    addscale(v[i],temp,weight,2,2);
    EMf->addJy(temp,ix,iy,0,ns);
    // add current density - Z
    eqValue(0.0,temp,2,2);
    addscale(w[i],temp,weight,2,2);
    EMf->addJz(temp,ix,iy,0,ns);
    //Pxx - add pressure tensor
    eqValue(0.0,temp,2,2);
    addscale(u[i]*u[i],temp,weight,2,2);
    EMf->addPxx(temp,ix,iy,0,ns);
    // Pxy - add pressure tensor
    eqValue(0.0,temp,2,2);
    addscale(u[i]*v[i],temp,weight,2,2);
    EMf->addPxy(temp,ix,iy,0,ns);
    // Pxz - add pressure tensor
    eqValue(0.0,temp,2,2);
    addscale(u[i]*w[i],temp,weight,2,2);
    EMf->addPxz(temp,ix,iy,0,ns);
    // Pyy - add pressure tensor
    eqValue(0.0,temp,2,2);
    addscale(v[i]*v[i],temp,weight,2,2);
    EMf->addPyy(temp,ix,iy,0,ns);
    // Pyz - add pressure tensor
    eqValue(0.0,temp,2,2);
    addscale(v[i]*w[i],temp,weight,2,2);
    EMf->addPyz(temp,ix,iy,0,ns);
    // Pzz - add pressure tensor
    eqValue(0.0,temp,2,2);
    addscale(w[i]*w[i],temp,weight,2,2);
    EMf->addPzz(temp,ix,iy,0,ns);
  }
  
  //storing the corner native values if OS operations are to be undertaken
  if (vct->getRefLevelAdj()==1 && grid->getLevel()>0)
    {
      storeCornerOsValues(EMf, vct);
    }


  //cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << "ns " << ns <<"norm2 bef commG " << sum<<endl;
  // communicate contribution from ghost cells     
  EMf->communicateGhostP2G(ns,0,0,0,0,vct);
  //delArr3(weight,2,2);
  freeArr3(&weight);
  //delArr3(temp,2,2);
  freeArr3(&temp);

  return;
}

/* AMR, ME: this communicate, without the Grid object, is kept for compiling issues
   but is NOT to be called by MLMD functions*/
//int Particles2Dcomm::communicate(VirtualTopology* ptVCT){
int Particles2Dcomm::communicate(VirtualTopology* ptVCT){
  cout << "\n\nCOMMUNICATE, NOT TO CALL BY MLMD FUNCTIONS...\nEXITING...\n\n";
  return -1;
}

/** communicate buffers */
// AMR, ME: this one, with the Grid object, is the one to be called be the MLMD functions
int Particles2Dcomm::communicate(VirtualTopology* ptVCT, Grid* grid, int BC_partCommunicate){

// particles also in the ghost cells also for the coarse grids
  // allocate buffers
  MPI_Status status;
  int new_buffer_size;
  int npExitingMax;
  // variable for memory availability of space for new particles
  int avail, availALL, avail1, avail2, avail3, avail4;
  // Feb4, substituted at the end by just putting to MIN_VAL the difference between the particles to really send and the total # of particles
  /*setToMINVAL_comm(b_XDX);
  setToMINVAL_comm(b_XSN);
  setToMINVAL_comm(b_YDX);
  setToMINVAL_comm(b_YSN);*/

  npExitXright =0, npExitXleft =0, npExitYright =0, npExitYleft =0, npExit=0, rightDomain = 0, rightDomainX=0, rightDomainY=0;
  npDeletedBoundary = 0;
  //int np_current = 0;  
 
  //cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << " start communicate nop " << nop<<endl;

  int np_current= 0;
  int nplast = nop-1;
  while (np_current < nplast+1){
    // BC on particles
    // values of BC_partCommunicate
    // 0: coarse grid       May4: no 0 option anymore                                                                                                                       
    // 1: refined grid, inner mover                                                                                                                
    // 2: refined grid, final position

    //cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "Bef ApplyParticleBC, ID " << ParticleID[np_current] << " x[np_current] " << x[np_current] << " y[np_current] " <<y[np_current] <<endl;
    bool Skip=applyParticleBC(BC_partCommunicate, np_current, &nplast, &npDeletedBoundary, grid, ptVCT);
    //cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "After ApplyParticleBC, ID " << ParticleID[np_current]<< " x[np_current] " <<x[np_current] << " y[np_current] " <<y[np_current]<<endl;
    if (!Skip)
      {
	double dx = grid->getDX();
	double dy = grid->getDY();
	// if the particle exits, apply the boundary conditions add the particle to communication buffer
	// enter here if you need to be communicated
	// the first two conditions legitimate, the others to catch errors and eventually to eliminate
	if ( (x[np_current] < xstart && ptVCT->getCoordinates(0) != 0 ) || (x[np_current] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1)) || (x[np_current] < Modified_xstart && ptVCT->getCoordinates(0) == 0 )|| (x[np_current] > Modified_xend && ptVCT->getCoordinates(0) == (ptVCT->getXLEN()-1) )  ){

	// communicate if they don't belong to the domain
	  if (x[np_current] < xstart && ptVCT->getCoordinates(0) != 0){
	    // check if there is enough space in the buffer before putting in the particle
	    if(((npExitXleft+1)*nVar)>=buffer_size){
	      cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicate resizing the sending buffer to " << (int) (buffer_size*2) << " buffer size" << endl;
	      if (!resize_buffers((int) buffer_size*2) )
		{
		  cout<<"communicate: increase MAX_BUFFER_SIZE\nExiting...";
		  return -1;
		}
	    }
	    // put it in the communication buffer
	    bufferXleft(b_XSN,np_current,ptVCT);
	    // delete the particle and pack the particle array, the value of nplast changes
	    del_pack(np_current,&nplast);
	    npExitXleft++;
	  } else if (x[np_current] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1)){
	    // check if there is enough space in the buffer before putting in the particle
	    if(((npExitXright+1)*nVar)>=buffer_size){
	      cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicate resizing the sending buffer " << (int) (buffer_size*2) << endl; 
	      if (!resize_buffers((int) (buffer_size*2)))
                {
                  cout<<"communicate: increase MAX_BUFFER_SIZE\nExiting...";
                  return -1;
                }
	    }
	    // put it in the communication buffer
	    bufferXright(b_XDX,np_current,ptVCT);
	    // delete the particle and pack the particle array, the value of nplast changes
	    del_pack(np_current,&nplast);
	    npExitXright++;
	  } else if (x[np_current] < Modified_xstart && ptVCT->getCoordinates(0) == 0 ){
	    cout <<"Communicate: problem with particle position, exiting";
	    return -1;
	  } else if (x[np_current] > Modified_xend && ptVCT->getCoordinates(0) == (ptVCT->getXLEN()-1)){
	    cout <<"Communicate: problem with particle position, exiting";
	    return -1;
	  }
	} //end condition on x
    else if ( (y[np_current] < ystart && ptVCT->getCoordinates(1) != 0 ) || (y[np_current] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1)) || (y[np_current] < Modified_ystart && ptVCT->getCoordinates(1) == 0 )|| (y[np_current] > Modified_yend && ptVCT->getCoordinates(1) == (ptVCT->getYLEN()-1))     ){

	  // communicate if they don't belong to the domain
	 if (y[np_current] < ystart && ptVCT->getCoordinates(1) != 0){
	 // check if there is enough space in the buffer before putting in the particle
	    if(((npExitYleft+1)*nVar)>=buffer_size){
	      cout << "R"<< ptVCT->getCartesian_rank_COMMTOTAL() << "communicate resizing the sending buffer " << (int) (buffer_size*2) << endl;
	      if (!resize_buffers((int) (buffer_size*2)))
                {
                  cout<<"communicate: increase MAX_BUFFER_SIZE\nExiting...";
                  return -1;
                }
	    }
	    // put it in the communication buffer
	    bufferYleft(b_YSN,np_current,ptVCT);  
	    // delete the particle and pack the particle array, the value of nplast changes
	    del_pack(np_current,&nplast);
	    npExitYleft++;
	 } else if (y[np_current] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1)){
	    // check if there is enough space in the buffer before putting in the particle
	    if(((npExitYright+1)*nVar)>=buffer_size){
	      cout << "R"<< ptVCT->getCartesian_rank_COMMTOTAL() << "communicate resizing the sending buffer " << (int) (buffer_size*2) << endl; 
	      if (!resize_buffers((int) (buffer_size*2)))
                {
                  cout<<"communicate: increase MAX_BUFFER_SIZE\nExiting...";
                  return -1;
                }
	    }
	    // put it in the communication buffer
	    bufferYright(b_YDX,np_current,ptVCT);
	    // delete the particle and pack the particle array, the value of nplast changes
	    del_pack(np_current,&nplast);
	    npExitYright++;
	 } else if (y[np_current] < Modified_ystart && ptVCT->getCoordinates(1) == 0 ){
	   cout <<"Communicate: problem with particle position, exiting";
	    return -1;
	 } else if (y[np_current] > Modified_yend && ptVCT->getCoordinates(1) == (ptVCT->getYLEN()-1)){
	    cout <<"Communicate: problem with particle position, exiting";
	    return -1;
	  }
	}// end condition on y
	else {
	  // AMR, ME
	  // particle is still in the domain, procede with the next particle
	  // PRA ops: check if this particle has to be communicated to finer grid
	  // this check does NOT have to be done during the inner iterations cycles,
	  // so check on LastCommunicate
	  
	  //BOHif (! ( grid->getLevel()==0 || (grid->getLevel()>0 && FinerLevels_PRAOps==1 ))  || PRAIntersection== false )

	  //if (grid->getLevel()==0 && LastCommunicate==1 && FinerLevels_PRAOps==1 && ( PRACollectionMethod== 0 && AlreadyAccumulated[np_current]== false))
	  if (LastCommunicate==1 && FinerLevels_PRAOps==1 && ( PRACollectionMethod== 0 && AlreadyAccumulated[np_current]== false)  && PRAIntersection)
	    {int res;

	      res=PRARepopulationAdd(np_current);
	      AlreadyAccumulated[np_current]=true;
	      if (res<0)
		{
		  cout << "Insufficient repopulation buffer size, exiting" << endl;
		  return -1;
		}
	    }
	  // end AMR, ME
	  np_current++;
	}// end particle in thr right proc
      }// end skip
  }// end while
    
  nop = nplast+1;
  npExitingMax = 0;
  // calculate the maximum number of particles exiting from this domain
  // use this value to check if communication is needed
  // and to  resize the buffer
  npExitingMax = maxNpExiting();
  // broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
  npExitingMax = reduceMaxNpExiting(ptVCT->getCART_COMM(),npExitingMax);
  /*****************************************************/
  /*           SEND AND RECEIVE MESSAGES               */
  /*****************************************************/
	
  
  new_buffer_size = buffer_size;
  while (npExitingMax*nVar + 1> new_buffer_size  )
    {
      new_buffer_size= (int)(new_buffer_size*2);
    }
  
  if (new_buffer_size!= buffer_size)
    {
      cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "ns " << ns << " communicate resizing the receiving buffer" << endl;
      cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "ns " << ns << " new_buffer_size " << new_buffer_size <<endl;
      if (!resize_buffers((int) new_buffer_size) )         
	{                   
	  cout<<"communicate: increase MAX_BUFFER_SIZE\nExiting...";   
	  return -1;                                                                      
	} 
    }

  // Feb 4: the stopping condition for the unbuffer
  b_XSN[nVar*npExitXleft]=MIN_VAL; 
  b_XDX[nVar*npExitXright]=MIN_VAL;
  b_YSN[nVar*npExitYleft]=MIN_VAL;
  b_YDX[nVar*npExitYright]=MIN_VAL;
  // end Feb 4: the stopping condition for the unbuffer  

  if (npExitingMax > 0){
    communicateParticles(npExitingMax*nVar + 1,b_XSN,b_XDX,b_YSN,b_YDX,ptVCT);
    
    // UNBUFFERING
    // message from XLEFT
    //cout << "R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "unbuffer in communicate" <<endl;
    avail1 = unbuffer(b_XDX, ptVCT);
    // message from XRIGHT
    avail2 = unbuffer(b_XSN, ptVCT);
    // message from XLEFT
    avail3 = unbuffer(b_YDX, ptVCT);
    // message from XRIGHT
    avail4 = unbuffer(b_YSN, ptVCT);
    // if one of these numbers is negative than there is not enough space for particles
    avail = avail1 + avail2 + avail3 + avail4;
    availALL = reduceNumberParticles(ptVCT->getCART_COMM(),avail);
    if (availALL < 0)
      return(-1);  // too many particles coming, save data nad stop simulation
    
  }

  //cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "end communicate nop " << nop<<endl;

  
  return(0); // everything was fine
	
}
/** resize the buffers */
bool Particles2Dcomm::resize_buffers(int new_buffer_size){
  cout << "RESIZING FROM " <<  buffer_size << " TO " << new_buffer_size << endl;

  if (new_buffer_size >MAX_BUFFER_SIZE)
    {
      cout << "new_buffer_size (nVar included): " << new_buffer_size <<endl;
      cout << "MAX_BUFFER_SIZE (nVar included): " << MAX_BUFFER_SIZE <<endl;
      cout << "Communication buffers cannot be resized to this dimension... Exiting...\n";
      return false;
    }
  // resize b_XSN
	double *temp = new double[buffer_size];
	
	for(int i=0; i < buffer_size; i++)
		temp[i] = b_XSN_ptr[i];
	//delete[] b_XSN_ptr;
	delete[] b_XSN;
	b_XSN = new double[new_buffer_size];
	for(int i=0; i < buffer_size; i++)
		b_XSN[i] = temp[i];
	for(int i= buffer_size; i < new_buffer_size; i++)
		b_XSN[i] = MIN_VAL;
		
	// resize b_XDX  
	for(int i=0; i < buffer_size; i++)
		temp[i] = b_XDX_ptr[i];
	//delete[] b_XDX_ptr;
	delete[] b_XDX;
	b_XDX = new double[new_buffer_size];
	for(int i=0; i < buffer_size; i++)
		b_XDX[i] = temp[i];
	for(int i=buffer_size; i < new_buffer_size; i++)
		b_XDX[i] = MIN_VAL;
	
	
	
	// resize b_YDX
	for(int i=0; i < buffer_size; i++)
		temp[i] = b_YDX_ptr[i];
	//delete[] b_YDX_ptr;
	delete[] b_YDX;
	b_YDX = new double[new_buffer_size];
	for(int i=0; i < buffer_size; i++)
		b_YDX[i] = temp[i];
	for(int i=buffer_size; i < new_buffer_size; i++)
		b_YDX[i] = MIN_VAL;
	
	// resize b_YSN
	for(int i=0; i < buffer_size; i++)
		temp[i] = b_YSN_ptr[i];
	//delete[] b_YSN_ptr;
	delete[] b_YSN;
	b_YSN = new double[new_buffer_size];
	for(int i=0; i < buffer_size; i++)
		b_YSN[i] = temp[i];
	for(int i=buffer_size; i < new_buffer_size; i++)
		b_YSN[i] = MIN_VAL;
	
	delete[] temp;
	
	b_XDX_ptr = b_XDX; 
	b_YDX_ptr = b_YDX;  
	b_YSN_ptr = b_YSN;
	b_XSN_ptr = b_XSN;
    
	
	buffer_size = new_buffer_size;
	
	return true;
}
/** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
void Particles2Dcomm::bufferXleft(double *b_, int np_current, VirtualTopology* vct){

  b_[npExitXleft*nVar]    = x[np_current];
  b_[npExitXleft*nVar +1] = y[np_current];
  b_[npExitXleft*nVar +2] = u[np_current];
  b_[npExitXleft*nVar +3] = v[np_current];
  b_[npExitXleft*nVar +4] = w[np_current];
  b_[npExitXleft*nVar +5] = q[np_current];
  b_[npExitXleft*nVar +6] = xptilde[np_current];
  b_[npExitXleft*nVar +7] = yptilde[np_current];
  b_[npExitXleft*nVar +8] = uptilde[np_current];
  b_[npExitXleft*nVar +9] = vptilde[np_current];
  b_[npExitXleft*nVar +10] = wptilde[np_current];
  if (TrackParticleID)
    b_[npExitXleft*nVar +11] = ParticleID[np_current];
  if (cVERBOSE)
    cout << "Particle exiting to Xleft: X=" << x[np_current] << " ("<< xstart<<"," << xend << ")"<< endl;
	
}

/** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
void Particles2Dcomm::bufferXright(double *b_, int np_current, VirtualTopology* vct){

  b_[npExitXright*nVar]    = x[np_current];
  b_[npExitXright*nVar +1] = y[np_current];
  b_[npExitXright*nVar +2] = u[np_current];
  b_[npExitXright*nVar +3] = v[np_current];
  b_[npExitXright*nVar +4] = w[np_current];
  b_[npExitXright*nVar +5] = q[np_current];
  b_[npExitXright*nVar +6] = xptilde[np_current];
  b_[npExitXright*nVar +7] = yptilde[np_current];
  b_[npExitXright*nVar +8] = uptilde[np_current];
  b_[npExitXright*nVar +9] = vptilde[np_current];
  b_[npExitXright*nVar +10] = wptilde[np_current];
  if (TrackParticleID)
    b_[npExitXright*nVar +11] = ParticleID[np_current];
  if(cVERBOSE)
    cout << "Particle exiting to Xright: X=" << x[np_current] << " ("<< xstart<< "," << xend << ")" << endl;
}

/** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles2Dcomm::bufferYleft(double *b_, int np, VirtualTopology* vct){
  b_[npExitYleft*nVar]    = x[np];
  b_[npExitYleft*nVar +1] = y[np]; 
  b_[npExitYleft*nVar +2] = u[np];
  b_[npExitYleft*nVar +3] = v[np];
  b_[npExitYleft*nVar +4] = w[np];
  b_[npExitYleft*nVar +5] = q[np];
  b_[npExitYleft*nVar +6] = xptilde[np];
  b_[npExitYleft*nVar +7] = yptilde[np];
  b_[npExitYleft*nVar +8] = uptilde[np];
  b_[npExitYleft*nVar +9] = vptilde[np];
  b_[npExitYleft*nVar +10] = wptilde[np];
  if (TrackParticleID)
    b_[npExitYleft*nVar +11] = ParticleID[np];
  
  if (cVERBOSE)
    cout << "Particle exiting to Yleft: Y=" << y[np] << " ("<< ystart<< "," << yend << ")" << endl;
}

/** put a particle exiting to Y-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles2Dcomm::bufferYright(double *b_, int np, VirtualTopology* vct){
  b_[npExitYright*nVar]    = x[np];
  b_[npExitYright*nVar +1] = y[np];
  b_[npExitYright*nVar +2] = u[np];
  b_[npExitYright*nVar +3] = v[np];
  b_[npExitYright*nVar +4] = w[np];
  b_[npExitYright*nVar +5] = q[np];
  b_[npExitYright*nVar +6] = xptilde[np];
  b_[npExitYright*nVar +7] = yptilde[np];
  b_[npExitYright*nVar +8] = uptilde[np];
  b_[npExitYright*nVar +9] = vptilde[np];
  b_[npExitYright*nVar +10] = wptilde[np];
  if (TrackParticleID)
    b_[npExitYright*nVar +11] = ParticleID[np];
  if (cVERBOSE)
    cout << "Particle exiting to Yright: Y=" << y[np] << " ("<< ystart<< "," << yend << ")" << endl;
}

/** Unpack the received buffer: 
* take the data from the buffer and add particles to the domain 
* check if it is the right domain:
* with implicit scheme particles can transverse more than one domain*/
int Particles2Dcomm::unbuffer(double *b_){
  cout << "USE UNBUFFER(buffer, vct), EXITING"<< endl;
  return -1;
}
int Particles2Dcomm::unbuffer(double *b_, VirtualTopology *ptVCT){
  int np_current =0;
  // put the new particles at the end of the array, and update the number of particles
  while(b_[np_current*nVar] != MIN_VAL){
    x[nop] = b_[nVar*np_current];
    y[nop] = b_[nVar*np_current+1];
    u[nop] = b_[nVar*np_current+2];
    v[nop] = b_[nVar*np_current+3];
    w[nop] = b_[nVar*np_current+4];
    q[nop] = b_[nVar*np_current+5];
    xptilde[nop] = b_[nVar*np_current+6];
    yptilde[nop] = b_[nVar*np_current+7];
    uptilde[nop] = b_[nVar*np_current+8];
    vptilde[nop] = b_[nVar*np_current+9];
    wptilde[nop] = b_[nVar*np_current+10];
    if (TrackParticleID)
      ParticleID[nop]=(unsigned long) b_[nVar*np_current+11];

    //Feb 3
    if(FinerLevels_PRAOps==1 && PRACollectionMethod ==0)
      AlreadyAccumulated[nop]= false;
    // end Feb 3
    np_current++;
    
    if (cVERBOSE)
      cout << "Receiving Particle: X=" << x[nop] << ",Y=" << y[nop] << " ("<< xstart<<"," << xend << ")"<< " x ("<< ystart<<"," << yend << ")"<<endl;

    bool XnotRightDom= (x[nop] < xstart && ptVCT->getCoordinates(0) != 0 ) || (x[nop] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1) );
    bool YnotRightDom=(y[nop] < ystart && ptVCT->getCoordinates(1) != 0 ) || (y[nop] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1));


    if(XnotRightDom || YnotRightDom)
      {
	rightDomain++; // the particle is not in the domain
	if(XnotRightDom) rightDomainX++;
	if(YnotRightDom) rightDomainY++;
      }
    else // the particle IS in the right domain, check for PRA ops
      {
	if (LastCommunicate==1 && FinerLevels_PRAOps==1 && PRACollectionMethod ==0)
	  // i.e.: this is the unbuffer after the definitive positions, not the one in the inner mover
	  // AND this level is coarser level to somebody else (needs to store particles for repop)
	  // AND I am collecting coarse particle for repop one by one
	  {int res;
	    //cout <<  "R" <<ptVCT->getCartesian_rank_COMMTOTAL() <<" PRARepopulationAdd in unbuffer\n";

	    /*if (ParticleID[nop]== 169310 && qom== -1)
	      cout <<"Particle 169310 -1 at PRARepopulationAdd in unbuffer\n";*/
	    
	    

	    res=PRARepopulationAdd(nop);
	    AlreadyAccumulated[nop]=true;
	    if (res<0)
	      {
		cout << "Insufficient repopulation buffer size, exiting" << endl;
		return -1;
	      }
	  }//end if LastCommunicate
      }// end particle in right domain
    // end AMR, ME
    nop++;
    if (nop > (npmax - (int) (.01*npmax) ) ){
      cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "Unbuffer exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
      return(-1); // end the simulation because you dont have enough space on the array
    }
  }// end while
  return(0); // everything was fine
  
}
/** Delete the a particle from the array and pack the the array, update the number of 
* particles that are exiting
* For deleting the particle from the array take the last particle and put it
* in the position of the particle you want to delete
* @param np = the index of the particle that must be deleted
* @param nplast = the index of the last particle in the array
*/
void Particles2Dcomm::del_pack(int np_current, int *nplast){
	x[np_current] = x[*nplast];
	y[np_current] = y[*nplast];
	u[np_current] = u[*nplast];
	v[np_current] = v[*nplast];
	w[np_current] = w[*nplast];
	q[np_current] = q[*nplast];
	xptilde[np_current] = xptilde[*nplast];
	yptilde[np_current] = yptilde[*nplast];
	uptilde[np_current] = uptilde[*nplast];
	vptilde[np_current] = vptilde[*nplast];
	wptilde[np_current] = wptilde[*nplast];
	if (TrackParticleID)
		ParticleID[np_current]=ParticleID[*nplast];

	// MLMD ops
	if (FinerLevels_PRAOps==1 && PRACollectionMethod== 0)
	  AlreadyAccumulated[np_current]= AlreadyAccumulated[*nplast];
	
	npExit++;
	(*nplast)--;
}
/** method to calculate how many particles are out of right domain */
int Particles2Dcomm::isMessagingDone(VirtualTopology* ptVCT){
	int result = 0;
	result = reduceNumberParticles(ptVCT->getCART_COMM(),rightDomain);
	if (result > 0 && cVERBOSE && ptVCT->getCartesian_rank()==0)
		cout << "Further Comunication: " << result << " particles not in the right domain" << endl;
	return(result);
	
}

/** method to calculate how many boundary particles are out of right domain */
int Particles2Dcomm::isMessagingDoneSP(VirtualTopology* ptVCT){
  /*int resultL = 0, resultR =0, resultB=0, resultT=0;
  if (ptVCT->getCOMM_B_BOTTOM()!= MPI_COMM_NULL)
    resultB = reduceNumberParticles(ptVCT->getCOMM_B_BOTTOM(),rightDomainX);
  if (ptVCT->getCOMM_B_TOP()!= MPI_COMM_NULL)
    resultT = reduceNumberParticles(ptVCT->getCOMM_B_TOP(),rightDomainX);
  if (ptVCT->getCOMM_B_LEFT()!= MPI_COMM_NULL)
    resultL = reduceNumberParticles(ptVCT->getCOMM_B_LEFT(),rightDomainY);
  if (ptVCT->getCOMM_B_RIGHT()!= MPI_COMM_NULL)
    resultR = reduceNumberParticles(ptVCT->getCOMM_B_RIGHT(),rightDomainY);
  if (resultB>0 or resultT>0 or resultL>0 or resultR>0)
    return(1) ; // if return is >0, communication is continued
  
    return (0);*/
  int result = 0;
  if (ptVCT->getCOMM_B_ALL()== MPI_COMM_NULL) return (result);
  
  result = reduceNumberParticles(ptVCT->getCOMM_B_ALL(),rightDomain);
  if (result > 0 && cVERBOSE && ptVCT->getCartesian_rank()==0)
    cout << "Further Comunication: " << result << " particles not in the right domain" << endl;
  return(result);
}

/** calculate the maximum number exiting from this domain */
int Particles2Dcomm::maxNpExiting(){
	int maxNp = 0;
	if (npExitXright > maxNp)
		maxNp = npExitXright;
	if (npExitXleft  > maxNp)
		maxNp = npExitXleft;
	if (npExitYright > maxNp)
		maxNp = npExitYright;
	if (npExitYleft  > maxNp)
		maxNp = npExitYleft;
	return(maxNp);
}
/** return X-coordinate of particle array */
double* Particles2Dcomm::getXall() const{ return(x);}
/** return Y-coordinate  of particle array */
double* Particles2Dcomm::getYall() const{ return(y);}
/** return Z-coordinate  of particle array*/
double* Particles2Dcomm::getZall() const{ 
	cout << "2D Particle in X-Y space. no need for calling Particles2Dcomm::getZall()" << endl;
	return(x);
}
/** get X-velocity of particle with label indexPart */
double* Particles2Dcomm::getUall() const{ return(u);}

/** get Y-velocity of particle with label indexPart */
double* Particles2Dcomm::getVall() const{ return(v);}

/**get Z-velocity of particle with label indexPart */
double* Particles2Dcomm::getWall() const{ return(w);}

/**get ID of particle with label indexPart */
unsigned long* Particles2Dcomm::getParticleIDall() const{return (ParticleID);}

/**get charge of particle with label indexPart */
double* Particles2Dcomm::getQall() const{ return(q);}

/** return X-coordinate of particle with index indexPart */
double Particles2Dcomm::getX(int indexPart) const{	return(x[indexPart]);}

/** return Y-coordinate  of particle with index indexPart */
double Particles2Dcomm::getY(int indexPart) const{ return(y[indexPart]);}

/** return Y-coordinate  of particle with index indexPart */
double Particles2Dcomm::getZ(int indexPart) const{ 
	cout << "2D Particle in X-Y space. no need for calling Particles2DcommXY::getZ(int indexPart) " << endl;
	return(x[0]);
}

/** get u (X-velocity) of particle with label indexPart */
double Particles2Dcomm::getU(int indexPart) const{ return(u[indexPart]);}


/** get v (Y-velocity) of particle with label indexPart */
double Particles2Dcomm::getV(int indexPart) const{ return(v[indexPart]);}

/**get w (Z-velocity) of particle with label indexPart */
double Particles2Dcomm::getW(int indexPart) const{ return(w[indexPart]);}


/**get ID of particle with label indexPart */
unsigned long Particles2Dcomm::getParticleID(int indexPart) const{ return(ParticleID[indexPart]);}


/**get charge of particle with label indexPart */
double Particles2Dcomm::getQ(int indexPart) const{ return(q[indexPart]);}


/** return the number of particles */
int Particles2Dcomm::getNOP() const{  return(nop);}

/** print particles info */
void Particles2Dcomm::Print(VirtualTopology* ptVCT)const{
    cout << endl;
    cout << "Number of Particles: " << nop << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ","  << ptVCT->getCoordinates(2) << ")"<< endl;
    cout <<  "Xin = " << xstart << "; Xfin = " << xend << endl;
    cout <<  "Yin = " << ystart << "; Yfin = " << yend << endl;
    cout <<  "Zin = " << 0 << "; Zfin = " << 0 << endl;
    cout <<  "Number of species = " << ns << endl;
    for (int i=0; i < nop; i++)
		cout << "Particles #" << i << " x=" << x[i] << " y=" << y[i] << " z=" << 0 << " u=" << u[i] << " v="<< v[i] << " w=" << w[i] << endl;
    cout << endl;
}
/** print just the number of particles */
void Particles2Dcomm::PrintNp(VirtualTopology* ptVCT)const{
    cout << endl;
    cout << "Number of Particles of species "<< ns << ": " << nop << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ","  << ptVCT->getCoordinates(2) << ")"<< endl;
    cout << endl;
}

/**AMR methods, ME*/
/**init operations connected with the repopulation of particles, sometimes fails on ALL_TARGETS== ALL_RECEIVERS*/
/*int Particles2Dcomm::initPRAVariables(int species, CollectiveIO* col,VirtualTopology* vct, Grid* grid, Field* EMf){
  // the buffer for communication of PRA variables are defined in allocate, 
  // together with the other particle communication buffers

  if (vct->getNgrids()<2)
    return 1;

  ratio= col->getRatio();

  //number of cells, ghost cell INCLUDED, for particle repopulation; x left
  PRA_Xleft = col->GetPRA_Xleft(); 
  //number of cells, ghost cell INCLUDED, for particle repopulation; x right
  PRA_Xright = col->GetPRA_Xright();
  //number of cells, ghost cell INCLUDED, for particle repopulation; y left
  PRA_Yleft = col->GetPRA_Yleft();
  //number of cells, ghost cell INCLUDED, for particle repopulation; y right
  PRA_Yright = col->GetPRA_Yright();
  
  if (grid->getLevel() >0){  //PRA area, native particles falling here are deleted and substituted with repopulated particles 
    PRA_oxStartLeft   = 0- grid->getDX();
    PRA_oxEndLeft     = 0+ (PRA_Xleft-1)*grid->getDX();
    PRA_oxStartRight  = Lx- (PRA_Xright-1)*grid->getDX();
    PRA_oxEndRight    = Lx+ grid->getDX();
    
    PRA_oyStartLeft   = 0- grid->getDY();
    PRA_oyEndLeft     = 0+ (PRA_Yleft-1)*grid->getDY();
    PRA_oyStartRight  = Ly- (PRA_Yright-1)*grid->getDY();
    PRA_oyEndRight    = Ly+ grid->getDY();
  }
  else{   //used only for the safety checks, NOT in the code; initialized anyhow
    PRA_oxStartLeft   = 0- grid->getDX();
    PRA_oxEndLeft     = 0;
    PRA_oxStartRight  = Lx;
    PRA_oxEndRight    = Lx+ grid->getDX();

    PRA_oyStartLeft   = 0- grid->getDY();
    PRA_oyEndLeft     = 0;
    PRA_oyStartRight  = Ly;
    PRA_oyEndRight    = Ly+ grid->getDY();
  }
      
  //if Level< Levels-1 (this grid is the coarser grid for some other grid), limits for thr PRA area of the child in local coords
  if ( grid->getLevel() < vct->getNgrids()-1){

    double Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid
    double Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid

    double finedx = grid->getDX()/col->getRatio();
    double finelx = col->getLx()/(double)pow(col->getRatio(),grid->getLevel()+1);
    double finedy = grid->getDY()/col->getRatio();
    double finely = col->getLy()/(double)pow(col->getRatio(),grid->getLevel()+1);
    
    // when parent particles enter this area, they have to be communicated;  
    // the parent dx or dy is already taken into account here

    PRA_CoxStartLeft   = Ox - finedx - grid->getDX() ;
    PRA_CoxEndLeft     = Ox + (PRA_Xleft-1)*finedx + grid->getDX();
    PRA_CoxStartRight  = Ox + finelx - (PRA_Xright-1)*finedx -  grid->getDX();
    PRA_CoxEndRight    = Ox + finelx + finedx + grid->getDX();

    PRA_CoyStartLeft   = Oy - finedy - grid->getDY();
    PRA_CoyEndLeft     = Oy + (PRA_Yleft-1)*finedy + grid->getDY();
    PRA_CoyStartRight  = Oy + finely - (PRA_Yright-1)*finedy - grid->getDY();
    PRA_CoyEndRight    = Oy + finely + finedy + grid->getDY();

    //if (Modified_xstart > PRA_CoxEndLeft or Modified_xend < PRA_CoxStartLeft or Modified_ystart > PRA_CoyEndLeft or Modified_yend < PRA_CoyStartLeft)
    if (Modified_xstart > PRA_CoxEndRight or Modified_xend < PRA_CoxStartLeft or Modified_ystart > PRA_CoyEndRight or Modified_yend < PRA_CoyStartLeft)
      {
	PRAIntersection= false;
      }	
    else
      {
	PRAIntersection= true;
      }

  }
  else{ //actually not used, initialized anyhow
    PRA_CoxStartLeft   = 0- grid->getDX();;
    PRA_CoxEndLeft     = 0- grid->getDX();;
    PRA_CoxStartRight  = Lx+ grid->getDX();
    PRA_CoxEndRight    = Lx+ grid->getDX();

    PRA_CoyStartLeft   = 0- grid->getDY();;
    PRA_CoyEndLeft     = 0- grid->getDY();;
    PRA_CoyStartRight  = Ly+ grid->getDY();
    PRA_CoyEndRight    = Ly+ grid->getDY();


    PRAIntersection= false;
  } 

  // for PRASend, modified from initWeightBC
  int i,j,nproc;
  double finedx, finedy,finelx,finely, xfirst, xlast,xfirstnext, Ox, Oy;
  double coarsedx, coarsedy,coarselx,coarsely;
  double finelxplusfinedx;
  double xshift, yshift;
  int xnnl, xnnu, ynnl, ynnu;
  // remember that there is always only 1 PRA ghost cell, the other in the active domain
  targetBC= new int[col->getXLEN()*col->getYLEN()];
  BCSide  = new int[col->getXLEN()*col->getYLEN()];
  
  targetBOTTOM=0;
  targetTOP=0;
  targetLEFT=0;
  targetRIGHT=0;
  nmessageBC = 0;
  
  //cout << "\nProc " << vct->getCartesian_rank_COMMTOTAL()<< " started building targets out of it\n\n";

  if ( grid->getLevel() < vct->getNgrids()-1) {

    //cout << "\nProc " << vct->getCartesian_rank_COMMTOTAL()<< " started building targets\n\n";

    finedx = grid->getDX()/col->getRatio();
    finelx = col->getLx()/(double)pow(col->getRatio(),grid->getLevel()+1);
    finedy = grid->getDY()/col->getRatio();
    finely = col->getLy()/(double)pow(col->getRatio(),grid->getLevel()+1);
    finelxplusfinedx = col->getLx()/(double)col->getNxc()*((double)col->getNxc()+1)/(double)pow((double)col->getRatio(),(double)grid->getLevel()+1.);
    Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid 
    Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid
    j=0;
    // last 2 conditions from initWeightBC
    if ( grid->getYend()>= PRA_CoyStartLeft && grid->getYstart()<= PRA_CoyEndLeft && grid->getXstart() <= PRA_CoxEndRight && grid->getXend() >= PRA_CoxStartLeft   ) { 
      xfirst = max(PRA_CoxStartLeft - Ox , ceil((grid->getXstart()-Ox)/finedx)*finedx);
      xfirstnext = ceil((grid->getXend()-Ox)/finedx)*finedx;
      xlast  = min(PRA_CoxEndRight-Ox , floor((grid->getXend()-xfirst-Ox)/finedx)*finedx+xfirst);

      //if(fabs(grid->getXend()-xlast-Ox) < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision 
      if(fabs(grid->getXend()-xlast-Ox)*0.99999 < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision; 0.99999 because the obvious choices of ratios often bring to overlapping that this alone does not resolve
	xlast = xlast - finedx;
      }


      xnnl = floor((xlast-xfirst)/finedx+0.5)+1; //floor(x+0.5) used to round to closest integer in case the division is not working properly                                     
      for (i =0;i<xnnl;i++) {
	xshift = xfirst + i * finedx ;
	nproc = col->getYLEN()*floor(xshift/((grid->getNXC()-2.)*finedx))+col->getXLEN()*col->getYLEN()*(grid->getLevel()+1); // rank of the proc on the fine grid receiving this point(in MPI_COMM_WORLD)           
	nproc = max(nproc, col->getXLEN()*col->getYLEN()*(grid->getLevel()+1));
	nproc = min(nproc, col->getYLEN()*((grid->getLevel()+2)*col->getXLEN()-1));
	if (i==0){
	  targetBC[0] = nproc;
	  nmessageBC++;                                                                                                                        
	  BCSide[0]=0;
	  targetBOTTOM++;                                                                                        
	}
	if(nproc != targetBC[j]){
	  j++;
	  nmessageBC++;                                                                                                                        
	  BCSide[j]=0;
	  targetBOTTOM++;                                                                                                                      
	  targetBC[j]=nproc;
	}
      }
    }
    
    if (grid->getYend()>= PRA_CoyStartRight && grid->getYstart()<=PRA_CoyEndRight && grid->getXstart() <= PRA_CoxEndRight && grid->getXend() >= PRA_CoxStartLeft) {
      xfirst = max(PRA_CoxStartLeft- Ox, ceil((grid->getXstart()-Ox)/finedx)*finedx);
      xfirstnext = ceil((grid->getXend()-Ox)/finedx)*finedx;
      xlast  = min(PRA_CoxEndRight - Ox, floor((grid->getXend()-xfirst-Ox)/finedx)*finedx+xfirst);
    
      //if(fabs(grid->getXend()-xlast-Ox) < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision 
      if(fabs(grid->getXend()-xlast-Ox)*0.99999 < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision; 0.99999 because the obvious choices of ratios often bring to overlapping that this alone does not resolve          
	xlast = xlast - finedx;
      }
      xnnu = floor((xlast-xfirst)/finedx+0.5)+1;
      for (i =0;i<xnnu;i++) {
	xshift = xfirst + i * finedx ;
	nproc = col->getYLEN()*(floor(xshift/((grid->getNXC()-2.)*finedx))+1)-1+col->getXLEN()*col->getYLEN()*(grid->getLevel()+1); // rank of the proc on the fine grid receiving this point(in MPI_COMM_WORLD)                  
	nproc = max(nproc, col->getYLEN()*(col->getXLEN()*(grid->getLevel()+1)+1)-1);
	nproc = min(nproc, col->getYLEN()*col->getXLEN()*(grid->getLevel()+2)-1);
	if (i==0){
	  if(nmessageBC>0){
	    j++;
	  }
	  targetBC[j] = nproc;
	  nmessageBC++;
	  BCSide[j]=1;
	  targetTOP++;                                  
	}
	if(nproc != targetBC[j]){
	  j++;
	  nmessageBC++;                                  
	  BCSide[j]=1;
	  targetTOP++;                      
	  targetBC[j]=nproc;
	}
      }
    }
    // left
    if (grid->getYstart() <= PRA_CoyStartRight && grid->getYend() >=  PRA_CoyEndLeft  && grid->getXstart() <= PRA_CoxEndLeft && grid->getXend() >= PRA_CoxStartLeft) {
      xfirst = max(PRA_CoyEndLeft -Oy , ceil((grid->getYstart()-Oy)/finedy)*finedy);
      xfirstnext = ceil((grid->getYend()-Oy)/finedy)*finedy;
      xlast  = min(  PRA_CoyStartRight - Oy, floor((grid->getYend()-xfirst-Oy)/finedy)*finedy+xfirst);
      
      if(fabs(grid->getYend()-xlast-Oy) < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision
      
	xlast = xlast - finedy;
      }
      ynnl = floor((xlast-xfirst)/finedy+0.5)+1;
      for (i =0;i<ynnl;i++) {
	yshift = xfirst + i * finedy ;                              
	nproc =floor(yshift/((grid->getNYC()-2.)*finedy))+col->getXLEN()*col->getYLEN()*(grid->getLevel()+1); // rank of the proc on the fine grid receiving this point(in MPI_\
	COMM_WORLD)                                                                                                                                                                         
	nproc = max(nproc, col->getYLEN()*col->getXLEN()*(grid->getLevel()+1));
	nproc = min(nproc, col->getYLEN()*(1+col->getXLEN()*(grid->getLevel()+1))-1);
	if (i==0){
	  if(nmessageBC>0){
	    j++;
	  }
	  targetBC[j] = nproc;                                                                                                                                            
	  BCSide[j]=2;
	  targetLEFT++;
	  nmessageBC++;
	}
	if(nproc != targetBC[j]){
	  j++;
	  nmessageBC++;
	  targetBC[j]=nproc;                                                                                                                                                        BCSide[j]=2;
	  targetLEFT++;                                                                                                                                            
	}
    } 
  }
  // right
  if (grid->getYstart() <=PRA_CoyStartRight && grid->getYend() >=PRA_CoyEndLeft  && grid->getXstart() <=PRA_CoxEndRight  && grid->getXend() >= PRA_CoxStartRight) {
    xfirst = max( PRA_CoyEndLeft - Oy, ceil((grid->getYstart()-Oy)/finedy)*finedy);
    xfirstnext = ceil((grid->getYend()-Oy)/finedy)*finedy;
    //    xlast  = min( PRA_CoyStartRight , floor((grid->getYend()-xfirst-Oy)/finedy)*finedy+xfirst);
    xlast  = min(  PRA_CoyStartRight - Oy, floor((grid->getYend()-xfirst-Oy)/finedy)*finedy+xfirst);  // fixed
    if(fabs(grid->getYend()-xlast-Oy) < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision 

      xlast = xlast - finedy;
    }
    ynnu = floor((xlast-xfirst)/finedy+0.5)+1;
    for (i =0;i<ynnu;i++) {
      yshift = xfirst + i * finedy ;                       
      nproc =floor(yshift/((grid->getNYC()-2.)*finedy))+col->getYLEN()*(col->getXLEN()*(grid->getLevel()+2)-1); // rank of the proc on the fine grid receiving this point(in MPI_COMM_WORLD)                                                                                                                                                                     
      nproc = max(nproc, col->getYLEN()*(col->getXLEN()*(grid->getLevel()+2)-1));
      nproc = min(nproc, col->getYLEN()*col->getXLEN()*(grid->getLevel()+2)-1); 
      if (i==0){
	if(nmessageBC>0){
	  j++;
	}
	targetBC[j] = nproc;                                                                                                                                                      
	BCSide[j]=3;
	targetRIGHT++;                                                                                                                                     
	nmessageBC++;
      }
      if(nproc != targetBC[j]){
	j++;
	nmessageBC++;
	targetBC[j]=nproc;                                                                                                                                                       
	BCSide[j]=3;
	targetRIGHT++;
      }
    } // end ynnu
  }   // end right

  //cout << "\nProc " << vct->getCartesian_rank_COMMTOTAL() << " finished building targets\n\n";

}  // end check on grid level



// end for PRASend, modified from initWeightBC   
// for PRAReceive, modified from initWeightBC
nmessagerecuBC=0;
// for debugging purposes, not actually used anywhere
int nmessagerecuBCLEFT=0;
int nmessagerecuBCRIGHT=0;
int nmessagerecuBCBOTTOM=0;
int nmessagerecuBCTOP=0;
// end for debugging purposes
if (grid->getLevel() > 0){
  fromBC= new int[col->getXLEN()*col->getYLEN()];  
  BCSidecu= new int[col->getXLEN()*col->getYLEN()];
  nmessagerecuBC=0;
  
  //If this grid is considered as fine by another grid
  coarsedx = grid->getDX()*col->getRatio();
  coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel()-1);
  coarsedy = grid->getDY()*col->getRatio();
  coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel()-1);
  Ox = grid->getOx(grid->getLevel()); //Origin x of the grid
  Oy = grid->getOy(grid->getLevel()); //Origin y of the grid
  j=0;
  if(vct->getCoordinates(1) == 0) {    // BOTTOM
    double xloc, yloc1, yloc2;
    if (vct->getCoordinates(0)== 0)
      {
	xfirst=PRA_oxStartLeft- coarsedx ; // coarsedx from def, to have the areas coincide in coarse and fine grid
      }
    else
      {
	xfirst= grid->getXstart();
      }
    if(vct->getCoordinates(0) == vct->getXLEN()-1) {
      xlast = PRA_oxEndRight + coarsedx;// coarsedx from def, to have the areas coincide in coarse and fine grid
    }else{
      xlast = grid->getXend();
      if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc
	{
	  xlast= xlast - grid->getDX();
	}
    }
    
    yloc1 = Oy + PRA_oyStartLeft- coarsedy;// coarsedy from def, to have the areas coincide in coarse and fine grid
    yloc2 = Oy + PRA_oyEndLeft + coarsedy; // coarsedy from def, to have the areas coincide in coarse and fine grid
    double YL= yloc1;
    int nproc;
    while (YL < yloc2 || fabs(YL - yloc2)< DBL_EPSILON ){
      xnnl = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
      for (i=0; i< xnnl; i++) {
	xloc = max(Ox + xfirst +i*grid->getDX(),0.);// Because when x < 0, it is on the same proc as if x=0  
	nproc =floor(YL/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1);
        bool found=false;
        for (int k=0; k< nmessagerecuBC; k++ )
          {
            if (BCSidecu[k]==0 && fromBC[k]==nproc)
	      { 
		found= true;
		break;
              }
          }
	if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	  fromBC[0] = nproc;
	  BCSidecu[0]=0;
	  nmessagerecuBC++;
	  nmessagerecuBCBOTTOM++;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new bottom from "<<fromBC[nmessagerecuBC-1] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	  }
	
	if(nproc != fromBC[j] && !found){ // first part to avoid repetitions with points at the same yloc, second with previous yloc
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCBOTTOM++;
	  fromBC[j]=nproc;
	  BCSidecu[j]=0;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new bottom from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC << " j " << j<<endl;
	}
      }// end xnnl
      YL += grid->getDY();
    }// end YL
  } // end BOTTOM
  
  if(vct->getCoordinates(1) == vct->getYLEN()-1) { // TOP
    double xloc, yloc1, yloc2;
    if (vct->getCoordinates(0)== 0)
      {
	xfirst=PRA_oxStartLeft - coarsedx; // coarsedx from def, to have the areas coincide in coarse and fine grid
      }
    else
      {
	xfirst= grid->getXstart();
      }
    if(vct->getCoordinates(0) == vct->getXLEN()-1) {
      xlast = PRA_oxEndRight + coarsedx;   // coarsedx from def, to have the areas coincide in coarse and fine grid
    }else{
      xlast = grid->getXend();
      if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc               
	{
	  xlast= xlast - grid->getDX();
	}
    }
    yloc1 = Oy +   PRA_oyStartRight- coarsedy;// coarsedy from def, to have the areas coincide in coarse and fine grid
    yloc2 = Oy +   PRA_oyEndRight +coarsedy; // coarsedy from def, to have the areas coincide in coarse and fine grid
    double YL= yloc1;
    int nproc;
    while (YL < yloc2 || fabs(YL - yloc2)< DBL_EPSILON ){
      xnnu = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
      for (i=0; i< xnnu; i++) {
	xloc = max(Ox + xfirst +i*grid->getDX(),0.);
	nproc =floor(YL/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	bool found=false;
        for (int k=0; k< nmessagerecuBC; k++ )
          {
            if (BCSidecu[k]==1 && fromBC[k]==nproc)
	      { 
		found= true;
		//cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<< " BCSidecu[k] " <<BCSidecu[k] << " fromBC[k] " << fromBC[k] << " nproc " <<nproc <<endl; 
		break;
              }
          }
        if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	  fromBC[0] = nproc;
	  BCSidecu[0]=1;
	  nmessagerecuBC++;
	  nmessagerecuBCTOP++;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new bottom from "<<fromBC[nmessagerecuBC-1] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	}
	// (nproc == fromBC[j] && BCSidecu[j]!= 1) othewise procs are not included if they are the last in a different side
	if((nproc != fromBC[j] || (nproc == fromBC[j] && BCSidecu[j]!= 1)) && !found){
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCTOP++;
	  fromBC[j]=nproc;
	  BCSidecu[j]=1;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new top from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	}
      }// end xnnu
      YL += grid->getDY();
    }// end YL
  } // end TOP
  
  if(vct->getCoordinates(0) == 0) {   // LEFT
    double yloc, xloc1, xloc2;      
    if (vct->getCoordinates(1)== 0)
      {
	xfirst=PRA_oyEndLeft + coarsedy;// this actually spans the y dir // coarsedy from def, to have the areas coincide in coarse and fine grid   
      }
    else
      {
	xfirst= grid->getYstart();
      }
    if(vct->getCoordinates(1) == vct->getYLEN()-1) {
      xlast = PRA_oyStartRight - coarsedy; // coarsedy from def, to have the areas coincide in coarse and fine grid
    }else{
      xlast = grid->getYend();
      if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc              
	{
	  xlast= xlast - grid->getDY();
	}
    }
    xloc1 = Ox + PRA_oxStartLeft -coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid
    xloc2 = Ox + PRA_oxEndLeft+ coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid
    double XL= xloc1;
    int nproc;
    while (XL < xloc2 || fabs(XL - xloc2)< DBL_EPSILON) {
      ynnl = floor((xlast-xfirst)/grid->getDY()+0.5)+1; 
      for (i=0; i< ynnl; i++) {
	yloc = Oy + xfirst +i*grid->getDY();
	nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(XL/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	bool found=false;
        for (int k=0; k< nmessagerecuBC; k++ )
          {
            if (BCSidecu[k]==2 && fromBC[k]==nproc)
	      { 
		found= true;
		break;
              }
          }
	if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	  fromBC[0] = nproc;
	  BCSidecu[0]=2;
	  nmessagerecuBC++;
	  nmessagerecuBCLEFT++;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new left from "<<fromBC[nmessagerecuBC-1] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	}
	if((nproc != fromBC[j] || (nproc == fromBC[j] && BCSidecu[j]!= 2)) && !found){
	  //if(nproc != fromBC[j] && !found){ // first part to avoid repetitions with points at the same yloc, second with previous yloc
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCLEFT++;
	  fromBC[j]=nproc;
	  BCSidecu[j]=2;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new left from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC << " j " << j<<endl;
	}
      } // end ynnl
      XL += grid->getDX();
    }// end XL
  } // end left
  
  if(vct->getCoordinates(0) == vct->getXLEN()-1) {  //RIGHT
    double yloc, xloc1, xloc2;
    if (vct->getCoordinates(1)== 0)
      {
	xfirst=PRA_oyEndLeft + coarsedy;// this actually spans the y dir // coarsedy from def, to have the areas coincide in coarse and fine grid
      }
    else
      {
	xfirst= grid->getYstart();
      }
    if(vct->getCoordinates(1) == vct->getYLEN()-1) {
      xlast = PRA_oyStartRight- coarsedy;// coarsedy from def, to have the areas coincide in coarse and fine grid 
    }else{
      xlast = grid->getYend();
      if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc
	{
	  xlast= xlast - grid->getDY();
	}
    }
    xloc1 = Ox + PRA_oxStartRight - coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid  
    xloc2 = Ox + PRA_oxEndRight + coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid
    double XL= xloc1;
    int nproc;
    ynnu = floor((xlast-xfirst)/grid->getDY()+0.5)+1;
    while (XL < xloc2 || fabs(XL - xloc2)< DBL_EPSILON ){
      for (i=0; i< ynnu; i++) {
	yloc = Oy + xfirst +i*grid->getDY(); 
	nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(XL/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	//cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << " right nproc1 " << nproc1 << " nproc2 "<< nproc2 << " xloc1 " <<xloc1 <<" xloc2 "<<xloc2 <<" yloc "<< yloc <<endl;
        bool found=false;
        for (int k=0; k< nmessagerecuBC; k++ )
          {
            if (BCSidecu[k]==3 && fromBC[k]==nproc)
	      { 
		found= true;
		break;
              }
          }
	if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	  fromBC[0] = nproc;
	  BCSidecu[0]=3;
	  nmessagerecuBC++;
	  nmessagerecuBCRIGHT++;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new right from "<<fromBC[0] << " nmessagerecuBC " << nmessagerecuBC << " nproc " << nproc<<endl;
	}
	if((nproc != fromBC[j] || (nproc == fromBC[j] && BCSidecu[j]!= 3)) && !found){
	  //if(nproc != fromBC[j] && !found){ // first part to avoid repetitions with points at the same yloc, second with previous yloc
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCRIGHT++;
	  fromBC[j]=nproc;
	  BCSidecu[j]=3;
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new right from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC << " j " << j<<endl;
	}
      }// end ynnu
      XL += grid->getDX();
    }// end XL
  } // end RIGHT
  
    
 }// end check on grid level 
  // end for PRAReceive, modified from initWeightBC 

  //some safety checks: for a grid (remind the PRA limits for the coarser grid), check that its finer grid's PRA does not fall into its PRA
  //it may be a problem with the modifications to the particle mover
  // NB: this part has never been tester
   if (grid->getLevel() < vct->getNgrids()-1){
    bool GoodCondition = (PRA_CoxStartLeft > PRA_oxEndLeft) && (PRA_CoxEndRight< PRA_oxStartRight) && (PRA_CoyStartLeft > PRA_oyEndLeft) && (PRA_CoyEndRight< PRA_oyStartRight);
    if (! GoodCondition){
      cout << "The child grid Particle Repopulation Area overlaps the current grid Particle Repopulation Area: recheck your init parameters...\n Some diagnostics then exiting..." <<endl;
     cout << "Lx: " << Lx << ", 0-dx: " <<0-grid->getDX() << ", Lx+dx: " <<Lx + grid->getDX() <<endl;
     cout << "PRA_CoxStartLeft: " <<PRA_CoxStartLeft <<", PRA_oxEndLeft: " <<PRA_oxEndLeft <<endl;
     if (! (PRA_CoxStartLeft > PRA_oxEndLeft))
       cout <<"Must be: PRA_CoxStartLeft > PRA_oxEndLeft!!!!" <<endl;
     cout << "PRA_CoxEndRight: " <<PRA_CoxEndRight <<", PRA_oxStartRight: " <<PRA_oxStartRight <<endl;
     if (!(PRA_CoxEndRight< PRA_oxStartRight))
       cout <<"Must be: PRA_CoxEndRight< PRA_oxStartRight!!!!" <<endl;
     cout << "Ly: " << Ly << ", 0-dy: " <<0-grid->getDY() << ", Ly+dy: " <<Ly + grid->getDY() <<endl;
     cout << "PRA_CoyStartLeft: " <<PRA_CoyStartLeft <<", PRA_oyEndLeft: " <<PRA_oyEndLeft <<endl;
     if(! (PRA_CoyStartLeft > PRA_oyEndLeft))
       cout <<"Must be: PRA_CoyStartLeft > PRA_oyEndLeft!!!!" <<endl;
     cout << "PRA_CoyEndRight: " <<PRA_CoyEndRight <<", PRA_oyStartRight: " <<PRA_oyStartRight <<endl;
     if (!(PRA_CoyEndRight< PRA_oyStartRight))
       cout<<"Must be: PRA_CoyEndRight< PRA_oyStartRight!!!!" <<endl;
     return -1;
    }
    // other safety check: the PRA should not fall into the parent grid's ghost area
    if (PRA_CoxStartLeft<0 || PRA_CoxEndLeft<0 || PRA_CoxStartRight> Lx || PRA_CoxEndRight> Lx || PRA_CoyStartLeft<0 || PRA_CoyEndLeft<0 || PRA_CoyStartRight> Ly || PRA_CoyEndRight> Ly )
      {
	cout << "The child Particle Repopulation Area is falling into the current grid  ghost area: recheck yout init parameters...\nSome diagnostics then exiting..." << endl;
	cout <<"PRA_CoxStartLeft: " << PRA_CoxStartLeft <<", PRA_CoxEndLeft: " <<PRA_CoxEndLeft <<endl;
	if (PRA_CoxEndLeft<0 || PRA_CoxEndLeft<0)
	  cout <<"PRA_CoxEndLeft and PRA_CoxEndLeft must be >0!!!"<<endl;
	cout <<"PRA_CoxStartRight: " <<PRA_CoxStartRight <<", PRA_CoxEndRight: " <<PRA_CoxEndRight <<endl;
	if (PRA_CoxStartRight> Lx || PRA_CoxEndRight> Lx)
	  cout <<"PRA_CoxStartRight and PRA_CoxEndRight must be <Lx!!!"<<endl;
	cout <<"PRA_CoyStartLeft: " << PRA_CoyStartLeft <<", PRA_CoyEndLeft: " <<PRA_CoyEndLeft<<endl;
        if (PRA_CoyEndLeft<0 || PRA_CoyEndLeft<0)
          cout <<"PRA_CoyEndLeft and PRA_CoyEndLeft must be >0!!!"<<endl;
	cout <<"PRA_CoyStartRight: " <<PRA_CoyStartRight <<", PRA_CoyEndRight: " <<PRA_CoyEndRight <<endl;
	if (PRA_CoyStartRight> Ly || PRA_CoyEndRight> Ly)
          cout <<"PRA_CoyStartRight and PRA_CoyEndRight must be <Ly!!!"<<endl;
	 
	return -1;
      }

  }
//for how the particle motion routines are modified, I need at least 1 PRA cell for side; exit if not//
    if (PRA_Xleft<1 || PRA_Xright<1 || PRA_Yleft<1 || PRA_Xright<1){
    cout << "At least 1 PRA cell per side is needed: modify your input file...\n Exiting...";
    return -1;
  }

  //the PRA area cannot be wider than a fine grid processor
if (grid->getLevel() >0)
  {// nxc is the number of cells in the processor, not the one from input file
    if ( PRA_Xleft-1> grid->getNXC()/2.0 || PRA_Xright-1> grid->getNXC()/2.0 || PRA_Yleft-1> grid->getNYC()/2.0  || PRA_Yright-1> grid->getNYC()/2.0   )
      {
	cout << "The Particle Repopulation Area is too wide compared with the number of cells per processor: revise input parameters.\nSome diagnostics then exiting...";
	cout << "Nxc: " << grid->getNXC() << ", PRA_Xleft: "<< PRA_Xleft <<", PRA_Xright: " <<PRA_Xright <<endl;
	cout << "Nyc: "<< grid->getNYC() << ", PRA_Yleft: "<< PRA_Yleft <<", PRA_Yright: " <<PRA_Yright <<endl;
	return -1;
      }
  }

// check that the number of sends match the number of receives
int ALL_TARGETS;
int ALL_RECEIVERS;

if (grid->getLevel() ==0)
  {
    MPI_Allreduce ( &nmessageBC, &ALL_TARGETS, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM());  
  }
 else
   {
     MPI_Allreduce ( &nmessagerecuBC, &ALL_RECEIVERS, 1,MPI_INT, MPI_SUM, vct->getCART_COMM());
   }
if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())) && grid->getLevel() ==0)
  {
    cout << "Level 0: ALL_TARGETS= " <<  ALL_TARGETS <<endl;
  }
if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())) && !grid->getLevel() ==0)
  {
    cout <<"Level 0: ALL_RECEIVERS= " <<  ALL_RECEIVERS <<endl;
  }
//appropriate iniializations done
MPI_Allreduce( &nmessageBC, &ALL_TARGETS, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce ( &nmessagerecuBC, &ALL_RECEIVERS, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

int ALLTARGETS_LEFT, ALLTARGETS_RIGHT, ALLTARGETS_BOTTOM, ALLTARGETS_TOP;
int ALLRECEIVERS_LEFT,  ALLRECEIVERS_RIGHT,  ALLRECEIVERS_BOTTOM,  ALLRECEIVERS_TOP;

MPI_Allreduce( &targetLEFT, &ALLTARGETS_LEFT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &targetRIGHT, &ALLTARGETS_RIGHT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &targetBOTTOM, &ALLTARGETS_BOTTOM, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &targetTOP, &ALLTARGETS_TOP, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

MPI_Allreduce( &nmessagerecuBCLEFT, &ALLRECEIVERS_LEFT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &nmessagerecuBCRIGHT, &ALLRECEIVERS_RIGHT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &nmessagerecuBCBOTTOM, &ALLRECEIVERS_BOTTOM, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &nmessagerecuBCTOP, &ALLRECEIVERS_TOP, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

if (ALL_TARGETS!=ALL_RECEIVERS)
  {

    
    if (vct->getCartesian_rank_COMMTOTAL()==0)
      {
	cout <<"ALL_TARGETS= " <<ALL_TARGETS <<"!= ALL_RECEIVERS= "<<ALL_RECEIVERS <<endl;
	cout <<"ALLTARGETS_LEFT= " << ALLTARGETS_LEFT << " ALLRECEIVERS_LEFT= " << ALLRECEIVERS_LEFT <<endl;
	cout <<"ALLTARGETS_RIGHT= " << ALLTARGETS_RIGHT << " ALLRECEIVERS_RIGHT= "<< ALLRECEIVERS_RIGHT <<endl;
	cout <<"ALLTARGETS_BOTTOM= " << ALLTARGETS_BOTTOM << " ALLRECEIVERS_BOTTOM= "<< ALLRECEIVERS_BOTTOM <<endl;
	cout <<"ALLTARGETS_TOP= " << ALLTARGETS_TOP << " ALLRECEIVERS_TOP= "<< ALLRECEIVERS_TOP <<endl;
	return -1;
      }
  }

return 1; // OK
}**/ // end old version of initPRAVariables

/**prints some info about the PRA area; only in verbose mode*/
void Particles2Dcomm::checkAfterInitPRAVariables(int species, CollectiveIO* col,VirtualTopology* vct, Grid* grid)
{

  /*if (!cVERBOSE)
    return;*/

  double lengthx=col->getLx()/(double)pow(col->getRatio(),grid->getLevel()); //Total length of the grid in x
  double lengthy=col->getLy()/(double)pow(col->getRatio(),grid->getLevel()); //Total length of the grid in y
  //if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())) and ns ==0) // print only from the lower rank from each level
  if (0)
      {

      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "Preliminary check on PRA variables, rank (local to the grid) " << vct->getCartesian_rank() <<", Level " << grid->getLevel() << ", species " << species <<endl;
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "Local dx: " << grid->getDX() <<", local dy: " << grid->getDY() <<endl;
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "PRA_Xleft: " << PRA_Xleft << ", PRA_Xright: " << PRA_Xright << ", PRA_Yleft: " << PRA_Yleft << ", PRA_Yright: " << PRA_Yright <<endl;
      cout << "Active grid limits (with respect to parent grid): x: "<< grid->getOx(grid->getLevel()) <<"-" << lengthx + grid->getOx(grid->getLevel());
      cout << ", y: " << grid->getOy(grid->getLevel()) << ", " <<lengthy + grid->getOy(grid->getLevel()) << endl;
  
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "PRA limits (with respect to local grid): \n";
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "x left: " <<PRA_oxStartLeft <<"-" << PRA_oxEndLeft << endl;
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "x right: " <<PRA_oxStartRight <<"-" << PRA_oxEndRight << endl;  
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "y left: " <<PRA_oyStartLeft <<"-" << PRA_oyEndLeft << endl;
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
      cout << "y right: " <<PRA_oyStartRight <<"-" << PRA_oyEndRight << endl;

      if (grid->getLevel() < vct->getNgrids()-1){
	double Clengthx=lengthx/col->getRatio();
	double Clengthy=lengthy/col->getRatio();

	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
	cout << "Limits of the active part of the child grid, with respect to the current grid:\n";
	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
	cout << "x: "<< grid->getOx(grid->getLevel()+1) <<"-" << Clengthx + grid->getOx(grid->getLevel()+1);
	cout << ", y: " << grid->getOy(grid->getLevel()+1) << ", " <<Clengthy + grid->getOy(grid->getLevel()+1) << endl;
	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L"  << grid->getLevel() <<":";
	cout << "Grid spacing of the child grid: dx: " << grid->getDX()/col->getRatio() <<", dy: " << grid->getDY()/col->getRatio() <<endl;
	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L"  << grid->getLevel() <<":";
	cout << "Child PRA limits (with respect to local grid): \n";
	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
	cout << "x left: " <<PRA_CoxStartLeft <<"-" << PRA_CoxEndLeft << endl;
	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
	cout << "x right: " <<PRA_CoxStartRight <<"-" << PRA_CoxEndRight << endl;
	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
	cout << "y left: " <<PRA_CoyStartLeft <<"-" << PRA_CoyEndLeft << endl;
	cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "L" << grid->getLevel() <<":";
	cout << "y right: " <<PRA_CoyStartRight <<"-" << PRA_CoyEndRight << endl;
      }
    }

  /*if ( grid->getLevel() < vct->getNgrids()-1)  // for coarse grids                                                                                                                                                 
    {  
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() << ": PRA nmessageBC: " <<nmessageBC <<", sum targets: " << targetBOTTOM+ targetTOP+ targetLEFT+ targetRIGHT <<endl;                                                                                                                      
      for (int i=0; i< nmessageBC; i++)                                                                                                                    
	{                                                                                                                                                  
	  cout << "R" << vct->getCartesian_rank_COMMTOTAL()<<"qom" <<qom << " i: " <<i << ", targetBC[i] "<<targetBC[i] <<", BCSide[i] " <<BCSide[i]<<endl;             
	}                                                                                                                                                  
      cout<< "R" << vct->getCartesian_rank_COMMTOTAL() << " targetBOTTOM " << targetBOTTOM << ", targetTOP " << targetTOP << ", targetLEFT " << targetLEFT << ", targetRIGHT " << targetRIGHT <<endl;   
    }
  if ( grid->getLevel() >0){  // for fine grids
    for (int i=0;i<nmessagerecuBC;i++){
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<",qom" << qom<<  ": receiving PRA BC from "<<fromBC[i] <<", " << i << " of "<< nmessagerecuBC<<endl;
    }
    }*/
}


void Particles2Dcomm::initPRABuffers(Grid* grid, VirtualTopology* vct)
{
  
  np_REPOP_b_BOTTOM =0;
  np_REPOP_b_TOP =0;
  np_REPOP_b_LEFT =0;
  np_REPOP_b_RIGHT =0;

  // Feb 3
  /*setToMINVAL(REPOP_b_BOTTOM);
  setToMINVAL(REPOP_b_TOP);
  setToMINVAL(REPOP_b_LEFT);
  setToMINVAL(REPOP_b_RIGHT);*/

  nop_OS=0;
}
int Particles2Dcomm::PRARepopulationAdd(int PN)
{

  //int prec= 6;
  // PN: index of the particle to examine in the species vector
  if ( y[PN]>= PRA_CoyStartLeft && y[PN]<= PRA_CoyEndLeft && x[PN]>= PRA_CoxStartLeft && x[PN] <= PRA_CoxEndRight ) // area 0, Bottom
    {

      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar]=    x[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +1]= y[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +2]= u[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +3]= v[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +4]= w[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +5]= q[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +6]= xptilde[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +7]= yptilde[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +8]= uptilde[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +9]= vptilde[PN];
      REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +10]= wptilde[PN];

      if (TrackParticleID)
        REPOP_b_BOTTOM[np_REPOP_b_BOTTOM*nVar +11]= ParticleID[PN];

      np_REPOP_b_BOTTOM++;
      if (np_REPOP_b_BOTTOM > (MAX_NP_REPOP_SIZE - (int) (.01*MAX_NP_REPOP_SIZE) ) ){
	cout << "PRARepopulationAdd bottom exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
	return(-1); // end the simulation because you dont have enough space on the array                      
      }
      
      return 1; 
    }
  else if (y[PN]>= PRA_CoyStartRight && y[PN]<= PRA_CoyEndRight && x[PN]>= PRA_CoxStartLeft && x[PN] <= PRA_CoxEndRight ) // area 1, top
    {
     
      REPOP_b_TOP[np_REPOP_b_TOP*nVar]=    x[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +1]= y[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +2]= u[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +3]= v[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +4]= w[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +5]= q[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +6]= xptilde[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +7]= yptilde[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +8]= uptilde[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +9]= vptilde[PN];
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +10]= wptilde[PN];
      if (TrackParticleID)
      REPOP_b_TOP[np_REPOP_b_TOP*nVar +11]= ParticleID[PN];

      np_REPOP_b_TOP++;
      if (np_REPOP_b_TOP > (MAX_NP_REPOP_SIZE - (int) (.01*MAX_NP_REPOP_SIZE) ) ){
        cout << "PRARepopulationAdd top exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
        return(-1); // end the simulation because you dont have enough space on the array                                      
      }
      
      return 1;
    }
  else if ( y[PN]>= PRA_CoyEndLeft && y[PN]<= PRA_CoyStartRight && x[PN]>= PRA_CoxStartLeft && x[PN] <= PRA_CoxEndLeft ) // area 2, left
    {

      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar]=    x[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +1]= y[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +2]= u[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +3]= v[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +4]= w[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +5]= q[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +6]= xptilde[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +7]= yptilde[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +8]= uptilde[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +9]= vptilde[PN];
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +10]= wptilde[PN];
      if (TrackParticleID)
      REPOP_b_LEFT[np_REPOP_b_LEFT*nVar +11]= ParticleID[PN];

      np_REPOP_b_LEFT++;
      if (np_REPOP_b_LEFT > (MAX_NP_REPOP_SIZE - (int) (.01*MAX_NP_REPOP_SIZE) ) ){
        cout << "PRARepopulationAdd left exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
        return(-1); // end the simulation because you dont have enough space on the array                                     \
                                                                                                                   
      } 
      return 1;
    }
  else if (y[PN]>= PRA_CoyEndLeft && y[PN]<= PRA_CoyStartRight && x[PN]>= PRA_CoxStartRight && x[PN] <= PRA_CoxEndRight ) // area 3, right
    {
      
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar]=    x[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +1]= y[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +2]= u[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +3]= v[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +4]= w[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +5]= q[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +6]= xptilde[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +7]= yptilde[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +8]= uptilde[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +9]= vptilde[PN];
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +10]= wptilde[PN];
      if (TrackParticleID)
      REPOP_b_RIGHT[np_REPOP_b_RIGHT*nVar +11]= ParticleID[PN];

      np_REPOP_b_RIGHT++;
      if (np_REPOP_b_RIGHT > (MAX_NP_REPOP_SIZE - (int) (.01*MAX_NP_REPOP_SIZE) ) ){
        cout << "PRARepopulationAdd right exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
        return(-1); // end the simulation because you dont have enough space on the array                                                                                                                         
      } 

      return 1;
    }
  // if it arrives here, particle not to communicate to finer
  return 1;
   
}

void Particles2Dcomm::resizePRAbuffers(int new_max_nop)
{

  double *tmp= new double[nVar* MAX_NP_REPOP_SIZE];


  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {tmp[i]= REPOP_b_BOTTOM[i];}
  delete []REPOP_b_BOTTOM;
  REPOP_b_BOTTOM= new double[nVar* new_max_nop];
  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {REPOP_b_BOTTOM[i]= tmp[i];}
  for (int i=nVar* MAX_NP_REPOP_SIZE; i< nVar* new_max_nop; i++)
    {REPOP_b_BOTTOM[i]= MIN_VAL;}

  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {tmp[i]= REPOP_b_TOP[i];}
  delete []REPOP_b_TOP;
  REPOP_b_TOP= new double[nVar* new_max_nop];
  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {REPOP_b_TOP[i]= tmp[i];}
  for (int i=nVar* MAX_NP_REPOP_SIZE; i< nVar* new_max_nop; i++)
    {REPOP_b_TOP[i]= MIN_VAL;}

  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {tmp[i]= REPOP_b_LEFT[i];}
  delete []REPOP_b_LEFT;
  REPOP_b_LEFT= new double[nVar* new_max_nop];
  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {REPOP_b_LEFT[i]= tmp[i];}
  for (int i=nVar* MAX_NP_REPOP_SIZE; i< nVar* new_max_nop; i++)
    {REPOP_b_LEFT[i]= MIN_VAL;}

  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {tmp[i]= REPOP_b_RIGHT[i];}
  delete []REPOP_b_RIGHT;
  REPOP_b_RIGHT= new double[nVar* new_max_nop];
  for (int i=0; i< nVar* MAX_NP_REPOP_SIZE; i++)
    {REPOP_b_RIGHT[i]= tmp[i];}
  for (int i=nVar* MAX_NP_REPOP_SIZE; i< nVar* new_max_nop; i++)
    {REPOP_b_RIGHT[i]= MIN_VAL;}

  delete[]tmp;

  REPOP_b_BOTTOM_ptr= REPOP_b_BOTTOM;
  REPOP_b_TOP_ptr= REPOP_b_TOP;
  REPOP_b_LEFT_ptr= REPOP_b_LEFT;
  REPOP_b_RIGHT_ptr= REPOP_b_RIGHT;

  MAX_NP_REPOP_SIZE = new_max_nop;
}

int Particles2Dcomm::PRASend(Grid* grid, VirtualTopology* vct)
{
  // if it returns <1, sim stopped gracefully
  int ierr;

  double *tmp= new double[MAX_NP_REPOP_SIZE*nVar];

  int sentBOTTOM=0;
  int sentTOP=0;
  int sentLEFT=0;
  int sentRIGHT=0;

  // number of particles to be sent to each proc in the selected direction
  // (divided equally, not according to the particle position)
  int PshareBOTTOM;
  int PshareTOP;
  int PshareLEFT;
  int PshareRIGHT;


  if (targetBOTTOM>1)
    {PshareBOTTOM=(int) floor(np_REPOP_b_BOTTOM/targetBOTTOM)+1;}
  else if (targetBOTTOM==1)
    {PshareBOTTOM=np_REPOP_b_BOTTOM;}
  else if (targetBOTTOM==0)
    {PshareBOTTOM=0;}

  if (targetTOP>1)
    {PshareTOP=(int) floor(np_REPOP_b_TOP/targetTOP)+1;}
  else if (targetTOP==1)
    {PshareTOP=np_REPOP_b_TOP;}
  else if (targetTOP==0)
    {PshareTOP=0;}

  if (targetLEFT>1)
    {PshareLEFT=(int) floor(np_REPOP_b_LEFT/targetLEFT)+1;}
  else if (targetLEFT==1)
    {PshareLEFT=np_REPOP_b_LEFT;}
  else if (targetLEFT==0)
    {PshareLEFT=0;}

  if (targetRIGHT>1)
    {PshareRIGHT=(int) floor(np_REPOP_b_RIGHT/targetRIGHT)+1;}
  else if (targetRIGHT==1)
    {PshareRIGHT=np_REPOP_b_RIGHT;}
  else if (targetRIGHT==0)
    {PshareRIGHT=0;}

  if (1)
    {
      if (nmessageBC!= targetTOP+targetBOTTOM+targetLEFT+targetRIGHT)
	{
	  cout << "Particles2Dcomm.cpp: problems with sending " <<endl;
	  return -1;
	}
    }
  
  
  bool test_without_part= np_REPOP_b_BOTTOM == 0 && np_REPOP_b_TOP ==0 && np_REPOP_b_LEFT==0 && np_REPOP_b_RIGHT==0;
  if (!test_without_part)
    {
      if (PshareBOTTOM >MAX_NP_REPOP_SIZE || PshareTOP >MAX_NP_REPOP_SIZE || PshareLEFT >MAX_NP_REPOP_SIZE || PshareRIGHT >MAX_NP_REPOP_SIZE)
	{
	  cout << "R" << vct->getCartesian_rank_COMMTOTAL()<< "qom"<< qom<< "Increase the size of the repopulation buffers... \nExiting";
	  fflush(stdout);
	  return -1;
	}
    }
  

  // count the number of sent particles
  if (0)// debug
    {
      int total_p_sent= np_REPOP_b_BOTTOM+ np_REPOP_b_TOP+ np_REPOP_b_LEFT+ np_REPOP_b_RIGHT;
      int total_p_sent_Clevel= 0;
      int level_nop=0;
      MPI_Barrier(vct->getCART_COMM());// at level level                                                               
      MPI_Allreduce(&total_p_sent, &total_p_sent_Clevel, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
      MPI_Allreduce(&nop, &level_nop, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
      if (vct->getCartesian_rank_COMMTOTAL() ==0)
	{
	  cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "Total number of particles sent by Coarse level:  " << total_p_sent_Clevel << ", qom " << qom<<" out of " << level_nop<< endl;
	}
    }
  //int debug_rank= vct->getCartesian_rank_COMMTOTAL(); //for debug
  
  //cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << " PshareBOTTOM " << PshareBOTTOM << " PshareTOP " << PshareTOP << " PshareLEFT " << PshareLEFT << " PshareRIGHT " << PshareRIGHT <<endl;

  for (int i=0; i< nmessageBC; i++)
    {
      if (BCSide[i]==0) //send bottom
	{
	  // arrange tmp for send
	  //setToMINVAL(tmp); //Feb3
	  int PshareBOTTOM_r= PshareBOTTOM;
	  if (sentBOTTOM== targetBOTTOM-1 && targetBOTTOM>1) // to deal with reminders
	    PshareBOTTOM_r= np_REPOP_b_BOTTOM - (PshareBOTTOM*(targetBOTTOM-1));

	  memcpy(tmp, REPOP_b_BOTTOM + sentBOTTOM*PshareBOTTOM*nVar, PshareBOTTOM_r*nVar*sizeof(double));
	  ierr=MPI_Send(tmp,PshareBOTTOM_r*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD);
	  //ierr=MPI_Ssend(tmp,PshareBOTTOM*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD);  
	  sentBOTTOM++;
	}
      else if (BCSide[i]==1) // send top 
	{
          //setToMINVAL(tmp); //Feb3
	  int PshareTOP_r= PshareTOP;
	  if (sentTOP== targetTOP-1 && targetTOP>1) // to deal with reminders            
            PshareTOP_r= np_REPOP_b_TOP -(PshareTOP*(targetTOP-1));
	  memcpy(tmp, REPOP_b_TOP + sentTOP*PshareTOP*nVar, PshareTOP_r*nVar*sizeof(double));
          ierr=MPI_Send(tmp,PshareTOP_r*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD);
	  //ierr=MPI_Ssend(tmp,PshareTOP*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD); 
          sentTOP++;
	}
      else if (BCSide[i]==2) // send left
	{
          //setToMINVAL(tmp); //Feb3
	  int PshareLEFT_r= PshareLEFT;
	  if (sentLEFT== targetLEFT-1 && targetLEFT>1) // to deal with reminders  
            PshareLEFT_r= np_REPOP_b_LEFT -(PshareLEFT*(targetLEFT-1));
	  memcpy(tmp, REPOP_b_LEFT + sentLEFT*PshareLEFT*nVar, PshareLEFT_r*nVar*sizeof(double));
          ierr=MPI_Send(tmp, PshareLEFT_r*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD);
	  //ierr=MPI_Ssend(tmp, PshareLEFT*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD);   
          sentLEFT++;
	}
      else if (BCSide[i]==3) // send right
	{
	  //setToMINVAL(tmp); //Feb3
	  int PshareRIGHT_r= PshareRIGHT;
	  if (sentRIGHT== targetRIGHT-1 && targetRIGHT>1) // to deal with reminders                         
            PshareRIGHT_r= np_REPOP_b_RIGHT -(PshareRIGHT*(targetRIGHT-1));
	  memcpy(tmp, REPOP_b_RIGHT + sentRIGHT*PshareRIGHT*nVar, PshareRIGHT_r*nVar*sizeof(double));
          ierr=MPI_Send(tmp,PshareRIGHT_r*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD);
	  //ierr=MPI_Ssend(tmp,PshareRIGHT*nVar, MPI_DOUBLE, targetBC[i], ns, MPI_COMM_WORLD);  
          sentRIGHT++;
	} 
      // comment away after debug
      else
	{
	  cout << "problem with sending PRA particles, exiting";
	  fflush(stdout);
	  return -1;
	}
    } // end sends
  
  // comment away after debug

  if ((sentBOTTOM != targetBOTTOM) || (sentTOP != targetTOP) || (sentLEFT != targetLEFT) || (sentRIGHT != targetRIGHT))
    {
      cout << "MESS WHEN SENDING PRA PARTICLES!!!, exiting" <<endl;
      fflush(stdout);
      return -1;

    }
  else if (sentBOTTOM + sentTOP + sentLEFT + sentRIGHT !=nmessageBC)
    {
      cout << "MESS WHEN SENDING PRA PARTICLES!!!, exiting" <<endl;
      fflush(stdout);
      return -1;
      }
  
  delete[] tmp;
  return 1;
}

void Particles2Dcomm::setToMINVAL(double *vec)
{// to be used just with vectors of size MAX_NP_REPOP_SIZE*nVar

  //cout << "IN SETTOMINVAL: MAX_NP_REPOP_SIZ " << MAX_NP_REPOP_SIZE <<" nVar " <<nVar << endl;
  memcpy(vec, MIN_VAL_VEC, MAX_NP_REPOP_SIZE*nVar*sizeof(double));
  return;
}

int Particles2Dcomm::PRAReceive(Grid* grid, VirtualTopology *vct, Field* EMf)
{

  
  // return <0 if problems
  MPI_Status status;
  int ierr;
  int out;
  LastCommunicate=0;

  /*if (0)
    {
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<<"qom" <<qom << "Nop before PRAReceive: " << nop << endl;
      }*/

  int nop_BeforeSplit = nop; // for communicateSP 
  n_ref_p_comm_BOTTOM=0;   // n refined grid particles generated from the received ones THAT NEED COMMUNICATION
  n_ref_p_comm_TOP=0;   // n refined grid particles generated from the received ones THAT NEED COMMUNICATION
  n_ref_p_comm_LEFT=0;   // n refined grid particles generated from the received ones THAT NEED COMMUNICATION
  n_ref_p_comm_RIGHT=0;   // n refined grid particles generated from the received ones THAT NEED COMMUNICATION
  GeneratedRP=0;    // how many refined particles are generated by splitting

  AcceptedRPAfterSplit=0; // how many of the generated particles are preliminarly accepted (fall in the PRA)
  AcceptedRPAfterComm=0; // how many of the generated particles are definitively accepted (after the communicateSP)
  RejectedRP=0;     // how many particles are rejected (do not fall in the PRA)          
  CommunicatedRP=0;  // particles generated by the local proc, accepetd and sent to other procs

  int received_p=0;

  // the tag for PRASend / PRAReceive is the n
  for (int i=0; i< nmessagerecuBC; i++)
    {
      int n_rec_p_eachRec=0;// index in REPOP_receive_b   
      //setToMINVAL(REPOP_receive_b);  //Feb3
      
      // receive
      ierr = MPI_Recv(REPOP_receive_b,MAX_NP_REPOP_SIZE*nVar,MPI_DOUBLE, MPI_ANY_SOURCE, ns,MPI_COMM_WORLD, &status);
      // count the number of doubles received
      //Feb3
      int Recv_Double;
      MPI_Get_count(&status, MPI_DOUBLE, &Recv_Double);
      if (Recv_Double>= (MAX_NP_REPOP_SIZE-1)*nVar )
	{
	  cout << "PRAReceive: insufficient buffer size" << endl;
          return -1;
	  }

      while ((n_rec_p_eachRec+1)*nVar<=Recv_Double ) //Feb3

	{
	  //cout << "Test:" << REPOP_receive_b[n_rec_p_eachRec*nVar] << " and " <<REPOP_receive_b[n_rec_p_eachRec*nVar+10]<< endl;
	  out= SplitCoarseParticle(vct, grid, n_rec_p_eachRec);    
	  
	  if (out<0)//in this case, the buffer hosting the generated refined particles is too small
	    {
	      cout << "PRAReceive: insufficient buffer size for particle split" << endl;
	      return -1;
	    }
	  n_rec_p_eachRec++;
	  received_p++;
	}     // end split
      
      
    }   // end receive + split
  
  if (0)// debug                                                                                                                                  
    {
      int total_received_p= 0;
      int total_accepted=0;
      int total_particle=0;
      MPI_Allreduce(&received_p, &total_received_p, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
      MPI_Allreduce( &AcceptedRPAfterSplit,&total_accepted, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
      MPI_Allreduce(&nop, &total_particle, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
      if (vct->getCartesian_rank() == 0) // any proc on this level
        {
          cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "Total number of particles received by Refined llevel:  " << total_received_p << ", qom " << qom<< endl;
	  cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "Total number of particles accepted by Refined llevel:  " << total_accepted << ", qom " << qom<< endl;
	  cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "Total number of particles on Refined llevel:  " << total_particle << ", qom " << qom<< endl;
        }

    }
  //cout << "R" << vct->getCartesian_rank_COMMTOTAL()  <<"received " << received_p << " qom " << qom << endl;
  
  int nop_afterSplit_beforeCommSP=nop;
  int timeCommunicateSP=0;
  // distribuite finer grid particles among the finer grid processors
  int avail = communicateSP(vct, grid);
  //int avail=2;
  timeCommunicateSP++;
  //cout << "R" << vct->getCartesian_rank_COMMTOTAL() << " ended first communicateSP" << endl;
  if (avail < 0)
    {
      cout << "Insufficient buffer size in communicating SP particles" << endl;
      return(-1);
    }
  
  while(isMessagingDoneSP(vct) >0){
    // COMMUNICATION
    avail = communicateSP(vct, grid);
    timeCommunicateSP++;
    if (avail < 0)
      {
	cout << "Insufficient buffer size in communicating SP particles" << endl;
	return(-1);     
      }
    //Nov8 MPI_Barrier(vct->getCART_COMM()); // do I really need this???
  }// end isMessagingDone
  
  if (grid->getLevel()>0 && vct->getRefLevelAdj()==1) // communication OS particles done under condition RefLevelAdj
    {

      cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << "nop_OS before communicateOS: " <<nop_OS <<endl;
      int availOS = communicateOS(vct, grid);
      if (availOS<0)
	{
	  cout << "Insufficient buffer size in communicating OS particles" << endl;
	  return(-1);
	}
      
      MPI_Barrier(vct->getCART_COMM());
      
      while(isMessagingDone(vct) >0){
	// COMMUNICATION                                                            
	availOS = communicateOS(vct, grid);
	if (availOS < 0)
	  {
	    cout << "Insufficient buffer size in communicating OS particles" << endl;
	    return(-1);
	  }
      }// end isMessagingDone   
      cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << "nop_OS after communicateOS: " <<nop_OS <<endl;
    }  // end communicate OS particles
  

  // debug                                                                                                                                               
  int nop_OS_all_after=0;
  // copying the newly received repopulated particles to buffers to set them up again during cycles when no particles from the coarse grid
  if (SubCycling)
    {
      RP_nop= nop- nop_BeforeSplit;

      if (RP_nop > SizeRP_Sub-1)
	{
	  cout << "Increase SizeRP_Sub, exiting" <<endl;
	  return -1;
	}

      memcpy(RP_x, x+nop_BeforeSplit, sizeof(double)*RP_nop);
      memcpy(RP_y, y+nop_BeforeSplit, sizeof(double)*RP_nop);
      memcpy(RP_u, u+nop_BeforeSplit, sizeof(double)*RP_nop);
      memcpy(RP_v, v+nop_BeforeSplit, sizeof(double)*RP_nop);
      memcpy(RP_w, w+nop_BeforeSplit, sizeof(double)*RP_nop);
      memcpy(RP_q, q+nop_BeforeSplit, sizeof(double)*RP_nop);
      if (TrackParticleID)
	memcpy(RP_ParticleID, ParticleID+nop_BeforeSplit, sizeof(unsigned long)*RP_nop);

  }
  if (0)// debug                                 
    {
      int total_particle=0;
      MPI_Allreduce(&nop, &total_particle, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
      if (vct->getCartesian_rank() == 0) // any proc on this level                          
        {
          cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "Total number of particles on Refined llevel after communicateSP:  " << total_particle << ", qom " << qom<< endl;
        }
    }
  return 1;
}
int Particles2Dcomm::SplitCoarseParticle(VirtualTopology* vct, Grid* grid, int n_rec_p)
// n_rec_p is the index of the particle in REPOP_receive_b
{

  // coarse particle values
  double x_C = REPOP_receive_b[n_rec_p * nVar + 0];
  double y_C = REPOP_receive_b[n_rec_p * nVar + 1];
  double u_C = REPOP_receive_b[n_rec_p * nVar + 2];        
  double v_C = REPOP_receive_b[n_rec_p * nVar + 3];                                            
  double w_C = REPOP_receive_b[n_rec_p * nVar + 4];
  double q_C = REPOP_receive_b[n_rec_p * nVar + 5];
  double xptilde_C = REPOP_receive_b[n_rec_p * nVar + 6];
  double yptilde_C = REPOP_receive_b[n_rec_p * nVar + 7];        
  double uptilde_C = REPOP_receive_b[n_rec_p * nVar + 8];        
  double vptilde_C = REPOP_receive_b[n_rec_p * nVar + 9];                                                                
  double wptilde_C = REPOP_receive_b[n_rec_p * nVar + 10];
  //ID
  unsigned long ParticleID_C;
  if (TrackParticleID)
    {            
      ParticleID_C = REPOP_receive_b[n_rec_p * nVar + 11];
      }
  

  // values for each of the finer particle
  double x_F;
  double y_F;
  // velocities                                                                                                                                                     
  unsigned long ParticleID_F;

  double CoarseDX= grid->getDX() * ratio;
  double CoarseDY= grid->getDY() * ratio;


  for (int spx=0; spx<ratio; spx ++)
    {
      // actual splitting
      
      x_F= x_C - CoarseDX/2.0 + grid->getDX()/2.0 + grid->getDX()* spx - grid->getOx(grid->getLevel());
      for (int spy=0; spy< ratio; spy++)
	{
	  y_F= y_C - CoarseDY/2.0 + grid->getDY()/2.0 +grid->getDY()* spy - grid->getOy(grid->getLevel());
	  GeneratedRP++;
	  // consider generated particle only if it falls in the repopulation area of the fine grid
	  //if (x_F< PRA_oxStartLeft || (x_F> PRA_oxEndLeft && x_F< PRA_oxStartRight) || x_F> PRA_oxEndRight || y_F< PRA_oyStartLeft || (y_F> PRA_oyEndLeft && y_F< PRA_oyStartRight) || y_F> PRA_oyEndRight )
	    
	  bool BOTTOM= (x_F> PRA_oxStartLeft && x_F< PRA_oxEndRight && y_F> PRA_oyStartLeft && y_F< PRA_oyEndLeft); 
	  bool TOP=    (x_F> PRA_oxStartLeft && x_F< PRA_oxEndRight && y_F> PRA_oyStartRight && y_F< PRA_oyEndRight);
	  bool LEFT=   (x_F> PRA_oxStartLeft && x_F< PRA_oxEndLeft && y_F> PRA_oyEndLeft && y_F< PRA_oyStartRight);
	  bool RIGHT=  (x_F> PRA_oxStartRight && x_F< PRA_oxEndRight && y_F> PRA_oyEndLeft && y_F< PRA_oyStartRight);
	  //cout << "BOTTOM " << BOTTOM << " TOP " << TOP <<" LEFT " << LEFT << " RIGHT " << RIGHT << endl;

	  if (! (BOTTOM || TOP || LEFT || RIGHT) )
	    {
	      // generated particle is outside of the repopulation area: reject
	      RejectedRP++;
	      // 0 for debugging OS part
	      if (grid->getLevel()>0 && vct->getRefLevelAdj()==1 )
		{
		  double dx= grid->getDX();
		  double dy= grid->getDY();
	      
		  // here, check if the particle is eligible for the OS buffers (OS= Other Side)
		  bool BOTTOM_OS=x_F> PRA_oxStartLeft- dx  && x_F< PRA_oxEndRight + dx && y_F> PRA_oyStartLeft - dy  && y_F< PRA_oyStartLeft;
		  bool TOP_OS=x_F> PRA_oxStartLeft- dx && x_F< PRA_oxEndRight + dx && y_F> PRA_oyEndRight && y_F< PRA_oyEndRight + dy;
		  bool LEFT_OS=x_F> PRA_oxStartLeft- dx && x_F<PRA_oxStartLeft && y_F> PRA_oyStartLeft && y_F< PRA_oyEndRight;
		  bool RIGHT_OS=x_F> PRA_oxEndRight && x_F< PRA_oxEndRight + dx && y_F> PRA_oyStartLeft && y_F< PRA_oyEndRight;
		  
		  if (BOTTOM_OS || TOP_OS || LEFT_OS || RIGHT_OS)
		    {
		      OS_x[nop_OS]=x_F;
		      OS_y[nop_OS]=y_F;
		      OS_u[nop_OS]=u_C;
		      OS_v[nop_OS]=v_C;
		      OS_w[nop_OS]=w_C;
		      OS_q[nop_OS]=q_C/ratio/ratio;
		      if (TrackParticleID){
			OS_ParticleID[nop_OS]=20000000;
		      }
		      
		      if (nop_OS > (npmax_OS - (int) (.01*npmax_OS) ) ){
			cout << "SplitCoarseParticle, OS bottom ops exceeding npmax_OS: Particles need to be resized Save Data and Stop the simulation" << endl;
			return(-1); // end the simulation because you dont have enough space on the array                                                      		  
		      }
		      nop_OS++;
		    }
		}// end OS particles
	      continue;
	    } // end if (! (BOTTOM || TOP || LEFT || RIGHT) )
	  AcceptedRPAfterSplit++; 
	  x[nop] = x_F;
	  y[nop] = y_F;
     
	  u[nop] = u_C;    // didn't even bother to create u_F
	  v[nop] = v_C;
	  w[nop] = w_C;
	  q[nop] = q_C/ratio/ratio;

	  xptilde[nop] = xptilde_C;
	  yptilde[nop] = yptilde_C;
	  uptilde[nop] = uptilde_C;
	  vptilde[nop] = vptilde_C;
	  wptilde[nop] = wptilde_C;
	  if (TrackParticleID)
	    {ParticleID[nop]= 10000000;}  // something better for this
	  
	  if (nop > (npmax - (int) (.01*npmax) ) ){
	    cout << "SplitCoarseParticle exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
	    return(-1); // end the simulation because you dont have enough space on the array                                  
	  }
	  nop++;
	  
	  // end actual splitting & storing
     
	}  // end spy cycle
    }// end spx cycle
  
  return 1;
}

void Particles2Dcomm::setToMINVAL_SP(double *vec, int np_SplitPartComm)
{// to be used only for vectors of size MAX_NP_SPLIPARTCOMM* nVar
  memcpy(vec, MIN_VAL_VEC_SP, np_SplitPartComm*nVar*sizeof(double));
  return;
}

void Particles2Dcomm::setToMINVAL_OS(double *vec, int np)
{// to be used only for vectors of size MAX_NP_OSPARTCOMM* nVarOS  
  memcpy(vec, MIN_VAL_VEC_OS, np*nVarOS*sizeof(double));
  return;
}

int Particles2Dcomm::pack_SP(double *vec, int part_index, double x_F, double y_F, double u_F, double v_F, double w_F, double q_F, unsigned long ParticleID_F)
{
  // dimension check
  if ((part_index > (MAX_NP_SPLIPARTCOMM - (int) (.01*MAX_NP_SPLIPARTCOMM) ) ) )
    {
      cout << "pack_SP: resize particle buffer, exiting..." << endl;
      return -1;
    }
  // pack
  vec[(part_index)*nVar+ 0] = x_F;
  vec[(part_index)*nVar+ 1] = y_F;
  vec[(part_index)*nVar+ 2] = u_F;
  vec[(part_index)*nVar+ 3] = v_F;
  vec[(part_index)*nVar+ 4] = w_F;
  vec[(part_index)*nVar+ 5] = q_F;
  vec[(part_index)*nVar+ 6] = 0.0; // xptilde
  vec[(part_index)*nVar+ 7] = 0.0; // yptilde
  vec[(part_index)*nVar+ 8] = 0.0; // uptilde
  vec[(part_index)*nVar+ 9] = 0.0; // vptilde
  vec[(part_index)*nVar+ 10] = 0.0; // wptilde
  if (TrackParticleID)
    {vec[(part_index)*nVar+ 11]= ParticleID_F;} // ParticleID
  // the counter in incremented outside

  return 1;
}

//int Particles2Dcomm::communicateSP(VirtualTopology* ptVCT){
int Particles2Dcomm::communicateSP(VirtualTopology* ptVCT, Grid* grid){ 
  // to distribute the splitted particles into adjacent procs; the accepted splitted particles were temporarely added to the x, .. of the current proc
  MPI_Status status;
  int new_np_SplitPartComm;
  int npExitingMax;
  // variable for memory availability of space for new particles
  int avail, availALL, avail1, avail2, avail3, avail4;


  //Feb 4
  /*setToMINVAL_SP(SplittedParticles_Comm_BOTTOM, max_np_SplitPartComm); 
  setToMINVAL_SP(SplittedParticles_Comm_TOP, max_np_SplitPartComm);                
  setToMINVAL_SP(SplittedParticles_Comm_LEFT, max_np_SplitPartComm);
  setToMINVAL_SP(SplittedParticles_Comm_RIGHT, max_np_SplitPartComm);  */
  // end Feb 4
  //cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "entered communicate SP" <<endl;
  npExitXright =0, npExitXleft =0, npExitYright =0, npExitYleft =0, npExit=0, rightDomain = 0, rightDomainX=0, rightDomainY=0;
  npDeletedBoundary = 0;


  /*MPI_Barrier(ptVCT->getCART_COMM());
  if (ptVCT->getCartesian_rank() ==0 || ptVCT->getCartesian_rank_COMMTOTAL() == 2000)
  cout <<"communictaeSP: before for on particles (just to check if everybody gets here)\n";*/
  

  int np_current=0;
  //int np_current=0;
  int nplast = nop-1;
  
  while (np_current < nplast+1){
    
    if ( (x[np_current] < xstart && ptVCT->getCoordinates(0) != 0 ) || (x[np_current] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1)) || (x[np_current] < Modified_xstart && ptVCT->getCoordinates(0) == 0 )|| (x[np_current] > Modified_xend && ptVCT->getCoordinates(0) == (ptVCT->getXLEN()-1))     ){
      // communicate if they don't belong to the domain
      if (x[np_current] < xstart && ptVCT->getCoordinates(0) != 0){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitXleft+1)>= max_np_SplitPartComm){ 
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size" << endl;
	  if (!resize_buffersSP((int) (max_np_SplitPartComm*2)))
	    {
	      cout<<"communicateSP: increase MAX_NP_SPLIPARTCOMM\nExiting...";
	      return -1;
	    }
	  cout << "End R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size, Xleft" << endl;
	}
	// put it in the communication buffer
	bufferXleft(SplittedParticles_Comm_LEFT,np_current,ptVCT);
	// delete the particle and pack the particle array, the value of nplast changes
	del_pack(np_current,&nplast);
	npExitXleft++;
      } else if (x[np_current] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1)){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitXright+1)>=max_np_SplitPartComm){
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size" << endl;
          if (!resize_buffersSP((int) (max_np_SplitPartComm*2)))
            {
              cout<<"communicateSP: increase MAX_NP_SPLIPARTCOMM\nExiting...";
              return -1;
            }
	  cout << "End R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size, Xright" << endl;
	}
	// put it in the communication buffer
	bufferXright(SplittedParticles_Comm_RIGHT,np_current,ptVCT); 
	// delete the particle and pack the particle array, the value of nplast changes
	del_pack(np_current,&nplast);	
	npExitXright++;
      } /*else if (x[np_current] < Modified_xstart && ptVCT->getCoordinates(0) == 0 ){
	cout <<"CommunicateSP: problem with particle position, exiting";
	cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << " x[np_current] " << x[np_current] << " Modified_xstart "<< Modified_xstart << " Lx " <<Lx<< endl;
	return -1;
      } else if (x[np_current] > Modified_xend && ptVCT->getCoordinates(0) == (ptVCT->getXLEN()-1)){
	cout <<"CommunicateSP: problem with particle position, exiting";
	cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << " x[np_current] " << x[np_current]<< " Modified_xend " << Modified_xend << " Lx " << Lx << " Lx + dx " << Lx + grid-> getDX ()<< endl;
	return -1;
	}*/
    } 
    else if ( (y[np_current] < ystart && ptVCT->getCoordinates(1) != 0 ) || (y[np_current] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1)) || (y[np_current] < Modified_ystart && ptVCT->getCoordinates(1) == 0 )|| (y[np_current] > Modified_yend && ptVCT->getCoordinates(1) == (ptVCT->getYLEN()-1))     ){
      // communicate if they don't belong to the domain
      if (y[np_current] < ystart && ptVCT->getCoordinates(1) != 0){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitYleft+1)>=max_np_SplitPartComm){
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size" << endl;
          if (!resize_buffersSP((int) (max_np_SplitPartComm*2)))
            {
              cout<<"communicateSP: increase MAX_NP_SPLIPARTCOMM\nExiting...";
              return -1;
            }
	  cout << "End R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size, yleft" << endl;
	}
	// put it in the communication buffer
	bufferYleft(SplittedParticles_Comm_BOTTOM,np_current,ptVCT); 
	// delete the particle and pack the particle array, the value of nplast changes
	del_pack(np_current,&nplast);
	npExitYleft++;
      } else if (y[np_current] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1)){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitYright+1)>=max_np_SplitPartComm){
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size" << endl;
          if (!resize_buffersSP((int) (max_np_SplitPartComm*2)))
            {
              cout<<"communicateSP: increase MAX_NP_SPLIPARTCOMM\nExiting...";
              return -1;
            }
	  cout << "End R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateSP doubling the sending buffer size, Yright" << endl;
	}
	// put it in the communication buffer
	bufferYright(SplittedParticles_Comm_TOP ,np_current,ptVCT); 
	// delete the particle and pack the particle array, the value of nplast changes
	del_pack(np_current,&nplast);
	npExitYright++;
      }
    }  
    else {
      // particle is still in the domain, procede with the next particle
      np_current++;
      AcceptedRPAfterComm++;
    }
  } // end while

  
  
  nop = nplast+1;
  if (nop > (npmax - (int) (.01*npmax) ) ){
    cout << "communicateSP exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
    return(-1); // end the simulation because you dont have enough space on the array                                                                          
  }
  npExitingMax = 0;
  // calculate the maximum number of particles exiting from this domain
  // use this value to check if communication is needed
  // and to  resize the buffer
  npExitingMax = maxNpExiting();
  // broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
  //npExitingMax = reduceMaxNpExiting(ptVCT->getCART_COMM(),npExitingMax);

  int local_npExitingMaxX= ((npExitXleft > npExitXright) ? npExitXleft : npExitXright);
  int local_npExitingMaxY= ((npExitYleft > npExitYright) ? npExitYleft : npExitYright);

  /** lots of checks, remove later **/
  if ((npExitXleft>0 or npExitXright>0 or npExitYleft>0 or npExitYright>0) and (nmessagerecuBC==0))
    {
      cout<< "There must be a bug, exiting\n";
      return -1;
    }
  if ((npExitXleft>0 or npExitXright>0) and (ptVCT->getCOMM_B_BOTTOM()== MPI_COMM_NULL and ptVCT->getCOMM_B_TOP()== MPI_COMM_NULL))
    {
      cout<< "There must be a bug, exiting\n";
      return -1;
    }
  if ((npExitYleft>0 or npExitYright>0) and (ptVCT->getCOMM_B_LEFT()== MPI_COMM_NULL and ptVCT->getCOMM_B_RIGHT()== MPI_COMM_NULL))
    {
      cout<< "There must be a bug, exiting\n";
      return -1;
    }  
  /** end checks **/
  int npExitingMaxX=0, npExitingMaxY=0;

  if (ptVCT->getCOMM_B_BOTTOM()!= MPI_COMM_NULL)
    npExitingMaxX= reduceMaxNpExiting(ptVCT->getCOMM_B_BOTTOM(),local_npExitingMaxX);   
  if (ptVCT->getCOMM_B_TOP()!= MPI_COMM_NULL)
    npExitingMaxX= reduceMaxNpExiting(ptVCT->getCOMM_B_TOP(),local_npExitingMaxX);
  if (ptVCT->getCOMM_B_LEFT()!= MPI_COMM_NULL)
    npExitingMaxY= reduceMaxNpExiting(ptVCT->getCOMM_B_LEFT(),local_npExitingMaxY);
  if (ptVCT->getCOMM_B_RIGHT()!= MPI_COMM_NULL)
    npExitingMaxY= reduceMaxNpExiting(ptVCT->getCOMM_B_RIGHT(),local_npExitingMaxY);

  npExitingMax=0;
  npExitingMax=((npExitingMaxX > npExitingMaxY) ? npExitingMaxX : npExitingMaxY);

  //int npExitingMaxX= reduceMaxNpExiting(ptVCT->getCART_COMM(),local_npExitingMaxX);
  //int npExitingMaxY= reduceMaxNpExiting(ptVCT->getCART_COMM(),local_npExitingMaxY); 

  /*****************************************************/
  /*           SEND AND RECEIVE MESSAGES               */
  /*****************************************************/
  
  
  //  new_np_SplitPartComm = npExitingMax + 1;
  new_np_SplitPartComm = max_np_SplitPartComm;
  while (npExitingMax +1 > new_np_SplitPartComm){        
    new_np_SplitPartComm= (int) (new_np_SplitPartComm*2);
  }
  if (new_np_SplitPartComm != max_np_SplitPartComm)
    {
      cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL()<<"ns " << ns << " communicateSP resizing the receiving buffer" << endl;                            
    if (!resize_buffersSP((int) new_np_SplitPartComm) )    
      {                                                                                                                                      
        cout<<"communicateSP: increase MAX_NP_SPLIPARTCOMM\nExiting...";                                                                     
        return -1;                                                                                                                           
      }
  }

  MPI_Comm Comm;
  int CommID; //commID:  0:left, 1:right, 2:bottom, 3:top
  if (npExitingMaxX>0 and (ptVCT->getCOMM_B_BOTTOM()!= MPI_COMM_NULL or ptVCT->getCOMM_B_TOP()!= MPI_COMM_NULL)){
    if (ptVCT->getCOMM_B_BOTTOM()!= MPI_COMM_NULL) 
      {
	Comm= ptVCT->getCOMM_B_BOTTOM();
	CommID=2;
      }
    else if (ptVCT->getCOMM_B_TOP()!= MPI_COMM_NULL)
      {
	Comm= ptVCT->getCOMM_B_TOP();
	CommID=3;
      }// if both are different from MPI_COMM_NULL, YLEN=1 and the communicators point to the same groups of procs; then, only BOTTOM will perform all the communication

    //Feb 4
    SplittedParticles_Comm_LEFT[nVar*npExitXleft]= MIN_VAL;
    SplittedParticles_Comm_RIGHT[nVar*npExitXright]=MIN_VAL;
    //end Feb 4

    if(!communicateParticles_BoundaryComm(CommID,(npExitingMaxX+1)*nVar,SplittedParticles_Comm_LEFT, SplittedParticles_Comm_RIGHT,ptVCT))
      return -1;
    //cout << "R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "unbuffer in communicateSP" <<endl;          
    // UNBUFFERING                                                                               
    avail1 = unbuffer(SplittedParticles_Comm_LEFT, ptVCT);                                             
    avail2 = unbuffer(SplittedParticles_Comm_RIGHT, ptVCT);                               
    // if one of these numbers is negative than there is not enough space for particles                        
    avail = avail1 + avail2;
    //cout << "R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "before reduceNumberParticles in communicateSP" <<endl;
    availALL = reduceNumberParticles(Comm,avail);                     
    //cout << "R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "after reduceNumberParticles in communicateSP" <<endl;
    if (availALL < 0)                                                                             
      {                                                                       
        cout << "CommunicateSP, top and bottom, too many particles coming: exiting " <<endl;                       
        return(-1);  // too many particles coming, save data nad stop simulation          
      }                                                                                             
  }

  if (npExitingMaxY>0 and (ptVCT->getCOMM_B_LEFT()!= MPI_COMM_NULL or ptVCT->getCOMM_B_RIGHT()!= MPI_COMM_NULL)){
    if (ptVCT->getCOMM_B_LEFT()!= MPI_COMM_NULL)
      {
        Comm= ptVCT->getCOMM_B_LEFT();
        CommID=0;
      }
    else if (ptVCT->getCOMM_B_RIGHT()!= MPI_COMM_NULL)
      {
        Comm= ptVCT->getCOMM_B_RIGHT();
        CommID=1;
      }// if both are different from MPI_COMM_NULL, YLEN=1 and the communicators point to the same groups of procs; then, only BOTTOM will perform all the communication       

    //Feb 4
    SplittedParticles_Comm_BOTTOM[nVar*npExitYleft]=MIN_VAL;
    SplittedParticles_Comm_TOP[nVar*npExitYright]=MIN_VAL;
    // end Feb 4

    if(!communicateParticles_BoundaryComm(CommID,(npExitingMaxY+1)*nVar,SplittedParticles_Comm_BOTTOM, SplittedParticles_Comm_TOP,ptVCT))
      return -1;
    //cout << "R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "unbuffer in communicateSP" <<endl;
    // UNBUFFERING                                                                                                                                                
    avail1 = unbuffer(SplittedParticles_Comm_BOTTOM, ptVCT);
    avail2 = unbuffer(SplittedParticles_Comm_TOP, ptVCT);
    // if one of these numbers is negative than there is not enough space for particles                                                                       
    avail = avail1 + avail2;
    availALL = reduceNumberParticles(Comm,avail);
    if (availALL < 0)
      {
	cout << "CommunicateSP, left and right, too many particles coming: exiting " <<endl;
	return(-1);  // too many particles coming, save data nad stop simulation                                                                                  
      }
  }

  return(0); // everything was fine
  
}


      
int Particles2Dcomm::CountPRAParticles(VirtualTopology* vct)
{
  int npTOP=0;
  int npBOTTOM=0;
  int npLEFT=0; 
  int npRIGHT=0;
  int nonPRA=0;

  for (int i=0; i<nop; i++)
  {
    // BOTTOM
    if (x[i]>= PRA_oxStartLeft && x[i]<= PRA_oxEndRight && y[i]>= PRA_oyStartLeft && y[i]<= PRA_oyEndLeft  )
    {
      npBOTTOM++;
    }
    else if (x[i]>= PRA_oxStartLeft && x[i]<= PRA_oxEndRight && y[i]>= PRA_oyStartRight && y[i]<= PRA_oyEndRight  )
    {
      npTOP++;
    }
    else if (x[i]>= PRA_oxStartLeft && x[i]<= PRA_oxEndLeft && y[i]>= PRA_oyEndLeft && y[i]<= PRA_oyStartRight   )
    {
      npLEFT++;
    }
    else if (x[i]>= PRA_oxStartRight && x[i]<= PRA_oxEndRight && y[i]>= PRA_oyEndLeft && y[i]<= PRA_oyStartRight   )
    {
      npRIGHT++;
    }
    else if (x[i]< PRA_oxStartLeft || x[i]> PRA_oxEndRight || y[i]< PRA_oyStartLeft || y[i]> PRA_oyEndRight ) 
    {
      cout <<"\n\nR" <<vct->getCartesian_rank_COMMTOTAL() << "There shouldn't be anybody here, exiting\n\n\n";
      return -1;
    }
    else
    {
      nonPRA++;
    }

  }

  if (npBOTTOM+ npTOP+ npLEFT+ npRIGHT + nonPRA!= nop)
    {
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << "Error during counting of particles\n";
      cout << "npBOTTOM+ npTOP+ npLEFT+ npRIGHT + nonPRA=  " << npBOTTOM+ npTOP+ npLEFT+ npRIGHT + nonPRA << ", nop= "<< nop <<endl;
      return -1;
    }
  //aggregate
  int allBOTTOM=-1, allTOP=-1, allLEFT=-1, allRIGHT=-1, allNonPRA=-1, allNop=-1;

  MPI_Barrier(vct->getCART_COMM());
  MPI_Allreduce(&npBOTTOM, &allBOTTOM, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
  MPI_Allreduce(&npTOP, &allTOP, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
  MPI_Allreduce(&npLEFT, &allLEFT, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
  MPI_Allreduce(&npRIGHT, &allRIGHT, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
  MPI_Allreduce(&nonPRA, &allNonPRA, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
  MPI_Allreduce(&nop, &allNop, 1, MPI_INT, MPI_SUM, vct->getCART_COMM());
  
  if (vct->getCartesian_rank_COMMTOTAL() == vct->getXLEN()*vct->getYLEN())
    {
      cout << "Refined level, species "<<ns <<": particles in BOTTOM: " << allBOTTOM << ", TOP: " << allTOP <<", LEFT: "<< allLEFT << ", RIGHT: "<<allRIGHT <<", nonPRA: " <<allNonPRA << ", nop: "<<allNop << endl; 
    }
  
  return 1;
}

int Particles2Dcomm::printPRAparticles(VirtualTopology* ptVCT, Grid* grid)
{
  if (grid->getLevel()==0)
  {
    if (np_REPOP_b_BOTTOM>0)
    {
      cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL()<< "np_REPOP_b_BOTTOM:" <<np_REPOP_b_BOTTOM<< endl;
      for (int i=0; i<np_REPOP_b_BOTTOM; i++)
      {
	cout<<"\t\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":x " << REPOP_b_BOTTOM[i*nVar+0] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":y " << REPOP_b_BOTTOM[i*nVar+1] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":u " << REPOP_b_BOTTOM[i*nVar+2] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":v " << REPOP_b_BOTTOM[i*nVar+3] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":w " << REPOP_b_BOTTOM[i*nVar+4] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":q " << REPOP_b_BOTTOM[i*nVar+5] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":xptilde " << REPOP_b_BOTTOM[i*nVar+6] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":yptilde " << REPOP_b_BOTTOM[i*nVar+7] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":uptilde " << REPOP_b_BOTTOM[i*nVar+8] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":vptilde " << REPOP_b_BOTTOM[i*nVar+9] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":wptilde " << REPOP_b_BOTTOM[i*nVar+10] ;

	if (TrackParticleID)
	  cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Bi " << i <<":ID " << REPOP_b_BOTTOM[i*nVar+11] <<endl;
      }
    } 
    if (np_REPOP_b_TOP>0)
    {
      cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL()<< "np_REPOP_b_TOP:" <<np_REPOP_b_TOP << endl;
      for (int i=0; i<np_REPOP_b_TOP; i++)
      {
	cout<<"\t\t\t\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":x " << REPOP_b_TOP[i*nVar+0] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":y " << REPOP_b_TOP[i*nVar+1] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":u " << REPOP_b_TOP[i*nVar+2] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":v " << REPOP_b_TOP[i*nVar+3] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":w " << REPOP_b_TOP[i*nVar+4] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":q " << REPOP_b_TOP[i*nVar+5] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":xptilde " << REPOP_b_TOP[i*nVar+6] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":yptilde " << REPOP_b_TOP[i*nVar+7] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":uptilde " << REPOP_b_TOP[i*nVar+8] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":vptilde " << REPOP_b_TOP[i*nVar+9] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":wptilde " << REPOP_b_TOP[i*nVar+10] ;
	
	if (TrackParticleID)
	  cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ti " << i <<":ID " << REPOP_b_TOP[i*nVar+11] <<endl;
      }
    }
    if (np_REPOP_b_LEFT>0)
    {
      cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL()<< "np_REPOP_b_LEFT:" <<np_REPOP_b_LEFT<< endl;
      for (int i=0; i<np_REPOP_b_LEFT; i++)
      {
	cout<<"\t\t\t\t\t\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":x " << REPOP_b_LEFT[i*nVar+0] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":y " << REPOP_b_LEFT[i*nVar+1] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":u " << REPOP_b_LEFT[i*nVar+2] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":v " << REPOP_b_LEFT[i*nVar+3] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":w " << REPOP_b_LEFT[i*nVar+4] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":q " << REPOP_b_LEFT[i*nVar+5] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":xptilde " << REPOP_b_LEFT[i*nVar+6] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":yptilde " << REPOP_b_LEFT[i*nVar+7] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":uptilde " << REPOP_b_LEFT[i*nVar+8] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":vptilde " << REPOP_b_LEFT[i*nVar+9] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":wptilde " << REPOP_b_LEFT[i*nVar+10] ;
	
	if (TrackParticleID)
	  cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Li " << i <<":ID " << REPOP_b_LEFT[i*nVar+11] <<endl;
      }
    }
    if (np_REPOP_b_RIGHT>0)
    {
      cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL()<< "np_REPOP_b_RIGHT:"<<np_REPOP_b_RIGHT << endl;
      for (int i=0; i<np_REPOP_b_RIGHT; i++)
      {
	cout<<"\t\t\t\t\t\t\t\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":x " << REPOP_b_RIGHT[i*nVar+0] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":y " << REPOP_b_RIGHT[i*nVar+1] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":u " << REPOP_b_RIGHT[i*nVar+2] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":v " << REPOP_b_RIGHT[i*nVar+3] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":w " << REPOP_b_RIGHT[i*nVar+4] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":q " << REPOP_b_RIGHT[i*nVar+5] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":xptilde " << REPOP_b_RIGHT[i*nVar+6] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":yptilde " << REPOP_b_RIGHT[i*nVar+7] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":uptilde " << REPOP_b_RIGHT[i*nVar+8] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":vptilde " << REPOP_b_RIGHT[i*nVar+9] ;
	cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":wptilde " << REPOP_b_RIGHT[i*nVar+10] ;
	
	if (TrackParticleID)
	  cout<<"\tR" <<ptVCT->getCartesian_rank_COMMTOTAL()<<" Ri " << i <<":ID " << REPOP_b_RIGHT[i*nVar+11] <<endl;
      }
    }
  }// end level
  return 1;
}

void Particles2Dcomm::setPRACollectionMethod(int which)
{
  PRACollectionMethod= which;
}
int Particles2Dcomm::getPRACollectionMethod()
{
  return PRACollectionMethod;
}
int Particles2Dcomm::CollectivePRARepopulationAdd(VirtualTopology* ptVCT, Grid* grid)
{
  int out;
  //if (! ( grid->getLevel()==0 || (grid->getLevel()>0 && FinerLevels_PRAOps==1 )) )
  if (! ( grid->getLevel()==0 || (grid->getLevel()>0 && FinerLevels_PRAOps==1 ))  || PRAIntersection== false ) 
  {
    //cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL()<<": no need for CollectivePRARepopulationAdd\n";
    return 1;
  }


  //cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL()<<": CollectivePRARepopulationAdd\n";
  for (int i=0; i<nop; i++)
  {
   out=PRARepopulationAdd(i);
   if (out<0)
   {
     cout << "Insufficient repopulation buffer size, exiting" << endl;
     return -1;
   }
  }


  return 1;
}

bool Particles2Dcomm::applyParticleBC(int BC_partCommunicate, int np_current, int *nplast, int *npDeletedBoundary, Grid * grid, VirtualTopology * ptVCT)
{
  // values of BC_partCommunicate
  // 0: coarse grid  
  // 1: refined grid, inner mover  
  // 2: refined grid, final position

  /*if (0 && ptVCT->getCartesian_rank_COMMTOTAL() <16 && ParticleID[np_current]== 13100 && qom== -256)
    {
      cout << "Particle " << ParticleID[np_current] <<" is in applyParticleBC, BC_partCommunicate: " <<BC_partCommunicate << endl;
      }*/
  
  if (BC_partCommunicate==0)
    {// coarse grid
      if (  (x[np_current] < 0- grid->getDX() && !ptVCT->getPERIODICX()) || (x[np_current] < 0 && ptVCT->getPERIODICX()) ){ 
	BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft, grid->getDX(), ptVCT, ptVCT->getPERIODICX()); 
      }
      else if ( ( x[np_current] > Lx+grid->getDX() && !ptVCT->getPERIODICX()) || ( x[np_current] > Lx && ptVCT->getPERIODICX())  ) { 
	BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft, grid->getDX(), ptVCT, ptVCT->getPERIODICX());
      }
      if ( (y[np_current] < 0- grid->getDY() && !ptVCT->getPERIODICY()) || (y[np_current] < 0 && ptVCT->getPERIODICY())  ){  
	BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft, grid->getDY(), ptVCT, ptVCT->getPERIODICY());  
      }
      else if ( (y[np_current] > Ly +grid->getDY() && !ptVCT->getPERIODICY()) || (y[np_current] > Ly  && ptVCT->getPERIODICY())  ){  
	BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft, grid->getDY(), ptVCT, ptVCT->getPERIODICY());
      }
      return false; //skip to false, do not skip the other chekcs on the particle
    }
  else if (BC_partCommunicate==1)
    {// refined grid, inner mover
      // if particles are outside the PRA, deleted
      if (x[np_current] < PRA_oxStartLeft){
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate
      }
      else if (x[np_current] > PRA_oxEndRight){
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate
      }
      if (y[np_current] < PRA_oyStartLeft ) {
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate
      }
      else if (y[np_current] > PRA_oyEndRight ) {
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate
      }
      return false; // particle passed the preliminary tests, check if buffering needed
    }
  else if (BC_partCommunicate==2)
    {// refined grid, final position
      if (x[np_current] < PRA_oxEndLeft){
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate
      }
      else if (x[np_current] > PRA_oxStartRight){
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate
      }
      if (y[np_current] < PRA_oyEndLeft ) {
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate  
      }
      else if (y[np_current] > PRA_oyStartRight ) {
	del_pack(np_current,&(*nplast));
	(*npDeletedBoundary)++;
	//continue;
	return true; //skip to true, skip the rest of communicate 
      }
      return false; // particle passed the preliminary tests, check if buffering needed 
    }// end check on BC
  cout << "\n\n\n\n\nPROBLEM: THERE SHOULD BE NOTHING THERE\n\n\n\n\n\n";
  return false;
}
// to communicate particles for OS operations (patching the "other side" of the moments on the refined grid) to the right processor   
int Particles2Dcomm::communicateOS(VirtualTopology* ptVCT, Grid* grid)
{

  // particles also in the ghost cells also for the coarse grids
  // allocate buffers
  MPI_Status status;
  int new_OSpart;
  int npExitingMax;
  // variable for memory availability of space for new particles
  int avail, availALL, avail1, avail2, avail3, avail4;
  setToMINVAL_OS(OSParticles_Comm_BOTTOM, max_np_OsPartComm);// particles going down
  setToMINVAL_OS(OSParticles_Comm_TOP, max_np_OsPartComm);  // particles going up
  setToMINVAL_OS(OSParticles_Comm_LEFT, max_np_OsPartComm); // particles going left
  setToMINVAL_OS(OSParticles_Comm_RIGHT, max_np_OsPartComm);   // particles going right

  //cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << " nop_OS bef sorting in communicateOS "<< nop_OS << endl;

  npExitXright =0, npExitXleft =0, npExitYright =0, npExitYleft =0, npExit=0, rightDomain = 0, rightDomainX=0, rightDomainY=0;
  npDeletedBoundary = 0; 

  // deal with the Bottom quantities (particles which will fix the Yleft side moments)
  int np_current = 0, nplast = nop_OS-1;
   while (np_current < nplast+1){
    if ( (OS_x[np_current] < xstart && ptVCT->getCoordinates(0) != 0 ) || (OS_x[np_current] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1)) ){
      // communicate if they don't belong to the domain
      if (OS_x[np_current] < xstart && ptVCT->getCoordinates(0) != 0){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitXleft+1)>=max_np_OsPartComm){ 
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateOS doubling the sending buffer size" << endl;
          if (!resize_buffersOS((int) (max_np_OsPartComm*2)))
            {
              cout<<"communicateOS: increase MAX_NP_OSPARTCOMM\nExiting...";
              return -1;
            }
	}
	// put it in the communication buffer
	bufferXleftOS(OSParticles_Comm_LEFT,np_current,ptVCT);
	// delete the particle and pack the particle array, the value of nplast changes
	del_packOS(np_current,&nplast);
	npExitXleft++;
	//	del_packOS(np_current,&nplast); // REMOVE
      } else if (OS_x[np_current] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1)){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitXright+1)>=max_np_OsPartComm){
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateOS doubling the sending buffer size" << endl;
          if (!resize_buffersOS((int) (max_np_OsPartComm*2)))
            {
              cout<<"communicateOS: increase MAX_NP_OSPARTCOMM\nExiting...";
              return -1;
            }
	}
	// put it in the communication buffer
	bufferXrightOS(OSParticles_Comm_RIGHT,np_current,ptVCT); 
	// delete the particle and pack the particle array, the value of nplast changes
	del_packOS(np_current,&nplast);
	npExitXright++;
	//del_packOS(np_current,&nplast);// REMOVE 
      }
    } 
    else if ( (OS_y[np_current] < ystart && ptVCT->getCoordinates(1) != 0 ) || (OS_y[np_current] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1)) ){
      // communicate if they don't belong to the domain
      if (OS_y[np_current] < ystart && ptVCT->getCoordinates(1) != 0){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitYleft+1)>= max_np_OsPartComm){
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateOS doubling the sending buffer size" << endl;
          if (!resize_buffersOS((int) (max_np_OsPartComm*2)))
            {
              cout<<"communicateOS: increase MAX_NP_OSPARTCOMM\nExiting...";
              return -1;
            }
	}
	// put it in the communication buffer
	bufferYleftOS(OSParticles_Comm_BOTTOM,np_current,ptVCT); 
	// delete the particle and pack the particle array, the value of nplast changes
	del_packOS(np_current,&nplast);
	npExitYleft++;
	//del_packOS(np_current,&nplast);// REMOVE 
      } else if (OS_y[np_current] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1)){
	// check if there is enough space in the buffer before putting in the particle
	if( (npExitYright+1)>=max_np_OsPartComm){
	  cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() << "communicateOS doubling the sending buffer size" << endl;
          if (!resize_buffersOS((int) (max_np_OsPartComm*2)))
            {
              cout<<"communicateOS: increase MAX_NP_OSPARTCOMM\nExiting...";
              return -1;
            }
	}
	// put it in the communication buffer
	bufferYrightOS(OSParticles_Comm_TOP ,np_current,ptVCT); 
	// delete the particle and pack the particle array, the value of nplast changes
	del_packOS(np_current,&nplast);
	npExitYright++;
	  //del_packOS(np_current,&nplast);// REMOVE 
      }
    }  
    else {
      // particle is still in the domain, procede with the next particle
      np_current++;
    }
    
  }// end while 
  
  nop_OS = nplast+1;
  if (nop_OS > (npmax_OS - (int) (.01*npmax_OS) ) ){
    cout << "communicateOS exceeding npmax_OS: Particles need to be resized Save Data and Stop the simulation" << endl;
    return(-1); // end the simulation because you dont have enough space on the array                                    \
                                                                                                                          
  }

  npExitingMax = 0;
  // calculate the maximum number of particles exiting from this domain
  // use this value to check if communication is needed
  // and to  resize the buffer
  npExitingMax = maxNpExiting();
  // broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
  npExitingMax = reduceMaxNpExiting(ptVCT->getCART_COMM(),npExitingMax);
  
  /*****************************************************/
  /*           SEND AND RECEIVE MESSAGES               */
  /*****************************************************/
  
  //  new_OSpart = npExitingMax + 1;
  new_OSpart = max_np_OsPartComm;
  while (npExitingMax +1 > new_OSpart){
    new_OSpart=(int) (new_OSpart*2);
  }
  if (new_OSpart!= max_np_OsPartComm)
    {
      cout << "R" << ptVCT->getCartesian_rank_COMMTOTAL() <<"ns " <<ns << " communicateOS resizing the receiving buffer" << endl; 
    if (!resize_buffersOS((int) (new_OSpart)))
      {                                                         
        cout<<"communicateOS: increase MAX_NP_OSPARTCOMM\nExiting...";    
        return -1;                                                  
      } 
  }

  /*if (ptVCT->getCartesian_rank_COMMTOTAL()== ptVCT->getXLEN()*ptVCT->getYLEN())
    {
      cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << " npExitingMax in communicateOS" << npExitingMax <<endl;
      }*/
  
  if (npExitingMax > 0){
    // the first parameter is the buffer size, not the number of particles
    //communicateParticles(new_OSpart*nVarOS,OSParticles_Comm_LEFT,OSParticles_Comm_RIGHT,OSParticles_Comm_BOTTOM, OSParticles_Comm_TOP,ptVCT);
    communicateParticles( (npExitingMax + 1)*nVarOS,OSParticles_Comm_LEFT,OSParticles_Comm_RIGHT,OSParticles_Comm_BOTTOM, OSParticles_Comm_TOP,ptVCT);     


    // UNBUFFERING
    // message from XLEFT
    //cout << "R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "unbuffer in communicateOS" <<endl;
    avail1 = unbufferOS(OSParticles_Comm_LEFT,ptVCT);
    // message from XRIGHT
    avail2 = unbufferOS(OSParticles_Comm_RIGHT,ptVCT);
    // message from YLEFT
    avail3 = unbufferOS(OSParticles_Comm_BOTTOM,ptVCT);
    // message from YRIGHT
    avail4 = unbufferOS(OSParticles_Comm_TOP,ptVCT);
    // if one of these numbers is negative than there is not enough space for particles
    avail = avail1 + avail2 + avail3 + avail4;
    availALL = reduceNumberParticles(ptVCT->getCART_COMM(),avail);
    if (availALL < 0)
      return(-1);  // too many particles coming, save data nad stop simulation
    
  }
  
  return(0); // everything was fine
}

/** put a OS particle exiting to X-LEFT in the bufferXLEFT for communication*/
void Particles2Dcomm::bufferXleftOS(double *b_, int np_current, VirtualTopology* vct){

  b_[npExitXleft*nVarOS]    = OS_x[np_current];
  b_[npExitXleft*nVarOS +1] = OS_y[np_current];
  b_[npExitXleft*nVarOS +2] = OS_u[np_current];
  b_[npExitXleft*nVarOS +3] = OS_v[np_current];
  b_[npExitXleft*nVarOS +4] = OS_w[np_current];
  b_[npExitXleft*nVarOS +5] = OS_q[np_current];

  if (TrackParticleID)
    b_[npExitXleft*nVarOS +6] = OS_ParticleID[np_current];
  if (cVERBOSE)
    cout << "Particle exiting to Xleft: X=" << OS_x[np_current] << " ("<< xstart<<"," << xend << ")"<< endl;
  
}
/** put a OS particle exiting to X-RIGHT in the bufferXRIGHT for communication*/
void Particles2Dcomm::bufferXrightOS(double *b_, int np_current, VirtualTopology* vct){

  b_[npExitXright*nVarOS]    = OS_x[np_current];
  b_[npExitXright*nVarOS +1] = OS_y[np_current];
  b_[npExitXright*nVarOS +2] = OS_u[np_current];
  b_[npExitXright*nVarOS +3] = OS_v[np_current];
  b_[npExitXright*nVarOS +4] = OS_w[np_current];
  b_[npExitXright*nVarOS +5] = OS_q[np_current];

  if (TrackParticleID)
    b_[npExitXright*nVarOS +6] = OS_ParticleID[np_current];
  if (cVERBOSE)
    cout << "Particle exiting to Xright: X=" << OS_x[np_current] << " ("<< xstart<<"," << xend << ")"<< endl;

}

/** put a OS particle exiting to Y-LEFT in the bufferYLEFT for communication*/
void Particles2Dcomm::bufferYleftOS(double *b_, int np_current, VirtualTopology* vct){

  b_[npExitYleft*nVarOS]    = OS_x[np_current];
  b_[npExitYleft*nVarOS +1] = OS_y[np_current];
  b_[npExitYleft*nVarOS +2] = OS_u[np_current];
  b_[npExitYleft*nVarOS +3] = OS_v[np_current];
  b_[npExitYleft*nVarOS +4] = OS_w[np_current];
  b_[npExitYleft*nVarOS +5] = OS_q[np_current];

  if (TrackParticleID)
    b_[npExitYleft*nVarOS +6] = OS_ParticleID[np_current];
  if (cVERBOSE)
    cout << "Particle exiting to Yleft: Y=" << OS_y[np_current] << " ("<< xstart<<"," << xend << ")"<< endl;

}

/** put a OS particle exiting to Y-RIGHT in the bufferYRIGHT for communication*/
void Particles2Dcomm::bufferYrightOS(double *b_, int np_current, VirtualTopology* vct){

  b_[npExitYright*nVarOS]    = OS_x[np_current];
  b_[npExitYright*nVarOS +1] = OS_y[np_current];
  b_[npExitYright*nVarOS +2] = OS_u[np_current];
  b_[npExitYright*nVarOS +3] = OS_v[np_current];
  b_[npExitYright*nVarOS +4] = OS_w[np_current];
  b_[npExitYright*nVarOS +5] = OS_q[np_current];

  if (TrackParticleID)
    b_[npExitYright*nVarOS +6] = OS_ParticleID[np_current];
  if (cVERBOSE)
    cout << "Particle exiting to Yright: Y=" << OS_y[np_current] << " ("<< xstart<<"," << xend << ")"<< endl;

}

void Particles2Dcomm::del_packOS(int np_current, int *nplast){
  OS_x[np_current] = OS_x[*nplast];
  OS_y[np_current] = OS_y[*nplast];
  OS_u[np_current] = OS_u[*nplast];
  OS_v[np_current] = OS_v[*nplast];
  OS_w[np_current] = OS_w[*nplast];
  OS_q[np_current] = OS_q[*nplast];
  
  if (TrackParticleID)
    OS_ParticleID[np_current]=OS_ParticleID[*nplast];

  npExit++;
  (*nplast)--;
}

int Particles2Dcomm::unbufferOS(double *b_, VirtualTopology *ptVCT){
  int np_current =0;
  // put the new particles at the end of the array, and update the number of particles
  while(b_[np_current*nVarOS] != MIN_VAL){
    OS_x[nop_OS] = b_[nVarOS*np_current];
    OS_y[nop_OS] = b_[nVarOS*np_current+1];
    OS_u[nop_OS] = b_[nVarOS*np_current+2];
    OS_v[nop_OS] = b_[nVarOS*np_current+3];
    OS_w[nop_OS] = b_[nVarOS*np_current+4];
    OS_q[nop_OS] = b_[nVarOS*np_current+5];
    if (TrackParticleID)
      OS_ParticleID[nop_OS]=(unsigned long) b_[nVarOS*np_current+6];
    np_current++;
    
    if (cVERBOSE)
      cout << "Receiving Particle: X=" << OS_x[nop_OS] << ",Y=" << OS_y[nop_OS] << " ("<< xstart<<"," << xend << ")"<< " x ("<< ystart<<"," << yend << ")"<<endl;


    bool XnotRightDom= (x[nop] < xstart && ptVCT->getCoordinates(0) != 0 ) || (x[nop] > xend && ptVCT->getCoordinates(0) != (ptVCT->getXLEN()-1) );
    bool YnotRightDom=(y[nop] < ystart && ptVCT->getCoordinates(1) != 0 ) || (y[nop] > yend && ptVCT->getCoordinates(1) != (ptVCT->getYLEN()-1));
    if(XnotRightDom || YnotRightDom)
      {rightDomain++;} // the particle is not in the domain
    
    nop_OS++;
    if (nop_OS > (npmax_OS - (int) (.01*npmax_OS) ) ){
      cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "UnbufferOS exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
      return(-1); // end the simulation because you dont have enough space on the array
    }
  }// end while
  return(0); // everything was fine
  
}

int Particles2Dcomm::interpP2G_OS(Field* EMf, Grid *grid, VirtualTopology* vct)
{

  if (! (vct->getRefLevelAdj()==1 && grid->getLevel()>0))
    {
      //cout<< "interpP2G_OS not needed, returning" <<endl;
      return 1;
    }
  /*else
    {
      cout<< "interpP2G_OS" <<endl;
      }*/
  /*// debug ops
  if (nop_OS>0 && !(vct->getYleft_neighbor()==MPI_PROC_NULL || vct->getYright_neighbor()==MPI_PROC_NULL || vct->getXleft_neighbor()==MPI_PROC_NULL || vct->getXright_neighbor()==MPI_PROC_NULL))
    {
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<< "OS in wrong processor, nop_OS_B: " <<"\n";
      return -1;
      }*/
  // Restoring the corner native values if OS operations are to be undertaken  
  // otherwise, on those points the communicateNode is done twice

  
  REstoreCornerOsValues(EMf, vct);

  /*//debug ops
  for (int i=0; i< nop_OS; i++)
    {
      if ((OS_x[i]< 0 && vct->getXleft_neighbor()!= MPI_PROC_NULL) || (OS_x[i]>Lx && vct->getXright_neighbor()!= MPI_PROC_NULL) || (OS_y[i]< 0 && vct->getYleft_neighbor()!= MPI_PROC_NULL) || (OS_y[i]>Ly && vct->getYright_neighbor()!= MPI_PROC_NULL))
	{
	  cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<< "OS in wrong processor" << endl;
	  cout << "OS x " << OS_x[i] <<" OS y " << OS_y[i] << " - dx "  <<0-grid->getDX() << " Lx + dx " << Lx + grid->getDX() << " - dy "  <<0-grid->getDY() << " Ly + dy " << Ly + grid->getDY() << " xend " << xend << " ystart " << ystart <<endl;
	  return -1;
	}
	}*/

  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << " InterpP2G OS" << endl;
  double*** weight;// = newArr3(double,2,2,1);
  allocArr3(&weight, 2, 2, 1);
  double*** temp;// = newArr3(double,2,2,1);
  allocArr3(&temp, 2, 2, 1);
  int ix,iy, temp2,temp1;
  double inv_dx, inv_dy;
  inv_dx = 1.0/dx;
  inv_dy = 1.0/dy;
  
  /*if (0 && vct->getCartesian_rank_COMMTOTAL() == 16 && ns ==0 )
    {
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << "rhons before interpP2G os \n";
      EMf->printRhons(ns, vct);
      }*/
  
  for (register int i=0; i < nop_OS; i++){
    ix = 2 +  int(floor((OS_x[i]-xstart)*inv_dx));
    iy = 2 +  int(floor((OS_y[i]-ystart)*inv_dy));
    
    weight[1][1][0] = ((OS_x[i] - grid->getModifiedXN(ix-1,iy-1,0))*inv_dx)*((OS_y[i] - grid->getModifiedYN(ix-1,iy-1,0))*inv_dy);  
    weight[1][0][0] = ((OS_x[i] - grid->getModifiedXN(ix-1,iy,0))*inv_dx)*((grid->getModifiedYN(ix-1,iy,0) - OS_y[i])*inv_dy);
    weight[0][1][0] = ((grid->getModifiedXN(ix,iy-1,0) - OS_x[i])*inv_dx)*((OS_y[i] - grid->getModifiedYN(ix,iy-1,0))*inv_dy);
    weight[0][0][0] = ((grid->getModifiedXN(ix,iy,0) - OS_x[i])*inv_dx)*((grid->getModifiedYN(ix,iy,0) - OS_y[i])*inv_dy);
    
    //cout << "ix " << ix << " iy " << iy << " x_OS[i] " << OS_x[i]  << " y_OS[i] " << OS_y[i] << " grid->getModifiedXN(ix,iy,0) " << grid->getModifiedXN(ix,iy,0) << " grid->getModifiedYN(ix,iy,0) " << grid->getModifiedYN(ix,iy,0) << " weight[0][0][0] " << weight[0][0][0] << " weight[0][1][0] " << weight[0][1][0]<< " weight[1][0][0] " << weight[1][0][0] << " weight[1][1][0] " << weight[1][1][0] << "-2 dx  " << -2* grid->getDX() << "- dx  " << grid->getDX()<<"Lx + dx" << Lx + grid->getDX() << "Lx + 2 dx" << Lx + 2*grid->getDX() <<endl;
    
    scale(weight,OS_q[i],2,2);
    // add charge density                                                                                                                   
    EMf->addRho_OS(weight,ix,iy,0,ns, grid);
    // add current density - X                                                                                                              
    eqValue(0.0,temp,2,2);
    addscale(OS_u[i],temp,weight,2,2);
    EMf->addJx_OS(temp,ix,iy,0,ns, grid);
    // add current density - Y                                                                                                              
    eqValue(0.0,temp,2,2);
    addscale(OS_v[i],temp,weight,2,2);
    EMf->addJy_OS(temp,ix,iy,0,ns, grid);
    // add current density - Z                                                                                                              
    eqValue(0.0,temp,2,2);
    addscale(OS_w[i],temp,weight,2,2);
    EMf->addJz_OS(temp,ix,iy,0,ns, grid);
    //Pxx - add pressure tensor                                                                                                             
    eqValue(0.0,temp,2,2);
    addscale(OS_u[i]*OS_u[i],temp,weight,2,2);
    EMf->addPxx_OS(temp,ix,iy,0,ns, grid);
    // Pxy - add pressure tensor                                                                                                            
    eqValue(0.0,temp,2,2);
    addscale(OS_u[i]*OS_v[i],temp,weight,2,2);
    EMf->addPxy_OS(temp,ix,iy,0,ns, grid);
    // Pxz - add pressure tensor                                                                                                            
    eqValue(0.0,temp,2,2);
    addscale(OS_u[i]*OS_w[i],temp,weight,2,2);
    EMf->addPxz_OS(temp,ix,iy,0,ns, grid);
    // Pyy - add pressure tensor                                                                                                            
    eqValue(0.0,temp,2,2);
    addscale(OS_v[i]*OS_v[i],temp,weight,2,2);
    EMf->addPyy_OS(temp,ix,iy,0,ns, grid);
    // Pyz - add pressure tensor                                                                                                            
    eqValue(0.0,temp,2,2);
    addscale(OS_v[i]*OS_w[i],temp,weight,2,2);
    EMf->addPyz_OS(temp,ix,iy,0,ns, grid);
    // Pzz - add pressure tensor                                                                                                            
    eqValue(0.0,temp,2,2);
    addscale(OS_w[i]*OS_w[i],temp,weight,2,2);
    EMf->addPzz_OS(temp,ix,iy,0,ns, grid);
    
  }
  
  // communicate contribution from ghost cells     
  // commented for debug
  EMf->communicateGhostP2GOS(ns,0,0,0,0,vct); // just to fix the ghost corners
  //delArr3(weight,2,2);
  freeArr3(&weight);
  //delArr3(temp,2,2);
  freeArr3(&temp);
  
  /*//debug ops
    if (0 && vct->getCartesian_rank_COMMTOTAL() == 16 && ns==0 )
    {
    cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << "rhons after interpP2G os \n";
    EMf->printRhons(ns, vct);
    }*/

  return 1;

}


void Particles2Dcomm::storeCornerOsValues(Field * EMf, VirtualTopology * vct)
{
  // upper left
  if(vct->getXleft_neighbor()==MPI_PROC_NULL)
    {
      OSfix_upperLeft_rhons = EMf->getRHOns(0,nyn-2,0,ns);

      OSfix_upperLeft_Jx = EMf->getJxs(0,nyn-2,0,ns);
      OSfix_upperLeft_Jy = EMf->getJys(0,nyn-2,0,ns);
      OSfix_upperLeft_Jz = EMf->getJxs(0,nyn-2,0,ns);

      OSfix_upperLeft_Pxx = EMf->getpXXns(0,nyn-2,0,ns);
      OSfix_upperLeft_Pxy = EMf->getpXYns(0,nyn-2,0,ns);
      OSfix_upperLeft_Pxz = EMf->getpXZns(0,nyn-2,0,ns);
      OSfix_upperLeft_Pyy = EMf->getpYYns(0,nyn-2,0,ns);
      OSfix_upperLeft_Pyz = EMf->getpYZns(0,nyn-2,0,ns);
      OSfix_upperLeft_Pzz = EMf->getpZZns(0,nyn-2,0,ns);
    }
  if(vct->getYright_neighbor()==MPI_PROC_NULL)
    {
      OSfix_upperLeft_rhons = EMf->getRHOns(1,nyn-1,0,ns);

      OSfix_upperLeft_Jx = EMf->getJxs(1,nyn-1,0,ns);
      OSfix_upperLeft_Jy = EMf->getJys(1,nyn-1,0,ns);
      OSfix_upperLeft_Jz = EMf->getJzs(1,nyn-1,0,ns);

      OSfix_upperLeft_Pxx = EMf->getpXXns(1,nyn-1,0,ns);
      OSfix_upperLeft_Pxy = EMf->getpXYns(1,nyn-1,0,ns);
      OSfix_upperLeft_Pxz = EMf->getpXZns(1,nyn-1,0,ns);
      OSfix_upperLeft_Pyy = EMf->getpYYns(1,nyn-1,0,ns);
      OSfix_upperLeft_Pyz = EMf->getpYZns(1,nyn-1,0,ns);
      OSfix_upperLeft_Pzz = EMf->getpZZns(1,nyn-1,0,ns);
    }

  // upper right
  if (vct->getXright_neighbor()==MPI_PROC_NULL)
    {
      OSfix_upperRight_rhons = EMf->getRHOns(nxn-1,nyn-2,0,ns);

      OSfix_upperRight_Jx = EMf->getJxs(nxn-1,nyn-2,0,ns);
      OSfix_upperRight_Jy = EMf->getJys(nxn-1,nyn-2,0,ns);
      OSfix_upperRight_Jz = EMf->getJzs(nxn-1,nyn-2,0,ns);

      OSfix_upperRight_Pxx = EMf->getpXXns(nxn-1,nyn-2,0,ns);
      OSfix_upperRight_Pxy = EMf->getpXYns(nxn-1,nyn-2,0,ns);
      OSfix_upperRight_Pxz = EMf->getpXZns(nxn-1,nyn-2,0,ns);
      OSfix_upperRight_Pyy = EMf->getpYYns(nxn-1,nyn-2,0,ns);
      OSfix_upperRight_Pyz = EMf->getpYZns(nxn-1,nyn-2,0,ns);
      OSfix_upperRight_Pzz = EMf->getpZZns(nxn-1,nyn-2,0,ns);
    }
  if (vct->getYright_neighbor()==MPI_PROC_NULL)
    {
      OSfix_upperRight_rhons = EMf->getRHOns(nxn-2,nyn-1,0,ns);

      OSfix_upperRight_Jx = EMf->getJxs(nxn-2,nyn-1,0,ns);
      OSfix_upperRight_Jy = EMf->getJys(nxn-2,nyn-1,0,ns);
      OSfix_upperRight_Jz = EMf->getJzs(nxn-2,nyn-1,0,ns);

      OSfix_upperRight_Pxx = EMf->getpXXns(nxn-2,nyn-1,0,ns);
      OSfix_upperRight_Pxy = EMf->getpXYns(nxn-2,nyn-1,0,ns);
      OSfix_upperRight_Pxz = EMf->getpXZns(nxn-2,nyn-1,0,ns);
      OSfix_upperRight_Pyy = EMf->getpYYns(nxn-2,nyn-1,0,ns);
      OSfix_upperRight_Pyz = EMf->getpYZns(nxn-2,nyn-1,0,ns);
      OSfix_upperRight_Pzz = EMf->getpZZns(nxn-2,nyn-1,0,ns);
    }

  // lower left              
  if(vct->getXleft_neighbor()==MPI_PROC_NULL)
    {   
      OSfix_lowerLeft_rhons = EMf->getRHOns(0,1,0,ns);

      OSfix_lowerLeft_Jx = EMf->getJxs(0,1,0,ns);
      OSfix_lowerLeft_Jy = EMf->getJys(0,1,0,ns);
      OSfix_lowerLeft_Jz = EMf->getJzs(0,1,0,ns);

      OSfix_lowerLeft_Pxx = EMf->getpXXns(0,1,0,ns);
      OSfix_lowerLeft_Pxy = EMf->getpXYns(0,1,0,ns);
      OSfix_lowerLeft_Pxz = EMf->getpXZns(0,1,0,ns);
      OSfix_lowerLeft_Pyy = EMf->getpYYns(0,1,0,ns);
      OSfix_lowerLeft_Pyz = EMf->getpYZns(0,1,0,ns);
      OSfix_lowerLeft_Pzz = EMf->getpZZns(0,1,0,ns);
    }
  if(vct->getYleft_neighbor()==MPI_PROC_NULL)
    {
      OSfix_lowerLeft_rhons = EMf->getRHOns(1,0,0,ns);

      OSfix_lowerLeft_Jx = EMf->getJxs(1,0,0,ns);
      OSfix_lowerLeft_Jy = EMf->getJys(1,0,0,ns);
      OSfix_lowerLeft_Jz = EMf->getJzs(1,0,0,ns);

      OSfix_lowerLeft_Pxx = EMf->getpXXns(1,0,0,ns);
      OSfix_lowerLeft_Pxy = EMf->getpXYns(1,0,0,ns);
      OSfix_lowerLeft_Pxz = EMf->getpXZns(1,0,0,ns);
      OSfix_lowerLeft_Pyy = EMf->getpYYns(1,0,0,ns);
      OSfix_lowerLeft_Pyz = EMf->getpYZns(1,0,0,ns);
      OSfix_lowerLeft_Pzz = EMf->getpZZns(1,0,0,ns);
    }

  // lower right
  if(vct->getXright_neighbor()==MPI_PROC_NULL)
    {
      OSfix_lowerRight_rhons = EMf->getRHOns(nxn-1,1,0,ns);

      OSfix_lowerRight_Jx = EMf->getJxs(nxn-1,1,0,ns);
      OSfix_lowerRight_Jy = EMf->getJys(nxn-1,1,0,ns);
      OSfix_lowerRight_Jz = EMf->getJzs(nxn-1,1,0,ns);

      OSfix_lowerRight_Pxx = EMf->getpXXns(nxn-1,1,0,ns);
      OSfix_lowerRight_Pxy = EMf->getpXYns(nxn-1,1,0,ns);
      OSfix_lowerRight_Pxz = EMf->getpXZns(nxn-1,1,0,ns);
      OSfix_lowerRight_Pyy = EMf->getpYYns(nxn-1,1,0,ns);
      OSfix_lowerRight_Pyz = EMf->getpYZns(nxn-1,1,0,ns);
      OSfix_lowerRight_Pzz = EMf->getpZZns(nxn-1,1,0,ns);

    }
  if(vct->getYleft_neighbor()==MPI_PROC_NULL)
    {
      OSfix_lowerRight_rhons = EMf->getRHOns(nxn-2,0,0,ns);

      OSfix_lowerRight_Jx = EMf->getJxs(nxn-2,0,0,ns);
      OSfix_lowerRight_Jy = EMf->getJys(nxn-2,0,0,ns);
      OSfix_lowerRight_Jz = EMf->getJzs(nxn-2,0,0,ns);

      OSfix_lowerRight_Pxx = EMf->getpXXns(nxn-2,0,0,ns);
      OSfix_lowerRight_Pxy = EMf->getpXYns(nxn-2,0,0,ns);
      OSfix_lowerRight_Pxz = EMf->getpXZns(nxn-2,0,0,ns);
      OSfix_lowerRight_Pyy = EMf->getpYYns(nxn-2,0,0,ns);
      OSfix_lowerRight_Pyz = EMf->getpYZns(nxn-2,0,0,ns);
      OSfix_lowerRight_Pzz = EMf->getpZZns(nxn-2,0,0,ns);
    }
  
}

void Particles2Dcomm::REstoreCornerOsValues(Field * EMf, VirtualTopology * vct)
{
  // upper left
  if(vct->getXleft_neighbor()==MPI_PROC_NULL)
    {
      EMf->setRHOns(OSfix_upperLeft_rhons,0,nyn-2,0,ns);

      EMf->setJxs(OSfix_upperLeft_Jx,0,nyn-2,0,ns);
      EMf->setJys(OSfix_upperLeft_Jy,0,nyn-2,0,ns);
      EMf->setJxs(OSfix_upperLeft_Jz,0,nyn-2,0,ns);

      EMf->setpXXns(OSfix_upperLeft_Pxx,0,nyn-2,0,ns);
      EMf->setpXYns(OSfix_upperLeft_Pxy,0,nyn-2,0,ns);
      EMf->setpXZns(OSfix_upperLeft_Pxz,0,nyn-2,0,ns);
      EMf->setpYYns(OSfix_upperLeft_Pyy,0,nyn-2,0,ns);
      EMf->setpYZns(OSfix_upperLeft_Pyz,0,nyn-2,0,ns);
      EMf->setpZZns(OSfix_upperLeft_Pzz,0,nyn-2,0,ns);
    }
  if(vct->getYright_neighbor()==MPI_PROC_NULL)
    {
      EMf->setRHOns(OSfix_upperLeft_rhons,1,nyn-1,0,ns);

      EMf->setJxs(OSfix_upperLeft_Jx,1,nyn-1,0,ns);
      EMf->setJys(OSfix_upperLeft_Jy,1,nyn-1,0,ns);
      EMf->setJzs(OSfix_upperLeft_Jz,1,nyn-1,0,ns);

      EMf->setpXXns(OSfix_upperLeft_Pxx,1,nyn-1,0,ns);
      EMf->setpXYns(OSfix_upperLeft_Pxy,1,nyn-1,0,ns);
      EMf->setpXZns(OSfix_upperLeft_Pxz,1,nyn-1,0,ns);
      EMf->setpYYns(OSfix_upperLeft_Pyy,1,nyn-1,0,ns);
      EMf->setpYZns(OSfix_upperLeft_Pyz,1,nyn-1,0,ns);
      EMf->setpZZns(OSfix_upperLeft_Pzz,1,nyn-1,0,ns);
    }

  // upper right
  if (vct->getXright_neighbor()==MPI_PROC_NULL)
    {
      EMf->setRHOns(OSfix_upperRight_rhons,nxn-1,nyn-2,0,ns);

      EMf->setJxs(OSfix_upperRight_Jx,nxn-1,nyn-2,0,ns);
      EMf->setJys(OSfix_upperRight_Jy,nxn-1,nyn-2,0,ns);
      EMf->setJzs(OSfix_upperRight_Jz,nxn-1,nyn-2,0,ns);

      EMf->setpXXns(OSfix_upperRight_Pxx,nxn-1,nyn-2,0,ns);
      EMf->setpXYns(OSfix_upperRight_Pxy,nxn-1,nyn-2,0,ns);
      EMf->setpXZns(OSfix_upperRight_Pxz,nxn-1,nyn-2,0,ns);
      EMf->setpYYns(OSfix_upperRight_Pyy,nxn-1,nyn-2,0,ns);
      EMf->setpYZns(OSfix_upperRight_Pyz,nxn-1,nyn-2,0,ns);
      EMf->setpZZns(OSfix_upperRight_Pzz,nxn-1,nyn-2,0,ns);
    }
  if (vct->getYright_neighbor()==MPI_PROC_NULL)
    {
      EMf->setRHOns(OSfix_upperRight_rhons,nxn-2,nyn-1,0,ns);

      EMf->setJxs(OSfix_upperRight_Jx,nxn-2,nyn-1,0,ns);
      EMf->setJys(OSfix_upperRight_Jy,nxn-2,nyn-1,0,ns);
      EMf->setJzs(OSfix_upperRight_Jz,nxn-2,nyn-1,0,ns);

      EMf->setpXXns(OSfix_upperRight_Pxx,nxn-2,nyn-1,0,ns);
      EMf->setpXYns(OSfix_upperRight_Pxy,nxn-2,nyn-1,0,ns);
      EMf->setpXZns(OSfix_upperRight_Pxz,nxn-2,nyn-1,0,ns);
      EMf->setpYYns(OSfix_upperRight_Pyy,nxn-2,nyn-1,0,ns);
      EMf->setpYZns(OSfix_upperRight_Pyz,nxn-2,nyn-1,0,ns);
      EMf->setpZZns(OSfix_upperRight_Pzz,nxn-2,nyn-1,0,ns);
    }

  // lower left              
  if(vct->getXleft_neighbor()==MPI_PROC_NULL)
    {   
      EMf->setRHOns(OSfix_lowerLeft_rhons,0,1,0,ns);

      EMf->setJxs(OSfix_lowerLeft_Jx,0,1,0,ns);
      EMf->setJys(OSfix_lowerLeft_Jy,0,1,0,ns);
      EMf->setJzs(OSfix_lowerLeft_Jz,0,1,0,ns);

      EMf->setpXXns(OSfix_lowerLeft_Pxx,0,1,0,ns);
      EMf->setpXYns(OSfix_lowerLeft_Pxy,0,1,0,ns);
      EMf->setpXZns(OSfix_lowerLeft_Pxz,0,1,0,ns);
      EMf->setpYYns(OSfix_lowerLeft_Pyy,0,1,0,ns);
      EMf->setpYZns(OSfix_lowerLeft_Pyz,0,1,0,ns);
      EMf->setpZZns(OSfix_lowerLeft_Pzz,0,1,0,ns);
    }
  if(vct->getYleft_neighbor()==MPI_PROC_NULL)
    {
      EMf->setRHOns(OSfix_lowerLeft_rhons,1,0,0,ns);

      EMf->setJxs(OSfix_lowerLeft_Jx,1,0,0,ns);
      EMf->setJys(OSfix_lowerLeft_Jy,1,0,0,ns);
      EMf->setJzs(OSfix_lowerLeft_Jz,1,0,0,ns);

      EMf->setpXXns(OSfix_lowerLeft_Pxx,1,0,0,ns);
      EMf->setpXYns(OSfix_lowerLeft_Pxy,1,0,0,ns);
      EMf->setpXZns(OSfix_lowerLeft_Pxz,1,0,0,ns);
      EMf->setpYYns(OSfix_lowerLeft_Pyy,1,0,0,ns);
      EMf->setpYZns(OSfix_lowerLeft_Pyz,1,0,0,ns);
      EMf->setpZZns(OSfix_lowerLeft_Pzz,1,0,0,ns);
    }

  // lower right
  if(vct->getXright_neighbor()==MPI_PROC_NULL)
    {
      EMf->setRHOns(OSfix_lowerRight_rhons,nxn-1,1,0,ns);

      EMf->setJxs(OSfix_lowerRight_Jx,nxn-1,1,0,ns);
      EMf->setJys(OSfix_lowerRight_Jy,nxn-1,1,0,ns);
      EMf->setJzs(OSfix_lowerRight_Jz,nxn-1,1,0,ns);

      EMf->setpXXns(OSfix_lowerRight_Pxx,nxn-1,1,0,ns);
      EMf->setpXYns(OSfix_lowerRight_Pxy,nxn-1,1,0,ns);
      EMf->setpXZns(OSfix_lowerRight_Pxz,nxn-1,1,0,ns);
      EMf->setpYYns(OSfix_lowerRight_Pyy,nxn-1,1,0,ns);
      EMf->setpYZns(OSfix_lowerRight_Pyz,nxn-1,1,0,ns);
      EMf->setpZZns(OSfix_lowerRight_Pzz,nxn-1,1,0,ns);

    }
  if(vct->getYleft_neighbor()==MPI_PROC_NULL)
    {
      EMf->setRHOns(OSfix_lowerRight_rhons,nxn-2,0,0,ns);

      EMf->setJxs(OSfix_lowerRight_Jx,nxn-2,0,0,ns);
      EMf->setJys(OSfix_lowerRight_Jy,nxn-2,0,0,ns);
      EMf->setJzs(OSfix_lowerRight_Jz,nxn-2,0,0,ns);

      EMf->setpXXns(OSfix_lowerRight_Pxx,nxn-2,0,0,ns);
      EMf->setpXYns(OSfix_lowerRight_Pxy,nxn-2,0,0,ns);
      EMf->setpXZns(OSfix_lowerRight_Pxz,nxn-2,0,0,ns);
      EMf->setpYYns(OSfix_lowerRight_Pyy,nxn-2,0,0,ns);
      EMf->setpYZns(OSfix_lowerRight_Pyz,nxn-2,0,0,ns);
      EMf->setpZZns(OSfix_lowerRight_Pzz,nxn-2,0,0,ns);
    }
}

bool Particles2Dcomm::resize_buffersSP(int new_particleNumber){
  cout << "RESIZING SP BUFFERS FROM " <<  max_np_SplitPartComm << "PARTICLES TO " << new_particleNumber << " PARTICLES"<< endl;
  
  // the resize is possible up to MAX_NP_SPLIPARTCOMM particles;
  // if needed, change this number in allocate
  if (new_particleNumber> MAX_NP_SPLIPARTCOMM)
    return false;

  double *temp = new double[max_np_SplitPartComm*nVar];

  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    temp[i] = SplittedParticles_Comm_BOTTOM_ptr[i];                                                                                       
  delete[] SplittedParticles_Comm_BOTTOM;
  SplittedParticles_Comm_BOTTOM = new double[new_particleNumber*nVar];
  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    SplittedParticles_Comm_BOTTOM[i] = temp[i];
  for(int i= max_np_SplitPartComm*nVar; i < new_particleNumber*nVar; i++)
    SplittedParticles_Comm_BOTTOM[i] = MIN_VAL;

  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    temp[i] = SplittedParticles_Comm_TOP_ptr[i];
  delete[] SplittedParticles_Comm_TOP;
  SplittedParticles_Comm_TOP = new double[new_particleNumber*nVar];
  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    SplittedParticles_Comm_TOP[i] = temp[i];
  for(int i= max_np_SplitPartComm*nVar; i < new_particleNumber*nVar; i++)
    SplittedParticles_Comm_TOP[i] = MIN_VAL;

  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    temp[i] = SplittedParticles_Comm_LEFT_ptr[i];
  delete[] SplittedParticles_Comm_LEFT;
  SplittedParticles_Comm_LEFT = new double[new_particleNumber*nVar];
  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    SplittedParticles_Comm_LEFT[i] = temp[i];
  for(int i= max_np_SplitPartComm*nVar; i < new_particleNumber*nVar; i++)
    SplittedParticles_Comm_LEFT[i] = MIN_VAL;

  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    temp[i] = SplittedParticles_Comm_RIGHT_ptr[i];
  delete[] SplittedParticles_Comm_RIGHT;
  SplittedParticles_Comm_RIGHT = new double[new_particleNumber*nVar];
  for(int i=0; i < max_np_SplitPartComm*nVar; i++)
    SplittedParticles_Comm_RIGHT[i] = temp[i];
  for(int i= max_np_SplitPartComm*nVar; i < new_particleNumber*nVar; i++)
    SplittedParticles_Comm_RIGHT[i] = MIN_VAL;

  SplittedParticles_Comm_BOTTOM_ptr = SplittedParticles_Comm_BOTTOM;
  SplittedParticles_Comm_TOP_ptr = SplittedParticles_Comm_TOP;
  SplittedParticles_Comm_LEFT_ptr = SplittedParticles_Comm_LEFT;
  SplittedParticles_Comm_RIGHT_ptr = SplittedParticles_Comm_RIGHT;

  max_np_SplitPartComm = new_particleNumber;

  cout << "END RESIZING SP BUFFERS FROM " <<  max_np_SplitPartComm << "PARTICLES TO " << new_particleNumber << " PARTICLES"<< endl;
  return true; //everything OK
}

bool Particles2Dcomm::resize_buffersOS(int new_particleNumber){
  cout << "RESIZING OS BUFFERS FROM " <<  max_np_OsPartComm << "PARTICLES TO " << new_particleNumber << " PARTICLES"<< endl;

  // the resize is possible up to MAX_NP_OSPARTCOMM particles;           
  // if needed, change this number in allocate                                                                                               
  if (new_particleNumber> MAX_NP_OSPARTCOMM)
    return false;

  double *temp = new double[max_np_OsPartComm*nVarOS];

  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    temp[i] = OSParticles_Comm_BOTTOM_ptr[i];
  delete[] OSParticles_Comm_BOTTOM;
  OSParticles_Comm_BOTTOM = new double[new_particleNumber*nVarOS];
  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    OSParticles_Comm_BOTTOM[i] = temp[i];
  for(int i= max_np_OsPartComm*nVarOS; i < new_particleNumber*nVarOS; i++)
    OSParticles_Comm_BOTTOM[i] = MIN_VAL;

  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    temp[i] = OSParticles_Comm_TOP_ptr[i];
  delete[] OSParticles_Comm_TOP;
  OSParticles_Comm_TOP = new double[new_particleNumber*nVarOS];
  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    OSParticles_Comm_TOP[i] = temp[i];
  for(int i= max_np_OsPartComm*nVarOS; i < new_particleNumber*nVarOS; i++)
    OSParticles_Comm_TOP[i] = MIN_VAL;

  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    temp[i] = OSParticles_Comm_LEFT_ptr[i];
  delete[] OSParticles_Comm_LEFT;
  OSParticles_Comm_LEFT = new double[new_particleNumber*nVarOS];
  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    OSParticles_Comm_LEFT[i] = temp[i];
  for(int i= max_np_OsPartComm*nVarOS; i < new_particleNumber*nVarOS; i++)
    OSParticles_Comm_LEFT[i] = MIN_VAL;

  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    temp[i] = OSParticles_Comm_RIGHT_ptr[i];
  delete[] OSParticles_Comm_RIGHT;
  OSParticles_Comm_RIGHT = new double[new_particleNumber*nVarOS];
  for(int i=0; i < max_np_OsPartComm*nVarOS; i++)
    OSParticles_Comm_RIGHT[i] = temp[i];
  for(int i= max_np_OsPartComm*nVarOS; i < new_particleNumber*nVarOS; i++)
    OSParticles_Comm_RIGHT[i] = MIN_VAL;

  OSParticles_Comm_BOTTOM_ptr = OSParticles_Comm_BOTTOM;
  OSParticles_Comm_TOP_ptr = OSParticles_Comm_TOP;
  OSParticles_Comm_LEFT_ptr = OSParticles_Comm_LEFT;
  OSParticles_Comm_RIGHT_ptr = OSParticles_Comm_RIGHT;

  max_np_OsPartComm = new_particleNumber;

  return true; //everything OK                                                                                                                
}

// used for internal ops in communicate; the first buffer_size elements of the vector to MINVAL                                           
void Particles2Dcomm::setToMINVAL_comm(double *vec)
{
  memcpy(vec, MIN_VAL_VEC_COMM, buffer_size*sizeof(double));
  // try if faster
  /*    for (int i=0; i<buffer_size; i++)
	vec[i]=MIN_VAL;*/
  return;
}
double Particles2Dcomm::roundPrec(double x, int prec)
{
  double power = 1.0;
  int i;

  if (prec > 0)
    for (i = 0; i < prec; i++)
      power *= 10.0;
  else if (prec < 0)
    for (i = 0; i < prec; i++)
      power /= 10.0;

  if (x > 0)
    x = floor(x * power + 0.5) / power;
  else if (x < 0)
    x = ceil(x * power - 0.5) / power;

  if (x == -0)
    x = 0;

  return x;
}
int Particles2Dcomm::RandomizePositionPRAParticles(int species, VirtualTopology* vct,  Grid* grid)
{
  // only for refined grids                                                                                                                                  
  if (grid->getLevel()==0) return 1;

  //only for boundary cores, since PRA cannot extent further than them                                                         
  if (! (vct->getXright_neighbor()==MPI_PROC_NULL or vct->getXleft_neighbor()==MPI_PROC_NULL or vct->getYright_neighbor()==MPI_PROC_NULL or vct->getYleft_neighbor()==MPI_PROC_NULL))
    return 1;

  srand((unsigned)time(NULL)+ species + vct->getCartesian_rank()+1); // try not to have the same sequence                                                        
  double  dx = grid->getDX(),dy = grid->getDY();
  int ix, iy;

  int ExtraCells=0;  //0 to Randomize only on repopulated particles


  if (nop> npmax-1)
    {
      cout << "nop problem, Exiting\n";
      return -1;
    }
  for (int i=0; i<nop; i++)
    {
      // X left boundary                                     
      if (vct->getXleft_neighbor()==MPI_PROC_NULL)
        {

          if (x[i]<= grid->getXN(PRA_Xleft, 1, 0) + ExtraCells*dx)
            {
              // fMin + (double)(rand()/(double)RAND_MAX * (fMax- fMin)          
              // x: random in the PRA, x dir                        
	      double range= dx*(PRA_Xleft + ExtraCells);
              x[i]=-dx + (double)(rand()/(double)RAND_MAX) * range;
	      while (! (x[i]> -dx and x[i] <Lx +dx))
		{
		  cout << "Randomize: x[i] " <<x[i] << " -dx: " << -dx << " Lx+dx: " << Lx+dx <<endl;
		  x[i]=-dx + (double)(rand()/(double)RAND_MAX) * range;
		  cout << "Randomize: this was a potential seg fault\n";
		}
	      //cout<< "lower: " << "limit " << grid->getXN(PRA_Xleft, 1, 1) + ExtraCells*dx  << " -dx+range " << -dx+range <<endl;
              // y: random in the cell                     
              iy = 1 +  int(floor((y[i]-ystart)/dy));
	      // diagn
	      if (iy <0 or iy> nyn-1)
		{
		  cout << "X left B, iy out of range; iy " <<iy <<", y[i]  " <<y[i] <<" Ly+ dy " <<Ly+dy<<endl; 
		  return -1;
		}
              y[i]= grid->getYN(1, iy, 0) + (double)(rand()/(double)RAND_MAX)* dy;
            }
        }
      // X right boundary   
      if (vct->getXright_neighbor()==MPI_PROC_NULL)
        {
          if (x[i]>= grid->getXN(nxn-1-PRA_Xright, 1, 0) - ExtraCells*dx)
            {
              // fMin + (double)(rand()/(double)RAND_MAX * (fMax- fMin)      
              // x: random in the PRA, x dir                                                  
	      double range= dx*(PRA_Xright + ExtraCells);
              x[i]=grid->getXN(nxn-1-PRA_Xright, 1, 0) - ExtraCells*dx  + (double)(rand()/(double)RAND_MAX) * range;
	      while (! (x[i]> -dx and x[i] <Lx +dx))
		{
		  cout << "Randomize: x[i] " <<x[i] << " -dx: "<< -dx << " Lx+dx: " <<Lx+dx <<endl;
		  x[i]=grid->getXN(nxn-1-PRA_Xright, 1, 0) - ExtraCells*dx  + (double)(rand()/(double)RAND_MAX) * range;
		  cout << "Randomize: this was a potential seg fault\n";
		}
              // y: random in the cell                            
	      //cout << "upper: " << " limit " << grid->getXN(nxn-1-PRA_Xright, 1, 1) - ExtraCells*dx << " limit + range " << grid->getXN(nxn-1-PRA_Xright, 1, 1) - ExtraCells*dx + range << " Lx +dx " << Lx+dx <<endl;
              iy = 1 +  int(floor((y[i]-ystart)/dy));
	      if (iy <0 or iy> nyn-1)
                {
                  cout << "X right B, iy out of range; iy " <<iy <<", y[i]  " <<y[i] <<" Ly+ dy " <<Ly+dy<<endl;
		  return -1;
                }
              y[i]= grid->getYN(1, iy, 0) + (double)(rand()/(double)RAND_MAX)* dy;
            }
        }
      // Y left boundary    
      if (vct->getYleft_neighbor()==MPI_PROC_NULL)
        {
          if (y[i]<= grid->getYN(1,PRA_Yleft, 0) + ExtraCells*dy)
            {
              // fMin + (double)(rand()/(double)RAND_MAX * (fMax- fMin)      
              // y: random in the PRA, y dir          
	      double range= dy*(PRA_Yleft + ExtraCells);
              y[i]=-dy + (double)(rand()/(double)RAND_MAX) * range;
	      while (! (y[i]> -dy and y[i] <Ly +dy))
		{
		  cout << "Randomize: y[i] " <<y[i] << " -dy: "<< -dy << " Ly+dy: " <<Ly+dy <<endl;
		  y[i]=-dy + (double)(rand()/(double)RAND_MAX) * range;
		  cout << "Randomize: this was a potential seg fault\n";
		}
              // x: random in the cell                                 
              ix = 1 +  int(floor((x[i]-xstart)/dx));
	      if (ix <0 or ix> nxn-1)
                {
                  cout << "Y left B, iy out of range; ix " <<ix <<", x[i]  " <<x[i] <<" Lx+ dx " <<Lx+dx<<endl;
		  return -1;
                }
              x[i]= grid->getXN(ix, 1, 0) + (double)(rand()/(double)RAND_MAX)* dx;
            }
        }
      // Y right boundary   
      if (vct->getYright_neighbor()==MPI_PROC_NULL)
        {
          if (y[i]>= grid->getYN(1, nyn-1-PRA_Yright, 0) - ExtraCells*dy)
            {
              // fMin + (double)(rand()/(double)RAND_MAX * (fMax- fMin)   
              // y: random in the PRA, y dir               
              double range = dy*(PRA_Yright + ExtraCells);
              y[i]=grid->getYN(1, nyn-1-PRA_Yright, 0)- ExtraCells*dy + (double)(rand()/(double)RAND_MAX) * range;
	      while (! (y[i]> -dy and y[i] <Ly +dy))
		{
		  cout << "Randomize: y[i] " <<y[i] << " -dy: "<< -dy << " Ly+dy: " <<Ly+dy <<endl;
		  y[i]=grid->getYN(1, nyn-1-PRA_Yright, 0)- ExtraCells*dy + (double)(rand()/(double)RAND_MAX) * range;
		  cout << "Randomize: this was a potential seg fault\n";
		}
              // x: random in the cell                                                
              ix = 1 +  int(floor((x[i]-xstart)/dx));
	      if (ix <0 or ix> nxn-1)
                {
                  cout << "Y right B, iy out of range; ix " <<ix <<", x[i]  " <<x[i] <<" Lx+ dx " <<Lx+dx<<endl;
		  return -1;
                }
              x[i]= grid->getXN(ix, 1, 0) + (double)(rand()/(double)RAND_MAX)* dx;
            }
        }

    }// end nop                                                                                                       
  return 1;
}

int Particles2Dcomm::SubCyclingParticles(Grid* grid, bool CoarseOp, VirtualTopology* ptVCT)
{
  if (SubCycling == false or grid->getLevel()==0 or CoarseOp== true)
    {
      cout << "You asked to execute SubCyclingParticles in suspicious conditions... Re-check what you are doing...\n";
      return -1;
    }

  memcpy (x+nop, RP_x, sizeof(double)*RP_nop);
  memcpy (y+nop, RP_y, sizeof(double)*RP_nop);
  memcpy (u+nop, RP_u, sizeof(double)*RP_nop);  
  memcpy (v+nop, RP_v, sizeof(double)*RP_nop);
  memcpy (w+nop, RP_w, sizeof(double)*RP_nop);
  memcpy (q+nop, RP_q, sizeof(double)*RP_nop);
  if (TrackParticleID)
    memcpy (ParticleID+nop, RP_ParticleID, sizeof(unsigned long)*RP_nop);

  nop=nop+RP_nop;

  if (nop > (npmax - (int) (.01*npmax) ) )
    {
    cout <<"R" <<ptVCT->getCartesian_rank_COMMTOTAL() << "SubCyclingParticles exceeding npmax: Particles \
need to be resized Save Data and Stop the simulation" << endl;
    return(-1); // end the simulation because you dont have enough space on the array          
  }


  //cout << "R" <<ptVCT->getCartesian_rank_COMMTOTAL() << " After subcycling particles, nop " << nop << " npmax " <<npmax <<endl;
  return 1;
}



/**AMR methods, ME*/
/**init operations connected with the repopulation of particles, new version with communiation map*/
int Particles2Dcomm::initPRAVariables(int species, CollectiveIO* col,VirtualTopology* vct, Grid* grid, Field* EMf){
  // the buffer for communication of PRA variables are defined in allocate, 
  // together with the other particle communication buffers

  if (vct->getNgrids()<2)
    return 1;

  ratio= col->getRatio();

  //number of cells, ghost cell INCLUDED, for particle repopulation; x left
  PRA_Xleft = col->GetPRA_Xleft(); 
  //number of cells, ghost cell INCLUDED, for particle repopulation; x right
  PRA_Xright = col->GetPRA_Xright();
  //number of cells, ghost cell INCLUDED, for particle repopulation; y left
  PRA_Yleft = col->GetPRA_Yleft();
  //number of cells, ghost cell INCLUDED, for particle repopulation; y right
  PRA_Yright = col->GetPRA_Yright();
  
  if (grid->getLevel() >0){  //PRA area, native particles falling here are deleted and substituted with repopulated particles 
    PRA_oxStartLeft   = 0- grid->getDX();
    PRA_oxEndLeft     = 0+ (PRA_Xleft-1)*grid->getDX();
    PRA_oxStartRight  = Lx- (PRA_Xright-1)*grid->getDX();
    PRA_oxEndRight    = Lx+ grid->getDX();
    
    PRA_oyStartLeft   = 0- grid->getDY();
    PRA_oyEndLeft     = 0+ (PRA_Yleft-1)*grid->getDY();
    PRA_oyStartRight  = Ly- (PRA_Yright-1)*grid->getDY();
    PRA_oyEndRight    = Ly+ grid->getDY();
  }
  else{   //used only for the safety checks, NOT in the code; initialized anyhow
    PRA_oxStartLeft   = 0- grid->getDX();
    PRA_oxEndLeft     = 0;
    PRA_oxStartRight  = Lx;
    PRA_oxEndRight    = Lx+ grid->getDX();

    PRA_oyStartLeft   = 0- grid->getDY();
    PRA_oyEndLeft     = 0;
    PRA_oyStartRight  = Ly;
    PRA_oyEndRight    = Ly+ grid->getDY();
  }
      
  //if Level< Levels-1 (this grid is the coarser grid for some other grid), limits for thr PRA area of the child in local coords
  if ( grid->getLevel() < vct->getNgrids()-1){

    double Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid
    double Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid

    double finedx = grid->getDX()/col->getRatio();
    double finelx = col->getLx()/(double)pow(col->getRatio(),grid->getLevel()+1);
    double finedy = grid->getDY()/col->getRatio();
    double finely = col->getLy()/(double)pow(col->getRatio(),grid->getLevel()+1);
    
    // when parent particles enter this area, they have to be communicated;  
    // the parent dx or dy is already taken into account here

    PRA_CoxStartLeft   = Ox - finedx - grid->getDX() ;
    PRA_CoxEndLeft     = Ox + (PRA_Xleft-1)*finedx + grid->getDX();
    PRA_CoxStartRight  = Ox + finelx - (PRA_Xright-1)*finedx -  grid->getDX();
    PRA_CoxEndRight    = Ox + finelx + finedx + grid->getDX();

    PRA_CoyStartLeft   = Oy - finedy - grid->getDY();
    PRA_CoyEndLeft     = Oy + (PRA_Yleft-1)*finedy + grid->getDY();
    PRA_CoyStartRight  = Oy + finely - (PRA_Yright-1)*finedy - grid->getDY();
    PRA_CoyEndRight    = Oy + finely + finedy + grid->getDY();

    //if (Modified_xstart > PRA_CoxEndLeft or Modified_xend < PRA_CoxStartLeft or Modified_ystart > PRA_CoyEndLeft or Modified_yend < PRA_CoyStartLeft)
    if (Modified_xstart > PRA_CoxEndRight or Modified_xend < PRA_CoxStartLeft or Modified_ystart > PRA_CoyEndRight or Modified_yend < PRA_CoyStartLeft)
      {
	PRAIntersection= false;
      }	
    else
      {
	PRAIntersection= true;
      }

  }
  else{ //actually not used, initialized anyhow
    PRA_CoxStartLeft   = 0- grid->getDX();;
    PRA_CoxEndLeft     = 0- grid->getDX();;
    PRA_CoxStartRight  = Lx+ grid->getDX();
    PRA_CoxEndRight    = Lx+ grid->getDX();

    PRA_CoyStartLeft   = 0- grid->getDY();;
    PRA_CoyEndLeft     = 0- grid->getDY();;
    PRA_CoyStartRight  = Ly+ grid->getDY();
    PRA_CoyEndRight    = Ly+ grid->getDY();


    PRAIntersection= false;
  } 

  // for PRASend, modified from initWeightBC
  int i,j,nproc;
  double finedx, finedy,finelx,finely, xfirst, xlast,xfirstnext, Ox, Oy;
  double coarsedx, coarsedy,coarselx,coarsely;
  double finelxplusfinedx;
  double xshift, yshift;
  int xnnl, xnnu, ynnl, ynnu;
  // remember that there is always only 1 PRA ghost cell, the other in the active domain

  int expSR;
  if (col->getXLEN()> col->getYLEN())
    {
      expSR= col->getXLEN()*2*4; // a occhio                                                                                                  
    }
  else
    {
      expSR= col->getYLEN()*2*4;
    }

  targetBC= new int[expSR];
  BCSide  = new int[expSR];
  
  targetBOTTOM=0;
  targetTOP=0;
  targetLEFT=0;
  targetRIGHT=0;
  nmessageBC = 0;
  
  nmessagerecuBC=0;
  // for debugging purposes, not actually used anywhere
  int nmessagerecuBCLEFT=0;
  int nmessagerecuBCRIGHT=0;
  int nmessagerecuBCBOTTOM=0;
  int nmessagerecuBCTOP=0;
  // end for debugging purposes
  
  if (grid->getLevel() > 0){
    fromBC= new int[expSR];//[col->getXLEN()*col->getYLEN()];  
    BCSidecu= new int[expSR];//[col->getXLEN()*col->getYLEN()];
    nmessagerecuBC=0;
    
    //If this grid is considered as fine by another grid
    coarsedx = grid->getDX()*col->getRatio();
    coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel()-1);
    coarsedy = grid->getDY()*col->getRatio();
    coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel()-1);
    Ox = grid->getOx(grid->getLevel()); //Origin x of the grid
    Oy = grid->getOy(grid->getLevel()); //Origin y of the grid
    j=0;
    if(vct->getCoordinates(1) == 0) {    // BOTTOM
      double xloc, yloc1, yloc2;
      if (vct->getCoordinates(0)== 0)
	{
	  xfirst=PRA_oxStartLeft- coarsedx ; // coarsedx from def, to have the areas coincide in coarse and fine grid
	}
      else
	{
	  xfirst= grid->getXstart();
	}
      if(vct->getCoordinates(0) == vct->getXLEN()-1) {
	xlast = PRA_oxEndRight + coarsedx;// coarsedx from def, to have the areas coincide in coarse and fine grid
      }else{
	xlast = grid->getXend();
	if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc
	  {
	    xlast= xlast - grid->getDX();
	  }
      }
      
      yloc1 = Oy + PRA_oyStartLeft- coarsedy;// coarsedy from def, to have the areas coincide in coarse and fine grid
      yloc2 = Oy + PRA_oyEndLeft + coarsedy; // coarsedy from def, to have the areas coincide in coarse and fine grid
      double YL= yloc1;
      int nproc;
      while (YL < yloc2 || fabs(YL - yloc2)< DBL_EPSILON ){
	xnnl = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
	for (i=0; i< xnnl; i++) {
	  xloc = max(Ox + xfirst +i*grid->getDX(),0.);// Because when x < 0, it is on the same proc as if x=0  
	  nproc =floor(YL/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1);
	  bool found=false;
	  for (int k=0; k< nmessagerecuBC; k++ )
	    {
	      if (BCSidecu[k]==0 && fromBC[k]==nproc)
		{ 
		  found= true;
		  break;
		}
	    }
	  if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	    fromBC[0] = nproc;
	    BCSidecu[0]=0;
	    nmessagerecuBC++;
	    nmessagerecuBCBOTTOM++;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new bottom from "<<fromBC[nmessagerecuBC-1] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	  }
	  
	  if(nproc != fromBC[j] && !found){ // first part to avoid repetitions with points at the same yloc, second with previous yloc
	    j++;
	    nmessagerecuBC++;
	    nmessagerecuBCBOTTOM++;
	    fromBC[j]=nproc;
	    BCSidecu[j]=0;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new bottom from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC << " j " << j<<endl;
	  }
	}// end xnnl
	YL += grid->getDY();
      }// end YL
    } // end BOTTOM
    
    if(vct->getCoordinates(1) == vct->getYLEN()-1) { // TOP
      double xloc, yloc1, yloc2;
      if (vct->getCoordinates(0)== 0)
	{
	  xfirst=PRA_oxStartLeft - coarsedx; // coarsedx from def, to have the areas coincide in coarse and fine grid
	}
      else
	{
	  xfirst= grid->getXstart();
	}
      if(vct->getCoordinates(0) == vct->getXLEN()-1) {
	xlast = PRA_oxEndRight + coarsedx;   // coarsedx from def, to have the areas coincide in coarse and fine grid
      }else{
	xlast = grid->getXend();
	if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc               
	  {
	    xlast= xlast - grid->getDX();
	  }
      }
      yloc1 = Oy +   PRA_oyStartRight- coarsedy;// coarsedy from def, to have the areas coincide in coarse and fine grid
      yloc2 = Oy +   PRA_oyEndRight +coarsedy; // coarsedy from def, to have the areas coincide in coarse and fine grid
      double YL= yloc1;
      int nproc;
      while (YL < yloc2 || fabs(YL - yloc2)< DBL_EPSILON ){
	xnnu = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
	for (i=0; i< xnnu; i++) {
	  xloc = max(Ox + xfirst +i*grid->getDX(),0.);
	  nproc =floor(YL/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	  bool found=false;
	  for (int k=0; k< nmessagerecuBC; k++ )
	    {
	      if (BCSidecu[k]==1 && fromBC[k]==nproc)
		{ 
		  found= true;
		  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<< " BCSidecu[k] " <<BCSidecu[k] << " fromBC[k] " << fromBC[k] << " nproc " <<nproc <<endl; 
		  break;
		}
	    }
	  if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	    fromBC[0] = nproc;
	    BCSidecu[0]=1;
	    nmessagerecuBC++;
	    nmessagerecuBCTOP++;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new bottom from "<<fromBC[nmessagerecuBC-1] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	  }
	  // (nproc == fromBC[j] && BCSidecu[j]!= 1) othewise procs are not included if they are the last in a different side
	  if((nproc != fromBC[j] || (nproc == fromBC[j] && BCSidecu[j]!= 1)) && !found){
	    j++;
	    nmessagerecuBC++;
	    nmessagerecuBCTOP++;
	    fromBC[j]=nproc;
	    BCSidecu[j]=1;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new top from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	  }
	}// end xnnu
	YL += grid->getDY();
      }// end YL
    } // end TOP
    
    if(vct->getCoordinates(0) == 0) {   // LEFT
      double yloc, xloc1, xloc2;      
      if (vct->getCoordinates(1)== 0)
	{
	  xfirst=PRA_oyEndLeft + coarsedy;// this actually spans the y dir // coarsedy from def, to have the areas coincide in coarse and fine grid   
	}
      else
	{
	  xfirst= grid->getYstart();
	}
      if(vct->getCoordinates(1) == vct->getYLEN()-1) {
	xlast = PRA_oyStartRight - coarsedy; // coarsedy from def, to have the areas coincide in coarse and fine grid
      }else{
	xlast = grid->getYend();
	if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc              
	  {
	    xlast= xlast - grid->getDY();
	  }
      }
      xloc1 = Ox + PRA_oxStartLeft -coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid
      xloc2 = Ox + PRA_oxEndLeft+ coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid
      double XL= xloc1;
      int nproc;
      while (XL < xloc2 || fabs(XL - xloc2)< DBL_EPSILON) {
	ynnl = floor((xlast-xfirst)/grid->getDY()+0.5)+1; 
	for (i=0; i< ynnl; i++) {
	  yloc = Oy + xfirst +i*grid->getDY();
	  nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(XL/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	  bool found=false;
	  for (int k=0; k< nmessagerecuBC; k++ )
	    {
	      if (BCSidecu[k]==2 && fromBC[k]==nproc)
		{ 
		  found= true;
		  break;
		}
	    }
	  if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	    fromBC[0] = nproc;
	    BCSidecu[0]=2;
	    nmessagerecuBC++;
	    nmessagerecuBCLEFT++;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new left from "<<fromBC[nmessagerecuBC-1] << " nmessagerecuBC " << nmessagerecuBC<<endl;
	  }
	  if((nproc != fromBC[j] || (nproc == fromBC[j] && BCSidecu[j]!= 2)) && !found){
	    //if(nproc != fromBC[j] && !found){ // first part to avoid repetitions with points at the same yloc, second with previous yloc
	    j++;
	    nmessagerecuBC++;
	    nmessagerecuBCLEFT++;
	    fromBC[j]=nproc;
	    BCSidecu[j]=2;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new left from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC << " j " << j<<endl;
	  }
	} // end ynnl
	XL += grid->getDX();
      }// end XL
    } // end left
    
    if(vct->getCoordinates(0) == vct->getXLEN()-1) {  //RIGHT
      double yloc, xloc1, xloc2;
      if (vct->getCoordinates(1)== 0)
	{
	  xfirst=PRA_oyEndLeft + coarsedy;// this actually spans the y dir // coarsedy from def, to have the areas coincide in coarse and fine grid
	}
      else
	{
	  xfirst= grid->getYstart();
	}
      if(vct->getCoordinates(1) == vct->getYLEN()-1) {
	xlast = PRA_oyStartRight- coarsedy;// coarsedy from def, to have the areas coincide in coarse and fine grid 
      }else{
	xlast = grid->getYend();
	if ( fabs(xlast-xfirst) > DBL_EPSILON )//If the fine subdivision overlap coarse subdivision, to avoid catching the next coarse proc
	  {
	    xlast= xlast - grid->getDY();
	  }
      }
      xloc1 = Ox + PRA_oxStartRight - coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid  
      xloc2 = Ox + PRA_oxEndRight + coarsedx;// coarsedy from def, to have the areas coincide in coarse and fine grid
      double XL= xloc1;
      int nproc;
      ynnu = floor((xlast-xfirst)/grid->getDY()+0.5)+1;
      while (XL < xloc2 || fabs(XL - xloc2)< DBL_EPSILON ){
	for (i=0; i< ynnu; i++) {
	  yloc = Oy + xfirst +i*grid->getDY(); 
	  nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(XL/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << " right nproc1 " << nproc1 << " nproc2 "<< nproc2 << " xloc1 " <<xloc1 <<" xloc2 "<<xloc2 <<" yloc "<< yloc <<endl;
	  bool found=false;
	  for (int k=0; k< nmessagerecuBC; k++ )
	    {
	      if (BCSidecu[k]==3 && fromBC[k]==nproc)
		{ 
		  found= true;
		  break;
		}
	    }
	  if (i==0 && nmessagerecuBC==0){// j must be nmessagerecuBC-1, so not updated the first time
	    fromBC[0] = nproc;
	    BCSidecu[0]=3;
	    nmessagerecuBC++;
	    nmessagerecuBCRIGHT++;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new right from "<<fromBC[0] << " nmessagerecuBC " << nmessagerecuBC << " nproc " << nproc<<endl;
	  }
	  if((nproc != fromBC[j] || (nproc == fromBC[j] && BCSidecu[j]!= 3)) && !found){
	    //if(nproc != fromBC[j] && !found){ // first part to avoid repetitions with points at the same yloc, second with previous yloc
	    j++;
	    nmessagerecuBC++;
	    nmessagerecuBCRIGHT++;
	    fromBC[j]=nproc;
	    BCSidecu[j]=3;
	    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<  " new right from "<<fromBC[j] << " nmessagerecuBC " << nmessagerecuBC << " j " << j<<endl;
	  }
	}// end ynnu
	XL += grid->getDX();
      }// end XL
    } // end RIGHT
    
    
  }// end check on grid level 
  
  
  // now, the refined grid distributes the info
  
  // broadcast part: refined grid broadcast the communication map to the coarse grid
  
  int TagMap= 666;
  MPI_Request requestISend;
  MPI_Status status;
  
  if (grid->getLevel() > 0 )
    {
      fromBC[nmessagerecuBC]= -1;// to signal the termination of the list                              
      for (int r= 0; r< vct->getXLEN()*vct->getYLEN(); r++)
	{
	  // +1 to send the -1 for termination of the list as well                                     
	  MPI_Isend (fromBC, nmessagerecuBC+1, MPI_INT, r, TagMap, vct->getCART_COMM_TOTAL(), &requestISend);
	  MPI_Wait(&requestISend, &status);
	}
    }
  
  int *buffer=new int[expSR];//[col->getXLEN()*col->getYLEN()*4];
  if (grid->getLevel() == 0 )
    {
      nmessageBC=0;
      for (int r=0; r< vct->getXLEN()*vct->getYLEN(); r++)
	{
	  MPI_Recv(buffer, expSR,MPI_INT, MPI_ANY_SOURCE, TagMap, vct->getCART_COMM_TOTAL(), &status);
	  for (int i=0; i< expSR; i++)
	    {
	      if (buffer[i]!=-1)
		{
		  if (buffer[i]== vct->getCartesian_rank_COMMTOTAL())
		    nmessageBC++;
		}
	      else
		break;
	    }
	}
    }
  
  
  // specific location part
  int NInfo= 2;
  int TagInfo=667; 
  int *info=new int[NInfo];
  if (grid->getLevel() > 0 )
    {
      for (int i=0; i< nmessagerecuBC; i++)
	{
	  info[0]= vct->getCartesian_rank_COMMTOTAL();
	  info[1]= BCSidecu[i];
	  
	  MPI_Isend (info, NInfo, MPI_DOUBLE, fromBC[i], TagInfo, vct->getCART_COMM_TOTAL(), &requestISend);
	  MPI_Wait(&requestISend, &status);
	}
    }
  if (grid->getLevel() == 0 )
    {
      targetTOP=0;
      targetBOTTOM=0;
      targetLEFT=0;
      targetRIGHT=0;
      for (int i=0; i<nmessageBC; i++)
	{
	  MPI_Recv(info, NInfo,MPI_DOUBLE, MPI_ANY_SOURCE, TagInfo, vct->getCART_COMM_TOTAL(), &status);
	  
	  targetBC[i]= info[0];
	  BCSide[i]= info[1];
	  switch(BCSide[i])
	    {
	    case 0:
	      {
		targetBOTTOM++;
		break;
	      }
	    case 1:
	      targetTOP++;
	      break;
	      
	    case 2:
	      targetLEFT++;
	      break;
	      
	    case 3:
	      targetRIGHT++;
	      break;
	      
	    default:
	      {
		cout << "Problems in initPRAVariables, exiting..." << endl;
		return -1;
	      }
	    }// end switch
	}// end for
    }// end coarse grid
  // end coarse grid receiving the info  
  
  //some safety checks: for a grid (remind the PRA limits for the coarser grid), check that its finer grid's PRA does not fall into its PRA
  //it may be a problem with the modifications to the particle mover
  // NB: this part has never been tester
   if (grid->getLevel() < vct->getNgrids()-1){
    bool GoodCondition = (PRA_CoxStartLeft > PRA_oxEndLeft) && (PRA_CoxEndRight< PRA_oxStartRight) && (PRA_CoyStartLeft > PRA_oyEndLeft) && (PRA_CoyEndRight< PRA_oyStartRight);
    if (! GoodCondition){
      cout << "The child grid Particle Repopulation Area overlaps the current grid Particle Repopulation Area: recheck your init parameters...\n Some diagnostics then exiting..." <<endl;
     cout << "Lx: " << Lx << ", 0-dx: " <<0-grid->getDX() << ", Lx+dx: " <<Lx + grid->getDX() <<endl;
     cout << "PRA_CoxStartLeft: " <<PRA_CoxStartLeft <<", PRA_oxEndLeft: " <<PRA_oxEndLeft <<endl;
     if (! (PRA_CoxStartLeft > PRA_oxEndLeft))
       cout <<"Must be: PRA_CoxStartLeft > PRA_oxEndLeft!!!!" <<endl;
     cout << "PRA_CoxEndRight: " <<PRA_CoxEndRight <<", PRA_oxStartRight: " <<PRA_oxStartRight <<endl;
     if (!(PRA_CoxEndRight< PRA_oxStartRight))
       cout <<"Must be: PRA_CoxEndRight< PRA_oxStartRight!!!!" <<endl;
     cout << "Ly: " << Ly << ", 0-dy: " <<0-grid->getDY() << ", Ly+dy: " <<Ly + grid->getDY() <<endl;
     cout << "PRA_CoyStartLeft: " <<PRA_CoyStartLeft <<", PRA_oyEndLeft: " <<PRA_oyEndLeft <<endl;
     if(! (PRA_CoyStartLeft > PRA_oyEndLeft))
       cout <<"Must be: PRA_CoyStartLeft > PRA_oyEndLeft!!!!" <<endl;
     cout << "PRA_CoyEndRight: " <<PRA_CoyEndRight <<", PRA_oyStartRight: " <<PRA_oyStartRight <<endl;
     if (!(PRA_CoyEndRight< PRA_oyStartRight))
       cout<<"Must be: PRA_CoyEndRight< PRA_oyStartRight!!!!" <<endl;
     return -1;
    }
    // other safety check: the PRA should not fall into the parent grid's ghost area
    if (PRA_CoxStartLeft<0 || PRA_CoxEndLeft<0 || PRA_CoxStartRight> Lx || PRA_CoxEndRight> Lx || PRA_CoyStartLeft<0 || PRA_CoyEndLeft<0 || PRA_CoyStartRight> Ly || PRA_CoyEndRight> Ly )
      {
	cout << "The child Particle Repopulation Area is falling into the current grid  ghost area: recheck yout init parameters...\nSome diagnostics then exiting..." << endl;
	cout <<"PRA_CoxStartLeft: " << PRA_CoxStartLeft <<", PRA_CoxEndLeft: " <<PRA_CoxEndLeft <<endl;
	if (PRA_CoxEndLeft<0 || PRA_CoxEndLeft<0)
	  cout <<"PRA_CoxEndLeft and PRA_CoxEndLeft must be >0!!!"<<endl;
	cout <<"PRA_CoxStartRight: " <<PRA_CoxStartRight <<", PRA_CoxEndRight: " <<PRA_CoxEndRight <<endl;
	if (PRA_CoxStartRight> Lx || PRA_CoxEndRight> Lx)
	  cout <<"PRA_CoxStartRight and PRA_CoxEndRight must be <Lx!!!"<<endl;
	cout <<"PRA_CoyStartLeft: " << PRA_CoyStartLeft <<", PRA_CoyEndLeft: " <<PRA_CoyEndLeft<<endl;
        if (PRA_CoyEndLeft<0 || PRA_CoyEndLeft<0)
          cout <<"PRA_CoyEndLeft and PRA_CoyEndLeft must be >0!!!"<<endl;
	cout <<"PRA_CoyStartRight: " <<PRA_CoyStartRight <<", PRA_CoyEndRight: " <<PRA_CoyEndRight <<endl;
	if (PRA_CoyStartRight> Ly || PRA_CoyEndRight> Ly)
          cout <<"PRA_CoyStartRight and PRA_CoyEndRight must be <Ly!!!"<<endl;
	 
	return -1;
      }

  }
//for how the particle motion routines are modified, I need at least 1 PRA cell for side; exit if not//
    if (PRA_Xleft<1 || PRA_Xright<1 || PRA_Yleft<1 || PRA_Xright<1){
    cout << "At least 1 PRA cell per side is needed: modify your input file...\n Exiting...";
    return -1;
  }

  //the PRA area cannot be wider than a fine grid processor
if (grid->getLevel() >0)
  {// nxc is the number of cells in the processor, not the one from input file
    if ( PRA_Xleft-1> grid->getNXC()/2.0 || PRA_Xright-1> grid->getNXC()/2.0 || PRA_Yleft-1> grid->getNYC()/2.0  || PRA_Yright-1> grid->getNYC()/2.0   )
      {
	cout << "The Particle Repopulation Area is too wide compared with the number of cells per processor: revise input parameters.\nSome diagnostics then exiting...";
	cout << "Nxc: " << grid->getNXC() << ", PRA_Xleft: "<< PRA_Xleft <<", PRA_Xright: " <<PRA_Xright <<endl;
	cout << "Nyc: "<< grid->getNYC() << ", PRA_Yleft: "<< PRA_Yleft <<", PRA_Yright: " <<PRA_Yright <<endl;
	return -1;
      }
  }

// check that the number of sends match the number of receives
int ALL_TARGETS;
int ALL_RECEIVERS;

if (grid->getLevel() ==0)
  {
    MPI_Allreduce ( &nmessageBC, &ALL_TARGETS, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM());  
  }
 else
   {
     MPI_Allreduce ( &nmessagerecuBC, &ALL_RECEIVERS, 1,MPI_INT, MPI_SUM, vct->getCART_COMM());
   }
if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())) && grid->getLevel() ==0)
  {
    cout << "Level 0: ALL_TARGETS= " <<  ALL_TARGETS <<endl;
  }
if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())) && !grid->getLevel() ==0)
  {
    cout <<"Level 0: ALL_RECEIVERS= " <<  ALL_RECEIVERS <<endl;
  }
//appropriate iniializations done
MPI_Allreduce( &nmessageBC, &ALL_TARGETS, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce ( &nmessagerecuBC, &ALL_RECEIVERS, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

int ALLTARGETS_LEFT, ALLTARGETS_RIGHT, ALLTARGETS_BOTTOM, ALLTARGETS_TOP;
int ALLRECEIVERS_LEFT,  ALLRECEIVERS_RIGHT,  ALLRECEIVERS_BOTTOM,  ALLRECEIVERS_TOP;

MPI_Allreduce( &targetLEFT, &ALLTARGETS_LEFT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &targetRIGHT, &ALLTARGETS_RIGHT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &targetBOTTOM, &ALLTARGETS_BOTTOM, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &targetTOP, &ALLTARGETS_TOP, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

MPI_Allreduce( &nmessagerecuBCLEFT, &ALLRECEIVERS_LEFT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &nmessagerecuBCRIGHT, &ALLRECEIVERS_RIGHT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &nmessagerecuBCBOTTOM, &ALLRECEIVERS_BOTTOM, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
MPI_Allreduce( &nmessagerecuBCTOP, &ALLRECEIVERS_TOP, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

if (ALL_TARGETS!=ALL_RECEIVERS)
  {
    if (vct->getCartesian_rank_COMMTOTAL()==0)
      {
	cout <<"ALL_TARGETS= " <<ALL_TARGETS <<"!= ALL_RECEIVERS= "<<ALL_RECEIVERS <<endl;
	cout <<"ALLTARGETS_LEFT= " << ALLTARGETS_LEFT << " ALLRECEIVERS_LEFT= " << ALLRECEIVERS_LEFT <<endl;
	cout <<"ALLTARGETS_RIGHT= " << ALLTARGETS_RIGHT << " ALLRECEIVERS_RIGHT= "<< ALLRECEIVERS_RIGHT <<endl;
	cout <<"ALLTARGETS_BOTTOM= " << ALLTARGETS_BOTTOM << " ALLRECEIVERS_BOTTOM= "<< ALLRECEIVERS_BOTTOM <<endl;
	cout <<"ALLTARGETS_TOP= " << ALLTARGETS_TOP << " ALLRECEIVERS_TOP= "<< ALLRECEIVERS_TOP <<endl;
	return -1;
      }
  }
 delete[] buffer;
  delete[] info;
return 1; // OK
} // end new version of initPRAVariables

/*******************************************************************************************
Particles1Dcomm.cpp  -  Class for particles of the same species, in a 1D space
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#include <iostream>
#include <math.h>
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
#include "../fields/ESfield.h"
#include "Particles1Dcomm.h"

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
Particles1Dcomm::Particles1Dcomm(){
   // see allocate(int species, CollectiveIO* col, VirtualTopology* vct, Grid* grid)

}
/** deallocate particles */
Particles1Dcomm::~Particles1Dcomm(){
    delete[] x;
    delete[] u;
    delete[] q;
    // deallocate buffers
    delete[] b_XDX;
    delete[] b_XSN;
}
/** constructors fo a single species*/
void Particles1Dcomm::allocate(int species, CollectiveIO* col, VirtualTopology* vct, Grid* grid){
    // info from collectiveIO
    ns = species;
    npcel  = col->getNpcel(species);
    npcelx = col->getNpcelx(species);
	nop   =  col->getNp(species)/(vct->getNprocs());
    npmax =  col->getNpMax(species)/(vct->getNprocs());
	 
    qom   = col->getQOM(species);
    uth   = col->getUth(species);
    u0    = col->getU0(species);
	
    dt    = col->getDt();
    Lx     = col->getLx();
    TrackParticleID =col->getTrackParticleID(species);
	// info from Grid
	xstart = grid->getXstart();
    xend   = grid->getXend();
	nxn = grid->getNXN();
    invVOL=grid->getInvVOL();
    // info from VirtualTopology
    cVERBOSE = vct->getcVERBOSE();
	// boundary condition for particles
    bcPfaceXright = col->getBcPfaceXright();
    bcPfaceXleft = col->getBcPfaceXleft();
    bcPfaceYright = col->getBcPfaceYright();
    bcPfaceYleft = col->getBcPfaceYleft();
    
	// CONSTANTS
	c = col->getC();
	pi = 3.1415; 
	mu0 = 4*pi*1.e-7;
	eps0 = 1.0/(mu0*c*c);
	////////////////////////////////////////////////////////////////
    ////////////////     ALLOCATE ARRAYS   /////////////////////////
    ////////////////////////////////////////////////////////////////
    // positions
    x = new double[npmax];
    // velocities
    u = new double[npmax];
    // charge
    q = new double[npmax];
    //ID
    if (TrackParticleID){
      ParticleID= new unsigned long[npmax];
      BirthRank[0]=vct->getCartesian_rank();
      if (vct->getNprocs()>1) 
	    BirthRank[1]= (int) ceil(log10((double) (vct->getNprocs()))); // Number of digits needed for # of process in ID
      else BirthRank[1]=1;
      if (BirthRank[1]+ (int) ceil(log10((double) (npmax)))>10 && BirthRank[0] == 0 ) {
	    cerr<< "Error: can't Track particles in Particles1Dcomm::allocate"<<endl;
	    cerr<< "Unsigned long 'ParticleID' cannot store all the particles"<<endl;
	    return ;
      }
     }
    // BUFFERS
    // the buffer size should be decided depending on number of particles
    // the buffer size should be decided depending on number of particles
    if (TrackParticleID)
	   nVar=3;
    else 
	   nVar=4;
   buffer_size = (int) (.05*nop*nVar+1); // max: 5% of the particles in the processors is going out
   //buffer_size = 10;
   
   b_XDX = new double[buffer_size];
   b_XDX_ptr = b_XDX; // alias to make the resize
   b_XSN = new double[buffer_size];
   b_XSN_ptr = b_XSN; // alias to make the resize
    // if RESTART is true initialize the particle in allocate method
    restart = col->getRestart_status();
    if (restart!=0){
         if (vct->getCartesian_rank()==0 && ns==0)
	    cout << "LOADING PARTICLES FROM RESTART FILE in " + col->getRestartDirName() + "/restart.hdf" << endl;
	 stringstream ss;
	 ss << vct->getCartesian_rank();
	 string name_file = col->getRestartDirName() + "/restart" + ss.str() + ".hdf";
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
           cout << "couldn't open file: " << name_file << endl;
	   cout << "RESTART NOT POSSIBLE" << endl;
	 }
         
	 stringstream species_name;
	 species_name << ns;
	 // the cycle of the last restart is set to 0
	 string name_dataset = "/particles/species_" + species_name.str() + "/x/cycle_0";
	 dataset_id = H5Dopen(file_id,name_dataset.c_str());
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
	 // get u
	 name_dataset = "/particles/species_" + species_name.str() + "/u/cycle_0";
	 dataset_id = H5Dopen(file_id, name_dataset.c_str());
	 status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,u);
	 status = H5Dclose(dataset_id);
	 // get q
	 name_dataset = "/particles/species_" + species_name.str() + "/q/cycle_0";
	 dataset_id = H5Dopen(file_id, name_dataset.c_str());
	 status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,q);
	 status = H5Dclose(dataset_id);
	 // ID	
	 if (TrackParticleID){
             herr_t (*old_func)(void*);
             void *old_client_data;
             H5Eget_auto(&old_func, &old_client_data);
            /* Turn off error handling */
            H5Eset_auto(NULL, NULL);
	    name_dataset = "/particles/species_" + species_name.str() + "/ID/cycle_0";
	    dataset_id = H5Dopen(file_id, name_dataset.c_str());

            H5Eset_auto(old_func, old_client_data);
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
void Particles1Dcomm::interpP2G(Field* EMf, Grid *grid, VirtualTopology* vct){
     double*** weight = newArr3(double,2,2,1);
     int ix;
     double inv_dx;
     inv_dx = 1.0/(grid->getDX());
     for (register int i=0; i < nop; i++){
         ix = 1 +  int(floor((x[i]-grid->getXstart())*inv_dx)); // check the cell
         weight[1][0][0] = (x[i] - grid->getXN(ix,0,0))*inv_dx;
         weight[0][0][0] = (grid->getXN(ix+1,0,0) - x[i])*inv_dx;
		 scale(weight,q[i],2,2);
         // add charge density
         EMf->addRho(weight,ix,0,0,ns);
         
     }
     // Here you do the summation on the ghost     
     EMf->communicateGhostP2G(ns,0,0,0,0,vct);
	 delArr3(weight,2,2);
    
}

     
/** communicate buffers */
void Particles1Dcomm::communicate(VirtualTopology* ptVCT){
   // allocate buffers
   MPI_Status status;
   int new_buffer_size;
   int npExitingMax;
   
   for (int i=0; i < buffer_size; i++){
     b_XDX[i] = MIN_VAL;
     b_XSN[i] = MIN_VAL;
   }
   npExitXright =0, npExitXleft =0, npExit=0, rightDomain = 0;
   int np_current = 0, nplast = nop-1;
   while (np_current < nplast+1){
      // BC on particles
      if (x[np_current] < 0 && ptVCT->getXleft_neighbor() == MPI_PROC_NULL)
			 BCpart(&x[np_current],&u[np_current],Lx,uth,bcPfaceXright,bcPfaceXleft);
      else if (x[np_current] > Lx && ptVCT->getXright_neighbor() == MPI_PROC_NULL)
			 BCpart(&x[np_current],&u[np_current],Lx,uth,bcPfaceXright,bcPfaceXleft); 
	  // if the particle exits, apply the boundary conditions add the particle to communication buffer
      if (x[np_current] < xstart || x[np_current] >xend){
        	// communicate if they don't belong to the domain
		if (x[np_current] < xstart && ptVCT->getXleft_neighbor() != MPI_PROC_NULL){
			// check if there is enough space in the buffer before putting in the particle
                        if(((npExitXleft+1)*nVar)>=buffer_size){
                           cout << "resizing the sending buffer to " << (int) (buffer_size*2) << " buffer size" << endl;
                           resize_buffers((int) (buffer_size*2)); 
                          
                           
                        }
                        // put it in the communication buffer
                        bufferXleft(b_XSN,np_current,ptVCT);
            		//cout << "after buffering" << endl;
                        // delete the particle and pack the particle array, the value of nplast changes
            		del_pack(np_current,&nplast);
			npExitXleft++;
        	} else if (x[np_current] > xend && ptVCT->getXright_neighbor() != MPI_PROC_NULL){
            		// check if there is enough space in the buffer before putting in the particle
                        if(((npExitXright+1)*nVar)>=buffer_size){
                           cout << "resizing the sending buffer " << (int) (buffer_size*2) << endl; 
                           resize_buffers((int) (buffer_size*2)); 
                            
                          
                        }
                        // put it in the communication buffer
                        bufferXright(b_XDX,np_current,ptVCT);
            		// delete the particle and pack the particle array, the value of nplast changes
            		del_pack(np_current,&nplast);
                       
			npExitXright++;
        	} 
        } else {
              // particle is still in the domain, procede with the next particle
	      np_current++;
        }
      
   }
     
   nop = nplast+1;
   npExitingMax = 0;
   // calculate the maximum number of particles exiting from this domain
   // use this value to check if communication is needed
   // and to  resize the buffer
   npExitingMax = maxNpExiting();
   // broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
   npExitingMax = reduceMaxNpExiting(npExitingMax);

   /*****************************************************/
   /*           SEND AND RECEIVE MESSAGES               */
   /*****************************************************/
   
   new_buffer_size = npExitingMax*nVar + 1;
   
   if (new_buffer_size > buffer_size){
	 cout << "resizing the receiving buffer" << endl;
         resize_buffers(new_buffer_size);         
   } 
   if (npExitingMax > 0){
     communicateParticles(new_buffer_size,b_XSN,b_XDX,ptVCT);
     // UNBUFFERING
	 // message from XLEFT
	 unbuffer(b_XDX);
     // message from XRIGHT
	 unbuffer(b_XSN);
   }
   
 }
 /** resize the buffers */
void Particles1Dcomm::resize_buffers(int new_buffer_size){
  cout << "RESIZING FROM " <<  buffer_size << " TO " << new_buffer_size << endl;
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
 
  delete[] temp;

  b_XDX_ptr = b_XDX; 
  b_XSN_ptr = b_XSN;
    

  buffer_size = new_buffer_size;
  
 
}
/** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
void Particles1Dcomm::bufferXleft(double *b_, int np_current, VirtualTopology* vct){
   if (x[np_current] < 0)
	   b_[npExitXleft*nVar]    = x[np_current] + Lx; // this applies to the the leftmost processor
   else
	   b_[npExitXleft*nVar]    = x[np_current];
   b_[npExitXleft*nVar +1] = u[np_current];
   b_[npExitXleft*nVar +2] = q[np_current];
   if (TrackParticleID)
	   b_[npExitXleft*nVar +3] = ParticleID[np_current];
   if (cVERBOSE)
        cout << "Particle exiting to Xleft: X=" << x[np_current] << " ("<< xstart<<"," << xend << ")"<< endl;
   
}
/** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
void Particles1Dcomm::bufferXright(double *b_, int np_current, VirtualTopology* vct){
  
  if (x[np_current] > Lx)
        b_[npExitXright*nVar]    = x[np_current] - Lx; // this applies to the right most processor
   else
	b_[npExitXright*nVar]    = x[np_current];
    b_[npExitXright*nVar + 1] = u[np_current];
    b_[npExitXright*nVar + 2] = q[np_current];
   if (TrackParticleID)
	   b_[npExitXright*nVar + 3] = ParticleID[np_current];
   if(cVERBOSE)
        cout << "Particle exiting to Xright: X=" << x[np_current] << " ("<< xstart<< "," << xend << ")" << endl;
   
}
/** Unpack the received buffer: 
  * take the data from the buffer and add particles to the domain 
  * check if it is the right domain:
  * with implicit scheme particles can transverse more than one domain*/
void Particles1Dcomm::unbuffer(double *b_){
   int np_current =0;
   // put the new particles at the end of the array, and update the number of particles
   while(b_[np_current*nVar] != MIN_VAL){
          x[nop] = b_[nVar*np_current];
	  u[nop] = b_[nVar*np_current+1];
	  q[nop] = b_[nVar*np_current+2];
	  if (TrackParticleID)
		    ParticleID[nop]=(unsigned long) b_[nVar*np_current+3];
          np_current++;
	  if (cVERBOSE)
              cout << "Receiving Particle: X=" << x[nop]  << " ("<< xstart<<"," << xend << ")" <<endl;
	  if (x[nop] < xstart || x[nop] > xend) 
	  	rightDomain++; // the particle is not in the domain
	  nop++;
          if (nop > (npmax - (int) (.01*npmax) ) )
             cout << "Exceeding npmax: PArticles need to be resized" << endl;
   }
  
}
/** Delete the a particle from the array and pack the the array, update the number of 
 * particles that are exiting
 * For deleting the particle from the array take the last particle and put it
 * in the position of the particle you want to delete
 * @param np = the index of the particle that must be deleted
 * @param nplast = the index of the last particle in the array
 */
void Particles1Dcomm::del_pack(int np_current, int *nplast){
  x[np_current] = x[*nplast];
  u[np_current] = u[*nplast];
  q[np_current] = q[*nplast];
  if (TrackParticleID)
	  ParticleID[np_current]=ParticleID[*nplast];
  npExit++;
  (*nplast)--;
}
/** method to calculate how many particles are out of right domain */
int Particles1Dcomm::isMessagingDone(VirtualTopology* ptVCT){
   int result = 0;
   result = reduceNumberParticles(rightDomain);
   if (result > 0 && cVERBOSE && ptVCT->getCartesian_rank()==0)
      cout << "Further Comunication: " << result << " particles not in the right domain" << endl;
   return(result);

}
/** calculate the maximum number exiting from this domain */
int Particles1Dcomm::maxNpExiting(){
   int maxNp = 0;
   if (npExitXright > maxNp)
    maxNp = npExitXright;
   if (npExitXleft  > maxNp)
    maxNp = npExitXleft;
   return(maxNp);
}
/** return X-coordinate of particle array */
double* Particles1Dcomm::getXall() const{ return(x);}

/** get X-velocity of particle with label indexPart */
double* Particles1Dcomm::getUall() const{ return(u);}

/**get ID of particle with label indexPart */
unsigned long* Particles1Dcomm::getParticleIDall() const{return (ParticleID);}

/**get charge of particle with label indexPart */
double* Particles1Dcomm::getQall() const{ return(q);}

/** return X-coordinate of particle with index indexPart */
double Particles1Dcomm::getX(int indexPart) const{	return(x[indexPart]);}


/**get ID of particle with label indexPart */
unsigned long Particles1Dcomm::getParticleID(int indexPart) const{

 return(ParticleID[indexPart]);}


/**get charge of particle with label indexPart */
double Particles1Dcomm::getQ(int indexPart) const{return(q[indexPart]);}


/** return the number of particles */
int Particles1Dcomm::getNOP() const{return(nop);}

/** print particles info */
void Particles1Dcomm::Print(VirtualTopology* ptVCT)const{
    cout << endl;
    cout << "Number of Particles: " << nop << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ","  << ptVCT->getCoordinates(2) << ")"<< endl;
    cout <<  "Xin = " << xstart << "; Xfin = " << xend << endl;
    cout <<  "Number of species = " << ns << endl;
    for (int i=0; i < nop; i++)
     cout << "Particles #" << i << " x=" << x[i] << " u=" << u[i] << endl;
    cout << endl;
}
/** print just the number of particles */
void Particles1Dcomm::PrintNp(VirtualTopology* ptVCT)const{
    cout << endl;
    cout << "Number of Particles of species "<< ns << ": " << nop << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ","  << ptVCT->getCoordinates(2) << ")"<< endl;
    cout << endl;
}



/** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles1Dcomm::bufferYleft(double *b_, int np, VirtualTopology* vct){
	
       cout << "Method not implemented in 1D particle" << endl;
  

}
/** put a particle exiting to Y-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
inline void Particles1Dcomm::bufferYright(double *b_, int np, VirtualTopology* vct){
	 cout << "Method not implemented in 1D particle" << endl;
}

/** return Y-coordinate  of particle with index indexPart */
double Particles1Dcomm::getY(int indexPart) const{ 
cout << "1D Particle in X space. no need for calling Particles1Dcomm::getY()" << endl;
return(0.0);}

/** return Y-coordinate  of particle with index indexPart */
double Particles1Dcomm::getZ(int indexPart) const{ 
  cout << "1D Particle in X space. no need for calling Particles1DcommX::getZ(int indexPart) " << endl;
  return(0.0);
}

/** get u (X-velocity) of particle with label indexPart */
double Particles1Dcomm::getU(int indexPart) const{ return(u[indexPart]);}


/** get v (Y-velocity) of particle with label indexPart */
double Particles1Dcomm::getV(int indexPart) const{ 
cout << "1D Particle in X space. no need for calling Particles1Dcomm::getZ(int indexPart) " << endl;
return(0.0);}

/**get w (Z-velocity) of particle with label indexPart */
double Particles1Dcomm::getW(int indexPart) const{
cout << "1D Particle in X space. no need for calling Particles1Dcomm::getW(int indexPart) " << endl;
 return(0.0);}

/** get Y-velocity of particle with label indexPart */
double* Particles1Dcomm::getVall() const{
 cout << "1D Particle in X space. no need for calling Particles1Dcomm::getVall()" << endl;
 return(u);}

/**get Z-velocity of particle with label indexPart */
double* Particles1Dcomm::getWall() const{ 
cout << "1D Particle in X space. no need for calling Particles1Dcomm::getWall()" << endl;
return(u);}

/** return Y-coordinate  of particle array */
double* Particles1Dcomm::getYall() const{ 
cout << "1D Particle in X space. no need for calling Particles1Dcomm::getYall()" << endl;
return(x);}
/** return Z-coordinate  of particle array*/
double* Particles1Dcomm::getZall() const{ 
  cout << "1D Particle in X space. no need for calling Particles1Dcomm::getZall()" << endl;
  return(x);
}








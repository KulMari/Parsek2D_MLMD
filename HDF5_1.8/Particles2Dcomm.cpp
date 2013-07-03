/*******************************************************************************************
Particles2Dcomm.cpp  -  Class for particles of the same species, in a 2D space and 3component velocity
-------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
Modified for HDF5 1.8 -- Pierre HENRI -- June 2012
********************************************************************************************/

//---HPM on HP6
// #include "/usr/lpp/ppe.hpct/include/libhpc.h"

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
}
/** constructors fo a single species*/
void Particles2Dcomm::allocate(int species, CollectiveIO* col, VirtualTopology* vct, Grid* grid){
    // info from collectiveIO
    ns = species;
    npcel  = col->getNpcel(species);
    npcelx = col->getNpcelx(species);
    npcely = col->getNpcely(species);
    nop   =  col->getNp(species)/(vct->getNprocs());
    npmax =  col->getNpMax(species)/(vct->getNprocs());
    qom   = col->getQOM(species);
    uth   = col->getUth(species);
    vth   = col->getVth(species);
    wth   = col->getWth(species);
    u0    = col->getU0(species);
    v0    = col->getV0(species);
    w0    = col->getW0(species);
    dt    = col->getDt();
    Lx     = col->getLx();
    Ly     = col->getLy();
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
    //ID
    if (TrackParticleID){
		ParticleID= new unsigned long long[npmax];
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
    // BUFFERS
    // the buffer size should be decided depending on number of particles
    // the buffer size should be decided depending on number of particles
    if (TrackParticleID)
		nVar=12;
    else
		nVar=11;
	
	buffer_size = (int) (.05*nop*nVar+1); // max: 5% of the particles in the processors is going out
										  //buffer_size = 10;
//cout << "INITIAL buffer size =   " << buffer_size            << endl;
//cout << "(buffer_size*2) =       " << (buffer_size*2)        << endl;
//cout << "(int) (buffer_size*2) = " << (int) (buffer_size*2)  << endl;
 
	b_XDX = new double[buffer_size];
	b_XDX_ptr = b_XDX; // alias to make the resize
	b_XSN = new double[buffer_size];
	b_XSN_ptr = b_XSN; // alias to make the resize
	b_YDX = new double[buffer_size];
	b_YDX_ptr = b_YDX; // alias to make the resize
	b_YSN = new double[buffer_size];
	b_YSN_ptr = b_YSN; // alias to make the resize
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
		dataset_id = H5Dopen1(file_id,name_dataset.c_str());  //--- HDF5 1.8
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
		dataset_id = H5Dopen1(file_id, name_dataset.c_str()); //--- HDF5 1.8
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,y);
		status = H5Dclose(dataset_id);

		// get u
		name_dataset = "/particles/species_" + species_name.str() + "/u/cycle_0";
		dataset_id = H5Dopen1(file_id, name_dataset.c_str()); //--- HDF5 1.8
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,u);
		status = H5Dclose(dataset_id);
		// get v
		name_dataset = "/particles/species_" + species_name.str() + "/v/cycle_0";
		dataset_id = H5Dopen1(file_id, name_dataset.c_str()); //--- HDF5 1.8
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,v);
		status = H5Dclose(dataset_id);
		// get w
		name_dataset = "/particles/species_" + species_name.str() + "/w/cycle_0";
		dataset_id = H5Dopen1(file_id, name_dataset.c_str()); //--- HDF5 1.8
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,w);
		status = H5Dclose(dataset_id);
		// get q
		name_dataset = "/particles/species_" + species_name.str() + "/q/cycle_0";
		dataset_id = H5Dopen1(file_id, name_dataset.c_str()); //--- HDF5 1.8
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,q);
		status = H5Dclose(dataset_id);
		// ID
		if (TrackParticleID){
			herr_t (*old_func)(void*);
			void *old_client_data;
			H5Eget_auto1(&old_func, &old_client_data); //--- HDF5 1.8
            /* Turn off error handling */
            H5Eset_auto1(NULL, NULL); //--- HDF5 1.8
			name_dataset = "/particles/species_" + species_name.str() + "/ID/cycle_0";
			dataset_id = H5Dopen1(file_id, name_dataset.c_str()); //--- HDF5 1.8

            H5Eset_auto1(old_func, old_client_data); //--- HDF5 1.8
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
	double*** weight = newArr3(double,2,2,1);
	double*** temp = newArr3(double,2,2,1);
	int ix,iy, temp2,temp1;
	double inv_dx, inv_dy;
	inv_dx = 1.0/dx;
	inv_dy = 1.0/dy;


//if (vct->getCartesian_rank() ==0){
//         cout << "entering  interpP2G" << endl;
//       	 cout << "nop = " << nop << endl; }


	for (register int i=0; i < nop; i++){

		ix = 2 +  int(floor((x[i]-xstart)*inv_dx));
		iy = 2 +  int(floor((y[i]-ystart)*inv_dy));
		weight[1][1][0] = ((x[i] - grid->getXN(ix-1,iy-1,0))*inv_dx)*((y[i] - grid->getYN(ix-1,iy-1,0))*inv_dy);
		weight[1][0][0] = ((x[i] - grid->getXN(ix-1,iy,0))*inv_dx)*((grid->getYN(ix-1,iy,0) - y[i])*inv_dy);
		weight[0][1][0] = ((grid->getXN(ix,iy-1,0) - x[i])*inv_dx)*((y[i] - grid->getYN(ix,iy-1,0))*inv_dy);
		weight[0][0][0] = ((grid->getXN(ix,iy,0) - x[i])*inv_dx)*((grid->getYN(ix,iy,0) - y[i])*inv_dy);
		scale(weight,q[i],2,2);

		// add charge density
		EMf->addRho(weight,ix,iy,0,ns);
                // rho tagged
                if (ParticleID[i] == 1)
                        EMf->addRhotag(weight,ix,iy,0,ns);
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

	// communicate contribution from ghost cells
	EMf->communicateGhostP2G(ns,0,0,0,0,vct);
	delArr3(weight,2,2);
	delArr3(temp,2,2);
}


/** communicate buffers */
int Particles2Dcomm::communicate(VirtualTopology* ptVCT){

	// allocate buffers
	MPI_Status status;
	int new_buffer_size;
	int npExitingMax;

	// variable for memory availability of space for new particles
	int avail, availALL, avail1, avail2, avail3, avail4;
	for (int i=0; i < buffer_size; i++){
		b_XDX[i] = MIN_VAL;
		b_XSN[i] = MIN_VAL;
		b_YDX[i] = MIN_VAL;
		b_YSN[i] = MIN_VAL;
	}
	npExitXright =0, npExitXleft =0, npExitYright =0, npExitYleft =0, npExit=0, rightDomain = 0;
	npDeletedBoundary = 0;
	int np_current = 0, nplast = nop-1;

	while (np_current < nplast+1){
		// BC on particles
		if (x[np_current] < 0 && ptVCT->getXleft_neighbor() == -1)
			BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft);
		else if (x[np_current] > Lx && ptVCT->getXright_neighbor() == -1)
			BCpart(&x[np_current],&u[np_current],&v[np_current],&w[np_current],Lx,uth,vth,wth,bcPfaceXright,bcPfaceXleft);
		if (y[np_current] < 0 && ptVCT->getYleft_neighbor() == -1)  // check it here
			BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft);
		else if (y[np_current] > Ly && ptVCT->getYright_neighbor() == -1) //check it here
			BCpart(&y[np_current],&v[np_current],&u[np_current],&w[np_current],Ly,vth,uth,wth,bcPfaceYright,bcPfaceYleft);
		// if the particle exits, apply the boundary conditions add the particle to communication buffer
		if (x[np_current] < xstart || x[np_current] >xend){
        	// communicate if they don't belong to the domain
			if (x[np_current] < xstart && ptVCT->getXleft_neighbor() != -1){
				// check if there is enough space in the buffer before putting in the particle
				if(((npExitXleft+1)*nVar)>=buffer_size){
//						cout << "resizing the sending buffer to " << (int) (buffer_size*2) << " buffer size" << endl;
//						resize_buffers((int) (buffer_size*2));
					new_buffer_size = (int) (1.5*((npExitXleft+1)*nVar));
                                      	cout << "resizing the sending buffer to " << new_buffer_size << " buffer size" << endl;
                                      	resize_buffers(buffer_size, new_buffer_size);
					buffer_size = new_buffer_size;
				}
				// put it in the communication buffer
				bufferXleft(b_XSN,np_current,ptVCT);
				//cout << "after buffering" << endl;
				// delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);
				npExitXleft++;
        	} else if (x[np_current] > xend && ptVCT->getXright_neighbor() != -1){
				// check if there is enough space in the buffer before putting in the particle
				if(((npExitXright+1)*nVar)>=buffer_size){
//						cout << "resizing the sending buffer " << (int) (buffer_size*2) << endl;
//						resize_buffers((int) (buffer_size*2));
                                        new_buffer_size = (int) (1.5*((npExitXright+1)*nVar));
                                        cout << "resizing the sending buffer to " << new_buffer_size << " buffer size" << endl;
                                        resize_buffers(buffer_size, new_buffer_size);
					buffer_size = new_buffer_size;
				}
				// put it in the communication buffer
				bufferXright(b_XDX,np_current,ptVCT);
				// delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);

				npExitXright++;
        	} else if (x[np_current] < xstart && ptVCT->getXleft_neighbor() == -1 && bcPfaceXleft == 0){  // exit from domain and delete
																										  // delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);
				npDeletedBoundary++;
			} else if (x[np_current] > xend && ptVCT->getXright_neighbor() == -1 && bcPfaceXright == 0){
				// delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);
				npDeletedBoundary++;

			}

		} else  if (y[np_current] < ystart || y[np_current] >yend){
        	// communicate if they don't belong to the domain
			if (y[np_current] < ystart && ptVCT->getYleft_neighbor() != -1){
				// check if there is enough space in the buffer before putting in the particle
				if(((npExitYleft+1)*nVar)>=buffer_size){
// 						cout << "resizing the sending buffer " << (int) (buffer_size*2) << endl;
// 						resize_buffers((int) (buffer_size*2));
                                        new_buffer_size = (int) (1.5*((npExitYleft+1)*nVar));
                                        cout << "resizing the sending buffer to " << new_buffer_size << " buffer size" << endl;
                                        resize_buffers(buffer_size, new_buffer_size);
                                        buffer_size = new_buffer_size;
				}
				// put it in the communication buffer
				bufferYleft(b_YSN,np_current,ptVCT);
				// delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);
				npExitYleft++;
        	} else if (y[np_current] > yend && ptVCT->getYright_neighbor() != -1){
				// check if there is enough space in the buffer before putting in the particle
				if(((npExitYright+1)*nVar)>=buffer_size){
//						cout << "resizing the sending buffer " << (int) (buffer_size*2) << endl;
//						resize_buffers((int) (buffer_size*2));
                                        new_buffer_size = (int) (1.5*((npExitYright+1)*nVar));
                                        cout << "resizing the sending buffer to " << new_buffer_size << " buffer size" << endl;
                                        resize_buffers(buffer_size, new_buffer_size);
                                        buffer_size = new_buffer_size;
				}
				// put it in the communication buffer
				bufferYright(b_YDX,np_current,ptVCT);
				// delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);

				npExitYright++;
		    } else if (y[np_current] < ystart && ptVCT->getYleft_neighbor() == -1 && bcPfaceYleft == 0 ){
				// delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);
				npDeletedBoundary++;
			} else if (y[np_current] > yend && ptVCT->getYright_neighbor() == -1 && bcPfaceYright == 0 ){
				// delete the particle and pack the particle array, the value of nplast changes
				del_pack(np_current,&nplast);
				npDeletedBoundary++;
			}
		}  else {
			// particle is still in the domain, procede with the next particle
			np_current++;
		}

	}

	nop = nplast+1;
	npExitingMax = 0;
	// calculate the maximum number of particles exiting from this domain
	// use this value to check if communication is needed
// no 	// and to  resize the buffer
	npExitingMax = maxNpExiting();
	// broadcast the maximum number of particles exiting for sizing the buffer and to check if communication is really needed
	npExitingMax = reduceMaxNpExiting(npExitingMax,ptVCT->getCommunicator());

	/*****************************************************/
	/*           SEND AND RECEIVE MESSAGES               */
	/*****************************************************/

//	new_buffer_size = npExitingMax*nVar + 1;
//	if (new_buffer_size > buffer_size){
//		cout << "resizing the receiving buffer " << endl;
//		resize_buffers(new_buffer_size);
//	}

	if (npExitingMax > 0){

		communicateParticles(new_buffer_size,b_XSN,b_XDX,b_YSN,b_YDX,ptVCT);

		// UNBUFFERING
		// message from XLEFT
		avail1 = unbuffer(b_XDX);
		// message from XRIGHT
		avail2 = unbuffer(b_XSN);
		// message from XLEFT
		avail3 = unbuffer(b_YDX);
		// message from XRIGHT
		avail4 = unbuffer(b_YSN);
		// if one of these numbers is negative than there is not enough space for particles
		avail = avail1 + avail2 + avail3 + avail4;
		availALL = reduceNumberParticles(avail,ptVCT->getCommunicator());
		if (availALL < 0)
			return(-1);  // too many particles coming, save data nad stop simulation

	}

	return(0); // everything was fine

}
/** resize the buffers */
void Particles2Dcomm::resize_buffers(int buffer_size, int new_buffer_size){
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


}
/** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
void Particles2Dcomm::bufferXleft(double *b_, int np_current, VirtualTopology* vct){
	if (x[np_current] < 0)
		b_[npExitXleft*nVar]    = x[np_current] + Lx*(floor(-x[np_current]/Lx) + 1); // this applies to the the leftmost processor
	else
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

	if (x[np_current] > Lx)
        b_[npExitXright*nVar]    = x[np_current] - Lx*(floor(x[np_current]/Lx)); // this applies to the right most processor
	else
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
	if (y[np] < 0){
		b_[npExitYleft*nVar +1] = y[np] + Ly*(floor(-y[np]/Ly) + 1);
		y[np] += Ly*(floor(-y[np]/Ly) + 1);
	} else {
		b_[npExitYleft*nVar +1] = y[np];
	}
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
	if (y[np]  > Ly){  b_[npExitYright*nVar +1] = y[np] - Ly*(floor(y[np]/Ly));
        y[np] -= Ly*(floor(y[np]/Ly));
	} else {
		b_[npExitYright*nVar +1] = y[np];
	}
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
		np_current++;

		if (cVERBOSE)
			cout << "Receiving Particle: X=" << x[nop] << ",Y=" << y[nop] << " ("<< xstart<<"," << xend << ")"<< " x ("<< ystart<<"," << yend << ")"<<endl;
		if (x[nop] < xstart || x[nop] > xend || y[nop] < ystart || y[nop] > yend)
			rightDomain++; // the particle is not in the domain
		nop++;
		if (nop > (npmax - (int) (.01*npmax) ) ){
			cout << "Exceeding npmax: Particles need to be resized Save Data and Stop the simulation" << endl;
			return(-1); // end the simulation because you dont have enough space on the array
		}
	}
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

	npExit++;
	(*nplast)--;
}
/** method to calculate how many particles are out of right domain */
int Particles2Dcomm::isMessagingDone(VirtualTopology* ptVCT){
	int result = 0;
	result = reduceNumberParticles(rightDomain, ptVCT->getCommunicator());
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
unsigned long long* Particles2Dcomm::getParticleIDall() const{return (ParticleID);}

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












/*******************************************************************************************
ESfield1D.h  -  Electrostatic Field for 1D geometry
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef ESfield1D_H
#define ESfield1D_H


#include <iostream>
#include <sstream>

#include <math.h>
#include <mpi.h>


#include "../utility/Alloc.h"
#include "../mathlib/Basic.h"
#include "../utility/TransArraySpace.h"
#include "../solvers/CG.h"
#include "../solvers/GMRES_new2.h"
#include "hdf5.h"

using std::cout;
using std::cerr;
using std::endl;



/**
* 
* Electrostatic Field in 1D Geometry
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/

class ESfield1D : public Field{
  public:
  /** constructor */
  ESfield1D(CollectiveIO *col,Grid *grid);
  /** destructor */
  ~ESfield1D();

  /** initialize ES field */
  void init(VirtualTopology *vct, Grid *grid);
 
  /** communicate ghost for densities and interp rho from node to center */
  void interpDensitiesN2C(VirtualTopology *vct, Grid *grid);
  /** Calculate Electric field solving Posson Equation */
  void calculateField(Grid *grid, VirtualTopology *vct);
  /** Image of Poisson Solver (for SOLVER)*/
  void PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct);
  /** Image of Maxwell Solver (for Solver) */
  void MaxwellImage(double *im, double *vector, Grid *grid,VirtualTopology *vct);     
  /** Maxwell source term (for SOLVER) */
  void MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct);
  /** Smoothing after the interpolation**/
  void smooth(double value ,double ***vector,bool type, Grid *grid, VirtualTopology *vct);
  /** SPECIES: Smoothing after the interpolation for species fields**/
  void smooth(double value ,double ****vector,int is, bool type, Grid *grid, VirtualTopology *vct);
  
  /** set to 0 all the densities fields */
  void setZeroDensities();
  /** Sum rhon over species */
  void sumOverSpecies(VirtualTopology *vct);
  /** Sum current over different species */
  void sumOverSpeciesJ();

  /** communicate ghost for grid -> Particles interpolation */
  void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct);
    /** add an amount of charge density to charge density field at node X,Y,Z */
  void addRho(double ***weight, int X, int Y,int Z, int ns);
  /** add an amount of current density - direction X to current density field at node X,Y,Z */
  void addJx(double ***weight, int X, int Y, int Z,int ns);
  /** add an amount of current density - direction Y to current density field at node X,Y,Z */
  void addJy(double ***weight, int X, int Y, int Z, int ns);
  /** add an amount of current density - direction Z to current density field at node X,Y,Z */
  void addJz(double ***weight, int X, int Y, int Z, int ns);
 
  /** add an amount of pressure density - direction XX to current density field at node X,Y,Z */
  void addPxx(double ***weight, int X, int Y, int Z, int ns);
  /** add an amount of pressure density - direction XY to current density field at node X,Y,Z */
  void addPxy(double ***weight, int X, int Y, int Z,int ns);
  /** add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
  void addPxz(double ***weight, int X, int Y, int Z,int ns);
  /** add an amount of pressure density - direction YY to current density field at node X,Y,Z */
  void addPyy(double ***weight, int X, int Y, int Z, int ns);
  /** add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
  void addPyz(double ***weight, int X, int Y, int Z, int ns);
  /** add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
  void addPzz(double ***weight, int X, int Y, int Z,int ns);
  
  

  
  /** get Potential array */
  double*** getPHI();
  /** get Electric Field component X defined on node(indexX,indexY,indexZ) */
  double &getEx(int indexX, int indexY, int indexZ) const;
  /** get Electric field X component array */
  double*** getEx() ;
   /** get Electric Field component Y defined on node(indexX,indexY,indexZ) */
  double &getEy(int indexX, int indexY, int indexZ) const;
  /** get Electric field Y component array */
  double*** getEy() ;
  /** get Electric Field component Z defined on node(indexX,indexY,indexZ) */
  double &getEz(int indexX, int indexY, int indexZ) const;
  /** get Electric field Z component array */
  double*** getEz() ;
  /** get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
  double &getBx(int indexX, int indexY, int indexZ) const;
  /** get Magnetic field X component array */
  double*** getBx();
  /** get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
  double &getBy(int indexX, int indexY, int indexZ) const;
  /** get Magnetic field Y component array */
  double*** getBy();
  /** get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
  double &getBz(int indexX, int indexY, int indexZ) const;
  /** get Magnetic field Z component array */
  double*** getBz() ;
  /** get density on cell(indexX,indexY,indexZ)  */
  double &getRHOc(int indexX, int indexY, int indexZ) const;
  /** get density array on center cell */
  double*** getRHOc() ;
  /** get density on nodes(indexX,indexY,indexZ)  */
  double &getRHOn(int indexX, int indexY, int indexZ) const;
  /** get density array on nodes  */
  double*** getRHOn() ;
  /** SPECIES: get density on nodes(indexX,indexY,indexZ)*/
  double &getRHOns(int indexX, int indexY,int indexZ,int ns) const;
  /** SPECIES: get density on center cell(indexX,indexY,indexZ)*/
  double &getRHOcs(int indexX, int indexY,int indexZ,int ns) const;
  /** SPECIES: get density array on nodes*/
  double**** getRHOns() ;



  /** get pressure tensor XX for species */
  double**** getpXXsn() ;
  /** get pressure tensor XY for species */
  double**** getpXYsn() ;
  /** get pressure tensor XZ for species*/
  double**** getpXZsn() ;
  /** get pressure tensor YY for species */
  double**** getpYYsn() ;
  /** get pressure tensor YZ for species */
  double**** getpYZsn() ;
  /** get pressure tensor ZZ for species */
  double**** getpZZsn() ;
    
  /** get Jx(X,Y,Z)  */
  double &getJx(int indexX, int indexY, int indexZ) const;
  /** get current -Direction X */
  double*** getJx();
  /** SPECIES: get current -Direction X */
  double**** getJxs();
  /** get Jy(X,Y,Z)  */
  double &getJy(int indexX, int indexY, int indexZ) const;
  /** get current -Direction Y */
  double*** getJy();
  /** SPECIES: get current -Direction Y */
  double**** getJys();
  /** get Jz(X,Y,Z)  */
  double &getJz(int indexX, int indexY, int indexZ) const;
  /** get current -Direction Z */
  double*** getJz();
  /** SPECIES: get current -Direction Z */
  double**** getJzs();
  
  
  /** print electromagnetic fields info */
  void print(void) const;
  


 private:
  /** light speed in vacuum*/
  double c;
  /** permeattivity in vacumm */
  double eps0;
  /** permeability in vacumm */
  double mu0;
  /** pi */
  double pi;
   /** Four PI, for normalization */
  double FourPI;
  /** time step */
  double dt;
  /** decentering parameter */
  double th;
  /** Smoothing value*/
  double Smooth;
  /** delt = c*th*dt */
  double delt;
  /** number of particles species */
  int ns;
  /** charge to mass ratio array for different species */
  double *qom;
 
 
  /** number of cells - X direction, including + 2 (guard cells) */
  int nxc;
  /** number of nodes - X direction, including + 2 extra nodes for guard cells */
  int nxn;
  /** grid spacing */
  double dx, invVOL;
  /** simulation box length - X direction   */
  double Lx;
  /** PHI: ES potential (indexX, indexY, indexZ), defined on central points between nodes*/
  double ***PHI;
  /** Ex: electric field X-component (indexX, indexY, indexZ), defined on nodes */
  double ***Ex;
  /** rhoc: charge density (indexX, indexY, indexZ), defined on central points between nodes */
  double ***rhoc;
  /** rhon: charge density (indexX, indexY, indexZ), defined on nodes */
  double ***rhon;
  /* rhons: charge density for each species (indexX, indexY, indexZ, #species), defined on nodes */
  double ****rhons;
  /* rhons: charge density for each species (indexX, indexY, indexZ, #species), defined on nodes */
  double ****rhocs;
  /** Field Boundary Condition
   0 =  Dirichlet Boundary Condition: specifies the value to take on the boundary of the domain
   1 =  Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain
   2 =  Periodic condition
  */
  /** Boundary Condition Electrostatic Potential: FaceXright */
  int bcPHIfaceXright;
  /** Boundary Condition Electrostatic Potential:FaceXleft */
  int bcPHIfaceXleft;
  /** Boundary Condition Electrostatic Potential:FaceYright */
  int bcPHIfaceYright;
  /** Boundary Condition Electrostatic Potential:FaceYleft */
  int bcPHIfaceYleft;
  /** RESTART BOOLEAN: useful for loading particles from restart file*/
  int restart1 ; 
  /** String with the directory for the restart file */
  string RestartDirName;
  
  /** CG tolerance for convergence*/
  double CGtol;

};
/** constructor */
inline ESfield1D::ESfield1D(CollectiveIO *col,Grid *grid){
	nxc = grid->getNXC();
	nxn = grid->getNXN();
	dx = grid->getDX();
	invVOL = grid->getInvVOL();
	Lx = col->getLx();
	ns  = col->getNs();
	c = col->getC();
	pi = 3.1415;
	mu0 = 4*pi*1.e-7;
	eps0 = 1.0/(mu0*c*c);
	dt = col->getDt();
	Smooth = col->getSmooth();
	th = col->getTh();
	delt = c*th*dt;
	FourPI =16*atan(1.0);
	CGtol = col->getCGtol();
	qom = new double[ns];
	for (int i=0; i < ns;i++)
		qom[i] = col->getQOM(i);
	// boundary conditions: PHI and EM fields
	bcPHIfaceXright  = col->getBcPHIfaceXright();
	bcPHIfaceXleft   = col->getBcPHIfaceXleft();
	// Restart 
	restart1 = col->getRestart_status();
	RestartDirName = col->getRestartDirName();
	////////////////// ARRAY ALLOCATION ////////////////
	// arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
	Ex = newArr3(double,nxn,1,1);
	rhon  = newArr3(double,nxn,1,1);
	// involving species
	rhons = newArr4(double,ns,nxn,1,1);
	rhocs = newArr4(double,ns,nxc,1,1);
	// arrays allocation: central points ---> the first central point has index 1, the last has index nxn-2!
	PHI    = newArr3(double,nxc,1,1);
	rhoc   = newArr3(double,nxc,1,1);

}
/** destructor: deallocate arrays*/
inline ESfield1D::~ESfield1D(){
        // nodes
	delArr3(Ex,nxn,1);	
	delArr3(rhon,nxn,1);
	// nodes and species
	delArr4(rhons,ns,nxn,1);
	// central points
	delArr3(PHI,nxc,1);
	delArr3(rhoc,nxc,1);
	// center and species 
	delArr4(rhocs,ns,nxc,1);
}
/** Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
inline void ESfield1D::calculateField(Grid *grid, VirtualTopology *vct){
  if (vct->getCartesian_rank() ==0)
    cout << "*** PHI CALCULATION ***" << endl;
  double *xkrylovPoisson = new double[nxc-2];
  double *bkrylovPoisson = new double[nxc-2];
  double ***gradPHI = newArr3(double,nxn,1,1);
  eqValue (0.0,gradPHI,nxn);
  // Calculate the laplacian(PHI) =   -4*PI*rho
  phys2solver(bkrylovPoisson,rhoc,nxc);    // bkrylov is place holder for source term
  scale(bkrylovPoisson,-FourPI,nxc-2);      
  // *******************************************************************
	// PUT THE POTENTIAL EQUAL TO ZERO IF YOU WANT TO USE PERIODIC WITH CG 
	// this is done to use CG even with the periodic
	// in image you put the value of PHI
	if (vct->getCartesian_rank()==0){
	   bkrylovPoisson[0] = 0.0;
	}
	if (vct->getCartesian_rank()==(vct->getXLEN()-1)){
	   bkrylovPoisson[nxc-3] = 0.0;
	}
	// *******************************************************************  
  eqValue(0.0,xkrylovPoisson,nxc-2); // in the initial guess xkrylov is 0.0 -> xkrylov is the unknown
  // CG does not converge for double periodic
  CG(xkrylovPoisson,(nxc-2),bkrylovPoisson, 100000,CGtol, &Field::PoissonImage, grid, vct, this);  // CG
  //GMRES(&Field::PoissonImage, xkrylovPoisson,(nxc-2),bkrylovPoisson,20,200,1E-4, grid, vct, this);  // GMRES
  solver2phys(PHI,xkrylovPoisson,nxc);
  communicateCenterBC(nxc,PHI,bcPHIfaceXright,bcPHIfaceXleft,vct);
  grid->gradC2N(Ex,Ex,PHI); // calculate the electric field on the nodes
  neg(Ex,nxn);
  communicateNode(nxn,Ex,vct);
  // smoothing
  smooth(Smooth,Ex,1,grid,vct);
  // deallocate 
  delete[] xkrylovPoisson;
  delete[] bkrylovPoisson;
  delArr3(gradPHI,nxn,1);
}
/** Image of Poisson Solver: it's the argument for the GMRES */
inline void ESfield1D::PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct){
    // allocate  2 two dimensional service vectors
    double ***temp = newArr3(double,nxc,1,1);
    double ***im   = newArr3(double,nxc,1,1);
    // initialize to zero, new vectors
    eqValue (0.0, image,nxc-2); // added by Enrico.
    eqValue (0.0, temp,nxc);
    eqValue (0.0, im,nxc);
    // move from krylov space to physical space and communicate ghost cells
    solver2phys(temp,vector,nxc);
    MPI_Barrier(MPI_COMM_WORLD);
    communicateCenterBC(nxc,temp,bcPHIfaceXright,bcPHIfaceXleft,vct);
    // calculate the laplacian
    grid->lapC2C(im,temp,vct);
    
	
	// move from physical space to krylov space, in this case ghost cells don't count   
    phys2solver(image,im,nxc);
	// *******************************************************************
	// PUT THE POTENTIAL EQUAL TO ZERO IF YOU WANT TO USE PERIODIC WITH CG 
	// this is done to use CG even with the periodic
	// in image you put the value of PHI
	if (vct->getCartesian_rank()==0){
	   image[0] = vector[0];
	}
	if (vct->getCartesian_rank()==(vct->getXLEN()-1)){
	   image[nxc-3] = vector[nxc-3];
	}
	// ******************************************************************* 
	
    // deallocate temporary array and objects
    delArr3(temp,nxc,1);
    delArr3(im,nxc,1);
}
/** set to 0 all the densities fields: this is need before the interpolation */
inline  void ESfield1D::setZeroDensities(){
 for (int i=0; i < nxn; i++)
 	rhon[i][0][0] = 0.0;
  for (int i=0; i < nxc; i++)
	rhoc[i][0][0] = 0.0;

 for (int kk=0; kk < ns; kk++)
    for (int i=0; i < nxn; i++)
	  rhons[kk][i][0][0] = 0.0;

 for (int kk=0; kk < ns; kk++)
    for (int i=0; i < nxc; i++)
	  rhocs[kk][i][0][0] = 0.0;
   

}
/** communicate ghost cells, for grid->Particles interpolation */
inline void ESfield1D::communicateGhostP2G(int ns,int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct){
         // communicate interpolation + BC for interpolated fields
	 communicateInterp(nxn,ns,rhons,0,0,dx,vct); // add the density on the same node
	 communicateNode(nxn,rhons,ns,vct);          // communicate ghost nodes node
}   

/** Sum rhon over species: here you can add background ions*/
inline void ESfield1D::sumOverSpecies(VirtualTopology *vct){
	// electrons
	for (int is=0; is<ns; is++) //{
 	   for (register int i=1; i < nxn-1; i++)
	     rhon[i][0][0]     += rhons[is][i][0][0]*invVOL;
	//******************************
	//******************************
	// THIS IS FOR THE BACKGROUND IONS
	for (register int i=1; i < nxn-1; i++)
	   rhon[i][0][0] += 0.07957747154595;
	//******************************
	//******************************
}
/** interpolate rho from node to center */
inline  void ESfield1D::interpDensitiesN2C(VirtualTopology *vct, Grid *grid){
    communicateNode(nxn,rhon,vct);
    grid->interpN2C(rhoc,rhon);
    communicateCenter(nxc,rhoc,vct);
}

/** add an amount of charge density to charge density field at node X,Y */
inline void ESfield1D::addRho(double ***weight, int X, int Y, int Z, int ns){
	 // check the INTERPOLATION!!!!! weight
	 rhons[ns][X+1][0][0]+= weight[1][0][0];
	 rhons[ns][X][0][0]  += weight[0][0][0];
}
/** 
Interpolation smoothing:
Smoothing (vector must already have ghost cells)
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector ; 
type = 1 --> node based vector   ;                     
*/
inline void  ESfield1D::smooth(double value,double ***vector, bool type, Grid *grid, VirtualTopology *vct){
 double alpha;
 int nx;
 if (value != 1.0){
   switch(type){ 
    case (0):
	nx=grid->getNXC();
	break;
     case(1):
	nx=grid->getNXN();
	break;
     }

   double **temp = newArr(double,nx,1);
   alpha=(1.0-value)/4;

   for (int i=1; i<nx-1; i++)
	temp[i][0]=value*vector[i][0][0]+alpha*(vector[i-1][0][0]+vector[i][0][0]+vector[i+1][0][0]+vector[i][0][0]);
 // Adjust boundary
 if (vct->getXleft_neighbor()==-1){
	alpha=(1.0-value)/3;
	temp[1][0]=value*vector[1][0][0]+alpha*(vector[1][0][0]+vector[2][0][0]+vector[1][0][0]);}
	
if (vct->getXright_neighbor()==-1){
	alpha=(1.0-value)/3;
	temp[nx-2][0]=value*vector[nx-2][0][0]+alpha*(vector[nx-2][0][0]+vector[nx-3][0][0]+vector[nx-2][0][0]);}


 //eq(vector,temp,nx,ny);
 // copy temp in vector
 for (int i=0; i < nx;i++)
   vector[i][0][0] = temp[i][0];
  // communicate
  if (type == 0)
   communicateCenter(nx,vector,vct);
  else 
   communicateNode(nx,vector,vct);
  delArr(temp,nx);
 } // end of if (value !=1)
}

/** 

SPECIES: Interpolation smoothing
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector  
type = 1 --> node based vector                        
*/
inline void  ESfield1D::smooth(double value,double ****vector,int is, bool type, Grid *grid, VirtualTopology *vct){
double alpha;
int nx;
if (value != 1.0){
switch(type){ 
case (0):
	nx=grid->getNXC();
	break;
case(1):
	nx=grid->getNXN();
	break;
}

double **temp = newArr(double,nx,1);
alpha=(1.0-value)/4;

for (int i=1; i<nx-1; i++)
	temp[i][0]=value*vector[is][i][0][0]+alpha*(vector[is][i-1][0][0]+vector[is][i][0][0]+vector[is][i+1][0][0]+vector[is][i][0][0]);
// Adjust boundary
if (vct->getXleft_neighbor()==-1){
	alpha=(1.0-value)/3;
	temp[1][0]=value*vector[is][1][0][0]+alpha*(vector[is][1][0][0]+vector[is][2][0][0]+vector[is][1][0][0]);}
	
if (vct->getXright_neighbor()==-1){
	alpha=(1.0-value)/3;
	temp[nx-2][0]=value*vector[is][nx-2][0][0]+alpha*(vector[is][nx-2][0][0]+vector[is][nx-3][0][0]+vector[is][nx-2][0][0]);
}
	

	

//eq(vector,temp,nx,ny,is);
// copy temp in vector
 for (int i=0; i < nx;i++)
     vector[is][i][0][0] = temp[i][0];
 // communicate
 if (type == 0)
   communicateCenter(nx,vector,is,vct);
 else 
  communicateNode(nx,vector,is,vct);
 delArr(temp,nx);
 } // end of if
}
/** Initialize the EM field with constants values or from restart*/
inline void ESfield1D::init(VirtualTopology *vct, Grid *grid){
// initialize E and rhos on nodes
if (restart1==0){ 
  if (vct->getCartesian_rank()==0){
    cout << "Initialized constant density" << endl;
  }
      
  for (int i=0; i < nxn; i++)
	   Ex[i][0][0] =  0.0;
  for (int kk=0; kk < ns; kk++)
    for (int i=0; i < nxn; i++){
	   rhons[kk][i][0][0] =  0.07957747154595;  // = ion rhons, for plasma neutrality
	     // ion plasma frequency = 1 (normalization)
    }
  for (int kk=0; kk < ns; kk++)
    for (int i=0; i < nxc; i++){
	   rhocs[kk][i][0][0] =  0.07957747154595;  // = ion rhons, for plasma neutrality
	    // ion plasma frequency = 1 (normalization)
    }
  
} else { // EM initialization from RESTART

	if (vct->getCartesian_rank()==0)
	    cout << "LOADING ES FIELD FROM RESTART FILE in " + RestartDirName + "/restart.hdf" << endl;
	stringstream ss;
	ss << vct->getCartesian_rank();
	string name_file = RestartDirName + "/restart" + ss.str() + ".hdf";
	// hdf stuff 
        hid_t    file_id, dataspace;
        hid_t    datatype, dataset_id;
        herr_t   status;
	size_t   size;
	hsize_t     dims_out[2];           /* dataset dimensions */
	int status_n;
	 
	 // open the hdf file
	 file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
         if (file_id < 0){
           cout << "couldn't open file: " << name_file << endl;
	   cout << "RESTART NOT POSSIBLE" << endl;
	 }
	 dataset_id = H5Dopen(file_id,"/fields/Ex/cycle_0");
	 datatype  = H5Dget_type(dataset_id);  
	 size  = H5Tget_size(datatype);
	 dataspace = H5Dget_space(dataset_id);    
	 status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
	 // check if the grid is consistent with the grid of the RESTART
	 if  (nxn != (dims_out[0] + 2)){
	    cout << "RESTART GRID NOT CONSISTENT!!!!!!!!!!!!!!!!!!" << endl;
	 }
	 // read Exn
	 double *temp_storage = new double[dims_out[0]*dims_out[1]];
	 dataset_id = H5Dopen(file_id,"/fields/Ex/cycle_0");
	 status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	 H5S_ALL,H5P_DEFAULT,temp_storage);
	 int k=0;
	 for (int i=1; i < nxn-1; i++)
	        Ex[i][0][0] = temp_storage[k++]; 
	communicateNode(nxn,Ex,vct);
	status = H5Dclose(dataset_id);
	// open the charge density for species
	stringstream *species_name = new stringstream[ns];
	for (int is=0; is < ns;is++){ 
	   species_name[is] << is;
	   string name_dataset = "/moments/species_" + species_name[is].str() + "/rho/cycle_0";
	   dataset_id = H5Dopen(file_id,name_dataset.c_str());
	   status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storage);
	   k=0;
	   for (int i=1; i < nxn-1; i++)
	        rhons[is][i][0][0] = temp_storage[k++]; 
	   communicateNode(nxn,rhons,is,vct);
	   status = H5Dclose(dataset_id);
	 }
	// close the hdf file
	status = H5Fclose(file_id);
	}
	// calculate densities on the center cells
	for (int is=0 ; is<ns; is++)
           grid->interpN2C(rhocs,is,rhons);
 }




/////////////////////////////////////////////////////////
//////////            RETURN METHODS                 ////
/////////////////////////////////////////////////////////
/** get Potential array ***/
inline double*** ESfield1D::getPHI() {return(PHI);}
/** get Electric Field component X defined on node(indexX,indexY,indexZ) */  
inline double &ESfield1D::getEx(int indexX, int indexY,  int indexZ) const{
  return(Ex[indexX][0][0]);}
/** get Electric Field X component array */
inline double*** ESfield1D::getEx(){return(Ex);}


/** get density on center cell(indexX,indexY,indexZ)  */
inline double &ESfield1D::getRHOc(int indexX, int indexY,int indexZ) const{
  return(rhoc[indexX][0][0]);}
 /** get density array defined on centers cells */
inline double*** ESfield1D::getRHOc() {return(rhoc);}
/** get density on node(indexX,indexY,indexZ)  */
inline double &ESfield1D::getRHOn(int indexX, int indexY,int indexZ) const{ 
return(rhon[indexX][0][0]);}
/** get density array defined on nodes*/
inline double*** ESfield1D::getRHOn() {return(rhon);}
/** SPECIES: get density array defined on nodes  */
inline double &ESfield1D::getRHOns(int indexX, int indexY, int indexZ,int ns) const{
  return(rhons[ns][indexX][0][0]);}
/** SPECIES: get density array defined on center cells  */
inline double &ESfield1D::getRHOcs(int indexX, int indexY,int indexZ, int ns) const{
  return(rhocs[ns][indexX][0][0]);}
/** SPECIES: get density array defined on nodes  */
inline double**** ESfield1D::getRHOns() {return(rhons);}


/** get Electric Field component Z defined on node(indexX,indexY,indexZ) */  
inline double &ESfield1D::getEz(int indexX, int indexY,int indexZ) const{
   cout << " no need to call getEz(): 1D Geometry" << endl;
  return(Ex[indexX][0][0]);}
/** get Electric Field Z component array */
inline double*** ESfield1D::getEz() {
   cout << " no need to call getEz(): 1D Geometry" << endl;
  return(Ex);}
/** get Electric Field component Y defined on node(indexX,indexY,indexZ) */  
inline double &ESfield1D::getEy(int indexX, int indexY,int indexZ) const{
  cout << " no need to call getEy(): 1D Geometry" << endl;
  return(Ex[indexX][0][0]);}
/** get Electric Field Y component array */
inline double*** ESfield1D::getEy() {
 cout << " no need to call getEy(): 1D Geometry" << endl;
return(Ex);}
/** Sum current over different species */
inline void ESfield1D::sumOverSpeciesJ(){
   cout << "ES field, no need to call ESfield1D::sumOverSpeciesJ() " << endl;

}
/** get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
inline double &ESfield1D::getBx(int indexX, int indexY,int indexZ) const{
   cout << "ES field, no magnetic field. You dont need to call getBx(), and save B" << endl; 
   return(Ex[0][0][0]);
}
/** get Magnetic Field component X array*/ 
inline double*** ESfield1D::getBx() {
cout << "ES field, no magnetic field. You dont need to call getBx(), and save B" << endl;
return(NULL);}
/** get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
inline double &ESfield1D::getBy(int indexX, int indexY,int indexZ) const{
  cout << "ES field, no magnetic field. You dont need to call getBy(), and save B" << endl; 
  return(Ex[0][0][0]);}
/** get Magnetic Field Y component array*/
inline double*** ESfield1D::getBy() {
   cout << "ES field, no magnetic field. You dont need to call getBy(), and save B" << endl;
   return(NULL);}
/** get Magnetic Field Z component defined on node(indexX,indexY,indexZ) */
inline double &ESfield1D::getBz(int indexX, int indexY,int indexZ) const{
   cout << "ES field, no magnetic field. You dont need to call getBz(), and save B" << endl; 
   return(Ex[0][0][0]);}
inline double*** ESfield1D::getBz() {cout << "ES field, no magnetic field. You dont need to call getBz(), and save B" << endl;
   return(NULL);}
/** SPECIES: get pressure tensor component XX defined on nodes */
inline double**** ESfield1D::getpXXsn() {
  cout << "ES field, no magnetic pressure. You dont need to call getpXXsn()" << endl; 
  return(NULL);}
/** SPECIES: get pressure tensor component XY defined on nodes */
inline double**** ESfield1D::getpXYsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpXYsn()" << endl;
return(NULL);}
 /** SPECIES: get pressure tensor component XZ defined on nodes */
inline double**** ESfield1D::getpXZsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpXZsn()" << endl;
return(NULL);}
/** SPECIES: get pressure tensor component YY defined on nodes */
inline double**** ESfield1D::getpYYsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpYYsn()" << endl;
return(NULL);}
/** SPECIES: get pressure tensor component YZ defined on nodes */
inline double**** ESfield1D::getpYZsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpYZsn()" << endl;
return(NULL);}
/** SPECIES: get pressure tensor component ZZ defined on nodes */
inline double**** ESfield1D::getpZZsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpZZsn()" << endl;
return(NULL);}
/** get current -Direction X */
inline double &ESfield1D::getJx(int indexX, int indexY,int indexZ) const{
 cout << "ES field, no current calculation" << endl; 
 return(Ex[0][0][0]);}
/** get current array X component **/ 
inline double*** ESfield1D::getJx() {
 cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
 return(NULL);}
/** get current -Direction Y */
inline double &ESfield1D::getJy(int indexX, int indexY,int indexZ) const{
  cout << "ES field, no current calculation" << endl; 
  return(Ex[0][0][0]);
}
/** get current array Y component */
inline double*** ESfield1D::getJy() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/** get current -Direction Z */
inline double &ESfield1D::getJz(int indexX, int indexY,int indexZ) const{
  cout << "ES field, no current calculation" << endl; 
  return(Ex[0][0][0]);}
/** get current array Z component */
inline double*** ESfield1D::getJz() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/**SPECIES: get current array X component */
inline double**** ESfield1D::getJxs() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/**SPECIES: get current array Y component */
inline double**** ESfield1D::getJys() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/**SPECIES: get current array Z component */
inline double**** ESfield1D::getJzs() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}

/** add an amount of charge density to current density - direction X to current density   */
inline void ESfield1D::addJx(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addJx!" << endl;	
}
/** add an amount of current density - direction Y to current density   */
inline  void ESfield1D::addJy(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addJx!" << endl;	
}
/** add an amount of current density - direction Z to current density   */
inline  void ESfield1D::addJz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ::addJz!" << endl;	
}
/** add an amount of pressure density - direction XX to current density  */
inline  void ESfield1D::addPxx(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addPxx!" << endl;	
}
/** add an amount of pressure density - direction XY to current density  */
inline  void ESfield1D::addPxy(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addPxy!" << endl;	
}
/** add an amount of pressure density - direction XZ to current density  */
inline  void ESfield1D::addPxz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addPxz!" << endl;	
}
/** add an amount of pressure density - direction YY to current density */
inline  void ESfield1D::addPyy(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addPyy!" << endl;	
}
/** add an amount of pressure density - direction YZ to current density  */
inline  void ESfield1D::addPyz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addPyz!" << endl;	
}
/** add an amount of pressure density - direction ZZ to current density  */
inline  void ESfield1D::addPzz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield1D::addPzz!" << endl;	
}




/** Image of Maxwell for SOLVER */
inline void ESfield1D::MaxwellImage(double *im, double *vector, Grid *grid,VirtualTopology *vct){
   cout << "No need to call ESfield1D::MaxwellImage" << endl;
}     
/** Maxwell Source for SOLVER */
inline void ESfield1D::MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct){
  cout << "No need to call ESfield1D::MaxwellSource" << endl;
}
/** Print info about electromagnetic field */
inline void ESfield1D::print(void) const{
   double sum_rho = 0.0;
   for (int i=1; i < nxn-2; i++)
      sum_rho += rhon[i][0][0];
   double result=0.0;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(&sum_rho, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
   cout << "Sum of net charge density on nodes: " << result << endl ;
   sum_rho = 0.0;
   for (int i=1; i < nxc-1; i++)
      sum_rho += rhoc[i][0][0];
   result=0.0;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(&sum_rho, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   
   cout << "Sum of net charge density on centers: " << result << endl ;

}

#endif

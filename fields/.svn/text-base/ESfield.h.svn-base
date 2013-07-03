/*******************************************************************************************
ESfield.h  -  Electrostatic Field with 3 components(x,y,z) defined on a 2-D grid
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef ESfield_H
#define ESfield_H


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
* Electrostatic Field with 3 components(x,y,z) defined on a 2-D grid
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/

class ESfield : public Field{
  public:
  /** constructor */
  ESfield(CollectiveIO *col,Grid *grid);
  /** destructor */
  ~ESfield();

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
  void sumOverSpecies();
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
  /** light speed */
  double c;
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
  /** number of cell - Y direction, including + 2 (guard cells) */
  int nyc;
  /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
  int nyn;
  /** grid spacing */
  double dx, dy, invVOL;
  /** simulation box length - X direction   */
  double Lx;
  /** simulation box length - Y direction   */
  double Ly;
  /** PHI: ES potential (indexX, indexY, indexZ), defined on central points between nodes*/
  double ***PHI;
  /** Ex: electric field X-component (indexX, indexY, indexZ), defined on nodes */
  double ***Ex;
  /** Ey: electric field Y-component (indexX, indexY, indexZ), defined on nodes */
  double ***Ey;
  /** Ez: electric field Z-component (indexX, indexY, indexZ, #species), defined on nodes */
  double ***Ez;
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

};

  /** Calculate Electric field with the implicit solver: the Maxwell solver method is called here */
inline void ESfield::calculateField(Grid *grid, VirtualTopology *vct){
  if (vct->getCartesian_rank() ==0)
    cout << "*** PHI CALCULATION ***" << endl;
  double *xkrylovPoisson = new double[(nxc-2)*(nyc-2)];
  double *bkrylovPoisson = new double[(nxc-2)*(nyc-2)];
  double ***gradPHIX = newArr3(double,nxn,nyn,1);
  double ***gradPHIY = newArr3(double,nxn,nyn,1);
  eqValue (0.0, gradPHIX, nxn,nyn);
  eqValue (0.0, gradPHIY, nxn,nyn);
  // Calculate the laplacian(PHI) =   -4*PI*rho
  phys2solver(bkrylovPoisson,rhoc,nxc,nyc);    // bkrylov is place holder for source term
  scale(bkrylovPoisson,-FourPI,(nxc-2)*(nyc-2));       // - rho*4pi
  eqValue(0.0,xkrylovPoisson,(nxc-2)*(nyc-2)); // in the initial guess xkrylov is 0.0 -> xkrylov is the unknown
  //CG(xkrylovPoisson,(nxc-2)*(nyc-2),bkrylovPoisson, 1000, 10E-7, &Field::PoissonImage, grid, vct, this);
  GMRES(&Field::PoissonImage, xkrylovPoisson, (nxc-2)*(nyc-2),bkrylovPoisson,20,200,1E-7, grid, vct, this);
  solver2phys(PHI,xkrylovPoisson,nxc,nyc);
  communicateCenterBC(nxc,nyc,PHI,bcPHIfaceXright,bcPHIfaceXleft,bcPHIfaceYright,bcPHIfaceYleft,vct);
  // calculate the electric field on the nodes
  grid->gradC2N(Ex,Ey,PHI);
  neg(Ex,nxn,nyn);
  neg(Ey,nxn,nyn);
  communicateNode(nxn,nyn,Ex,vct);
  communicateNode(nxn,nyn,Ey,vct);
  // smoothing
  smooth(Smooth,Ex,1,grid,vct);
  smooth(Smooth,Ey,1,grid,vct);
  // deallocate 
  delete[] xkrylovPoisson;
  delete[] bkrylovPoisson;
  delArr3(gradPHIX,nxn,nyn);
  delArr3(gradPHIY,nxn,nyn);
}
/** Image of Poisson Solver: it's the argument for the GMRES */
inline void ESfield::PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct){
    // allocate  2 two dimensional service vectors
    double ***temp = newArr3(double,nxc,nyc,1);
    double ***im   = newArr3(double,nxc,nyc,1);
    // initialize to zero, new vectors
    eqValue (0.0, image,(nxc-2)*(nyc-2)); // added by Enrico.
    eqValue (0.0, temp,nxc,nyc);
    eqValue (0.0, im,nxc,nyc);
    // move from krylov space to physical space and communicate ghost cells
    solver2phys(temp,vector,nxc,nyc);
    MPI_Barrier(MPI_COMM_WORLD);
    communicateCenterBC(nxc,nyc,temp,bcPHIfaceXright,bcPHIfaceXleft,bcPHIfaceYright,bcPHIfaceYleft,vct);
    // calculate the laplacian
    grid->lapC2C(im,temp,vct);
    // move from physical space to krylov space, in this case ghost cells don't count
    phys2solver(image,im,nxc,nyc);
    // deallocate temporary array and objects
    delArr3(temp,nxc,nyc);
    delArr3(im,nxc,nyc);
}
/** set to 0 all the densities fields: this is need before the interpolation */
inline  void ESfield::setZeroDensities(){
for (int i=0; i < nxn; i++)
	for (int j=0; j < nyn; j++){
	  rhon[i][j][0] = 0.0;
}
for (int i=0; i < nxc; i++)
	for (int j=0; j < nyc; j++){
	 rhoc[i][j][0] = 0.0;
	 
        }
for (int kk=0; kk < ns; kk++)
    for (int i=0; i < nxn; i++)
	for (int j=0; j < nyn; j++){
		  rhons[kk][i][j][0] = 0.0;
	}     

}
/** Sum rhon over species*/
inline void ESfield::sumOverSpecies(){
	for (int is=0; is<ns; is++){
 	 for (register int i=0; i <nxn; i++)
  	  for (register int j=0; j <nyn; j++)
	   rhon[i][j][0]     +=rhons[is][i][j][0];}
}
/** interpolate rho from node to center */
inline  void ESfield::interpDensitiesN2C(VirtualTopology *vct, Grid *grid){
 grid->interpN2C(rhoc,rhon);
}
/** communicate ghost cells, for grid->Particles interpolation */
inline void ESfield::communicateGhostP2G(int ns,int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct){
	// communicate interpolation + BC for interpolated fields
	communicateInterp(nxn,nyn,ns,rhons,0,0,0,0,dx,dy,vct);
	communicateNode(nxn, nyn, rhons, ns,vct);
}   
/** add an amount of charge density to charge density field at node X,Y */
inline void ESfield::addRho(double ***weight, int X, int Y, int Z, int ns){
  	rhons[X-1][Y-1][0][ns]+= weight[0][0][0]*invVOL;
  	rhons[X-1][Y][0][ns]  += weight[0][1][0]*invVOL;
	rhons[X][Y-1][0][ns]  += weight[1][0][0]*invVOL;
	rhons[X][Y][0][ns]    += weight[1][1][0]*invVOL;
}
/** 
Interpolation smoothing:
Smoothing (vector must already have ghost cells)
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector ; 
type = 1 --> node based vector   ;                     
*/
inline void  ESfield::smooth(double value,double ***vector, bool type, Grid *grid, VirtualTopology *vct){
 double alpha;
 int nx, ny;
 if (value != 1.0){
   switch(type){ 
    case (0):
	nx=grid->getNXC();
	ny=grid->getNYC();
	break;
     case(1):
	nx=grid->getNXN();
	ny=grid->getNYN();
	break;
     }

   double **temp = newArr(double,nx,ny);
   alpha=(1.0-value)/4;

   for (int i=1; i<nx-1; i++)
     for (int j=1; j<ny-1; j++)
	temp[i][j]=value*vector[i][j][0]+alpha*(vector[i-1][j][0]+vector[i][j-1][0]+vector[i+1][j][0]+vector[i][j+1][0]);
 // Adjust boundary
 if (vct->getXleft_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int j=1; j<ny-1;j++)
	temp[1][j]=value*vector[1][j][0]+alpha*(vector[1][j-1][0]+vector[2][j][0]+vector[1][j+1][0]);}
	
if (vct->getXright_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int j=1; j<ny-1;j++)
	temp[nx-2][j]=value*vector[nx-2][j][0]+alpha*(vector[nx-2][j-1][0]+vector[nx-3][j][0]+vector[nx-2][j+1][0]);}
	
if (vct->getYleft_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int i=1; i<nx-1;i++)
	temp[i][1]=value*vector[i][1][0]+alpha*(vector[i-1][1][0]+vector[i+1][1][0]+vector[i][2][0]);}
	
if (vct->getYright_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int i=1; i<nx-1;i++)
	temp[i][ny-2]=value*vector[i][ny-2][0]+alpha*(vector[i-1][ny-2][0]+vector[i][ny-3][0]+vector[i+1][ny-2][0]);}
	
if (vct->getXleft_neighbor()==-1 && vct->getYleft_neighbor()==-1){
  alpha=(1.0-value)/2;
  temp[1][1]=value*vector[1][1][0]+alpha*(vector[2][1][0]+vector[1][2][0]);
}
if (vct->getXleft_neighbor()==-1 && vct->getYright_neighbor()==-1){
  alpha=(1.0-value)/2;
  temp[1][ny-2]=value*vector[1][ny-2][0]+alpha*(vector[2][ny-1][0]+vector[1][ny-3][0]);
}
if (vct->getXright_neighbor()==-1 && vct->getYleft_neighbor()==-1){
  alpha=(1.0-value)/2; 
  temp[nx-2][1]=value*vector[nx-2][1][0]+alpha*(vector[nx-3][1][0]+vector[nx-2][2][0]);
}
if (vct->getXright_neighbor()==-1 && vct->getYright_neighbor()==-1){
  alpha=(1.0-value)/2;
  temp[nx-2][ny-2]=value*vector[nx-2][ny-2][0]+alpha*(vector[nx-3][ny-2][0]+vector[nx-2][ny-3][0]);
}

 //eq(vector,temp,nx,ny);
 // copy temp in vector
 for (int i=0; i < nx;i++)
   for (int j=0; j < ny;j++)
     vector[i][j][0] = temp[i][j];
 // communicate
  if (type == 0)
  communicateCenter(nx, ny, vector, vct);
 else 
  communicateNode(nx, ny, vector, vct);

  delArr(temp,nx);
 } // end of if (value !=1)
}

/** 

SPECIES: Interpolation smoothing
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector  
type = 1 --> node based vector                        
*/
inline void  ESfield::smooth(double value,double ****vector,int is, bool type, Grid *grid, VirtualTopology *vct){
double alpha;
int nx, ny;
if (value != 1.0){
switch(type){ 
case (0):
	nx=grid->getNXC();
	ny=grid->getNYC();
	break;
case(1):
	nx=grid->getNXN();
	ny=grid->getNYN();
	break;
}

double **temp = newArr(double,nx,ny);
alpha=(1.0-value)/4;

for (int i=1; i<nx-1; i++)
 for (int j=1; j<ny-1; j++)
	temp[i][j]=value*vector[is][i][j][0]+alpha*(vector[is][i-1][j][0]+vector[is][i][j-1][0]+vector[is][i+1][j][0]+vector[is][i][j+1][0]);
// Adjust boundary
if (vct->getXleft_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int j=1; j<ny-1;j++)
	temp[1][j]=value*vector[is][1][j][0]+alpha*(vector[is][1][j-1][0]+vector[is][2][j][0]+vector[is][1][j+1][0]);}
	
if (vct->getXright_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int j=1; j<ny-1;j++)
	temp[nx-2][j]=value*vector[is][nx-2][j][0]+alpha*(vector[is][nx-2][j-1][0]+vector[is][nx-3][j][0]+vector[is][nx-2][j+1][0]);}
	
if (vct->getYleft_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int i=1; i<nx-1;i++)
	temp[i][1]=value*vector[is][i][1][0]+alpha*(vector[is][i-1][1][0]+vector[is][i+1][1][0]+vector[is][i][2][0]);}
	
if (vct->getYright_neighbor()==-1){
	alpha=(1.0-value)/3;
	for (int i=1; i<nx-1;i++)
	temp[i][ny-2]=value*vector[is][i][ny-2][0]+alpha*(vector[is][i-1][ny-2][0]+vector[is][i][ny-3][0]+vector[is][i+1][ny-2][0]);}
	
if (vct->getXleft_neighbor()==-1 && vct->getYleft_neighbor()==-1){
   alpha=(1.0-value)/2;
   temp[1][1]=value*vector[is][1][1][0]+alpha*(vector[is][2][1][0]+vector[is][1][2][0]);}
if (vct->getXleft_neighbor()==-1 && vct->getYright_neighbor()==-1){
alpha=(1.0-value)/2;
temp[1][ny-2]=value*vector[is][1][ny-2][0]+alpha*(vector[is][2][ny-1][0]+vector[is][1][ny-3][0]);}
if (vct->getXright_neighbor()==-1 && vct->getYleft_neighbor()==-1){
alpha=(1.0-value)/2;
temp[nx-2][1]=value*vector[is][nx-2][1][0]+alpha*(vector[is][nx-3][1][0]+vector[is][nx-2][2][0]);}
if (vct->getXright_neighbor()==-1 && vct->getYright_neighbor()==-1){
alpha=(1.0-value)/2;
temp[nx-2][ny-2]=value*vector[is][nx-2][ny-2][0]+alpha*(vector[is][nx-3][ny-2][0]+vector[is][nx-2][ny-3][0]);}

//eq(vector,temp,nx,ny,is);
// copy temp in vector
 for (int i=0; i < nx;i++)
   for (int j=0; j < ny;j++)
     vector[is][i][j][0] = temp[i][j];
 // communicate
 if (type == 0)
   communicateCenter(nx, ny, vector, is, vct);
 else 
  communicateNode(nx, ny, vector, is, vct);

 delArr(temp,nx);
 } // end of if
}
/** Initialize the EM field with constants values or from restart*/
inline void ESfield::init(VirtualTopology *vct, Grid *grid){
// initialize E and rhos on nodes
if (restart1==0){  
  for (int i=0; i < nxn; i++)
	for (int j=0; j <nyn; j++){
	   Ex[i][j][0] =  0.0;
	   Ey[i][j][0] =  0.0;
	   Ez[i][j][0] =  0.0;
	   rhons[0][i][j][0] =  0.07957747154595;  // = ion rhons, for plasma neutrality
	   rhons[1][i][j][0] =  0.07957747154595;  // ion plasma frequency = 1 (normalization)
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
	 if  (nxn != (dims_out[0] + 2) || nyn != (dims_out[1] + 2)){
	    cout << "RESTART GRID NOT CONSISTENT!!!!!!!!!!!!!!!!!!" << endl;
	 }
	 // read Exn
	double *temp_storage = new double[dims_out[0]*dims_out[1]];
	dataset_id = H5Dopen(file_id,"/fields/Ex/cycle_0");
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	int k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
	        Ex[i][j][0] = temp_storage[k++]; 
	communicateNode(nxn, nyn, Ex,vct);
	status = H5Dclose(dataset_id);
	// read Eyn
	dataset_id = H5Dopen(file_id,"/fields/Ey/cycle_0");
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
	        Ey[i][j][0] = temp_storage[k++]; 
	communicateNode(nxn, nyn, Ey,vct);
	status = H5Dclose(dataset_id);
        // read Ezn 
	dataset_id = H5Dopen(file_id,"/fields/Ez/cycle_0");
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
	H5S_ALL,H5P_DEFAULT,temp_storage);
	k=0;
	for (int i=1; i < nxn-1; i++)
	      for (int j=1; j <nyn-1; j++)
	        Ez[i][j][0] = temp_storage[k++]; 
	communicateNode(nxn, nyn, Ez,vct);
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
	      for (int j=1; j <nyn-1; j++)
	        rhons[is][i][j][0] = temp_storage[k++]; 
	   communicateNode(nxn, nyn, rhons,is,vct);
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
inline double*** ESfield::getPHI() {return(PHI);}
/** get Electric Field component X defined on node(indexX,indexY,indexZ) */  
inline double &ESfield::getEx(int indexX, int indexY,  int indexZ) const{
  return(Ex[indexX][indexY][0]);}
/** get Electric Field X component array */
inline double*** ESfield::getEx(){return(Ex);}
/** get Electric Field component Y defined on node(indexX,indexY,indexZ) */  
inline double &ESfield::getEy(int indexX, int indexY,int indexZ) const{
  return(Ey[indexX][indexY][0]);}
/** get Electric Field Y component array */
inline double*** ESfield::getEy() {return(Ey);}
/** get Electric Field component Z defined on node(indexX,indexY,indexZ) */  
inline double &ESfield::getEz(int indexX, int indexY,int indexZ) const{
  return(Ez[indexX][indexY][0]);}
/** get Electric Field Z component array */
inline double*** ESfield::getEz() {return(Ez);}
/** get density on center cell(indexX,indexY,indexZ)  */
inline double &ESfield::getRHOc(int indexX, int indexY,int indexZ) const{
  return(rhoc[indexX][indexY][0]);}
 /** get density array defined on centers cells */
inline double*** ESfield::getRHOc() {return(rhoc);}
/** get density on node(indexX,indexY,indexZ)  */
inline double &ESfield::getRHOn(int indexX, int indexY,int indexZ) const{ 
return(rhon[indexX][indexY][0]);}
/** get density array defined on nodes*/
inline double*** ESfield::getRHOn() {return(rhon);}
/** SPECIES: get density array defined on nodes  */
inline double &ESfield::getRHOns(int indexX, int indexY, int indexZ,int ns) const{
  return(rhons[ns][indexX][indexY][0]);}
/** SPECIES: get density array defined on center cells  */
inline double &ESfield::getRHOcs(int indexX, int indexY,int indexZ, int ns) const{
  return(rhocs[ns][indexX][indexY][0]);}
/** SPECIES: get density array defined on nodes  */
inline double**** ESfield::getRHOns() {return(rhons);}



/** Sum current over different species */
inline void ESfield::sumOverSpeciesJ(){
   cout << "ES field, no need to call ESfield::sumOverSpeciesJ() " << endl;

}
/** get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
inline double &ESfield::getBx(int indexX, int indexY,int indexZ) const{
   cout << "ES field, no magnetic field. You dont need to call getBx(), and save B" << endl; 
   return(Ex[0][0][0]);
}
/** get Magnetic Field component X array*/ 
inline double*** ESfield::getBx() {
cout << "ES field, no magnetic field. You dont need to call getBx(), and save B" << endl;
return(NULL);}
/** get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
inline double &ESfield::getBy(int indexX, int indexY,int indexZ) const{
  cout << "ES field, no magnetic field. You dont need to call getBy(), and save B" << endl; 
  return(Ex[0][0][0]);}
/** get Magnetic Field Y component array*/
inline double*** ESfield::getBy() {
   cout << "ES field, no magnetic field. You dont need to call getBy(), and save B" << endl;
   return(NULL);}
/** get Magnetic Field Z component defined on node(indexX,indexY,indexZ) */
inline double &ESfield::getBz(int indexX, int indexY,int indexZ) const{
   cout << "ES field, no magnetic field. You dont need to call getBz(), and save B" << endl; 
   return(Ex[0][0][0]);}
inline double*** ESfield::getBz() {cout << "ES field, no magnetic field. You dont need to call getBz(), and save B" << endl;
   return(NULL);}
/** SPECIES: get pressure tensor component XX defined on nodes */
inline double**** ESfield::getpXXsn() {
  cout << "ES field, no magnetic pressure. You dont need to call getpXXsn()" << endl; 
  return(NULL);}
/** SPECIES: get pressure tensor component XY defined on nodes */
inline double**** ESfield::getpXYsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpXYsn()" << endl;
return(NULL);}
 /** SPECIES: get pressure tensor component XZ defined on nodes */
inline double**** ESfield::getpXZsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpXZsn()" << endl;
return(NULL);}
/** SPECIES: get pressure tensor component YY defined on nodes */
inline double**** ESfield::getpYYsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpYYsn()" << endl;
return(NULL);}
/** SPECIES: get pressure tensor component YZ defined on nodes */
inline double**** ESfield::getpYZsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpYZsn()" << endl;
return(NULL);}
/** SPECIES: get pressure tensor component ZZ defined on nodes */
inline double**** ESfield::getpZZsn() {
cout << "ES field, no magnetic pressure. You dont need to call getpZZsn()" << endl;
return(NULL);}
/** get current -Direction X */
inline double &ESfield::getJx(int indexX, int indexY,int indexZ) const{
 cout << "ES field, no current calculation" << endl; 
 return(Ex[0][0][0]);}
/** get current array X component **/ 
inline double*** ESfield::getJx() {
 cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
 return(NULL);}
/** get current -Direction Y */
inline double &ESfield::getJy(int indexX, int indexY,int indexZ) const{
  cout << "ES field, no current calculation" << endl; 
  return(Ex[0][0][0]);
}
/** get current array Y component */
inline double*** ESfield::getJy() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/** get current -Direction Z */
inline double &ESfield::getJz(int indexX, int indexY,int indexZ) const{
  cout << "ES field, no current calculation" << endl; 
  return(Ex[0][0][0]);}
/** get current array Z component */
inline double*** ESfield::getJz() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/**SPECIES: get current array X component */
inline double**** ESfield::getJxs() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/**SPECIES: get current array Y component */
inline double**** ESfield::getJys() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}
/**SPECIES: get current array Z component */
inline double**** ESfield::getJzs() {
cout << "ES field, you are trying to save current density, that is not being calculated" << endl;
return(NULL);}

/** add an amount of charge density to current density - direction X to current density   */
inline void ESfield::addJx(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addJx!" << endl;	
}
/** add an amount of current density - direction Y to current density   */
inline  void ESfield::addJy(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addJx!" << endl;	
}
/** add an amount of current density - direction Z to current density   */
inline  void ESfield::addJz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ::addJz!" << endl;	
}
/** add an amount of pressure density - direction XX to current density  */
inline  void ESfield::addPxx(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addPxx!" << endl;	
}
/** add an amount of pressure density - direction XY to current density  */
inline  void ESfield::addPxy(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addPxy!" << endl;	
}
/** add an amount of pressure density - direction XZ to current density  */
inline  void ESfield::addPxz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addPxz!" << endl;	
}
/** add an amount of pressure density - direction YY to current density */
inline  void ESfield::addPyy(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addPyy!" << endl;	
}
/** add an amount of pressure density - direction YZ to current density  */
inline  void ESfield::addPyz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addPyz!" << endl;	
}
/** add an amount of pressure density - direction ZZ to current density  */
inline  void ESfield::addPzz(double ***weight, int X, int Y, int Z, int ns){
      cout << "No need to call ESfield::addPzz!" << endl;	
}




/** Image of Maxwell for SOLVER */
inline void ESfield::MaxwellImage(double *im, double *vector, Grid *grid,VirtualTopology *vct){
   cout << "No need to call ESfield::MaxwellImage" << endl;
}     
/** Maxwell Source for SOLVER */
inline void ESfield::MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct){
  cout << "No need to call ESfield::MaxwellSource" << endl;
}
/** Print info about electromagnetic field */
inline void ESfield::print(void) const{


}
/** constructor */
inline ESfield::ESfield(CollectiveIO *col,Grid *grid){
	nxc = grid->getNXC();
	nxn = grid->getNXN();
	nyc = grid->getNYC();
	nyn = grid->getNYN();
	dx = grid->getDX();
	dy = grid ->getDY();
	invVOL = grid->getInvVOL();
	Lx = col->getLx();
	Ly = col->getLy();
	ns  = col->getNs();
	c = col->getC();
	dt = col->getDt();
	Smooth = col->getSmooth();
	th = col->getTh();
	delt = c*th*dt;
	FourPI =16*atan(1.0);
	qom = new double[ns];
	for (int i=0; i < ns;i++)
		qom[i] = col->getQOM(i);
        // boundary conditions: PHI and EM fields
	bcPHIfaceXright  = col->getBcPHIfaceXright();
	bcPHIfaceXleft   = col->getBcPHIfaceXleft();
	bcPHIfaceYright  = col->getBcPHIfaceYright();
	bcPHIfaceYleft   = col->getBcPHIfaceYleft();
        // Restart 
	restart1 = col->getRestart_status();
        RestartDirName = col->getRestartDirName();
        ////////////////// ARRAY ALLOCATION ////////////////
        // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
	Ex = newArr3(double,nxn,nyn,1);
	Ey = newArr3(double,nxn,nyn,1);
	Ez = newArr3(double,nxn,nyn,1);
	rhon  = newArr3(double,nxn,nyn,1);
	// involving species
	rhons = newArr4(double,ns,nxn,nyn,1);
	rhocs = newArr4(double,ns,nxc,nyc,1);
	// arrays allocation: central points ---> the first central point has index 1, the last has index nxn-2!
	PHI    = newArr3(double,nxc,nyc,1);
	rhoc   = newArr3(double,nxc,nyc,1);

}
/** destructor: deallocate arrays*/
inline ESfield::~ESfield(){
        // nodes
	delArr3(Ex,nxn,nyn);
	delArr3(Ey,nxn,nyn);
	delArr3(Ez,nxn,nyn);	
	delArr3(rhon,nxn,nyn);
        // nodes and species
	delArr4(rhons,ns,nxn,nyn);
	// central points
	delArr3(PHI,nxc,nyc);
	delArr3(rhoc,nxc,nxc);
	// center and species 
	delArr4(rhocs,ns,nxc,nyc);
}
#endif

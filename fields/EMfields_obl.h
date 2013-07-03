/*******************************************************************************************
EMfields.h  -  Electromagnetic Field with 3 components(x,y,z) defined on a 2-D grid. Solved
using the implicit Maxwell solver.
-------------------
developers: Stefano Markidis,  Giovanni Lapenta
	********************************************************************************************/

#ifndef EMfields_H
#define EMfields_H


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
*  Electromagnetic fields and sources defined for each local grid, and for an implicit Maxwell's solver
 *
 * @date Fri Jun 4 2007 KUL
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */


class EMfields : public Field{
public:
	/** constructor */
	EMfields(CollectiveIO *col,Grid *grid);
	/** destructor */
	~EMfields();

	/** initialize the electromagnetic fields with constant values */
	void init(VirtualTopology *vct, Grid *grid);
	/** initialize with reconnection system */
	void initGEM(VirtualTopology *vct, Grid *grid);
	/** initialize with uniform distribution*/
	void initUniform(VirtualTopology *vct, Grid *grid);
	/** initialize LHDI */
	void initLHDI(VirtualTopology *vct, Grid *grid);
	/** initialize with Alfven wave*/
	void initAlfven(VirtualTopology *vct, Grid *grid);
	/** initialized with rotated magnetic field */
	void initEM_rotate(VirtualTopology *vct, Grid *grid, double B, double theta);
	/** add magnetic dipole in position (x0,y0) with magnetic moment mo in direction Y*/
	void addDipole(double B0, double x0, double y0,Grid *grid);
	/** add IMF field*/
	void addIMF(double B_IMF, double theta, Grid *grid);
	/** add a perturbattion to charge density */
	void AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid *grid);
	/** add a perturbattion to the EM field */
	void AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid *grid,VirtualTopology *vct);


	/** Calculate Electric field using the implicit Maxwell solver */
	void calculateField(Grid *grid, VirtualTopology *vct);
	/** Image of Poisson Solver (for SOLVER)*/
	void PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct);
	/** Image of Maxwell Solver (for Solver) */
	void MaxwellImage(double *im, double *vector, Grid *grid,VirtualTopology *vct);
	/** Maxwell source term (for SOLVER) */
	void MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct);
	/** Calculate Magnetic field with the implicit solver: calculate B defined on nodes
		With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
	void calculateB(Grid *grid, VirtualTopology *vct);

	/** Calculate the three components of Pi(implicit pressure) cross image vector */
	void PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid *grid);
	/** Calculate the three components of mu (implicit permeattivity) cross image vector */
	void MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid *grid);
	/* Calculate susceptibility tensor on the boundaries */
	void sustensorLeft(double *susxx, double *susxy, double *susxz, double* susyx, double* susyy, double* susyz, double* suszx, double* suszy, double* suszz,  int dir);
	/** calculate the suceptibility tensor on right boundary */
	void sustensorRight(double *susxx, double *susxy, double *susxz, double* susyx, double* susyy, double* susyz, double* suszx, double* suszy, double* suszz,  int dir);
	/** Calculate rho hat, Jx hat, Jy hat, Jz hat */
	void calculateHatFunctions(Grid *grid, VirtualTopology *vct);


	/** communicate ghost for densities and interp rho from node to center */
	void interpDensitiesN2C(VirtualTopology *vct, Grid *grid);
	/** set to 0 all the densities fields */
	void setZeroDensities();
	/** Sum rhon over species */
	void sumOverSpecies(VirtualTopology *vct);
	/** Sum current over different species */
	void sumOverSpeciesJ();
	/** Smoothing after the interpolation**/
	void smooth(int nvolte, double value ,double ***vector,bool type, Grid *grid, VirtualTopology *vct);
	/** Smoothing for the electric field**/
	void smoothE(int nvolte, double value ,double ***vector,bool type, Grid *grid, VirtualTopology *vct);
	/** SPECIES: Smoothing after the interpolation for species fields**/
	void smooth(int nvolte, double value ,double ****vector,int is, bool type, Grid *grid, VirtualTopology *vct);

	/** communicate ghost for grid -> Particles interpolation */
	void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct);
    /** add an amount of charge density to charge density field at node X,Y,Z */
	void addRho(double ***weight, int X, int Y,int Z, int is);
	/** add an amount of current density - direction X to current density field at node X,Y,Z */
	void addJx(double ***weight, int X, int Y, int Z,int is);
	/** add an amount of current density - direction Y to current density field at node X,Y,Z */
	void addJy(double ***weight, int X, int Y, int Z, int is);
	/** add an amount of current density - direction Z to current density field at node X,Y,Z */
	void addJz(double ***weight, int X, int Y, int Z, int is);

	/** add an amount of pressure density - direction XX to current density field at node X,Y,Z */
	void addPxx(double ***weight, int X, int Y, int Z, int is);
	/** add an amount of pressure density - direction XY to current density field at node X,Y,Z */
	void addPxy(double ***weight, int X, int Y, int Z,int is);
	/** add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
	void addPxz(double ***weight, int X, int Y, int Z,int is);
	/** add an amount of pressure density - direction YY to current density field at node X,Y,Z */
	void addPyy(double ***weight, int X, int Y, int Z, int is);
	/** add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
	void addPyz(double ***weight, int X, int Y, int Z, int is);
	/** add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
	void addPzz(double ***weight, int X, int Y, int Z,int is);

	/** adjust densities on boundaries that are not periodic */
	void adjustNonPeriodicDensities(VirtualTopology *vct);

	/** Perfect conductor boundary conditions LEFT wall */
	void perfectConductorLeftImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Perfect conductor boundary conditions RIGHT wall */
	void perfectConductorRightImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid);
	/** Perfect conductor boundary conditions for source LEFT wall*/
    void perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Perfect conductor boundary conditions for source RIGHT wall*/
    void perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Magnetic mirror boundary conditions for source LEFT wall */
    void magneticMirrorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Magnetic mirror boundary conditions for source RIGHT wall */
    void magneticMirrorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Magnetic mirror boundary conditions for source LEFT wall */
    void openLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Magnetic mirror boundary conditions for source RIGHT wall */
    void openRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);


	/** Perfect conductor boundary conditions LEFT wall */
	void perfectConductorLeft(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Perfect conductor boundary conditions RIGHT wall */
	void perfectConductorRight(double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid);
	/** Perfect conductor boundary conditions LEFT wall */
	void BperfectConductorLeft(int dir,Grid *grid,VirtualTopology *vct);
	/** Perfect conductor boundary conditions RIGHT wall */
	void BperfectConductorRight(int dir,Grid *grid,VirtualTopology *vct);
	/** Magnetic mirror boundary conditions LEFT wall*/
	void magneticMirrorLeft(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Magnetic mirror boundary conditions RIGHT wall*/
	void magneticMirrorRight(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Open Boundary BC LEFT wall*/
	void openLeft(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Opne Boundary BC RIGHT wall*/
	void openRight(double ***vectorX,double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Open Boundary condition with no plasma injection LEFT wall */
	void openNoPlasmaLeft(double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid);
	/** Open Boundary condition with no plasma injection RIGHT wall */
	void openNoPlasmaRight(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** BC for PHI: divergence cleaning */
	void bcPHI_Left(double ***im,int dir);
	/** BC for PHI: divergence cleaning */
	void bcPHI_Right(double ***im,int dir);


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
	double &getRHOns(int indexX, int indexY,int indexZ,int is) const;
	/** SPECIES: get density on center cell(indexX,indexY,indexZ)*/
	double &getRHOcs(int indexX, int indexY,int indexZ,int is) const;
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
	/* 4*PI for normalization */
	double FourPI;
	/** time step */
	double dt;
	/** decentering parameter */
	double th;
	/** Number of smoothing steps */
	int Nvolte;
	/** Smoothing value*/
	double Smooth;
	/** delt = c*th*dt */
	double delt;
	/** number of particles species */
	int ns;
	/** GEM challenge parameters */
	double B0x, B0y, B0z, delta;
	/** charge to mass ratio array for different species */
	double *qom;


	// KEEP IN MEMORY GUARD CELLS ARE INCLUDED
	/** number of cells - X direction, including + 2 (guard cells) */
	int nxc;
	/** number of nodes - X direction, including + 2 extra nodes for guard cells */
	int nxn;
	/** number of cell - Y direction, including + 2 (guard cells) */
	int nyc;
	/** number of nodes - Y direction, including + 2 extra nodes for guard cells */
	int nyn;
	/** local grid boundaries coordinate  */
	double xStart, xEnd, yStart, yEnd;
	/** grid spacing dx*/
	double dx;
	/** grid spacing dy*/
	double dy;
	/** inverse of the cell volume */
	double invVOL;
	/** simulation box length - X direction   */
	double Lx;
	/** simulation box length - Y direction   */
	double Ly;

	/** Electric potential, defined on central points of the cell*/
	double***  PHI;
	/** Electric field X-component, defined on nodes */
	double***  Ex;
	/** Implicit electric field X-component, defined on nodes */
	double***  Exth;
	/** Electric field Y-component, defined on nodes */
	double***  Ey;
	/** Implicit electric field Y-component, defined on nodes */
	double***  Eyth;
	/** Electric field Z-component, defined on nodes */
	double***  Ez;
	/** Implicit electric field Z-component, defined on nodes */
	double***  Ezth;
	/** Magnetic field X-component, defined on central points of the cell*/
	double***  Bxc;
	/** Magnetic field Y-component, defined on central points of the cell*/
	double***  Byc;
	/** Magnetic field Z-component, defined on central points of the cell*/
	double***  Bzc;
	/** Magnetic field X-component, defined on nodes*/
	double***  Bxn;
	/** Magnetic field Y-component, defined on nodes*/
	double***  Byn;
	/** Magnetic field Z-component, defined on nodes*/
	double***  Bzn;


	/**some temporary arrays (for calculate hat functions)*/
	double*** tempXC;
	double*** tempYC;
	double*** tempZC;
	double*** tempXN;
	double*** tempYN;
	double*** tempZN;
	/** other temporary arrays (in MaxwellSource) */
	double*** tempC;
	double*** tempX;
	double*** tempY;
	double*** tempZ;
	double*** temp2X;
	double*** temp2Y;
	double*** temp2Z;
	/** and some for MaxwellImage */
	double*** imageX;
	double*** imageY;
	double*** imageZ;
	double*** Dx;
	double*** Dy;
	double*** Dz;
	double*** vectX;
	double*** vectY;
	double*** vectZ;
	double*** divC;

	/*********************************************************************************
		/*************                SOURCES                                           **
		/********************************************************************************/

		/** Charge density, defined on central points of the cell */
		double*** rhoc;
		/** Charge density, defined on nodes */
		double*** rhon;
		/** Implicit charge density, defined on central points of the cell */
		double***  rhoh;
		/** SPECIES: charge density for each species, defined on nodes */
		double****  rhons;
		/** SPECIES: charge density for each species, defined on central points of the cell */
		double**** rhocs;
		/** Current density component-X, defined on nodes */
		double***  Jx;
		/** Current density component-Y, defined on nodes */
		double***  Jy;
		/** Current density component-Z, defined on nodes */
		double***  Jz;
		/** Implicit current density X-component, defined on nodes */
		double***  Jxh;
		/** Implicit current density Y-component, defined on nodes */
		double***  Jyh;
		/** Implicit current density Z-component, defined on nodes */
		double***  Jzh;
		/** SPECIES: current density component-X for species, defined on nodes */
		double**** Jxs;
		/** SPECIES: current density component-Y for species, defined on nodes */
		double****  Jys;
		/** SPECIES: current density component-Z for species, defined on nodes */
		double****  Jzs;

		/** SPECIES: pressure tensor component-XX, defined on nodes */
		double**** pXXsn;
		/** SPECIES: pressure tensor component-XY, defined on nodes */
		double ****pXYsn;
		/** SPECIES: pressure tensor component-XZ, defined on nodes */
		double ****pXZsn;
		/** SPECIES: pressure tensor component-XZ, defined on nodes */
		double ****pYYsn;
		/** SPECIES: pressure tensor component-YZ, defined on nodes */
		double ****pYZsn;
		/** SPECIES: pressure tensor component-ZZ, defined on nodes */
		double ****pZZsn;


		/** Field Boundary Condition
			0 =  Dirichlet Boundary Condition: specifies the value to take on the boundary of the domain
			1 =  Neumann Boundary Condition: specifies the value of derivative to take on the boundary of the domain
			2 =  Periodic condition
			*/

		/** Boundary Condition Electrostatic Potential: X right wall */
		int bcPHIfaceXright;
		/** Boundary Condition Electrostatic Potential: X left wall */
		int bcPHIfaceXleft;
		/** Boundary Condition Electrostatic Potential: Y right wall */
		int bcPHIfaceYright;
		/** Boundary Condition Electrostatic Potential: Y left wall */
		int bcPHIfaceYleft;

		/** Boundary condition for electric field
			0 = perfect conductor
			1 = magnetic mirror
			*/
		/** Boundary Condition EM Field: FaceXright */
		int bcEMfaceXright;
		/** Boundary Condition EM Field: FaceXleft */
		int bcEMfaceXleft;
		/** Boundary Condition EM Field: FaceYright */
		int bcEMfaceYright;
		/** Boundary Condition EM Field: FaceYleft */
		int bcEMfaceYleft;


		/** boolean for divergence cleaning */
		bool PoissonCorrection;

		/** GEM Challenge background ion */
		double *rhoINIT;
		/** Drift of the species */
		bool *DriftSpecies;
		/** RESTART BOOLEAN */
		int restart1 ;
		/** String with the directory for the restart file */
		string RestartDirName;

		/** CG tolerance criterium for stopping iterations */
		double CGtol;
		/**  GMRES tolerance criterium for stopping iterations */
		double GMREStol;

		/** Injection Velocity from the wall */
		double Vinj;



};


/** Calculate Electric field with the implicit Maxwell solver */
inline void EMfields::calculateField(Grid *grid, VirtualTopology *vct){
	if (vct->getCartesian_rank() ==0)
		cout << "*** E, B CALCULATION ***" << endl;
	double ***divE = newArr3(double,nxc,nyc,1);
	double ***gradPHIX = newArr3(double,nxn,nyn,1);
	double ***gradPHIY = newArr3(double,nxn,nyn,1);
	double ***gradPHIZ = newArr3(double,nxn,nyn,1);

	double *xkrylov = new double[3*(nxn-2)*(nyn-2)];
	double *bkrylov = new double[3*(nxn-2)*(nyn-2)];

	double *xkrylovPoisson = new double[(nxc-2)*(nyc-2)];
	double *bkrylovPoisson = new double[(nxc-2)*(nyc-2)];

	eqValue (0.0, xkrylov, 3*(nxn-2)*(nyn-2));
	eqValue (0.0, bkrylov, 3*(nxn-2)*(nyn-2));
	eqValue (0.0, divE,nxc,nyc);
	eqValue (0.0, tempC,nxc,nyc);
	eqValue (0.0, gradPHIX,nxn,nyn);
	eqValue (0.0, gradPHIY,nxn,nyn);
	// Adjust E calculating laplacian(PHI) = div(E)  -4*PI*rho
	if (PoissonCorrection){
		if (vct->getCartesian_rank() ==0)
			cout << "*** POISSON CORRECTION ***" << endl;
		grid->divN2C(divE,Ex,Ey);
		scale(tempC,rhoc,-FourPI,nxc,nyc);
		sum(divE,tempC,nxc,nyc);
		phys2solver(bkrylovPoisson,divE,nxc,nyc); // bkrylovPoisson is defined only in active centers
		eqValue(0.0,xkrylovPoisson,(nxc-2)*(nyc-2));
		if (vct->getPERIODICX() && vct->getPERIODICY())
			GMRES(&Field::PoissonImage, xkrylovPoisson, (nxc-2)*(nyc-2),bkrylovPoisson,20,200,GMREStol, grid, vct, this);
		else
			CG(xkrylovPoisson,(nxc-2)*(nyc-2),bkrylovPoisson, 10000, CGtol, &Field::PoissonImage, grid, vct, this);
		solver2phys(PHI,xkrylovPoisson,nxc,nyc);
		communicateCenter(nxc,nyc,PHI,vct);
		// adjust boundaries before calculating the lap
		if(vct->getXleft_neighbor()==-1 && (bcPHIfaceXleft ==0 ||  bcPHIfaceXleft ==1)) // perfect conductor
			bcPHI_Left(PHI,0);
		// boundary condition: Xright
		if(vct->getXright_neighbor()==-1 && (bcPHIfaceXright ==0 ||  bcPHIfaceXright ==1)) // perfect conductor
			bcPHI_Right(PHI,0);
		// boundary condition: Yleft
		if(vct->getYleft_neighbor()==-1 && (bcPHIfaceYleft ==0 ||  bcPHIfaceYleft ==1)) // perfect conductor
			bcPHI_Left(PHI,1);
		// boundary condition: Yright
		if(vct->getYright_neighbor()==-1 && (bcPHIfaceYright ==0 ||  bcPHIfaceYright ==1)) // perfect conductor
			bcPHI_Right(PHI,1);
		grid->gradC2N(gradPHIX,gradPHIY,PHI,vct);
		sub(Ex,gradPHIX,nxn,nyn);
		sub(Ey,gradPHIY,nxn,nyn);

    }
    // Advance timestep of E fron n to n+theta
    // find the solution with GMRES: the solution is the in the krylov space
    if (vct->getCartesian_rank() ==0)
		cout << "*** MAXWELL SOLVER ***" << endl;
    MaxwellSource(bkrylov,grid,vct);
    phys2solver(xkrylov,Ex,Ey,Ez,nxn,nyn);
    // GMRES
    GMRES(&Field::MaxwellImage, xkrylov, 3*(nxn-2)*(nyn-2),bkrylov,20, 200,GMREStol, grid, vct, this);
    // move from krylov space to physical space
    solver2phys(Exth,Eyth,Ezth,xkrylov,nxn,nyn);
	communicateNode(nxn,nyn,Exth,vct);
    communicateNode(nxn,nyn,Eyth,vct);
    communicateNode(nxn,nyn,Ezth,vct);
	// here you have to put the boundaries conditions
    // boundary condition: Xleft
    if(vct->getXleft_neighbor()==-1 && (bcEMfaceXleft ==0 || bcEMfaceXleft ==3)){
		perfectConductorLeft(Exth,Eyth,Ezth,0,grid);
	}
    // boundary condition: Xright
    if(vct->getXright_neighbor()==-1 && (bcEMfaceXright==0 || bcEMfaceXright==3)){
		perfectConductorRight(Exth,Eyth,Ezth,0,grid);
	}
    // boundary condition: Yleft
    if(vct->getYleft_neighbor()==-1 && (bcEMfaceYleft ==0 || bcEMfaceYleft ==3)){
		perfectConductorLeft(Exth,Eyth,Ezth,1,grid);
	}
    // boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && (bcEMfaceYright==0 ||  bcEMfaceYright==3)){
		perfectConductorRight(Exth,Eyth,Ezth,1,grid);
	}
	// Advance all the fields one full timestep
    addscale(1/th,-(1.0-th)/th,Ex,Exth,nxn,nyn);
    addscale(1/th,-(1.0-th)/th,Ey,Eyth,nxn,nyn);
    addscale(1/th,-(1.0-th)/th,Ez,Ezth,nxn,nyn);
	// SMOOTHING
    smooth(Nvolte,Smooth,Ex,1,grid,vct);
    smooth(Nvolte,Smooth,Ey,1,grid,vct);
    smooth(Nvolte,Smooth,Ez,1,grid,vct);
    // calculate B
    calculateB(grid,vct);
    // deallocate
	delete[] xkrylov;
	delete[] bkrylov;
	delete[] xkrylovPoisson;
	delete[] bkrylovPoisson;
	delArr3(divE,nxc,nyc);
	delArr3(gradPHIX,nxn,nyn);
	delArr3(gradPHIY,nxn,nyn);
	delArr3(gradPHIZ,nxn,nyn);

}



/** Image of Poisson Solver */
inline void EMfields::PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct){
    // allocate  2 two dimensional service vectors
    double ***temp = newArr3(double,nxc,nyc,1);
    double ***im  = newArr3(double,nxc,nyc,1);
    eqValue (0.0, image,(nxc-2)*(nyc-2));
    eqValue (0.0, temp,nxc,nyc);
    eqValue (0.0, im,nxc,nyc);
    // move from krylov space to physical space and communicate ghost cells
    solver2phys(temp,vector,nxc,nyc);
    MPI_Barrier(MPI_COMM_WORLD);
    communicateCenter(nxc,nyc,temp,vct);
	// adjust boundaried before calculating the lap
	if(vct->getXleft_neighbor()==-1 && (bcPHIfaceXleft ==0 ||  bcPHIfaceXleft ==1)) // perfect conductor
		bcPHI_Left(temp,0);
	// boundary condition: Xright
	if(vct->getXright_neighbor()==-1 && (bcPHIfaceXright ==0 ||  bcPHIfaceXright ==1)) // perfect conductor
		bcPHI_Right(temp,0);
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==-1 && (bcPHIfaceYleft ==0 ||  bcPHIfaceYleft ==1)) // perfect conductor
		bcPHI_Left(temp,1);
	// boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && (bcPHIfaceYright ==0 ||  bcPHIfaceYright ==1)) // perfect conductor
		bcPHI_Right(temp,1);
	// calculate the laplacian
    grid->lapC2C(im,temp,vct);
    phys2solver(image,im,nxc,nyc);
    // deallocate temporary array and objects
    delArr3(temp,nxc,nyc);
    delArr3(im,nxc,nyc);


}


/** Calculate sorgent for Maxwell solver */
inline void EMfields::MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct){
	eqValue (0.0, tempC,nxc,nyc);
	eqValue (0.0, tempX,nxn,nyn);
	eqValue (0.0, tempY,nxn,nyn);
	eqValue (0.0, tempZ,nxn,nyn);
	eqValue (0.0, tempXN,nxn,nyn);
	eqValue (0.0, tempYN,nxn,nyn);
	eqValue (0.0, tempZN,nxn,nyn);
	eqValue (0.0, temp2X,nxn,nyn);
	eqValue (0.0, temp2Y,nxn,nyn);
	eqValue (0.0, temp2Z,nxn,nyn);

	//prepare curl of B for known term of Maxwell solver: for the source term
	// communicate ghost
	communicateCenter(nxc,nyc,Bxc,vct);
	communicateCenter(nxc,nyc,Byc,vct);
	communicateCenter(nxc,nyc,Bzc,vct);
	// apply boundary conditions on B
	if(vct->getXleft_neighbor()==-1 && (bcEMfaceXleft ==0 ||  bcEMfaceXleft ==3)) // perfect conductor
		BperfectConductorLeft(0,grid,vct);
	// boundary condition: Xright
	if(vct->getXright_neighbor()==-1 && (bcEMfaceXright==0 || bcEMfaceXright==3)) // perfect conductor
		BperfectConductorRight(0,grid,vct);
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==-1 && (bcEMfaceYleft ==0 || bcEMfaceYleft ==3)) // perfect conductor
		BperfectConductorLeft(1,grid,vct);
	// boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && (bcEMfaceYright==0 || bcEMfaceYright==3)) // perfect conductor
		BperfectConductorRight(1,grid,vct);
	grid->curlC2N(tempXN,tempYN,tempZN,Bxc,Byc,Bzc,vct);
	scale(temp2X,Jxh,-FourPI/c,nxn,nyn);
	scale(temp2Y,Jyh,-FourPI/c,nxn,nyn);
	scale(temp2Z,Jzh,-FourPI/c,nxn,nyn);

	sum(temp2X,tempXN,nxn,nyn);
	sum(temp2Y,tempYN,nxn,nyn);
	sum(temp2Z,tempZN,nxn,nyn);
	scale(temp2X,delt,nxn,nyn);
	scale(temp2Y,delt,nxn,nyn);
	scale(temp2Z,delt,nxn,nyn);
	// gradient of rho hat: rhohat already communicate on by calculate Hat function
	grid->gradC2N_onesided_derBC(tempX,tempY,rhoh,vct);
	scale(tempX,-delt*delt*FourPI,nxn,nyn);
	scale(tempY,-delt*delt*FourPI,nxn,nyn);
	// sum E, past values
	sum(tempX,Ex,nxn,nyn);
	sum(tempY,Ey,nxn,nyn);
	sum(tempZ,Ez,nxn,nyn);
	// sum curl(B) + jhat part
	sum(tempX,temp2X,nxn,nyn);
	sum(tempY,temp2Y,nxn,nyn);
	sum(tempZ,temp2Z,nxn,nyn);
	// Boundary condition in the known term
	// boundary condition: Xleft
	if(vct->getXleft_neighbor()==-1 && (bcEMfaceXleft ==0 ||  bcEMfaceXleft ==3)) // perfect conductor
		perfectConductorLeftS(tempX,tempY,tempZ,0);
	else if (vct->getXleft_neighbor()==-1 && bcEMfaceXleft ==1) // magnetic mirror
		magneticMirrorLeftS(tempX,tempY,tempZ,0);
	else if (vct->getXleft_neighbor()==-1 && (bcEMfaceXleft == 2 || bcEMfaceXleft == 4)) // open boundary
		openLeftS(tempX,tempY,tempZ,0);
	// boundary condition: Xright
	if(vct->getXright_neighbor()==-1 && (bcEMfaceXright==0 || bcEMfaceXright==3)) // perfect conductor
		perfectConductorRightS(tempX,tempY,tempZ,0);
	else if (vct->getXright_neighbor()==-1 && bcEMfaceXright==1 ) // magnetic mirror
		magneticMirrorRightS(tempX,tempY,tempZ,0);
	else if (vct->getXright_neighbor()==-1 && (bcEMfaceXright==2 || bcEMfaceXright==4)) // open boundary
		openRightS(tempX,tempY,tempZ,0);
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==-1 && (bcEMfaceYleft ==0 || bcEMfaceYleft ==3)) // perfect conductor
		perfectConductorLeftS(tempX,tempY,tempZ,1);
	else if (vct->getYleft_neighbor()==-1 && bcEMfaceYleft ==1 ) // magnetic mirror
		magneticMirrorLeftS(tempX,tempY,tempZ,1);
	else if (vct->getYleft_neighbor()==-1 && (bcEMfaceYleft ==2 || bcEMfaceYleft ==4)) // open boundary
		openLeftS(tempX,tempY,tempZ,1);
	// boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && (bcEMfaceYright==0 || bcEMfaceYright==3)) // perfect conductor
		perfectConductorRightS(tempX,tempY,tempZ,1);
	else if (vct->getYright_neighbor()==-1 && bcEMfaceYright==1 ) // magnetic mirror
		magneticMirrorRightS(tempX,tempY,tempZ,1);
	else if (vct->getYright_neighbor()==-1 && (bcEMfaceYright==2 || bcEMfaceYright==4)) // open boundary
		openRightS(tempX,tempY,tempZ,1);
	// physical space -> Krylov space
	phys2solver(bkrylov,tempX,tempY,tempZ,nxn,nyn);

}
/** Mapping of Maxwell image to give to solver */
inline void EMfields::MaxwellImage(double *im, double *vector, Grid *grid, VirtualTopology *vct){
	eqValue (0.0, im, 3*(nxn-2)*(nyn-2));
	eqValue (0.0, imageX,nxn,nyn);
	eqValue (0.0, imageY,nxn,nyn);
	eqValue (0.0, imageZ,nxn,nyn);
	eqValue (0.0, tempZ,nxn,nyn);
	// move from krylov space to physical space
	solver2phys(vectX,vectY,vectZ,vector,nxn,nyn);
	// communicate
	communicateNode(nxn,nyn,vectX,vct);
	communicateNode(nxn,nyn,vectY,vct);
	communicateNode(nxn,nyn,vectZ,vct);
	// boundary condition: Xleft
	if(vct->getXleft_neighbor()==-1 && (bcEMfaceXleft ==0 || bcEMfaceXleft ==3))
		perfectConductorLeft(vectX,vectY,vectZ,0,grid);
	// boundary condition: Xright
	if(vct->getXright_neighbor()==-1 && (bcEMfaceXright==0 || bcEMfaceXright==3))
		perfectConductorRight(vectX,vectY,vectZ,0,grid);
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==-1 && (bcEMfaceYleft ==0 || bcEMfaceYleft ==3))
		perfectConductorLeft(vectX,vectY,vectZ,1,grid);
	// boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && (bcEMfaceYright==0 ||  bcEMfaceYright==3))
		perfectConductorRight(vectX,vectY,vectZ,1,grid);
	// end of the boundaries

	grid->lapN2N(imageX,vectX,vct);
	grid->lapN2N(imageY,vectY,vct);
	grid->lapN2N(imageZ,vectZ,vct);

	neg(imageX,nxn,nyn);
	neg(imageY,nxn,nyn);
	neg(imageZ,nxn,nyn);

	// grad(div(mu dot E(n + theta))       mu dot E(n + theta) = D
	MUdot(Dx, Dy, Dz, vectX, vectY, vectZ, grid);
	grid->divN2C(divC,Dx,Dy);
	// communicate
	communicateCenter(nxc,nyc,divC,vct);
	// here too you should put the BC
	if(vct->getXleft_neighbor()==-1 && (bcPHIfaceXleft ==0 ||  bcPHIfaceXleft ==1)) // perfect conductor
		bcPHI_Left(divC,0);
	// boundary condition: Xright
	if(vct->getXright_neighbor()==-1 && (bcPHIfaceXright ==0 ||  bcPHIfaceXright ==1)) // perfect conductor
		bcPHI_Right(divC,0);
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==-1 && (bcPHIfaceYleft ==0 ||  bcPHIfaceYleft ==1)) // perfect conductor
		bcPHI_Left(divC,1);
	// boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && (bcPHIfaceYright ==0 ||  bcPHIfaceYright ==1)) // perfect conductor
		bcPHI_Right(divC,1);
    grid->gradC2N(tempX,tempY,divC,vct);

	// -lap(E(n +theta)) -  grad(div(mu dot E(n + theta))
	sub(imageX,tempX,nxn,nyn);
	sub(imageY,tempY,nxn,nyn);
	sub(imageZ,tempZ,nxn,nyn);

	//scale delt*delt
	scale(imageX,delt*delt,nxn,nyn);
	scale(imageY,delt*delt,nxn,nyn);
	scale(imageZ,delt*delt,nxn,nyn);

	//  -lap(E(n +theta)) -  grad(div(mu dot E(n + theta)) + eps dot E(n + theta)
	sum(imageX,Dx,nxn,nyn);
	sum(imageY,Dy,nxn,nyn);
	sum(imageZ,Dz,nxn,nyn);
	sum(imageX,vectX,nxn,nyn);
	sum(imageY,vectY,nxn,nyn);
	sum(imageZ,vectZ,nxn,nyn);

	sum(Dx,vectX,nxn,nyn);
	sum(Dy,vectY,nxn,nyn);
	sum(Dz,vectZ,nxn,nyn);
	// boundary condition: Xleft
	if(vct->getXleft_neighbor()==-1 && bcEMfaceXleft ==0) {
		// perfect conductor
		perfectConductorLeftImage(imageX,imageY,imageZ,vectX,vectY,vectZ,0,grid);
	}
	// boundary condition: Xright
	if(vct->getXright_neighbor()==-1 && bcEMfaceXright==0) {
		// perfect conductor
		perfectConductorRightImage(imageX,imageY,imageZ,vectX,vectY,vectZ,0,grid);
	}
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==-1 && bcEMfaceYleft ==0){ // perfect conductor
		perfectConductorLeftImage(imageX,imageY,imageZ,vectX,vectY,vectZ,1,grid);
	}
	// boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && bcEMfaceYright==0){ // perfect conductor
		perfectConductorRightImage(imageX,imageY,imageZ,vectX,vectY,vectZ,1,grid);
  	}
	// move from physical space to krylov space
	phys2solver(im,imageX,imageY,imageZ,nxn,nyn);
}

/** Calculate Magnetic field with the implicit solver: calculate B defined on nodes
With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
inline void EMfields::calculateB(Grid *grid, VirtualTopology *vct){
	eqValue (0.0,tempXC,nxc,nyc);
	eqValue (0.0,tempYC,nxc,nyc);
	eqValue (0.0,tempZC,nxc,nyc);
	// calculate the curl of Eth
	grid->curlN2C(tempXC,tempYC,tempZC,Exth,Eyth,Ezth);
	// update the magnetic field
	addscale_no_BC(-c*dt,1,Bxc,tempXC,nxc,nyc);
	addscale_no_BC(-c*dt,1,Byc,tempYC,nxc,nyc);
	addscale_no_BC(-c*dt,1,Bzc,tempZC,nxc,nyc);
	// communicate ghost
	communicateCenter(nxc,nyc,Bxc,vct);
	communicateCenter(nxc,nyc,Byc,vct);
	communicateCenter(nxc,nyc,Bzc,vct);
	// apply boundary conditions on B
	if(vct->getXleft_neighbor()==-1 && (bcEMfaceXleft ==0 ||  bcEMfaceXleft ==3)) // perfect conductor
		BperfectConductorLeft(0,grid,vct);
	// boundary condition: Xright
	if(vct->getXright_neighbor()==-1 && (bcEMfaceXright==0 || bcEMfaceXright==3)) // perfect conductor
		BperfectConductorRight(0,grid,vct);
	// boundary condition: Yleft
	if(vct->getYleft_neighbor()==-1 && (bcEMfaceYleft ==0 || bcEMfaceYleft ==3)) // perfect conductor
		BperfectConductorLeft(1,grid,vct);
	// boundary condition: Yright
	if(vct->getYright_neighbor()==-1 && (bcEMfaceYright==0 || bcEMfaceYright==3)) // perfect conductor
		BperfectConductorRight(1,grid,vct);
	// interpolate C2N: only after you have applied the BC
	grid->interpC2N(Bxn,Bxc,vct);
	grid->interpC2N(Byn,Byc,vct);
	grid->interpC2N(Bzn,Bzc,vct);



}
/** initialize EM field with transverse electric waves 1D and rotate
anticlockwise (theta degrees)*/
inline void EMfields::initEM_rotate(VirtualTopology *vct, Grid *grid, double B, double theta){
	// initialize E and rhos on nodes
	for (int i=0; i < nxn; i++)
		for (int j=0; j <nyn; j++){
			Ex[i][j][0] =  0.0;
			Ey[i][j][0] =  0.0;
			Ez[i][j][0] =  0.0;
			Bxn[i][j][0] =  B*cos(theta*M_PI/180);
			Byn[i][j][0] =  B*sin(theta*M_PI/180);
			Bzn[i][j][0] = 0.0;
			rhons[0][i][j][0] =  0.07957747154595; // electrons: species is now first index
			rhons[1][i][j][0] =  0.07957747154595; // protons: species is now first index
		}
			// initialize B on centers
			grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);
	// communicate ghost
	communicateNode(nxn,nyn,Bxn,vct);
	communicateNode(nxn,nyn,Byn,vct);
	communicateNode(nxn,nyn,Bzn,vct);
	communicateNode(nxn,nyn,Ex,vct);
	communicateNode(nxn,nyn,Ey,vct);
	communicateNode(nxn,nyn,Ez,vct);

	for (int is=0 ; is<ns; is++)
		grid->interpN2C(rhocs,is,rhons);

}
/** Initialize the EM field with constants values or from restart*/
inline void EMfields::init(VirtualTopology *vct, Grid *grid){
	// initialize E and rhos on nodes
	if (restart1==0){
		for (int i=0; i < nxn; i++)
			for (int j=0; j <nyn; j++){
				Ex[i][j][0] =  0.0;
				Ey[i][j][0] =  0.0;
				Ez[i][j][0] =  0.0;
				Bxn[i][j][0] =  B0x;
				Byn[i][j][0] =  B0y;
				Bzn[i][j][0] =  B0z;
				for (int is=0; is < ns; is++){
					rhons[is][i][j][0] = rhoINIT[is]/FourPI;  // initialize with constant density
				}
			}
				// initialize B on centers
				grid->interpN2C(Bxc,Bxn);
		grid->interpN2C(Byc,Byn);
		grid->interpN2C(Bzc,Bzn);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons); // calculate density on the center
	} else { // EM initialization from RESTART
		MPI_Barrier(MPI_COMM_WORLD);
		if (vct->getCartesian_rank()==0)
			cout << "LOADING EM FIELD FROM RESTART FILE in " + RestartDirName + "/restart.hdf" << endl;
		string name_file;
        stringstream ss;
	    ss << vct->getCartesian_rank();
		name_file = RestartDirName + "/restart" + ss.str() + ".hdf";
		// hdf stuff
		hid_t    file_id, dataspace;
		hid_t    datatype, dataset_id;
		herr_t   status;
		size_t   size;
		hsize_t     dims_out[2];           /* dataset dimensions */
		int status_n;
		MPI_Barrier(MPI_COMM_WORLD);
		// open the hdf file
		file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_id < 0){
			cout << "couldn't open file: " << name_file << endl;
			cout << "RESTART NOT POSSIBLE" << endl;
		}

		dataset_id = H5Dopen(file_id,"/fields/Bx/cycle_0");
		datatype  = H5Dget_type(dataset_id);
		size  = H5Tget_size(datatype);
		dataspace = H5Dget_space(dataset_id);
		status_n  = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);


		// check if the grid is consistent with the grid of the RESTART
		if  (nxn != (dims_out[0] + 2) || nyn != (dims_out[1] + 2)){
			cout << "RESTART GRID NOT CONSISTENT!!!!!!!!!!!!!!!!!!" << endl;
		}
		// Bxn
		double *temp_storage = new double[dims_out[0]*dims_out[1]];
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		int k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Bxn[i][j][0] = temp_storage[k++];

		communicateNode(nxn, nyn, Bxn,vct);
		status = H5Dclose(dataset_id);

		// Byn
		dataset_id = H5Dopen(file_id,"/fields/By/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Byn[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Byn,vct);
		status = H5Dclose(dataset_id);


		// Bzn
		dataset_id = H5Dopen(file_id,"/fields/Bz/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Bzn[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Bzn,vct);
		status = H5Dclose(dataset_id);

		// Ex
		dataset_id = H5Dopen(file_id,"/fields/Ex/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Ex[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Ex,vct);
		status = H5Dclose(dataset_id);


		// Ey
		dataset_id = H5Dopen(file_id,"/fields/Ey/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Ey[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Ey,vct);
		status = H5Dclose(dataset_id);

		// Ez
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


	// initialize B on centers
	grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);
	// communicate ghost
	communicateNode(nxn,nyn,Bxn,vct);
	communicateNode(nxn,nyn,Byn,vct);
	communicateNode(nxn,nyn,Bzn,vct);
	communicateNode(nxn,nyn,Ex,vct);
	communicateNode(nxn,nyn,Ey,vct);
	communicateNode(nxn,nyn,Ez,vct);
	// communicate ghost
	communicateCenter(nxc,nyc,Bxc,vct);
	communicateCenter(nxc,nyc,Byc,vct);
	communicateCenter(nxc,nyc,Bzc,vct);

	for (int is=0 ; is<ns; is++)
		grid->interpN2C(rhocs,is,rhons);
}


/**Add a periodic perturbation in rho exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void EMfields::AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid *grid){

	double alpha;
	alpha=deltaBoB*B0/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);

	ne_mod *= alpha;
	ni_mod *= alpha;
	//cout<<" ne="<<ne_mod<<" ni="<<ni_mod<<" alpha="<<alpha<<endl;
	for (int i=0; i < nxn; i++)
		for (int j=0; j <nyn; j++){
			rhons[0][i][j][0] +=  ne_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + ne_phase);
			rhons[1][i][j][0] +=  ni_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + ni_phase);
		}

			for (int is=0 ; is<ns; is++)
				grid->interpN2C(rhocs,is,rhons);
}


/**Add a periodic perturbation exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void EMfields::AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid *grid, VirtualTopology *vct){

	double alpha;

	alpha=deltaBoB*B0/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);

	Ex_mod *= alpha;
	Ey_mod *= alpha;
	Ez_mod *= alpha;
	Bx_mod *= alpha;
	By_mod *= alpha;
	Bz_mod *= alpha;

	for (int i=0; i < nxn; i++)
		for (int j=0; j <nyn; j++){
			Ex[i][j][0] +=  Ex_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Ex_phase);
			Ey[i][j][0] +=  Ey_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Ey_phase);
			Ez[i][j][0] +=  Ez_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Ez_phase);
			Bxn[i][j][0] +=  Bx_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Bx_phase);
			Byn[i][j][0] +=  By_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + By_phase);
			Bzn[i][j][0] +=  Bz_mod * cos(kx*grid->getXN(i,j,0) + ky*grid->getYN(i,j,0) + Bz_phase);

		}

			// initialize B on centers
			grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);
	// communicate ghost
	communicateNode(nxn,nyn,Bxn,vct);
	communicateNode(nxn,nyn,Byn,vct);
	communicateNode(nxn,nyn,Bzn,vct);
	communicateNode(nxn,nyn,Ex,vct);
	communicateNode(nxn,nyn,Ey,vct);
	communicateNode(nxn,nyn,Ez,vct);



}


/**
*
 *
 * Initialize the EM for the GEM Challenge. The equilibrium chosen for the reconnection challange problem is a Harris equilibrium with
 * a floor in the density outside the current layer. The magnetic field is given by
 *   Bx(y) = B0x tanh(y/l)
 *
 * The charge is set to 1/(4 pi) in order to satisfy the omega_pi = 1 . The 2 species have same charge
 * density to guarantee plasma neutrality
 */
inline void EMfields::initGEM(VirtualTopology *vct, Grid *grid){
        double pertGEM;
	    pertGEM=1.0/10.0*0.0;
	    double pertX;
	    pertX=0.4;
	    double xpert;
	    xpert=0.0;
	    double ypert;
	    ypert=0.0;
	    double exp_pert;
	    exp_pert=1.0;
	if (restart1 ==0){

		// initialize
		if (vct->getCartesian_rank() ==0){
			cout << "-------------------------" << endl;
			cout << "Initialize GEM Challenge " << endl;
			cout << "-------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			cout << "Delta (current sheet thickness) = " << delta << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i];
				if (DriftSpecies[i])
					cout << " DRIFTING " << endl;
				else
					cout << " BACKGROUND " << endl;
			}
			cout << "-------------------------" << endl;
		}

		for (int i=0; i <nxn; i++)
			for (int j=0; j <nyn; j++){
				// initialize the density for species
				for (int is=0; is < ns; is++){
					//if (DriftSpecies[is])
					if (is == 0 || is == 1) // at the beginning no drift is added at the beginning
						rhons[is][i][j][0] = ((rhoINIT[is]/(cosh((grid->getYN(i,j,0)-Ly/2)/delta)*cosh((grid->getYN(i,j,0)-Ly/2)/delta))))/FourPI;
					else
						rhons[is][i][j][0] = rhoINIT[is]/FourPI;
				}
				// electric field
				Ex[i][j][0] =  0.0;
				Ey[i][j][0] =  0.0;
				Ez[i][j][0] =  0.0;
				// Magnetic field
				Bxn[i][j][0] = B0x*tanh((grid->getYN(i,j,0) - Ly/2)/delta);
				// add the initial GEM perturbation
				Bxn[i][j][0] +=(B0x*pertGEM)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,0)/Lx)*sin(M_PI*(grid->getYN(i,j,0)- Ly/2)/Ly  );
				Byn[i][j][0] = B0y -(B0x*pertGEM)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,0)/Lx)*cos(M_PI*(grid->getYN(i,j,0)- Ly/2)/Ly);
				// add the initial X perturbation
				xpert = grid->getXN(i,j,0)- Lx/2;
				ypert = grid->getYN(i,j,0)- Ly/2;
				exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

				Bxn[i][j][0] +=(B0x*pertX)*exp_pert*(
				           -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
				           -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
				           );

				Byn[i][j][0] +=(B0x*pertX)*exp_pert*(
				           cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
				           +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
				           );
				// guide field
				Bzn[i][j][0] = B0z*tanh((grid->getYN(i,j,0) - Ly/2)/delta);
			}
				// communicate ghost
				communicateNode(nxn,nyn,Bxn,vct);
		communicateNode(nxn,nyn,Byn,vct);
		communicateNode(nxn,nyn,Bzn,vct);
		for (int i=0; i <nxc; i++)
			for (int j=0; j <nyc; j++){
				Bxc[i][j][0] = B0x*tanh((grid->getYC(i,j,0) - Ly/2)/delta);
				// add the initial GEM perturbation
				Bxc[i][j][0] +=(B0x*pertGEM)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,0)/Lx)*sin(M_PI*(grid->getYC(i,j,0)- Ly/2)/Ly  );
				Byc[i][j][0] = B0y -(B0x*pertGEM)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,0)/Lx)*cos(M_PI*(grid->getYC(i,j,0)- Ly/2)/Ly);
				// add the initial X perturbation
				xpert = grid->getXC(i,j,0)- Lx/2;
				ypert = grid->getYC(i,j,0)- Ly/2;
				exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

				Bxc[i][j][0] +=(B0x*pertX)*exp_pert*(
				           -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
				           -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
				           );

				Byc[i][j][0] +=(B0x*pertX)*exp_pert*(
				           cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
				           +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
				           );
				// guide field
				Bzc[i][j][0] = B0z*tanh((grid->getYC(i,j,0) - Ly/2)/delta);
			}
				communicateCenter(nxc, nyc,Bxc,vct);
		communicateCenter(nxc, nyc,Byc,vct);
		communicateCenter(nxc, nyc,Bzc,vct);
		communicateNode(nxn,nyn,Ex,vct);
		communicateNode(nxn,nyn,Ey,vct);
		communicateNode(nxn,nyn,Ez,vct);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
	} else {
		init(vct,grid);  // use the fields from restart file
	}
}
/** initialize with uniform distribution */
inline void EMfields::initUniform(VirtualTopology *vct, Grid *grid){
	if (restart1 ==0){

		// initialize
		if (vct->getCartesian_rank() ==0){
			cout << "-------------------------" << endl;
			cout << "Initialize Uniform Distribution " << endl;
			cout << "-------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i] << endl;
			}
			cout << "-------------------------" << endl;
		}
		for (int i=0; i <nxn; i++)
			for (int j=0; j <nyn; j++){
				// initialize the density for species
				for (int is=0; is < ns; is++){
					rhons[is][i][j][0] = rhoINIT[is]/FourPI;
				}
				// electric field
				Ex[i][j][0] =  0.0;
				Ey[i][j][0] =  0.0;
				Ez[i][j][0] =  0.0;
				// Magnetic field
				Bxn[i][j][0] = B0x;
				Byn[i][j][0] = B0y;
				Bzn[i][j][0] = B0z;
			}
				// initialize B on centers
				grid->interpN2C(Bxc,Bxn);
		grid->interpN2C(Byc,Byn);
		grid->interpN2C(Bzc,Bzn);
		// communicate ghost
		communicateNode(nxn,nyn,Bxn,vct);
		communicateNode(nxn,nyn,Byn,vct);
		communicateNode(nxn,nyn,Bzn,vct);
		communicateNode(nxn,nyn,Ex,vct);
		communicateNode(nxn,nyn,Ey,vct);
		communicateNode(nxn,nyn,Ez,vct);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
	} else {
		init(vct,grid);  // use the fields from restart file
	}
}
/** Initialize the field to study the lower hybrid drift instability */
inline void EMfields::initLHDI(VirtualTopology *vct, Grid *grid){
	if (restart1 ==0){

		// initialize
		if (vct->getCartesian_rank() ==0){
			cout << "-------------------------" << endl;
			cout << "Initialize field for LHDI" << endl;
			cout << "-------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			cout << "Delta (current sheet thickness) = " << delta << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i];
				if (DriftSpecies[i])
					cout << " DRIFTING " << endl;
				else
					cout << " BACKGROUND " << endl;
			}
			cout << "-------------------------" << endl;
		}
		for (int i=0; i <nxn; i++)
			for (int j=0; j <nyn; j++){
				// initialize the density for species
				for (int is=0; is < ns; is++){
					if (DriftSpecies[is])
						rhons[is][i][j][0] = ((rhoINIT[is]/(cosh((grid->getYN(i,j,0)-Ly/2)/delta)*cosh((grid->getYN(i,j,0)-Ly/2)/delta))))/FourPI;
					else
						rhons[is][i][j][0] = rhoINIT[is]/FourPI;
				}
				// electric field
				Ex[i][j][0] =  0.0;
				Ey[i][j][0] =  0.0;
				Ez[i][j][0] =  0.0;
				// Bx
				Bxn[i][j][0] = B0x;
				// By
				Byn[i][j][0] = B0y;
				// Bz
				Bzn[i][j][0] = B0z*tanh((grid->getYN(i,j,0)-Ly/2)/delta);
			}
				// initialize B on centers
				grid->interpN2C(Bxc,Bxn);
		grid->interpN2C(Byc,Byn);
		grid->interpN2C(Bzc,Bzn);
		// communicate ghost
		communicateNode(nxn,nyn,Bxn,vct);
		communicateNode(nxn,nyn,Byn,vct);
		communicateNode(nxn,nyn,Bzn,vct);
		communicateNode(nxn,nyn,Ex,vct);
		communicateNode(nxn,nyn,Ey,vct);
		communicateNode(nxn,nyn,Ez,vct);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
	} else {
		init(vct,grid);  // use the fields from restart file
	}
}
/** initialize with uniform distribution */
inline void EMfields::initAlfven(VirtualTopology *vct, Grid *grid){
        double pertGEM;
	pertGEM=1.0/10.0;
	if (restart1 ==0){

		// initialize
		if (vct->getCartesian_rank() ==0){
			cout << "-------------------------" << endl;
			cout << "Initialize GEM Challenge " << endl;
			cout << "-------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			cout << "Delta (current sheet thickness) = " << delta << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i];
				if (DriftSpecies[i])
					cout << " DRIFTING " << endl;
				else
					cout << " BACKGROUND " << endl;
			}
			cout << "-------------------------" << endl;
		}

		for (int i=0; i <nxn; i++)
			for (int j=0; j <nyn; j++){
				// initialize the density for species
				for (int is=0; is < ns; is++){
				rhons[is][i][j][0] = rhoINIT[is]/FourPI;
				}
				// electric field
				Ex[i][j][0] =  0.0;
				Ey[i][j][0] =  0.0;
				Ez[i][j][0] =  0.0;
				// Magnetic field
				Bxn[i][j][0] = B0x;
				Byn[i][j][0] = B0y *sin(2*M_PI*grid->getXN(i,j,0)/Lx);
				Bzn[i][j][0] = B0z;
			}
				// communicate ghost
				communicateNode(nxn,nyn,Bxn,vct);
		communicateNode(nxn,nyn,Byn,vct);
		communicateNode(nxn,nyn,Bzn,vct);
		for (int i=0; i <nxc; i++)
			for (int j=0; j <nyc; j++){
				Bxc[i][j][0] = B0x;
				Byc[i][j][0] = B0y *sin(2*M_PI*grid->getXC(i,j,0)/Lx);
				Bzc[i][j][0] = B0z;
			}
				communicateCenter(nxc, nyc,Bxc,vct);
		communicateCenter(nxc, nyc,Byc,vct);
		communicateCenter(nxc, nyc,Bzc,vct);
		communicateNode(nxn,nyn,Ex,vct);
		communicateNode(nxn,nyn,Ey,vct);
		communicateNode(nxn,nyn,Ez,vct);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
	} else {
		init(vct,grid);  // use the fields from restart file
	}
}
/** add magnetic dipole in position (x0,y0) with magnetic field ampl B0 in direction Y*/
inline  void EMfields::addDipole(double B0, double x0, double y0,Grid *grid){
	double x_displ, y_displ;
	for (int i=0; i <nxn; i++)
		for (int j=0; j <nyn; j++){
			x_displ = grid->getXN(i,j,0) - x0;
			y_displ = grid->getYN(i,j,0) - y0;
			if (x_displ!= 0 && y_displ !=0){
				Bxn[i][j][0] += B0*(x_displ*y_displ)/pow(x_displ*x_displ + y_displ*y_displ,2.5);
				Byn[i][j][0] += B0*(-(1/3)*(x_displ*x_displ+ y_displ*y_displ))/pow(x_displ*x_displ + y_displ*y_displ,2.5);
				Bzn[i][j][0] += 0.0;
			} else { // eliminate the singularity
				Bxn[i][j][0] += B0;
				Byn[i][j][0] += 0.0;
				Bzn[i][j][0] += 0.0;
			}

		}
			// calculate on the centers
			grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);

}
/** add IMF field*/
inline  void EMfields::addIMF(double B_IMF, double theta, Grid *grid){
	for (int i=0; i <nxn; i++)
		for (int j=0; j <nyn; j++){
			Bxn[i][j][0] += B_IMF*cos((theta/90)*M_PI);
			Byn[i][j][0] += B_IMF*sin((theta/90)*M_PI);
			Bzn[i][j][0] += 0.0;
		}
			// calculate on the centers
			grid->interpN2C(Bxc,Bxn);
	grid->interpN2C(Byc,Byn);
	grid->interpN2C(Bzc,Bzn);

}


/** set to 0 all the densities fields */
inline  void EMfields::setZeroDensities(){
	for (register int i=0; i < nxn; i++)
		for (register int j=0; j < nyn; j++){
			Jx[i][j][0]   = 0.0;
			Jxh[i][j][0]  = 0.0;
			Jy[i][j][0]   = 0.0;
			Jyh[i][j][0]  = 0.0;
			Jz[i][j][0]   = 0.0;
			Jzh[i][j][0]  = 0.0;
			rhon[i][j][0] = 0.0;
        }
			for (register int i=0; i < nxc; i++)
				for (register int j=0; j < nyc; j++){
					rhoc[i][j][0] = 0.0;
					rhoh[i][j][0] = 0.0;
				}
					for (register int kk=0; kk < ns; kk++)
						for (register int i=0; i < nxn; i++)
							for (register int j=0; j < nyn; j++){
								rhons[kk][i][j][0] = 0.0;
								Jxs[kk][i][j][0]   = 0.0;
								Jys[kk][i][j][0]   = 0.0;
								Jzs[kk][i][j][0]   = 0.0;
								pXXsn[kk][i][j][0] = 0.0;
								pXYsn[kk][i][j][0] = 0.0;
								pXZsn[kk][i][j][0] = 0.0;
								pYYsn[kk][i][j][0] = 0.0;
								pYZsn[kk][i][j][0] = 0.0;
								pZZsn[kk][i][j][0] = 0.0;
							}
}


/** communicate ghost for grid -> Particles interpolation */
inline void EMfields::communicateGhostP2G(int ns,int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct){

	communicateInterp(nxn,nyn,ns,rhons,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,Jxs,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,Jys,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,Jzs,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,pXXsn,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,pXYsn,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,pXZsn,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,pYYsn,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,pYZsn,0,0,0,0,dx,dy,vct);
	communicateInterp(nxn,nyn,ns,pZZsn,0,0,0,0,dx,dy,vct);

	communicateNode(nxn, nyn, rhons, ns,vct);
	communicateNode(nxn, nyn, Jxs, ns,vct);
	communicateNode(nxn, nyn, Jys, ns,vct);
	communicateNode(nxn, nyn, Jzs, ns,vct);
	communicateNode(nxn, nyn, pXXsn, ns,vct);
	communicateNode(nxn, nyn, pXYsn, ns,vct);
	communicateNode(nxn, nyn, pXZsn, ns,vct);
	communicateNode(nxn, nyn, pYYsn, ns,vct);
	communicateNode(nxn, nyn, pYZsn, ns,vct);
	communicateNode(nxn, nyn, pZZsn, ns,vct);

}
/** adjust densities on boundaries that are not periodic */
inline  void EMfields::adjustNonPeriodicDensities(VirtualTopology *vct){
	// correct Xleft side if not periodic
	if (vct->getXleft_neighbor()==-1){
		for (int is=0; is < ns; is++)
			for (int i=2; i < nyn-2;i++){
				rhons[is][1][i][0]+=   rhons[is][1][i][0];
				Jxs[is][1][i][0]  +=   Jxs[is][1][i][0];
				Jys[is][1][i][0]  +=   Jys[is][1][i][0];
				Jzs[is][1][i][0]  +=   Jzs[is][1][i][0];
				pXXsn[is][1][i][0]  += pXXsn[is][1][i][0];
				pXYsn[is][1][i][0]  += pXYsn[is][1][i][0];
				pXZsn[is][1][i][0]  += pXZsn[is][1][i][0];
				pYYsn[is][1][i][0]  += pYYsn[is][1][i][0];
				pYZsn[is][1][i][0]  += pYZsn[is][1][i][0];
				pZZsn[is][1][i][0]  += pZZsn[is][1][i][0];
			}
	}
		// correct Xright side if not periodic
		if (vct->getXright_neighbor()==-1){
			for (int is=0; is < ns; is++)
				for (int i=2; i < nyn-2;i++){
					rhons[is][nxn-2][i][0]+=   rhons[is][nxn-2][i][0];
					Jxs[is][nxn-2][i][0]  +=   Jxs[is][nxn-2][i][0];
					Jys[is][nxn-2][i][0]  +=   Jys[is][nxn-2][i][0];
					Jzs[is][nxn-2][i][0]  +=   Jzs[is][nxn-2][i][0];
					pXXsn[is][nxn-2][i][0]  += pXXsn[is][nxn-2][i][0];
					pXYsn[is][nxn-2][i][0]  += pXYsn[is][nxn-2][i][0];
					pXZsn[is][nxn-2][i][0]  += pXZsn[is][nxn-2][i][0];
					pYYsn[is][nxn-2][i][0]  += pYYsn[is][nxn-2][i][0];
					pYZsn[is][nxn-2][i][0]  += pYZsn[is][nxn-2][i][0];
					pZZsn[is][nxn-2][i][0]  += pZZsn[is][nxn-2][i][0];
				}
		}
			// correct Yleft side if not periodic
			if (vct->getYleft_neighbor()==-1){
				for (int is=0; is < ns; is++)
					for (int i=2; i < nxn-2;i++){
						rhons[is][i][1][0]+=   rhons[is][i][1][0];
						Jxs[is][i][1][0]  +=   Jxs[is][i][1][0];
						Jys[is][i][1][0]  +=   Jys[is][i][1][0];
						Jzs[is][i][1][0]  +=   Jzs[is][i][1][0];
						pXXsn[is][i][1][0]  += pXXsn[is][i][1][0];
						pXYsn[is][i][1][0]  += pXYsn[is][i][1][0];
						pXZsn[is][i][1][0]  += pXZsn[is][i][1][0];
						pYYsn[is][i][1][0]  += pYYsn[is][i][1][0];
						pYZsn[is][i][1][0]  += pYZsn[is][i][1][0];
						pZZsn[is][i][1][0]  += pZZsn[is][i][1][0];
					}
			}
				// correct Yright side if not periodic
				if (vct->getYright_neighbor()==-1){
					for (int is=0; is < ns; is++)
						for (int i=2; i < nxn-2;i++){
							rhons[is][i][nyn-2][0]+=   rhons[is][i][nyn-2][0];
							Jxs[is][i][nyn-2][0]  +=   Jxs[is][i][nyn-2][0];
							Jys[is][i][nyn-2][0]  +=   Jys[is][i][nyn-2][0];
							Jzs[is][i][nyn-2][0]  +=   Jzs[is][i][nyn-2][0];
							pXXsn[is][i][nyn-2][0]  += pXXsn[is][i][nyn-2][0];
							pXYsn[is][i][nyn-2][0]  += pXYsn[is][i][nyn-2][0];
							pXZsn[is][i][nyn-2][0]  += pXZsn[is][i][nyn-2][0];
							pYYsn[is][i][nyn-2][0]  += pYYsn[is][i][nyn-2][0];
							pYZsn[is][i][nyn-2][0]  += pYZsn[is][i][nyn-2][0];
							pZZsn[is][i][nyn-2][0]  += pZZsn[is][i][nyn-2][0];
						}
				}
					//correct XleftYleft
					if (vct->getXleftYleft_neighbor()==-1){
						for (int is=0; is < ns; is++){
							rhons[is][1][1][0]+=   rhons[is][1][1][0];
							Jxs[is][1][1][0]  +=   Jxs[is][1][1][0];
							Jys[is][1][1][0]  +=   Jys[is][1][1][0];
							Jzs[is][1][1][0]  +=   Jzs[is][1][1][0];
							pXXsn[is][1][1][0]  += pXXsn[is][1][1][0];
							pXYsn[is][1][1][0]  += pXYsn[is][1][1][0];
							pXZsn[is][1][1][0]  += pXZsn[is][1][1][0];
							pYYsn[is][1][1][0]  += pYYsn[is][1][1][0];
							pYZsn[is][1][1][0]  += pYZsn[is][1][1][0];
							pZZsn[is][1][1][0]  += pZZsn[is][1][1][0];
						}
					}
					//correct XleftYright
					if (vct->getXleftYright_neighbor()==-1){
						for (int is=0; is < ns; is++){
							rhons[is][1][nyn-2][0]+=   rhons[is][1][nyn-2][0];
							Jxs[is][1][nyn-2][0]  +=   Jxs[is][1][nyn-2][0];
							Jys[is][1][nyn-2][0]  +=   Jys[is][1][nyn-2][0];
							Jzs[is][1][nyn-2][0]  +=   Jzs[is][1][nyn-2][0];
							pXXsn[is][1][nyn-2][0]  += pXXsn[is][1][nyn-2][0];
							pXYsn[is][1][nyn-2][0]  += pXYsn[is][1][nyn-2][0];
							pXZsn[is][1][nyn-2][0]  += pXZsn[is][1][nyn-2][0];
							pYYsn[is][1][nyn-2][0]  += pYYsn[is][1][nyn-2][0];
							pYZsn[is][1][nyn-2][0]  += pYZsn[is][1][nyn-2][0];
							pZZsn[is][1][nyn-2][0]  += pZZsn[is][1][nyn-2][0];
						}
					}
					//correct XrightYleft
					if (vct->getXrightYleft_neighbor()==-1){
						for (int is=0; is < ns; is++){
							rhons[is][nxn-2][1][0]+=   rhons[is][nxn-2][1][0];
							Jxs[is][nxn-2][1][0]  +=   Jxs[is][nxn-2][1][0];
							Jys[is][nxn-2][1][0]  +=   Jys[is][nxn-2][1][0];
							Jzs[is][nxn-2][1][0]  +=   Jzs[is][nxn-2][1][0];
							pXXsn[is][nxn-2][1][0]  += pXXsn[is][nxn-2][1][0];
							pXYsn[is][nxn-2][1][0]  += pXYsn[is][nxn-2][1][0];
							pXZsn[is][nxn-2][1][0]  += pXZsn[is][nxn-2][1][0];
							pYYsn[is][nxn-2][1][0]  += pYYsn[is][nxn-2][1][0];
							pYZsn[is][nxn-2][1][0]  += pYZsn[is][nxn-2][1][0];
							pZZsn[is][nxn-2][1][0]  += pZZsn[is][nxn-2][1][0];
						}
					}
					//correct XrightYright
					if (vct->getXrightYright_neighbor()==-1){
						for (int is=0; is < ns; is++){
							rhons[is][nxn-2][nyn-2][0]+=   rhons[is][nxn-2][nyn-2][0];
							Jxs[is][nxn-2][nyn-2][0]  +=   Jxs[is][nxn-2][nyn-2][0];
							Jys[is][nxn-2][nyn-2][0]  +=   Jys[is][nxn-2][nyn-2][0];
							Jzs[is][nxn-2][nyn-2][0]  +=   Jzs[is][nxn-2][nyn-2][0];
							pXXsn[is][nxn-2][nyn-2][0]  += pXXsn[is][nxn-2][nyn-2][0];
							pXYsn[is][nxn-2][nyn-2][0]  += pXYsn[is][nxn-2][nyn-2][0];
							pXZsn[is][nxn-2][nyn-2][0]  += pXZsn[is][nxn-2][nyn-2][0];
							pYYsn[is][nxn-2][nyn-2][0]  += pYYsn[is][nxn-2][nyn-2][0];
							pYZsn[is][nxn-2][nyn-2][0]  += pYZsn[is][nxn-2][nyn-2][0];
							pZZsn[is][nxn-2][nyn-2][0]  += pZZsn[is][nxn-2][nyn-2][0];
						}
					}


}

/** add an amount of charge density to charge density field at node X,Y */
inline void EMfields::addRho(double ***weight, int X, int Y, int Z, int is){
  	rhons[is][X-1][Y-1][0] += weight[0][0][0]*invVOL;
  	rhons[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
	rhons[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
	rhons[is][X][Y][0]    += weight[1][1][0]*invVOL;
}
/** add an amount of charge density to current density - direction X to current density field on the node*/
inline void EMfields::addJx(double ***weight, int X, int Y, int Z, int is){
	Jxs[is][X-1][Y-1][0] += weight[0][0][0]*invVOL;
	Jxs[is][X-1][Y][0]   += weight[0][1][0]*invVOL;
	Jxs[is][X][Y-1][0]   += weight[1][0][0]*invVOL;
	Jxs[is][X][Y][0]     += weight[1][1][0]*invVOL;}
/** add an amount of current density - direction Y to current density field on the node */
inline  void EMfields::addJy(double ***weight, int X, int Y, int Z, int is){
	Jys[is][X-1][Y-1][0]  += weight[0][0][0]*invVOL;
	Jys[is][X-1][Y][0]    += weight[0][1][0]*invVOL;
	Jys[is][X][Y-1][0]    += weight[1][0][0]*invVOL;
	Jys[is][X][Y][0]      += weight[1][1][0]*invVOL;
}
/** add an amount of current density - direction Z to current density field on the node */
inline  void EMfields::addJz(double ***weight, int X, int Y, int Z, int is){
	Jzs[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	Jzs[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	Jzs[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	Jzs[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XX to current density field on the node */
inline  void EMfields::addPxx(double ***weight, int X, int Y, int Z,int is){
	pXXsn[is][X-1][Y-1][0]   += weight[0][0][0]*invVOL;
	pXXsn[is][X-1][Y][0]     += weight[0][1][0]*invVOL;
	pXXsn[is][X][Y-1][0]     += weight[1][0][0]*invVOL;
	pXXsn[is][X][Y][0]       += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XY to current density field on the node*/
inline  void EMfields::addPxy(double ***weight, int X, int Y, int Z,int is){
	pXYsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	pXYsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	pXYsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	pXYsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XZ to current density field on the node */
inline  void EMfields::addPxz(double ***weight, int X, int Y,int Z, int is){
	pXZsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	pXZsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	pXZsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	pXZsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction YY to current density field on the node*/
inline  void EMfields::addPyy(double ***weight, int X, int Y, int Z, int is){
	pYYsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	pYYsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	pYYsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	pYYsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction YZ to current density field on the node */
inline  void EMfields::addPyz(double ***weight, int X, int Y, int Z, int is){
	pYZsn[is][X-1][Y-1][0]+= weight[0][0][0]*invVOL;
	pYZsn[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
	pYZsn[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
	pYZsn[is][X][Y][0]    += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction ZZ to current density field on the node */
inline  void EMfields::addPzz(double ***weight, int X, int Y, int Z, int is){
	pZZsn[is][X-1][Y-1][0]+= weight[0][0][0]*invVOL;
	pZZsn[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
	pZZsn[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
	pZZsn[is][X][Y][0]    += weight[1][1][0]*invVOL;
}


/**SPECIES: Sum the charge density of different species on NODES*/
inline void EMfields::sumOverSpecies(VirtualTopology *vct){
	// calculate the correct densities on the boundaries
	adjustNonPeriodicDensities(vct);

	// sum to get the total density
	for (int is=0; is<ns; is++)
		for (register int i=0; i <nxn; i++)
			for (register int j=0; j <nyn; j++)
				rhon[i][j][0]     +=rhons[is][i][j][0];
	communicateNode(nxn,nyn,rhon,vct);  // this added by STEF in KU Leuven

}

/**SPECIES: Sum current density for different species */
inline void EMfields::sumOverSpeciesJ(){
	for (int is=0; is < ns; is++)
		for (register int i=0; i <nxn; i++)
			for (register int j=0; j <nyn; j++) {
				Jx[i][j][0]     +=Jxs[is][i][j][0];
				Jy[i][j][0]     +=Jys[is][i][j][0];
				Jz[i][j][0]     +=Jzs[is][i][j][0];
			}
}

/** interpolate rho from node to center */
inline  void EMfields::interpDensitiesN2C(VirtualTopology *vct, Grid *grid){
	grid->interpN2C(rhoc,rhon);
	communicateCenter(nxc,nyc,rhoc,vct); // this added by stef in KU Leuven
}
/** Calculate hat rho hat, Jx hat, Jy hat, Jz hat */
inline void EMfields::calculateHatFunctions(Grid *grid, VirtualTopology *vct){
	eqValue (0.0,tempXC,nxc,nyc);
	eqValue (0.0,tempYC,nxc,nyc);
	eqValue (0.0,tempZC,nxc,nyc);

	//smoothing
	smooth(Nvolte,Smooth,rhoc,0,grid,vct);

	// calculate j hat
	for (int is=0; is < ns; is++){
		grid->divSymmTensorN2C(tempXC,tempYC,tempZC,pXXsn,pXYsn,pXZsn,pYYsn,pYZsn,pZZsn,is);
		scale(tempXC,-dt/2.0,nxc,nyc);
		scale(tempYC,-dt/2.0,nxc,nyc);
		scale(tempZC,-dt/2.0,nxc,nyc);
		// communicate
		communicateCenter(nxc,nyc,tempXC,vct);
		communicateCenter(nxc,nyc,tempYC,vct);
		communicateCenter(nxc,nyc,tempZC,vct);
		grid->interpC2N_BC(tempXN,tempXC,vct);
		grid->interpC2N_BC(tempYN,tempYC,vct);
		grid->interpC2N_BC(tempZN,tempZC,vct);

		sum(tempXN,Jxs,nxn,nyn,is);
		sum(tempYN,Jys,nxn,nyn,is);
		sum(tempZN,Jzs,nxn,nyn,is);
		// PIDOT
		PIdot(Jxh,Jyh,Jzh,tempXN,tempYN,tempZN,is,grid);}
	// Jxh, Jyh, Jzh are here correctly defined only on  active nodes
	communicateNode(nxn, nyn, Jxh, vct);
	communicateNode(nxn, nyn, Jyh, vct);
	communicateNode(nxn, nyn, Jzh, vct);

	smooth(Nvolte,Smooth,Jxh,1,grid,vct);
	smooth(Nvolte,Smooth,Jyh,1,grid,vct);
	smooth(Nvolte,Smooth,Jzh,1,grid,vct);

	// calculate rhoh
	grid->divN2C(tempXC,Jxh,Jyh);
	scale(tempXC,-dt*th,nxc,nyc);
	sum(tempXC,rhoc,nxc,nyc);
	eq(rhoh,tempXC,nxc,nyc);
	// at this point is communicated
	communicateCenter(nxc, nyc, rhoh, vct);



}

/** Calculate PI dot (vectX, vectY, vectZ)*/
inline void EMfields::PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid *grid){
	double beta, edotb, omcx, omcy, omcz,denom;
	beta = .5*qom[ns]*dt/c;
	for(int i=1; i <nxn-1;i++)
		for(int j=1; j <nyn-1;j++){
			omcx = beta*Bxn[i][j][0];
			omcy = beta*Byn[i][j][0];
			omcz = beta*Bzn[i][j][0];
			edotb =  vectX[i][j][0]*omcx +  vectY[i][j][0]*omcy + vectZ[i][j][0]*omcz;
			denom  = 1/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
			PIdotX[i][j][0] += (vectX[i][j][0] +(vectY[i][j][0]*omcz - vectZ[i][j][0]*omcy + edotb*omcx))*denom;
			PIdotY[i][j][0] += (vectY[i][j][0] +(vectZ[i][j][0]*omcx - vectX[i][j][0]*omcz + edotb*omcy))*denom;
			PIdotZ[i][j][0] += (vectZ[i][j][0] +(vectX[i][j][0]*omcy - vectY[i][j][0]*omcx + edotb*omcz))*denom;
		}
}



/** Calculate MU dot (vectX, vectY, vectZ)*/
inline void EMfields::MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid *grid){
	double beta, edotb, omcx, omcy, omcz, denom;
	for(int i=1; i < nxn-1;i++)
		for(int j=1; j < nyn-1;j++){
			MUdotX[i][j][0] = 0.0;
			MUdotY[i][j][0] = 0.0;
			MUdotZ[i][j][0] = 0.0;
		}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int i=1; i < nxn-1;i++)
					for(int j=1; j < nyn-1;j++){
						omcx = beta*Bxn[i][j][0];
						omcy = beta*Byn[i][j][0];
						omcz = beta*Bzn[i][j][0];
						edotb =  vectX[i][j][0]*omcx +  vectY[i][j][0]*omcy + vectZ[i][j][0]*omcz;
						denom  = FourPI/2*delt*dt/c*qom[is]*rhons[is][i][j][0]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);

						MUdotX[i][j][0] += (vectX[i][j][0] +(vectY[i][j][0]*omcz - vectZ[i][j][0]*omcy + edotb*omcx))*denom;
						MUdotY[i][j][0] += (vectY[i][j][0] +(vectZ[i][j][0]*omcx - vectX[i][j][0]*omcz + edotb*omcy))*denom;
						MUdotZ[i][j][0] += (vectZ[i][j][0] +(vectX[i][j][0]*omcy - vectY[i][j][0]*omcx + edotb*omcz))*denom;
					}
			}
}
/** calculate the suceptibility tensor on left boundary */
inline void EMfields::sustensorLeft(double *susxx, double *susxy, double *susxz, double* susyx, double* susyy, double* susyz, double* suszx, double* suszy, double* suszz,  int dir){
	double beta, omcx, omcy, omcz, denom;
	switch(dir){
		case 0: // direction X
			for(int i=0; i < nyn;i++){
				susxx[i] = 1.0;
				susxy[i] = 0.0;
				susxz[i] = 0.0;
				susyx[i] = 0.0;
				susyy[i] = 1.0;
				susyz[i] = 0.0;
				suszx[i] = 0.0;
				suszy[i] = 0.0;
				suszz[i] = 1.0;
			}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int j=0; j < nyn;j++){
					omcx = beta*Bxn[1][j][0];
					omcy = beta*Byn[1][j][0];
					omcz = beta*Bzn[1][j][0];
					denom  = FourPI/2*delt*dt/c*qom[is]*rhons[is][1][j][0]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
					susxx[j] += (1.0+omcx*omcx)*denom;
					susxy[j] += (omcz+omcx*omcy)*denom;
					susxz[j] += (omcx*omcz-omcy)*denom;
					susyx[j] += (omcx*omcy-omcz)*denom;
					susyy[j] += (1.0+omcy*omcy)*denom;
					susyz[j] += (omcx+omcy*omcz)*denom;
					suszx[j] += (omcy+omcx*omcz)*denom;
					suszy[j] += (omcy*omcz-omcx)*denom;
					suszz[j] += (1.0+omcz*omcz)*denom;

				}
			}
			break;
		case 1: // direction Y

			for(int i=0; i < nxn;i++){
				susxx[i] = 1.0;
				susxy[i] = 0.0;
				susxz[i] = 0.0;
				susyx[i] = 0.0;
				susyy[i] = 1.0;
				susyz[i] = 0.0;
				suszx[i] = 0.0;
				suszy[i] = 0.0;
				suszz[i] = 1.0;
			}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int i=0; i < nxn;i++){
					omcx = beta*Bxn[i][1][0];
					omcy = beta*Byn[i][1][0];
					omcz = beta*Bzn[i][1][0];
					denom  = FourPI/2*delt*dt/c*qom[is]*rhons[is][i][1][0]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
					susxx[i] += (1.0+omcx*omcx)*denom;
					susxy[i] += (omcz+omcx*omcy)*denom;
					susxz[i] += (omcx*omcz-omcy)*denom;
					susyx[i] += (omcx*omcy-omcz)*denom;
					susyy[i] += (1.0+omcy*omcy)*denom;
					susyz[i] += (omcx+omcy*omcz)*denom;
					suszx[i] += (omcy+omcx*omcz)*denom;
					suszy[i] += (omcy*omcz-omcx)*denom;
					suszz[i] += (1.0+omcz*omcz)*denom;
				}
			}
			break;

	}
}
/** calculate the suceptibility tensor on right boundary */
inline void EMfields::sustensorRight(double *susxx, double *susxy, double *susxz, double* susyx, double* susyy, double* susyz, double* suszx, double* suszy, double* suszz,  int dir){
	double beta, omcx, omcy, omcz, denom;
	switch(dir){
		case 0: // direction X
			for(int i=0; i < nyn;i++){
				susxx[i] = 1.0;
				susxy[i] = 0.0;
				susxz[i] = 0.0;
				susyx[i] = 0.0;
				susyy[i] = 1.0;
				susyz[i] = 0.0;
				suszx[i] = 0.0;
				suszy[i] = 0.0;
				suszz[i] = 1.0;
			}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int j=0; j < nyn;j++){
					omcx = beta*Bxn[nxn-2][j][0];
					omcy = beta*Byn[nxn-2][j][0];
					omcz = beta*Bzn[nxn-2][j][0];
					denom  = FourPI/2*delt*dt/c*qom[is]*rhons[is][nxn-2][j][0]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
					susxx[j] += (1.0+omcx*omcx)*denom;
					susxy[j] += (omcz+omcx*omcy)*denom;
					susxz[j] += (omcx*omcz-omcy)*denom;
					susyx[j] += (omcx*omcy-omcz)*denom;
					susyy[j] += (1.0+omcy*omcy)*denom;
					susyz[j] += (omcx+omcy*omcz)*denom;
					suszx[j] += (omcy+omcx*omcz)*denom;
					suszy[j] += (omcy*omcz-omcx)*denom;
					suszz[j] += (1.0+omcz*omcz)*denom;

				}
			}
			break;
		case 1: // direction Y
			for(int i=0; i < nxn;i++){
				susxx[i] = 1.0;
				susxy[i] = 0.0;
				susxz[i] = 0.0;
				susyx[i] = 0.0;
				susyy[i] = 1.0;
				susyz[i] = 0.0;
				suszx[i] = 0.0;
				suszy[i] = 0.0;
				suszz[i] = 1.0;
			}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int i=0; i < nxn;i++){
					omcx = beta*Bxn[i][nyn-2][0];
					omcy = beta*Byn[i][nyn-2][0];
					omcz = beta*Bzn[i][nyn-2][0];
					denom  = FourPI/2*delt*dt/c*qom[is]*rhons[is][i][nyn-2][0]/(1.0 + omcx*omcx + omcy*omcy + omcz*omcz);
					susxx[i] += (1.0+omcx*omcx)*denom;
					susxy[i] += (omcz+omcx*omcy)*denom;
					susxz[i] += (omcx*omcz-omcy)*denom;
					susyx[i] += (omcx*omcy-omcz)*denom;
					susyy[i] += (1.0+omcy*omcy)*denom;
					susyz[i] += (omcx+omcy*omcz)*denom;
					suszx[i] += (omcy+omcx*omcz)*denom;
					suszy[i] += (omcy*omcz-omcx)*denom;
					suszz[i] += (1.0+omcz*omcz)*denom;
				}
			}
			break;

	}
}

/**
Interpolation smoothing:
 Smoothing (vector must already have ghost cells)
 TO MAKE SMOOTH value as to be different from 1.0
 type = 0 --> center based vector ;
 type = 1 --> node based vector   ;
 */
inline void  EMfields::smooth(int nvolte, double value,double ***vector, bool type, Grid *grid, VirtualTopology *vct){
	double alpha;
	int nx, ny;

	for (int icount=1; icount<nvolte+1; icount++)
	{
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
			for (int i=1; i < nx-1;i++)
				for (int j=1; j < ny-1;j++)
					vector[i][j][0] = temp[i][j];
			// communicate
			if (type == 0)
				communicateCenter(nx, ny, vector, vct);
			else
				communicateNode(nx, ny, vector, vct);

			delArr(temp,nx);
		} // end of if (value !=1)
	}
}

/**
Interpolation smoothing:
 Smoothing (vector must already have ghost cells)
 TO MAKE SMOOTH value as to be different from 1.0
 type = 0 --> center based vector ;
 type = 1 --> node based vector   ;
 */
inline void  EMfields::smoothE(int nvolte, double value,double ***vector, bool type, Grid *grid, VirtualTopology *vct){
	double alpha;
	int nx, ny;

	for (int icount=1; icount<nvolte+1; icount++)
	{
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

			//eq(vector,temp,nx,ny);
			// copy temp in vector
			for (int i=1; i < nx-1;i++)
				for (int j=1; j < ny-1;j++)
					vector[i][j][0] = temp[i][j];
			// communicate
			if (type == 0)
				communicateCenter(nx, ny, vector, vct);
			else
				communicateNode(nx, ny, vector, vct);

			delArr(temp,nx);
		} // end of if (value !=1)
	}
}
/**

SPECIES: Interpolation smoothing
TO MAKE SMOOTH value as to be different from 1.0
type = 0 --> center based vector
type = 1 --> node based vector
*/
inline void  EMfields::smooth(int nvolte, double value,double ****vector,int is, bool type, Grid *grid, VirtualTopology *vct){
	double alpha;
	int nx, ny;

	for (int icount=1; icount<nvolte+1; icount++)
	{
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
			// corners
			if (vct->getXleftYleft_neighbor()==-1){
				alpha=(1.0-value)/2;
				temp[1][1]=value*vector[is][1][1][0]+alpha*(vector[is][2][1][0]+vector[is][1][2][0]);
			}
			if (vct->getXleftYright_neighbor()==-1){
				alpha=(1.0-value)/2;
				temp[1][ny-2]=value*vector[is][1][ny-2][0]+alpha*(vector[is][2][ny-1][0]+vector[is][1][ny-3][0]);
			}
			if (vct->getXrightYleft_neighbor()==-1){
				alpha=(1.0-value)/2;
				temp[nx-2][1]=value*vector[is][nx-2][1][0]+alpha*(vector[is][nx-3][1][0]+vector[is][nx-2][2][0]);
			}
			if (vct->getXrightYright_neighbor()==-1){
				alpha=(1.0-value)/2;
				temp[nx-2][ny-2]=value*vector[is][nx-2][ny-2][0]+alpha*(vector[is][nx-3][ny-2][0]+vector[is][nx-2][ny-3][0]);
			}
			// copy temp in vector
			for (int i=1; i < nx-1;i++)
				for (int j=1; j < ny-1;j++)
					vector[is][i][j][0] = temp[i][j];
			// communicate
			if (type == 0)
				communicateCenter(nx, ny, vector, is, vct);
			else
				communicateNode(nx, ny, vector, is, vct);

			delArr(temp,nx);
		} // end of if
	}
}



/** Perfect conductor boundary conditions for magnetic field */
inline  void EMfields::BperfectConductorLeft(int dir, Grid *grid, VirtualTopology *vct){
	switch(dir){
		case 0: // boundary condition on X-DIRECTION left
			for (int i=0; i <  nyc;i++){
				Bxc[0][i][0] = 0.0;
				Byc[0][i][0] = Byc[1][i][0];
				Bzc[0][i][0] = Bzc[1][i][0];
			}
			break;
		case 1: // boundary condition on Y-DIRECTION left
			for (int i=0; i < nxc;i++){
				Bxc[i][0][0] = B0x*tanh((grid->getYC(i,0,0) - Ly/2)/delta);
				Byc[i][0][0] = B0y;
				Bzc[i][0][0] = B0z*tanh((grid->getYC(i,0,0) - Ly/2)/delta);;
				Bzc[i][1][0] = B0z*tanh((grid->getYC(i,0,0) - Ly/2)/delta);
				Bzc[i][2][0] = B0z*tanh((grid->getYC(i,0,0) - Ly/2)/delta);
			}
			break;
	}
}
/** Perfect conductor boundary conditions for magnetic field  */
inline  void EMfields::BperfectConductorRight(int dir,Grid *grid,VirtualTopology *vct){
	switch(dir){
		case 0: // boundary condition on X-DIRECTION left
			for (int i=0; i <  nyc;i++){
				Bxc[nxc-1][i][0] = 0.0;
				Byc[nxc-1][i][0] = Byc[nxc-2][i][0];
				Bzc[nxc-1][i][0] = Bzc[nxc-2][i][0];
			}
		case 1: // boundary condition on Y-DIRECTION left
			for (int i=0; i < nxc;i++){
				Bxc[i][nyc-1][0] = B0x*tanh((grid->getYC(i,nyc-1,0) - Ly/2)/delta);
				Byc[i][nyc-1][0] = B0y;
				Bzc[i][nyc-1][0] = B0z*tanh((grid->getYC(i,nyc-1,0) - Ly/2)/delta);
				Bzc[i][nyc-2][0] = B0z*tanh((grid->getYC(i,nyc-1,0) - Ly/2)/delta);
				Bzc[i][nyc-3][0] = B0z*tanh((grid->getYC(i,nyc-1,0) - Ly/2)/delta);
			}
			break;
	}
}


/** Perfect conductor boundary conditions: LEFT wall */
inline  void EMfields::perfectConductorLeft(double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid){
	   double* susxx;
	   double* susxy;
	   double* susxz;
	   double* susyx;
	   double* susyy;
	   double* susyz;
	   double* suszx;
	   double* suszy;
	   double* suszz;
	   double bavg2;
	   switch(dir){
		   case 0: // boundary condition on X-DIRECTION left
			   for (int i=0; i <  nyn;i++){
				   vectorX[0][i][0] =  rhon[1][i][0]*dy;
				   vectorY[0][i][0] =  0.0;
				   vectorZ[0][i][0] =  0.0;
			   }
			   break;
		   case 1: // boundary condition on Y-DIRECTION left
			   susxx = new double[nxn];
			   susxy = new double[nxn];
			   susxz = new double[nxn];
			   susyx = new double[nxn];
			   susyy = new double[nxn];
			   susyz = new double[nxn];
			   suszx = new double[nxn];
			   suszy = new double[nxn];
			   suszz = new double[nxn];
			   sustensorLeft(susxx,susxy,susxz,susyx,susyy,susyz,suszx,suszy,suszz,dir);
			   for (int i=0; i < nxn;i++){
				   bavg2 = Bxn[i][1][0]*Bxn[i][1][0] + Byn[i][1][0]*Byn[i][1][0] + Bzn[i][1][0]*Bzn[i][1][0];
				   vectorX[i][0][0] = 0.0;
				   vectorZ[i][0][0] = 0.0;
				   vectorY[i][0][0] = (Ey[i][1][0] - susxy[i]*vectorX[i][0][0] -susyz[i]*vectorZ[i][0][0] -Jyh[i][1][0]*dt*th)/susyy[i];

			   }
				   break;
	   }
	   delete[] susxx;
	   delete[] susxy;
	   delete[] susxz;
	   delete[] susyx;
	   delete[] susyy;
	   delete[] susyz;
	   delete[] suszx;
	   delete[] suszy;
	   delete[] suszz;

}
/** Open boundary conditions: LEFT wall */
inline  void EMfields::openLeft( double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid){

}
/** Open boundary conditions: LEFT wall */
inline  void EMfields::openNoPlasmaLeft(double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid){


}
/** Perfect conductor boundary conditions: RIGHT wall */
inline  void EMfields::perfectConductorRight(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){
	   double* susxx;
	   double* susxy;
	   double* susxz;
	   double* susyx;
	   double* susyy;
	   double* susyz;
	   double* suszx;
	   double* suszy;
	   double* suszz;
	   double bavg2;
	   switch(dir){
		   case 0: // boundary condition on X-DIRECTION RIGHT
			   for (int i=0; i <  nyn;i++){
				   vectorX[nxn-1][i][0] = -rhon[nxn-2][i][0]*dy;
				   vectorY[nxn-1][i][0] = 0.0;
				   vectorZ[nxn-1][i][0] = 0.0;
			   }
			   break;
		   case 1: // boundary condition on Y-DIRECTION RIGHT
			   susxx = new double[nxn];
			   susxy = new double[nxn];
			   susxz = new double[nxn];
			   susyx = new double[nxn];
			   susyy = new double[nxn];
			   susyz = new double[nxn];
			   suszx = new double[nxn];
			   suszy = new double[nxn];
			   suszz = new double[nxn];
			   sustensorRight(susxx,susxy,susxz,susyx,susyy,susyz,suszx,suszy,suszz,dir);
			   for (int i=0; i < nxn;i++){
				   bavg2 = Bxn[i][nyn-2][0]*Bxn[i][nyn-2][0] + Byn[i][nyn-2][0]*Byn[i][nyn-2][0] + Bzn[i][nyn-2][0]*Bzn[i][nyn-2][0];
				   vectorX[i][nyn-1][0] = 0.0;
				   vectorZ[i][nyn-1][0] = 0.0;
				   vectorY[i][nyn-1][0] = (Ey[i][nyn-2][0] - susxy[i]*vectorX[i][nyn-1][0] -susyz[i]*vectorZ[i][nyn-1][0] -Jyh[i][nyn-2][0]*dt*th)/susyy[i];

			   }
				   break;
	   }
	   delete[] susxx;
	   delete[] susxy;
	   delete[] susxz;
	   delete[] susyx;
	   delete[] susyy;
	   delete[] susyz;
	   delete[] suszx;
	   delete[] suszy;
	   delete[] suszz;
}
inline  void EMfields::perfectConductorLeftImage(double ***imageX, double ***imageY, double ***imageZ, double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid){
	double* susxx;
	double* susxy;
	double* susxz;
	double* susyx;
	double* susyy;
	double* susyz;
	double* suszx;
	double* suszy;
	double* suszz;
	double bavg2;
	   switch(dir){
		   case 0: // boundary condition on X-DIRECTION left
			   susxx = new double[nyn];
			   susxy = new double[nyn];
			   susxz = new double[nyn];
			   susyx = new double[nyn];
			   susyy = new double[nyn];
			   susyz = new double[nyn];
			   suszx = new double[nyn];
			   suszy = new double[nyn];
			   suszz = new double[nyn];

			   sustensorLeft(susxx,susxy,susxz,susyx,susyy,susyz,suszx,suszy,suszz,dir);
			   for (int i=1; i <  nyn-1;i++){
				   bavg2 = Bxn[1][i][0]*Bxn[1][i][0] + Byn[1][i][0]*Byn[1][i][0] + Bzn[1][i][0]*Bzn[1][i][0];
				   imageX[1][i][0] = (susxx[i]*vectorX[1][i][0] + rhoh[1][i][0]*dx)/1000.0;
				   imageY[1][i][0] = vectorY[1][i][0];
				   imageZ[1][i][0] = vectorZ[1][i][0] - Vinj*bavg2/(fabs(Byn[1][i][0]) + 1e-10);
			   }
				   break;
		   case 1: // boundary condition on Y-DIRECTION left
			   susxx = new double[nxn];
			   susxy = new double[nxn];
			   susxz = new double[nxn];
			   susyx = new double[nxn];
			   susyy = new double[nxn];
			   susyz = new double[nxn];
			   suszx = new double[nxn];
			   suszy = new double[nxn];
			   suszz = new double[nxn];
			   sustensorLeft(susxx,susxy,susxz,susyx,susyy,susyz,suszx,suszy,suszz,dir);
			   for (int i=1; i < nxn-1;i++){
				   bavg2 = Bxn[i][1][0]*Bxn[i][1][0] + Byn[i][1][0]*Byn[i][1][0] + Bzn[i][1][0]*Bzn[i][1][0];
				   imageX[i][1][0] = vectorX[i][1][0];
				   //imageY[i][1][0] = (susyy[i]*vectorY[i][1][0] + rhoh[i][1][0]*dy)/1000.0;
				   imageZ[i][1][0] = vectorZ[i][1][0] - Vinj*bavg2/(fabs(Bxn[i][1][0]) + 1e-10);
				   imageY[i][1][0] = vectorY[i][1][0] - (Ey[i][1][0] - susxy[i]*vectorX[i][1][0] -susyz[i]*vectorZ[i][1][0] -Jyh[i][1][0]*dt*th*FourPI)/susyy[i];

			   }

				   break;
	   }
	   delete[] susxx;
	   delete[] susxy;
	   delete[] susxz;
	   delete[] susyx;
	   delete[] susyy;
	   delete[] susyz;
	   delete[] suszx;
	   delete[] suszy;
	   delete[] suszz;
}

/** Perfect conductor boundary conditions in Image Maxwell */
inline  void EMfields::perfectConductorRightImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){
	double* susxx;
	double* susxy;
	double* susxz;
	double* susyx;
	double* susyy;
	double* susyz;
	double* suszx;
	double* suszy;
	double* suszz;
	double bavg2;
	   switch(dir){
		   case 0: // boundary condition on X-DIRECTION left
			   susxx = new double[nyn];
			   susxy = new double[nyn];
			   susxz = new double[nyn];
			   susyx = new double[nyn];
			   susyy = new double[nyn];
			   susyz = new double[nyn];
			   suszx = new double[nyn];
			   suszy = new double[nyn];
			   suszz = new double[nyn];
			   sustensorRight(susxx,susxy,susxz,susyx,susyy,susyz,suszx,suszy,suszz,dir);
			   for (int i=1; i <  nyn-1;i++){
				   bavg2 = Bxn[nxn-2][i][0]*Bxn[nxn-2][i][0] + Byn[nxn-2][i][0]*Byn[nxn-2][i][0] + Bzn[nxn-2][i][0]*Bzn[nxn-2][i][0];
				   imageX[nxn-2][i][0] = (susxx[i]*vectorX[nxn-2][i][0] + rhoh[nxc-2][i][0]*dx)/1000.0;
				   imageY[nxn-2][i][0] = vectorY[nxn-2][i][0];
				   imageZ[nxn-2][i][0] = vectorZ[nxn-2][i][0] - Vinj*bavg2/(fabs(Byn[nxn-2][i][0]) + 1e-10);
			   }
				   break;
		   case 1: // boundary condition on Y-DIRECTION RIGHT
			   susxx = new double[nxn];
			   susxy = new double[nxn];
			   susxz = new double[nxn];
			   susyx = new double[nxn];
			   susyy = new double[nxn];
			   susyz = new double[nxn];
			   suszx = new double[nxn];
			   suszy = new double[nxn];
			   suszz = new double[nxn];
			   sustensorRight(susxx,susxy,susxz,susyx,susyy,susyz,suszx,suszy,suszz,dir);
			   for (int i=1; i < nxn-1;i++){
				   bavg2 = Bxn[i][nyn-2][0]*Bxn[i][nyn-2][0] + Byn[i][nyn-2][0]*Byn[i][nyn-2][0] + Bzn[i][nyn-2][0]*Bzn[i][nyn-2][0];
				   imageX[i][nyn-2][0] = vectorX[i][nyn-2][0];
				   //imageY[i][nyn-2][0] = (susyy[i]*vectorY[i][nyn-2][0] + rhoh[i][nyc-2][0]*dy)/1000.0;
				   imageY[i][nyn-2][0] = vectorY[i][nyn-2][0] - (Ey[i][nyn-2][0] - susxy[i]*vectorX[i][nyn-2][0] -susyz[i]*vectorZ[i][nyn-2][0] -Jyh[i][nyn-2][0]*dt*th*FourPI)/susyy[i];
				   imageZ[i][nyn-2][0] = vectorZ[i][nyn-2][0] - Vinj*bavg2/(fabs(Bxn[i][nyn-2][0]) + 1e-10);
			   }

				   break;
	   }
	   delete[] susxx;
	   delete[] susxy;
	   delete[] susxz;
	   delete[] susyx;
	   delete[] susyy;
	   delete[] susyz;
	   delete[] suszx;
	   delete[] suszy;
	   delete[] suszz;


}

/** Open boundary conditions: RIGHT wall */
inline  void EMfields::openRight(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){

}
/** Open boundary conditions with no plasma injection: RIGHT wall */
inline  void EMfields::openNoPlasmaRight(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){

}

/** Magnetic mirror boundary conditions: LEFT WALL */
inline  void EMfields::magneticMirrorLeft(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){


}
/** Magnetic mirror boundary conditions: RIGHT WALL */
inline  void EMfields::magneticMirrorRight(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){

}


/**
* Perfect conductor boundary conditions for source: LEFT WALL
 *
 *
 */
inline  void EMfields::perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
	switch(dir){
		case 0: // boundary condition on X-DIRECTION LEFT
			for (int i=1; i <  nyn-1;i++){
				vectorX[1][i][0] = 0.0;
				vectorY[1][i][0] = 0.0;
				vectorZ[1][i][0] = 0.0;
			}
			break;
		case 1: // boundary condition on Y-DIRECTION LEFT
			for (int i=1; i <  nxn-1 ;i++){
				vectorX[i][1][0] = 0.0;
				vectorY[i][1][0] = 0.0;
				vectorZ[i][1][0] = 0.0;
			}

			break;
	}
}
/**
* open conductor boundary conditions for source: LEFT WALL
 *
 *
 */
inline  void EMfields::openLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
	switch(dir){
		case 0: // boundary condition on X-DIRECTION LEFT
			for (int i=1; i <  nyn-1;i++){
				vectorX[1][i][0] = 0.0;
				vectorY[1][i][0] = 0.0;
				vectorZ[1][i][0] = 0.0;
			}
			break;
		case 1: // boundary condition on Y-DIRECTION LEFT
			for (int i=1; i <  nxn-1 ;i++){
				vectorX[i][1][0] = 0.0;
				vectorY[i][1][0] = 0.0;
				vectorZ[i][1][0] = 0.0;
			}

			break;
	}
}

/**
* Perfect conductor boundary conditions for source: LEFT WALL
 *
 *
 */
inline  void EMfields::perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
	switch(dir){
		case 0:  // boundary condition on X-DIRECTION RIGHT
			for (int i=1; i < nyn-1;i++){
				vectorX[nxn-2][i][0] = 0.0;
				vectorY[nxn-2][i][0] = 0.0;
				vectorZ[nxn-2][i][0] = 0.0;
			}
			break;
		case 1: // boundary condition on Y-DIRECTION RIGHT
			for (int i=1; i <  nxn-1;i++){
				vectorX[i][nyn-2][0] = 0.0;
				vectorY[i][nyn-2][0] = 0.0;
				vectorZ[i][nyn-2][0] = 0.0;
			}
			break;
	}

}
inline  void EMfields::openRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
	switch(dir){
		case 0:  // boundary condition on X-DIRECTION RIGHT
			for (int i=1; i < nyn-1;i++){
				vectorX[nxn-2][i][0] = 0.0;
				vectorY[nxn-2][i][0] = 0.0;
				vectorZ[nxn-2][i][0] = 0.0;
			}
			break;
		case 1: // boundary condition on Y-DIRECTION RIGHT
			for (int i=1; i <  nxn-1;i++){
				vectorX[i][nyn-2][0] = 0.0;
				vectorY[i][nyn-2][0] = 0.0;
				vectorZ[i][nyn-2][0] = 0.0;
			}
			break;
	}

}
/** Magnetic mirror boundary conditions for source: LEFT wall */
inline  void EMfields::magneticMirrorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
	switch(dir){
		case 0:
			for (int i=1; i <  nyn-1;i++){
				vectorX[1][i][0] = 0.0;
				vectorY[1][i][0] = 0.0;
				vectorZ[1][i][0] = 0.0;
			}
			break;
		case 1:
			for (int i=1; i < nxn-1;i++){
				vectorX[i][1][0] = 0.0;
				vectorY[i][1][0] = 0.0;
				vectorZ[i][1][0] = 0.0;
			}
			break;
	}
}
/** Magnetic mirror boundary conditions for source: RIGHT wall */
inline  void EMfields::magneticMirrorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir){
	switch(dir){
		case 0:
			for (int i=1; i <  nyn-1;i++){
				vectorX[nxn-2][i][0] = 0.0;
				vectorY[nxn-2][i][0] = 0.0;
				vectorZ[nxn-2][i][0] = 0.0;
			}
			break;
		case 1:
			for (int i=1; i <  nxn-1;i++){
				vectorX[i][nyn-2][0] = 0.0;
				vectorY[i][nyn-2][0] = 0.0;
				vectorZ[i][nyn-2][0] = 0.0;
			}
			break;
    }
}




/** BC for PHI  left side*/
inline void EMfields::bcPHI_Left(double ***im,int dir){
	switch(dir){
		case 0:
			for (int i=0; i <  nyc;i++){
				im[0][i][0] = 0.0;

			}
			break;
		case 1:
			for (int i=1; i < nxc;i++){
				im[i][0][0] = 0.0;

			}
			break;
	}
}
/** BC for PHI  right side*/
inline void EMfields::bcPHI_Right(double ***im,int dir){
	switch(dir){
		case 0:
			for (int i=0; i <  nyc;i++){
				im[nxc-1][i][0] = 0.0;

			}
			break;
		case 1:
			for (int i=0; i < nxc;i++){
				im[i][nyc-1][0] = 0.0;
			}
			break;
	}
}


/** get Potential array ***/
inline double*** EMfields::getPHI() {return(PHI);}
/** get Ex(X,Y,Z)  */
inline double &EMfields::getEx(int indexX, int indexY, int indexZ) const{
	return(Ex[indexX][indexY][0]);}
/** get Electric field  component X array*/
inline double*** EMfields::getEx() {return(Ex);}
/** get Ey(X,Y,Z)  */
inline double &EMfields::getEy(int indexX, int indexY, int indexZ) const{
	return(Ey[indexX][indexY][0]);}
/** get Electric field  component Y array*/
inline double*** EMfields::getEy() {return(Ey);}
/** get Ez(X,Y,Z)  */
inline double &EMfields::getEz(int indexX, int indexY, int indexZ) const{
	return(Ez[indexX][indexY][0]);}
/** get Electric field  component Z array*/
inline double*** EMfields::getEz() {return(Ez);}
/** get Bx(X,Y,Z)  */
inline double &EMfields::getBx(int indexX, int indexY, int indexZ) const{
	return(Bxn[indexX][indexY][0]);}
/** get Magnetic Field component X array*/
inline double*** EMfields::getBx() {return(Bxn);}
/**  get By(X,Y,Z) */
inline double &EMfields::getBy(int indexX, int indexY, int indexZ) const{
	return(Byn[indexX][indexY][0]);}
/** get Magnetic Field component Y array*/
inline double*** EMfields::getBy() {return(Byn);}
/**  get Bz(X,Y,Z) */
inline double &EMfields::getBz(int indexX, int indexY, int indexZ) const{
	return(Bzn[indexX][indexY][0]);}
/** get Magnetic Field component Z array*/
inline double*** EMfields::getBz() {return(Bzn);}
/** get rhoc(X,Y,Z) */
inline double &EMfields::getRHOc(int indexX, int indexY, int indexZ) const{
	return(rhoc[indexX][indexY][0]);}
inline double*** EMfields::getRHOc() {return(rhoc);}
/** get density on node(indexX,indexY,indexZ)  */
inline double &EMfields::getRHOn(int indexX, int indexY, int indexZ) const{
	return(rhon[indexX][indexY][0]);}
/** get density array defined on nodes*/
inline double*** EMfields::getRHOn() {return(rhon);}
/** get rhos(X,Y,Z) : density for species*/
inline double &EMfields::getRHOns(int indexX, int indexY, int indexZ, int is) const{
	return(rhons[is][indexX][indexY][0]);}
/** SPECIES: get density array defined on center cells  */
inline double &EMfields::getRHOcs(int indexX, int indexY, int indexZ,int is) const{
	return(rhocs[is][indexX][indexY][0]);}
/** get density array defined on nodes*/
inline double**** EMfields::getRHOns() {return(rhons);}
/** SPECIES: get pressure tensor component XX defined on nodes */
inline double**** EMfields::getpXXsn() {return(pXXsn);}
/** SPECIES: get pressure tensor component XY defined on nodes */
inline double**** EMfields::getpXYsn() {return(pXYsn);}
/** SPECIES: get pressure tensor component XZ defined on nodes */
inline double**** EMfields::getpXZsn() {return(pXZsn);}
/** SPECIES: get pressure tensor component YY defined on nodes */
inline double**** EMfields::getpYYsn() {return(pYYsn);}
/** SPECIES: get pressure tensor component YZ defined on nodes */
inline double**** EMfields::getpYZsn() {return(pYZsn);}
/** SPECIES: get pressure tensor component ZZ defined on nodes */
inline double**** EMfields::getpZZsn() {return(pZZsn);}
/** get current -Direction X */
inline double &EMfields::getJx(int indexX, int indexY,int indexZ) const{
	return(Jx[indexX][indexY][0]);}
/** get current array X component **/
inline double*** EMfields::getJx() {return(Jx);}
/** get current -Direction Y */
inline double &EMfields::getJy(int indexX, int indexY, int indexZ) const{
	return(Jy[indexX][indexY][0]);}
/** get current array Y component **/
inline double*** EMfields::getJy() {return(Jy);}
/** get current -Direction Z */
inline double &EMfields::getJz(int indexX, int indexY, int indexZ) const{
	return(Jz[indexX][indexY][0]);}
/** get current array Z component **/
inline double*** EMfields::getJz() {return(Jz);}
/**SPECIES: get current array X component */
inline double**** EMfields::getJxs() {return(Jxs);}
/**SPECIES: get current array Y component */
inline double**** EMfields::getJys() {return(Jys);}
/**SPECIES: get current array Z component */
inline double**** EMfields::getJzs() {return(Jzs);}

/** Print info about electromagnetic field */
inline void EMfields::print(void) const{


}
/** constructor */
inline EMfields::EMfields(CollectiveIO *col,Grid *grid){
	nxc = grid->getNXC();
	nxn = grid->getNXN();
	nyc = grid->getNYC();
	nyn = grid->getNYN();
	dx = grid->getDX();
	dy = grid ->getDY();
	invVOL = grid->getInvVOL();
	xStart = grid->getXstart();
	xEnd   = grid->getXend();
	yStart = grid->getYstart();
	yEnd   = grid->getYend();
	Lx = col->getLx();
	Ly = col->getLy();
	ns  = col->getNs();
	c = col->getC();
	dt = col->getDt();
	Smooth = col->getSmooth();
	Nvolte = col ->getNvolte();
	th = col->getTh();
	delt = c*th*dt;

    Vinj = col->getVinj();
	// FLAG ON POISSON CORRECTION
	// make always the divergence cleaning
	PoissonCorrection = true;


	qom = new double[ns];
	for (int i=0; i < ns;i++)
		qom[i] = col->getQOM(i);
	// GEM challenge parameters
	B0x = col->getB0x();
	B0y = col->getB0y();
	B0z = col->getB0z();
	delta = col->getDelta();
	// get the density background for the gem Challange
	rhoINIT = new double[ns];
	DriftSpecies = new bool[ns];
	for (int i=0; i < ns;i++){
	    rhoINIT[i] = col->getRHOinit(i);
	    if ((fabs(col->getW0(i))!=0) || (fabs(col->getU0(i))!=0)) // GEM and LHDI
			DriftSpecies[i] = true;
		else
			DriftSpecies[i] = false;
	}
	// boundary conditions: PHI and EM fields
	bcPHIfaceXright  = col->getBcPHIfaceXright();
	bcPHIfaceXleft   = col->getBcPHIfaceXleft();
	bcPHIfaceYright  = col->getBcPHIfaceYright();
	bcPHIfaceYleft   = col->getBcPHIfaceYleft();

	bcEMfaceXright   = col->getBcEMfaceXright();
	bcEMfaceXleft    = col->getBcEMfaceXleft();
	bcEMfaceYright   = col->getBcEMfaceYright();
	bcEMfaceYleft    = col->getBcEMfaceYleft();
	// set fourpi
	FourPI =16*atan(1.0);
	// Restart
	restart1 = col->getRestart_status();
	RestartDirName = col->getRestartDirName();
	/** get solvers tolerance */
	CGtol = col->getCGtol();
	GMREStol = col->getGMREStol();
	// arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
	Ex = newArr3(double,nxn,nyn,1);
	Ey = newArr3(double,nxn,nyn,1);
	Ez = newArr3(double,nxn,nyn,1);
	Exth = newArr3(double,nxn,nyn,1);
	Eyth = newArr3(double,nxn,nyn,1);
	Ezth = newArr3(double,nxn,nyn,1);
	Bxn   = newArr3(double,nxn,nyn,1);
	Byn   = newArr3(double,nxn,nyn,1);
	Bzn   = newArr3(double,nxn,nyn,1);
	rhon  = newArr3(double,nxn,nyn,1);
	Jx    = newArr3(double,nxn,nyn,1);
	Jy    = newArr3(double,nxn,nyn,1);
	Jz    = newArr3(double,nxn,nyn,1);
	Jxh   = newArr3(double,nxn,nyn,1);
	Jyh   = newArr3(double,nxn,nyn,1);
	Jzh   = newArr3(double,nxn,nyn,1);
	// involving species
	rhons = newArr4(double,ns,nxn,nyn,1);
	rhocs = newArr4(double,ns,nxc,nyc,1);
	Jxs   = newArr4(double,ns,nxn,nyn,1);
	Jys   = newArr4(double,ns,nxn,nyn,1);
	Jzs   = newArr4(double,ns,nxn,nyn,1);
	pXXsn  = newArr4(double,ns,nxn,nyn,1);
	pXYsn  = newArr4(double,ns,nxn,nyn,1);
	pXZsn  = newArr4(double,ns,nxn,nyn,1);
	pYYsn  = newArr4(double,ns,nxn,nyn,1);
	pYZsn  = newArr4(double,ns,nxn,nyn,1);
	pZZsn  = newArr4(double,ns,nxn,nyn,1);
	// arrays allocation: central points ---> the first central point has index 1, the last has index nxn-2!
	PHI    = newArr3(double,nxc,nyc,1);
	Bxc    = newArr3(double,nxc,nyc,1);
	Byc    = newArr3(double,nxc,nyc,1);
	Bzc    = newArr3(double,nxc,nyc,1);
	rhoc   = newArr3(double,nxc,nyc,1);
	rhoh   = newArr3(double,nxc,nyc,1);

	// temporary arrays
	tempXC = newArr3(double,nxc,nyc,1);
	tempYC = newArr3(double,nxc,nyc,1);
	tempZC = newArr3(double,nxc,nyc,1);
	tempXN = newArr3(double,nxn,nyn,1);
	tempYN = newArr3(double,nxn,nyn,1);
	tempZN = newArr3(double,nxn,nyn,1);
	tempC  = newArr3(double,nxc,nyc,1);
	tempX  = newArr3(double,nxn,nyn,1);
	tempY  = newArr3(double,nxn,nyn,1);
	tempZ  = newArr3(double,nxn,nyn,1);
	temp2X = newArr3(double,nxn,nyn,1);
	temp2Y = newArr3(double,nxn,nyn,1);
	temp2Z = newArr3(double,nxn,nyn,1);
	imageX = newArr3(double,nxn,nyn,1);
	imageY = newArr3(double,nxn,nyn,1);
	imageZ = newArr3(double,nxn,nyn,1);
	Dx = newArr3(double,nxn,nyn,1);
	Dy = newArr3(double,nxn,nyn,1);
	Dz = newArr3(double,nxn,nyn,1);
	vectX  = newArr3(double,nxn,nyn,1);
	vectY  = newArr3(double,nxn,nyn,1);
	vectZ  = newArr3(double,nxn,nyn,1);
	divC  = newArr3(double,nxc,nyc,1);

}
/** destructor: deallocate arrays*/
inline EMfields::~EMfields(){
	// nodes
	delArr3(Ex,nxn,nyn);
	delArr3(Ey,nxn,nyn);
	delArr3(Ez,nxn,nyn);
	delArr3(Exth,nxn,nyn);
	delArr3(Eyth,nxn,nyn);
	delArr3(Ezth,nxn,nyn);
	delArr3(Bxn,nxn,nyn);
	delArr3(Byn,nxn,nyn);
	delArr3(Bzn,nxn,nyn);
	delArr3(rhon,nxn,nyn);
	delArr3(Jx,nxn,nyn);
	delArr3(Jy,nxn,nyn);
	delArr3(Jz,nxn,nyn);
	delArr3(Jxh,nxn,nyn);
	delArr3(Jyh,nxn,nyn);
	delArr3(Jzh,nxn,nyn);
	// nodes and species
	delArr4(rhons,ns,nxn,nyn);
	delArr4(rhocs,ns,nxc,nyc);
	delArr4(Jxs,ns,nxn,nyn);
	delArr4(Jys,ns,nxn,nyn);
	delArr4(Jzs,ns,nxn,nyn);
	delArr4(pXXsn,ns,nxn,nyn);
	delArr4(pXYsn,ns,nxn,nyn);
	delArr4(pXZsn,ns,nxn,nyn);
	delArr4(pYYsn,ns,nxn,nyn);
	delArr4(pYZsn,ns,nxn,nyn);
	delArr4(pZZsn,ns,nxn,nyn);
	// central points
	delArr3(PHI,nxc,nyc);
	delArr3(Bxc,nxc,nyc);
	delArr3(Byc,nxc,nyc);
	delArr3(Bzc,nxc,nyc);
	delArr3(rhoc,nxc,nyc);
	delArr3(rhoh,nxc,nyc);
	// various stuff needs to be deallocated too
	delArr3(tempXC,nxc,nyc);
	delArr3(tempYC,nxc,nyc);
	delArr3(tempZC,nxc,nyc);
	delArr3(tempXN,nxn,nyn);
	delArr3(tempYN,nxn,nyn);
	delArr3(tempZN,nxn,nyn);
	delArr3(tempC,nxc,nyc);
	delArr3(tempX,nxn,nyn);
	delArr3(tempY,nxn,nyn);
	delArr3(tempZ,nxn,nyn);
	delArr3(temp2X,nxn,nyn);
	delArr3(temp2Y,nxn,nyn);
	delArr3(temp2Z,nxn,nyn);
	delArr3(imageX,nxn,nyn);
	delArr3(imageY,nxn,nyn);
	delArr3(imageZ,nxn,nyn);
	delArr3(Dx,nxn,nyn);
	delArr3(Dy,nxn,nyn);
	delArr3(Dz,nxn,nyn);
	delArr3(vectX,nxn,nyn);
	delArr3(vectY,nxn,nyn);
	delArr3(vectZ,nxn,nyn);
	delArr3(divC,nxc,nyc);


}
#endif

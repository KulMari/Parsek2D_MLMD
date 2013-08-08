/******************************************************************************************
EMfields.h  -  Electromagnetic Field with 3 components(x,y,z) defined on a 2-D grid. Solved
using the implicit Maxwell solver.
-------------------
developers: Stefano Markidis,  Giovanni Lapenta, Arnaud Beck, Maria Elena Innocenti
********************************************************************************************/

#ifndef EMfields_H
#define EMfields_H


#include <iostream>
#include <sstream>
#include <fstream>

#include <math.h>
#include <mpi.h>
#include <float.h>


#include "../utility/Alloc.h"
#include "../mathlib/Basic.h"
#include "../utility/TransArraySpace.h"
#include "../solvers/CG.h"
#include "../solvers/GMRES_new2.h"
#include "hdf5.h"

#define PI 4*atan(1)

using std::cout;
using std::cerr;
using std::endl;



/**
*  Electromagnetic fields and sources defined for each local grid, and for an implicit Maxwell's solver
 *
 * @date Fri Jun 4 2007 KUL
 * @author Stefano Markidis, Giovanni Lapenta, Arnaud Beck
 * @version 2.0
 *
 */


class EMfields : public Field{
public:
	/** constructor */
	EMfields(CollectiveIO *col,Grid *grid, VCtopology *vct);
	/** destructor */
	~EMfields();

	/** initialize the electromagnetic fields with constant values */
	void init(VirtualTopology *vct, Grid *grid);
	/** initialize the magnetic field as a dipole */
	void initLightwave(VirtualTopology *vct, Grid *grid);
	/** initialize with reconnection system */
	void initGEM(VirtualTopology *vct, Grid *grid);
	/** initialize with double Harris system */
	void initDoubleHarris(VirtualTopology *vct, Grid *grid);
	/** initialize with Force Free system (JxB=0) */
	void initForceFree(VirtualTopology *vct, Grid *grid);
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
        /** Initialize the weights for interpolation of BC between grids **/
        int initWeightBC(VirtualTopology *vct, Grid *grid, CollectiveIO* col );
        /** Initialize the weights for projections of refined fields between grids **/
        int initWeightProj(VirtualTopology *vct, Grid *grid, CollectiveIO* col );
        /** Output the value of ghost nodes **/
        void outputghost(VirtualTopology *vct, CollectiveIO *col, int cycle);

	/** Calculate Electric field using the implicit Maxwell solver */
	void calculateField(Grid *grid, VirtualTopology *vct);
        /** Apply electric field boundary conditions to finer grids*/
        void electricfieldBC(double ***componentX,double ***componentY, double ***componentZ, Grid *grid, VirtualTopology *vct); 
        /** Apply magnetic field boundary conditions to finer grids*/
        void magneticfieldBC(double ***component, Grid *grid, double ***from, VirtualTopology *vct);
	/** Send the boundary conditions to finer level */
	void sendBC(Grid *grid, VirtualTopology *vct);
	/** Send the refined fields to coarser level */
	void sendProjection(Grid *grid, VirtualTopology *vct);
	/** Receive the boundary conditions from coarser level*/
	void receiveBC(Grid *grid, VirtualTopology *vct, CollectiveIO *col);
	/** Receive the refined fields from finer level*/
	void receiveProjection(CollectiveIO *col, Grid *grid, VirtualTopology *vct);
	/** Image of Poisson Solver (for SOLVER)*/
	void PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct);
	/** Image of Maxwell Solver (for Solver) */
	void MaxwellImage(double *im, double *vector, Grid *grid,VirtualTopology *vct);
	/** Maxwell source term (for SOLVER) */
	void MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct);
	/** Calculate Magnetic field with the implicit solver: calculate B defined on nodes
		With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
	void calculateB(Grid *grid, VirtualTopology *vct);
	void calculateB_afterProj(Grid *grid, VirtualTopology *vct);
	/** Calculate the three components of Pi(implicit pressure) cross image vector */
	void PIdot(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid *grid);
	// PIdot also on ghost node
	void PIdot_alsoGN(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid *grid);
	/** Calculate the three components of mu (implicit permeattivity) cross image vector */
	void MUdot(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid *grid);
	/** Calculate the three components of mu (implicit permeattivity) cross image vector including ghost cells*/
	void MUdot_plusghost(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid *grid);
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
	/** Smoothing for the electric field (never called)**/
	void smoothE(int nvolte, double value ,double ***vector,bool type, Grid *grid, VirtualTopology *vct);
	/** SPECIES: Smoothing after the interpolation for species fields**/
	void smooth(int nvolte, double value ,double ****vector,int is, bool type, Grid *grid, VirtualTopology *vct);

	/** communicate ghost for grid -> Particles interpolation */
	void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct);
	/** communicate ghost for grid -> Particles interpolation for OS */
        void communicateGhostP2GOS(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct);
	/** add an amount to the specified moment at node X,Y,Z, species is**/  
	void addMoment(double **** moment, double ***weight, int X, int Y,int Z, int is);
	void addMoment_OS(double **** moment, double ***weight, int X, int Y, int Z, int is, Grid *g);  
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

	/** add an amount of charge density to charge density field at node X,Y,Z */
        void addRho_OS(double ***weight, int X, int Y,int Z, int is, Grid *g);
        /** add an amount of current density - direction X to current density field at node X,Y,Z */
        void addJx_OS(double ***weight, int X, int Y, int Z,int is, Grid *g);
        /** add an amount of current density - direction Y to current density field at node X,Y,Z */
        void addJy_OS(double ***weight, int X, int Y, int Z, int is, Grid *g);
	/** add an amount of current density - direction Z to current density field at node X,Y,Z */
        void addJz_OS(double ***weight, int X, int Y, int Z, int is, Grid *g);

	/** add an amount of pressure density - direction XX to current density field at node X,Y,Z */
        void addPxx_OS(double ***weight, int X, int Y, int Z, int is, Grid *g);
        /** add an amount of pressure density - direction XY to current density field at node X,Y,Z */
        void addPxy_OS(double ***weight, int X, int Y, int Z,int is, Grid *g);
	/** add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
	void addPxz_OS(double ***weight, int X, int Y, int Z,int is, Grid *g);
        /** add an amount of pressure density - direction YY to current density field at node X,Y,Z */
        void addPyy_OS(double ***weight, int X, int Y, int Z, int is, Grid *g);
        /** add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
	void addPyz_OS(double ***weight, int X, int Y, int Z, int is, Grid *g);
	/** add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
        void addPzz_OS(double ***weight, int X, int Y, int Z,int is, Grid *g);
	/** adjust densities on boundaries that are not periodic */
	void adjustNonPeriodicDensities(VirtualTopology *vct, int cycle);

	/** Perfect conductor boundary conditions LEFT wall */
	void perfectConductorLeftImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Perfect conductor boundary conditions RIGHT wall */
	void perfectConductorRightImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid);
	/** Magnetic Mirror boundary conditions LEFT wall */
	void magneticMirrorLeftImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid);
	/** Magnetic Mirror boundary conditions RIGHT wall */
	void magneticMirrorRightImage(double ***imageX, double ***imageY, double ***imageZ,double ***vectorX, double ***vectorY, double ***vectorZ, int dir,Grid *grid);
	/** Perfect conductor boundary conditions for source LEFT wall*/
    void perfectConductorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Perfect conductor boundary conditions for source RIGHT wall*/
    void perfectConductorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Magnetic mirror boundary conditions for source LEFT wall */
    void magneticMirrorLeftS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);
    /** Magnetic mirror boundary conditions for source RIGHT wall */
    void magneticMirrorRightS(double ***vectorX, double ***vectorY, double ***vectorZ,int dir);


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
	/** Magnetic mirror boundary conditions LEFT wall */
	void BmagneticMirrorLeft(int dir,Grid *grid,VirtualTopology *vct);
	/** Magnetic mirror boundary conditions RIGHT wall */
	void BmagneticMirrorRight(int dir,Grid *grid,VirtualTopology *vct);
	/** BC for PHI: divergence cleaning */
	void bcPHI_Left(double ***im,int dir);
	/** BC for PHI: divergence cleaning */
	void bcPHI_Right(double ***im,int dir);
	/** Continuous BC for PHI: divergence cleaning and divC*/
	void bcPHI_Left_continuous(double ***im,int dir);
	/** Continuous BC for PHI: divergence cleaning and divC*/
	void bcPHI_Right_continuous(double ***im,int dir);


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

	//ME
	/** get pressure tensor XX for species */
        double &getpXXns(int indexX, int indexY, int indexZ, int is) const;
        /** get pressure tensor XY for species */
        double &getpXYns(int indexX, int indexY, int indexZ, int is) const;
        /** get pressure tensor XZ for species*/
        double &getpXZns(int indexX, int indexY, int indexZ, int is) const;
        /** get pressure tensor YY for species */
        double &getpYYns(int indexX, int indexY, int indexZ, int is) const;
        /** get pressure tensor YZ for species */
        double &getpYZns(int indexX, int indexY, int indexZ, int is) const;
        /** get pressure tensor ZZ for species */
        double &getpZZns(int indexX, int indexY, int indexZ, int is) const;
	//end ME

	/** get Jx(X,Y,Z)  */
	double &getJx(int indexX, int indexY, int indexZ) const;
	/** get current -Direction X */
	double*** getJx();
	/** SPECIES: get current -Direction X */
	double**** getJxs();
	/** SPECIES: get current on center cell(indexX,indexY,indexZ)*/
	double &getJxs(int indexX, int indexY,int indexZ,int is) const;
	/** get Jy(X,Y,Z)  */
	double &getJy(int indexX, int indexY, int indexZ) const;
	/** get current -Direction Y */
	double*** getJy();
	/** SPECIES: get current -Direction Y */
	double**** getJys();
	/** SPECIES: get current on center cell(indexX,indexY,indexZ)*/
	double &getJys(int indexX, int indexY,int indexZ,int is) const;
	/** get Jz(X,Y,Z)  */
	double &getJz(int indexX, int indexY, int indexZ) const;
	/** get current -Direction Z */
	double*** getJz();
	/** SPECIES: get current -Direction Z */
	double**** getJzs();
	/** SPECIES: get current on center cell(indexX,indexY,indexZ)*/
	double &getJzs(int indexX, int indexY,int indexZ,int is) const;

	//sets
	void setRHOns(double value, int indexX, int indexY,int indexZ,int is);
	void setJxs(double value, int indexX, int indexY,int indexZ,int is);
	void setJys(double value, int indexX, int indexY,int indexZ,int is);
	void setJzs(double value, int indexX, int indexY,int indexZ,int is);
	void setpXXns(double value, int indexX, int indexY,int indexZ,int is);
	void setpXYns(double value, int indexX, int indexY,int indexZ,int is);
	void setpXZns(double value, int indexX, int indexY,int indexZ,int is);
	void setpYYns(double value, int indexX, int indexY,int indexZ,int is);
	void setpYZns(double value, int indexX, int indexY,int indexZ,int is);
	void setpZZns(double value, int indexX, int indexY,int indexZ,int is);
		
	/** print electromagnetic fields info */
	void print(void) const;

	// AMR, ME
	// for level < last
	int getNmessageBC();
	int *getTargetBC();
	int *getBCSide();
	int getTargetBOTTOM();
	int getTargetTOP();
	int getTargetLEFT();
	int getTargetRIGHT();
	// for level >0
	int getNmessagerecuBC();
	// end AMR, ME

	// ME
	void printRhons(int ns, VirtualTopology *vct);
	void printRhoc(VirtualTopology *vct);
	void setRhoc(double val);  // for debugging purposes
	void printVecN(double ***vec, VirtualTopology *vct);
	void printVecC(double ***vec, VirtualTopology *vct);
	void printVecN(double ****vec, VirtualTopology *vct, int is);
        void printVecC(double ****vec, VirtualTopology *vct, int is);
	void printJxs(int ns, VirtualTopology *vct);
	void printPxxsn(int ns, VirtualTopology *vct);
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
        /** Frequency of the incoming light wave **/
        double omega;


	// KEEP IN MEMORY GUARD CELLS ARE INCLUDED
        /** Ratio between grids **/
        double ratio;
	/** number of cells - X direction, including + 2 (guard cells) */
	int nxc;
	/** number of nodes - X direction, including + 2 extra nodes for guard cells */
	int nxn;
	/** Maximum number of projected nodes - X direction */
	int nxnproj;
	/** number of cell - Y direction, including + 2 (guard cells) */
	int nyc;
	/** number of nodes - Y direction, including + 2 extra nodes for guard cells */
	int nyn;
	/** Maximum number of projected nodes - Y direction */
	int nynproj;
    /** Indice of the last projected node */
    int lastindicex;
    int lastindicey;
        /** Number of points projected to coarser grid in x and y direction before/after possible separation**/
        int nxmsend,nxpsend,nymsend,nypsend;
        /** number of nodes of finer grid intersecting the grid in the x direction at y = Oy - dy**/
        int xnnl;
        /** number of nodes of finer grid intersecting the grid in the x direction at y = Oy +ly + dy**/
        int xnnu;
        /** number of nodes of finer grid intersecting the grid in the y direction at x = Ox - dx**/
        int ynnl;
        /** number of nodes of finer grid intersecting the grid in the y direction at x = Ox +lx +dx**/
        int ynnu;
        /** number of message the process has to send to coarser grid process for the projection of refined fields **/
        int nmessageProj;
        /** number of message the process has to receive from finer grid process for the projection of refined fields **/
        int nmessagerecvProj;
        /** number of message the process has to send to finer grid process for the projection of BC **/
        int nmessageBC;
        /** number of message the process has to receive from coarser grid process for the projection of BC **/
        int nmessagerecuBC;
	int nmessagerecuBCBOTTOM;
	int nmessagerecuBCTOP;
	int nmessagerecuBCLEFT;
	int nmessagerecuBCRIGHT;
        /** number of points received in each direction by each senders **/
        int *xmrecv;
        int *xprecv;
        int *ymrecv;
        int *yprecv;
        int *nxrecvProj;
        int *nyrecvProj;
        /** first index of the ghost nodes received from BC**/
        int *ixmrecvfirst;
        int *ixprecvfirst;
        int *iymrecvfirst;
        int *iyprecvfirst;
        /** first index of the nodes received from Projection **/
        int *ixrecvfirstProj;
        int *iyrecvfirstProj;
        int ixrecvfirstProjglobal;
        int iyrecvfirstProjglobal;
        int ixrecvlastProjglobal;
        int iyrecvlastProjglobal;
        int nxmrecv; // Number of points received in the x direction during projection
        int nymrecv; // Number of points received in the y direction during projection
        /** list of the rank of the receiving process for projection of BC **/
        int *targetBC;
	// AMR, ME
	// I need to know on which side of the BC the procs are, for dealing with particles
	int *BCSide;
	int *BCSidecu;
	// to help the coarse grid locate the first point it has to send to the refined grid
	double *xfirst_COARSE;
	double *yfirst_COARSE;
	// how many target BC processors per side
	int targetBOTTOM;
	int targetTOP;
	int targetLEFT;
	int targetRIGHT;
	int pointssentBOTTOM;
	int pointssentTOP;
	int pointssentLEFT;
	int pointssentRIGHT;
	int npointsreceivedBOTTOM;
	int npointsreceivedTOP;
	int npointsreceivedLEFT;	
	int npointsreceivedRIGHT; 
	int npointssentBOTTOM; 	
	int npointssentTOP; 	
	int npointssentLEFT;	
	int npointssentRIGHT;  
	// end AMR, ME
        /** list of the rank of the receiving process for projection of refined fields **/
        int *targetProj;
        /** list of the rank of the sending process for projection of BC **/
        int *fromBC;
        /** list of the rank of the sending process for projection of refined fields **/
        int *fromProj;
        /** total number of points that must be sent to each process for the projection of BC **/
        int *npointssent;
        /** Number of points received for all 4 ghost borders (xlow, xup,ylow,yup) for all sending process **/
        int *npointsreceived;
        /** Number of points received for the refined fields for all sending process **/
        int *npointsreceivedProj;
        /** list of the x indices of the grid that must be used during interpolation for projection of BC **/
        int *ixsent;
        /** list of the y indices of the grid that must be used during interpolation for projection of BC **/
        int *iysent;
        /** list of the x and y indices that must be used during projection of fields to coarser grid **/
        int* ixsentProj;
        int* iysentProj;
	/** local grid boundaries coordinate  */
	double xStart, xEnd, yStart, yEnd;
        /** Limits for exchange buffers after receiving projection **/
        int ixend, iyend;
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
        /** Origin of a grid **/
        double Ox;
        /** Origin of a grid **/
        double Oy;

	/** Electric potential, defined on central points of the cell*/
	double***  PHI;
	/** Electric field X-component, defined on nodes */
	double***  Ex;
	double***  Ex_old;
	double***  Ex_recvbufferproj;
	/** Implicit electric field X-component, defined on nodes */
	double***  Exth;
	/** Electric field Y-component, defined on nodes */
	double***  Ey;
	double***  Ey_old;
	double***  Ey_recvbufferproj;
	/** Implicit electric field Y-component, defined on nodes */
	double***  Eyth;
	/** Electric field Z-component, defined on nodes */
	double***  Ez;
	double***  Ez_old;
	double***  Ez_recvbufferproj;
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
	double***  Bxc_old;
	double***  Bxn_new;
	double***  Bxn_recvbufferproj;
	/** Magnetic field Y-component, defined on nodes*/
	double***  Byn;
	double***  Byc_old;
	double***  Byn_new;
	double***  Byn_recvbufferproj;
	/** Magnetic field Z-component, defined on nodes*/
	double***  Bzn;
	double***  Bzc_old;
	double***  Bzn_new;
	double***  Bzn_recvbufferproj;
        /** Weights for interpolation between grids (boundary conditions) **/
        double***  weightBC;
        /** Weights for projection between grids **/
        double****  weightProj;
        double***  normalizeProj;
        double***  reducedEx;
        double***  reducedEy;
        double***  reducedEz;
        double***  reducedBxn;
        double***  reducedByn;
        double***  reducedBzn;
        double***  normalizerecvProj;
        /** Buffer for electric field being projected to finer grid **/
        double*  bufferProj;
        double*  bufferBC;
        double*  bufferProjsend;
        double*  bufferProjrecv;
	// ME: these two buffer are bufferProjsend and bufferProjrecv with an header of 6 doubles
	// to build the projection map
	double*  INFObufferProjsend;
        double*  INFObufferProj;
	// end ME
        double*  bufferBCExxm ;  
        double*  bufferBCExxp ;
        double*  bufferBCExym ;
        double*  bufferBCExyp ;
        double*  bufferBCEyxm ;
        double*  bufferBCEyxp ;
        double*  bufferBCEyym ;
        double*  bufferBCEyyp ;
        double*  bufferBCEzxm ;
        double*  bufferBCEzxp ;
        double*  bufferBCEzym ;
        double*  bufferBCEzyp ;
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

	/*********************************************************************************/
		/*************                SOURCES                                           **/
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
        int i,j;
	if (vct->getCartesian_rank() ==0){
	  cout << "Level " << grid->getLevel() <<": *** E, B CALCULATION ***" << endl;
        }
	double ***divE = newArr3(double,nxc,nyc,1);
	double ***gradPHIX = newArr3(double,nxn,nyn,1);
	double ***gradPHIY = newArr3(double,nxn,nyn,1);
	double ***gradPHIZ = newArr3(double,nxn,nyn,1);

	double *xkrylov = new double[3*(nxn-2)*(nyn-2)];
	double *bkrylov = new double[3*(nxn-2)*(nyn-2)];

	double *bkrylovPoisson;
	double *xkrylovPoisson;

	//if (grid->getLevel()==0)//to be used in Poisson correction
	if (grid->getLevel()>-1)//to be used in Poisson correction
	  { bkrylovPoisson = new double[(nxc-2)*(nyc-2)];
	    xkrylovPoisson = new double[(nxc-2)*(nyc-2)];
    	eqValue (0.0, xkrylovPoisson, (nxc-2)*(nyc-2));
    	eqValue (0.0, bkrylovPoisson, (nxc-2)*(nyc-2));
	  }
	else
	  { bkrylovPoisson = new double[(nxc)*(nyc)];
	    xkrylovPoisson = new double[(nxc)*(nyc)];
    	eqValue (0.0, xkrylovPoisson, (nxc)*(nyc));
    	eqValue (0.0, bkrylovPoisson, (nxc)*(nyc));
	  }
	  
	eqValue (0.0, xkrylov, 3*(nxn-2)*(nyn-2));
	eqValue (0.0, bkrylov, 3*(nxn-2)*(nyn-2));
	eqValue (0.0, divE,nxc,nyc);
	eqValue (0.0, tempC,nxc,nyc);
	eqValue (0.0, gradPHIX,nxn,nyn);
	eqValue (0.0, gradPHIY,nxn,nyn);
          //Store old E and B before it is changed. It will be needed after receiving projection.
          if (grid->getLevel() < vct->getNgrids()-1){
              for (i=0;i<nxn;i++){
                  for (j=0;j<nyn;j++){
                      Ex_old[i][j][0] = Ex[i][j][0];
                      Ey_old[i][j][0] = Ey[i][j][0];
                      Ez_old[i][j][0] = Ez[i][j][0];
                  }
              }        
              for (i=0;i<nxc;i++){
                  for (j=0;j<nyc;j++){
                      Bxc_old[i][j][0] = Bxc[i][j][0];
                      Byc_old[i][j][0] = Byc[i][j][0];
                      Bzc_old[i][j][0] = Bzc[i][j][0];
                  }
              }        
          }

	// Adjust E calculating laplacian(PHI) = div(E)  -4*PI*rho
        PoissonCorrection=true;
	if (PoissonCorrection){
		if (vct->getCartesian_rank() ==0)
		  cout << "Level " << grid->getLevel() << ": *** POISSON CORRECTION ***" << endl;
                //Compute div(E)
		//if (grid->getLevel()==0)
		if (grid->getLevel()>-1)
		  grid->divN2C(divE,Ex,Ey);   //Ok in interior centers
		else
		  grid->divN2C_plusghost(divE,Ex,Ey);   //Ok also in ghost centers               
		//Compute -4*PI*rho [at the centers] 
		scale(tempC,rhoc,-FourPI,nxc,nyc);
                //Store div(E) - 4*PI*rho in the divE variable
		sum(divE,tempC,nxc,nyc);
                //Do the solve. Right hand side is in bkrylovPoisson, unknowns in xkrylovPoisson.
		if (grid->getLevel()==0)
		  {
		    phys2solver(bkrylovPoisson,divE,nxc,nyc); // bkrylovPoisson is defined only in interior centers
		    eqValue(0.0,xkrylovPoisson,(nxc-2)*(nyc-2));
		    if (vct->getPERIODICX() && vct->getPERIODICY())
		      GMRES(&Field::PoissonImage, xkrylovPoisson, (nxc-2)*(nyc-2),bkrylovPoisson,20,400,GMREStol, grid, vct, this);
		    else
		      CG(xkrylovPoisson,(nxc-2)*(nyc-2),bkrylovPoisson, 10000, CGtol, &Field::PoissonImage, grid, vct, this);
		    //Get PHI at the interior centers
		    solver2phys(PHI,xkrylovPoisson,nxc,nyc);
		  }
		else
		  {
		    //phys2solver_plusghost(bkrylovPoisson,divE,nxc,nyc); // bkrylovPoisson is defined also in the ghost centers
            //eqValue(0.0,xkrylovPoisson,(nxc)*(nyc));
		    phys2solver(bkrylovPoisson,divE,nxc,nyc); // bkrylovPoisson is defined only in interior centers
		    eqValue(0.0,xkrylovPoisson,(nxc-2)*(nyc-2));
		    
		    //CG does not work: with the current BC on the lapc2c, not symmetric anymore (maybe) 
		    //GMRES(&Field::PoissonImage, xkrylovPoisson, (nxc)*(nyc),bkrylovPoisson,20,400,GMREStol, grid, vct, this);
		    GMRES(&Field::PoissonImage, xkrylovPoisson, (nxc-2)*(nyc-2),bkrylovPoisson,20,400,GMREStol, grid, vct, this);
		    
		    //Get PHI at interior centers + ghost centers                                           
            //solver2phys_plusghost(PHI,xkrylovPoisson,nxc,nyc);
		    solver2phys(PHI,xkrylovPoisson,nxc,nyc);
		    
		  }
                //Communicate PHI to neighbours
		communicateCenter(nxc,nyc,PHI,vct);
               
                //Apply boundary conditions for PHI 
                // BC for PHI are needed only at level 0, ghost centers already calculated in the CG for refined levels
                if (grid->getLevel() == 0) {
		// adjust boundaries before calculating the lap
		  if(vct->getXleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceXleft ==0 )) 
		    bcPHI_Left(PHI,0);
		  if(vct->getXleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceXleft ==1 )){ 
		    bcPHI_Left_continuous(PHI,0);
          }
		  // boundary condition: Xright
		  if(vct->getXright_neighbor()==MPI_PROC_NULL && (bcPHIfaceXright ==0 )) 
		    bcPHI_Right(PHI,0);
		  if(vct->getXright_neighbor()==MPI_PROC_NULL && (bcPHIfaceXright ==1 )) 
		    bcPHI_Right_continuous(PHI,0);
		  // boundary condition: Yleft
		  if(vct->getYleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceYleft ==0 )) 
		    bcPHI_Left(PHI,1);
		  if(vct->getYleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceYleft ==1 )) 
		    bcPHI_Left_continuous(PHI,1);
		// boundary condition: Yright
		  if(vct->getYright_neighbor()==MPI_PROC_NULL && (bcPHIfaceYright ==0 )) 
		    bcPHI_Right(PHI,1);
		  if(vct->getYright_neighbor()==MPI_PROC_NULL && (bcPHIfaceYright ==1 )) 
		    bcPHI_Right_continuous(PHI,1);
                } else {
		  if(vct->getXleft_neighbor()==MPI_PROC_NULL ) 
		    bcPHI_Left_continuous(PHI,0);
		  if(vct->getXright_neighbor()==MPI_PROC_NULL) 
		    bcPHI_Right_continuous(PHI,0);
		  if(vct->getYleft_neighbor()==MPI_PROC_NULL ) 
		    bcPHI_Left_continuous(PHI,1);
		  if(vct->getYright_neighbor()==MPI_PROC_NULL) 
		    bcPHI_Right_continuous(PHI,1);



          }

                //Compute the gradient of PHI
		grid->gradC2N(gradPHIX,gradPHIY,PHI,vct); //ghost nodes not affected
                //Correct the electric field
		sub(Ex,gradPHIX,nxn,nyn);
		sub(Ey,gradPHIY,nxn,nyn);

    }
    // Advance timestep of E fron n to n+theta
    // find the solution with GMRES: the solution is the in the krylov space
    if (vct->getCartesian_rank() ==0){
      cout << "Level " << grid->getLevel()  << ": *** MAXWELL SOLVER ***" << endl;
    }
    MaxwellSource(bkrylov,grid,vct);
    phys2solver(xkrylov,Ex,Ey,Ez,nxn,nyn);
    // GMRES
    GMRES(&Field::MaxwellImage, xkrylov, 3*(nxn-2)*(nyn-2),bkrylov,20, 400,GMREStol, grid, vct, this);
    // move from krylov space to physical space
    solver2phys(Exth,Eyth,Ezth,xkrylov,nxn,nyn);
    communicateNode(nxn,nyn,Exth,vct);
    communicateNode(nxn,nyn,Eyth,vct);
    communicateNode(nxn,nyn,Ezth,vct);
	// here you have to put the boundaries conditions

    //BC are applied once again. This is usefull because all nodes of Eth and Ex are added later to get new Ex
    if (grid->getLevel() == 0) {
    // boundary condition: Xleft
    if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==0 ){
		perfectConductorLeft(Exth,Eyth,Ezth,0,grid);
	}else if (vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==1 ){
		magneticMirrorLeft(Exth,Eyth,Ezth,0,grid);
        }
    // boundary condition: Xright
    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==0 ){
		perfectConductorRight(Exth,Eyth,Ezth,0,grid);
	}else if (vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==1 ){
		magneticMirrorRight(Exth,Eyth,Ezth,0,grid);
        }
    // boundary condition: Yleft
    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==0 ){
		perfectConductorLeft(Exth,Eyth,Ezth,1,grid);
	}else if (vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==1 ){
		magneticMirrorLeft(Exth,Eyth,Ezth,1,grid);
        }
    // boundary condition: Yright
	if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==0 ){
		perfectConductorRight(Exth,Eyth,Ezth,1,grid);
	}else if (vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==1 ){
		magneticMirrorRight(Exth,Eyth,Ezth,1,grid);
        }  
     } else {
         electricfieldBC(Exth,Eyth,Ezth,grid,vct);
     } 

	// Advance all the fields one full timestep
    addscale(1/th,-(1.0-th)/th,Ex,Exth,nxn,nyn);
    addscale(1/th,-(1.0-th)/th,Ey,Eyth,nxn,nyn);
    addscale(1/th,-(1.0-th)/th,Ez,Ezth,nxn,nyn);
	// SMOOTHING
    //if (grid->getLevel() == 1) {
    smooth(Nvolte,Smooth,Ex,1,grid,vct);
    smooth(Nvolte,Smooth,Ey,1,grid,vct);
    smooth(Nvolte,Smooth,Ez,1,grid,vct);
    //}
    
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
    //if (grid->getLevel()==0)
    if (grid->getLevel()>-1)
      eqValue (0.0, image,(nxc-2)*(nyc-2));
    else
      eqValue (0.0, image,(nxc)*(nyc));

    eqValue (0.0, temp,nxc,nyc);
    eqValue (0.0, im,nxc,nyc);
    // move from krylov space to physical space and communicate ghost cells
    //if (grid->getLevel()==0)
    if (grid->getLevel()>-1)
      solver2phys(temp,vector,nxc,nyc);
    else
      solver2phys_plusghost(temp,vector,nxc,nyc);
    MPI_Barrier(vct->getCART_COMM());
    communicateCenter(nxc,nyc,temp,vct);
    //BC are applied only on the coarse grid
    // no need for that on the refined grid
    if (grid->getLevel() == 0) {
		// adjust boundaries before calculating the lap
		  if(vct->getXleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceXleft ==0 )) // perfect conductor
		    bcPHI_Left(PHI,0);
		  if(vct->getXleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceXleft ==1 )) // perfect conductor
		    bcPHI_Left_continuous(PHI,0);
		  // boundary condition: Xright
		  if(vct->getXright_neighbor()==MPI_PROC_NULL && (bcPHIfaceXright ==0 )) // perfect conductor
		    bcPHI_Right(PHI,0);
		  if(vct->getXright_neighbor()==MPI_PROC_NULL && (bcPHIfaceXright ==1 )) // perfect conductor
		    bcPHI_Right_continuous(PHI,0);
		  // boundary condition: Yleft
		  if(vct->getYleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceYleft ==0 )) // perfect conductor
		    bcPHI_Left(PHI,1);
		  if(vct->getYleft_neighbor()==MPI_PROC_NULL && (bcPHIfaceYleft ==1 )) // perfect conductor
		    bcPHI_Left_continuous(PHI,1);
		// boundary condition: Yright
		  if(vct->getYright_neighbor()==MPI_PROC_NULL && (bcPHIfaceYright ==0 )) // perfect conductor
		    bcPHI_Right(PHI,1);
		  if(vct->getYright_neighbor()==MPI_PROC_NULL && (bcPHIfaceYright ==1 )) // perfect conductor
		    bcPHI_Right_continuous(PHI,1);
                } else {

		  if(vct->getXleft_neighbor()==MPI_PROC_NULL ) 
		    bcPHI_Left_continuous(PHI,0);
		  if(vct->getXright_neighbor()==MPI_PROC_NULL) 
		    bcPHI_Right_continuous(PHI,0);
		  if(vct->getYleft_neighbor()==MPI_PROC_NULL ) 
		    bcPHI_Left_continuous(PHI,1);
		  if(vct->getYright_neighbor()==MPI_PROC_NULL) 
		    bcPHI_Right_continuous(PHI,1);
          }

	// calculate the laplacian
    //if (grid->getLevel()==0)
    if (grid->getLevel()>-1)
    {
	grid->lapC2C(im,temp,vct);
	phys2solver(image,im,nxc,nyc);
	}
    else
      {
	grid->lapC2C_plusghost(im,temp,vct);  // check this on boundaries
	phys2solver_plusghost(image,im,nxc,nyc);
	}
    // deallocate temporary array and objects
    delArr3(temp,nxc,nyc);
    delArr3(im,nxc,nyc);


}


/** Calculate sorgent for Maxwell solver */
inline void EMfields::MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct){
        int i,j;
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
        //Fine grid should NOT use the _onesided_derBC version but the regular version allowed by the fact that it receives boundary conditions
	grid->gradC2N(tempX,tempY,rhoh,vct);
	
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

        //BC are applied only on coarsest grid
        // apply BC on Ex, Ey, Ez
        if (grid->getLevel() == 0) {
	    // boundary condition: Xleft
           if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==0 ){
	    	perfectConductorLeft(vectX,vectY,vectZ,0,grid);
            } else if (vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==1 ){
	    	magneticMirrorLeft(vectX,vectY,vectZ,0,grid);
            }
	    // boundary condition: Xright
	    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==0 ){
	    	perfectConductorRight(vectX,vectY,vectZ,0,grid);
            } else if (vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==1 ){
	    	magneticMirrorRight(vectX,vectY,vectZ,0,grid);
            }

	    // boundary condition: Yleft
	    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==0 ){
	    	perfectConductorLeft(vectX,vectY,vectZ,1,grid);
            } else if (vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==1 ){
	    	magneticMirrorLeft(vectX,vectY,vectZ,1,grid);
            }

	    // boundary condition: Yright
	    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==0 ){
	    	perfectConductorRight(vectX,vectY,vectZ,1,grid);
            } else if (vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==1 ){
	    	magneticMirrorRight(vectX,vectY,vectZ,1,grid);
            }
	    // end of the boundaries  
        } else {
        // apply BC from coarser grid to finer grid
            electricfieldBC(vectX, vectY, vectZ, grid,vct);
        }
        
	grid->lapN2N_plusghost(imageX,vectX,vct);
	grid->lapN2N_plusghost(imageY,vectY,vct);
	grid->lapN2N_plusghost(imageZ,vectZ,vct);
       
	neg(imageX,nxn,nyn);
	neg(imageY,nxn,nyn);
	neg(imageZ,nxn,nyn);

	// grad(div(mu dot E(n + theta))       mu dot E(n + theta) = D

        //Here vectX, vectY and vectZ are storing E(n + theta) and BC have been applied to their ghost nodes.
        //These ghost nodes can be used to determine ghost nodes of MUdot which can be use later to dertermine
        //ghost centers of divC instead of having to use strange BC on ghost centers of divC.
	//MUdot(Dx, Dy, Dz, vectX, vectY, vectZ, grid);
	MUdot_plusghost(Dx, Dy, Dz, vectX, vectY, vectZ, grid);//Now acts on whole domain including ghost cells (A. Beck) This requires correct values on the ghost nodes of Bxyzn, vectxyz and rhons. for the moment it is correct on the coarse grid but not on the fine grid where rhons might be different from 0 at the ghost nodes.
        /*if (grid->getLevel() == 1 && vct->getCartesian_rank()==0) {
            cout << "apres MUdot Dx = " << Dx[8][1][0];
        }*/
	//grid->divN2C(divC,Dx,Dy);
	grid->divN2C_plusghost(divC,Dx,Dy);//Now ghost centers of divC are calculated using ghost nodes of MUdot
	// communicate
	communicateCenter(nxc,nyc,divC,vct);
        // No need for boundary conditions anymore here.

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

        //What is the use of that ??
	sum(Dx,vectX,nxn,nyn);
	sum(Dy,vectY,nxn,nyn);
	sum(Dz,vectZ,nxn,nyn);

	// move from physical space to krylov space
	phys2solver(im,imageX,imageY,imageZ,nxn,nyn);
}

/** Calculate Magnetic field with the implicit solver: calculate B defined on nodes
With E(n+ theta) computed, the magnetic field is evaluated from Faraday's law */
/*inline void EMfields::calculateB(Grid *grid, VirtualTopology *vct){
        int i,j;
	eqValue (0.0,tempXC,nxc,nyc);
	eqValue (0.0,tempYC,nxc,nyc);
	eqValue (0.0,tempZC,nxc,nyc);
	// calculate the curl of Eth
	grid->curlN2C_withghost(tempXC,tempYC,tempZC,Exth,Eyth,Ezth);
	// update the magnetic field
	addscale(-c*dt,1,Bxc,tempXC,nxc,nyc);
	addscale(-c*dt,1,Byc,tempYC,nxc,nyc);
	addscale(-c*dt,1,Bzc,tempZC,nxc,nyc);
	// communicate ghost
	communicateCenter(nxc,nyc,Bxc,vct);
	communicateCenter(nxc,nyc,Byc,vct);
	communicateCenter(nxc,nyc,Bzc,vct);
        //BC are applied only on coarsest grid
        if (grid->getLevel() == 0) {
	    // apply boundary conditions on B
	    if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==0 ){ // perfect conductor
	    	BperfectConductorLeft(0,grid,vct);
            } else if (vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==1) {
	    	BmagneticMirrorLeft(0,grid,vct);
            }
	    // boundary condition: Xright
	    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==0 ){ // perfect conductor
	    	BperfectConductorRight(0,grid,vct);
            } else if (vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==1){
	    	BmagneticMirrorRight(0,grid,vct);
            } 
	    // boundary condition: Yleft
	    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==0){ // perfect conductor
	    	BperfectConductorLeft(1,grid,vct);
            } else if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==1) {
	    	BmagneticMirrorLeft(1,grid,vct);
            }
	    // boundary condition: Yright
	    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==0){ // perfect conductor
	    	BperfectConductorRight(1,grid,vct);
            } else if (vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==1){
	    	BmagneticMirrorRight(1,grid,vct);
            }
        }
        if (grid->getLevel() > 0) {
            if (nmessagerecuBC > 0){
                magneticfieldBC(Bxn,grid,Bxn_old,vct); //Reminder: Bxn_old on a fine grid is storing the updated B field from the coarse one.
                magneticfieldBC(Byn,grid,Byn_old,vct);
                magneticfieldBC(Bzn,grid,Bzn_old,vct);
            }

	    // prova
	    if (0)
	      {
	    for (int i=0; i< nxc ; i++)
	      {
		if (vct->getYleft_neighbor()==MPI_PROC_NULL)
		  {
		     Bxc[i][0][0]=Bxc[i][5][0];
		    Byc[i][0][0]=Byc[i][5][0];
		    Bzc[i][0][0]=Bzc[i][5][0];
		  }
		if (vct->getYright_neighbor()==MPI_PROC_NULL)
		  {
		    Bxc[i][nyc-1][0]=Byc[i][nyc-1-5][0];
		    Byc[i][nyc-1][0]=Byc[i][nyc-1-5][0];
		    Bzc[i][nyc-1][0]=Bzc[i][nyc-1-5][0];
		  }
	      }
	    for (int j=0; j<nyc; j++)
	      {
		if (vct->getXleft_neighbor()==MPI_PROC_NULL)
                  {
		    Bxc[0][j][0]=Bxc[5][j][0];
		    Byc[0][j][0]=Byc[5][j][0];
		    Bzc[0][j][0]=Bzc[5][j][0];
		  }
		if (vct->getXright_neighbor()==MPI_PROC_NULL)
		  {
		    Bxc[nxc-1][j][0]=Bxc[nxc-1-5][j][0];
		    Byc[nxc-1][j][0]=Byc[nxc-1-5][j][0];
		    Bzc[nxc-1][j][0]=Bzc[nxc-1-5][j][0];
		  }
	      }
	    for (int i=0; i<nxn; i++)
	      {
		if (vct->getYleft_neighbor()==MPI_PROC_NULL)
		  {
		    Bxn[i][0][0]=Bxn[i][5][0];
		    Byn[i][0][0]=Byn[i][5][0];
		    Bzn[i][0][0]=Bzn[i][5][0];

Bxn[i][1][0]=Bxn[i][6][0];
		    Byn[i][1][0]=Byn[i][6][0];
		    Bzn[i][1][0]=Bzn[i][6][0];
		  }
		if (vct->getYright_neighbor()==MPI_PROC_NULL)
		  {
		    Bxn[i][nyn-1][0]=Bxn[i][nyn-1-5][0];
		    Byn[i][nyn-1][0]=Byn[i][nyn-1-5][0];
		    Bzn[i][nyn-1][0]=Bzn[i][nyn-1-5][0];
		    
		    Bxn[i][nyn-2][0]=Bxn[i][nyn-2-5][0];
		    Byn[i][nyn-2][0]=Byn[i][nyn-2-5][0];
		    Bzn[i][nyn-2][0]=Bzn[i][nyn-2-5][0];
		  }
	      }
	    for (int j=0; j<nyn; j++)
	      {
		if (vct->getXleft_neighbor()==MPI_PROC_NULL)
		  {
		    Bxn[0][j][0]=Bxn[5][j][0];
		    Byn[0][j][0]=Byn[5][j][0];
		    Bzn[0][j][0]=Bzn[5][j][0];
		    		    
		    Bxn[1][j][0]=Bxn[6][j][0];
		    Byn[1][j][0]=Byn[6][j][0];
		    Bzn[1][j][0]=Bzn[6][j][0];
		    
		  }
		if (vct->getXright_neighbor()==MPI_PROC_NULL)
                  {
		    //Bxn[nxn-1][j][0]=Bxn[nxn-1-5][j][0];
		    //Byn[nxn-1][j][0]=Byn[nxn-1-5][j][0];
		    //Bzn[nxn-1][j][0]=Bzn[nxn-1-5][j][0];

		    //Bxn[nxn-2][j][0]=Bxn[nxn-2-5][j][0];
		    //Byn[nxn-2][j][0]=Byn[nxn-2-5][j][0];
		    //Bzn[nxn-2][j][0]=Bzn[nxn-2-5][j][0];
                  }
	      }   
	    // prova
	      } // end prova
        }
	grid->interpC2N(Bxn,Bxc,vct);
        grid->interpC2N(Byn,Byc,vct);
        grid->interpC2N(Bzn,Bzc,vct);

	if (vct->getYleft_neighbor()==MPI_PROC_NULL)
	  {
	    for (int i=0; i<nxc; i++)
	      {
		Bxc[i][0][0] = .25*(Bxn[i+1][0][0] + Bxn[i][0][0] + Bxn[i+1][1][0] + Bxn[i][1][0]);
		Byc[i][0][0] = .25*(Byn[i+1][0][0] + Byn[i][0][0] + Byn[i+1][1][0] + Byn[i][1][0]);
		Bzc[i][0][0] = .25*(Bzn[i+1][0][0] + Bzn[i][0][0] + Bzn[i+1][1][0] + Bzn[i][1][0]);
	      }
	  }
	if (vct->getYright_neighbor()==MPI_PROC_NULL)
	  {
	    for(int i=0; i<nxc; i++)
	      { 
		Bxc[i][nyc-1][0] = .25*(Bxn[i+1][nyc-1][0] + Bxn[i][nyc-1][0] + Bxn[i+1][nyc][0] + Bxn[i][nyc][0]);
		Byc[i][nyc-1][0] = .25*(Byn[i+1][nyc-1][0] + Byn[i][nyc-1][0] + Byn[i+1][nyc][0] + Byn[i][nyc][0]);
		Bzc[i][nyc-1][0] = .25*(Bzn[i+1][nyc-1][0] + Bzn[i][nyc-1][0] + Bzn[i+1][nyc][0] + Bzn[i][nyc][0]);
              }
	  }
	if (vct->getXleft_neighbor()==MPI_PROC_NULL)
          {
            for(int j=0; j<nyc; j++)
	      { 
		Bxc[0][j][0] = .25*(Bxn[1][j][0] + Bxn[0][j][0] + Bxn[1][j+1][0] + Bxn[0][j+1][0]);
		Byc[0][j][0] = .25*(Byn[1][j][0] + Byn[0][j][0] + Byn[1][j+1][0] + Byn[0][j+1][0]);
		Bzc[0][j][0] = .25*(Bzn[1][j][0] + Bzn[0][j][0] + Bzn[1][j+1][0] + Bzn[0][j+1][0]);
              }
          }
	if (vct->getXright_neighbor()==MPI_PROC_NULL)
          {
            for(int i=0; i<nxc; i++)
              {
		Bxc[nxc-1][j][0] = .25*(Bxn[nxc-1+1][j][0] + Bxn[nxc-1][j][0] + Bxn[nxc-1+1][j+1][0] + Bxn[nxc-1][j+1][0]);
		Byc[nxc-1][j][0] = .25*(Byn[nxc-1+1][j][0] + Byn[nxc-1][j][0] + Byn[nxc-1+1][j+1][0] + Byn[nxc-1][j+1][0]);
		Bzc[nxc-1][j][0] = .25*(Bzn[nxc-1+1][j][0] + Bzn[nxc-1][j][0] + Bzn[nxc-1+1][j+1][0] + Bzn[nxc-1][j+1][0]);
              }
	  }
	// interpolate C2N to update the active nodes. Ghost nodes have just been updated by the BCs.                                          
	
}*/

inline void EMfields::calculateB(Grid *grid, VirtualTopology *vct){
  int i,j;
  eqValue (0.0,tempXC,nxc,nyc);
  eqValue (0.0,tempYC,nxc,nyc);
  eqValue (0.0,tempZC,nxc,nyc);

   // calculate the curl of Eth
  if (grid->getLevel() > 0) {
    grid->curlN2C(tempXC,tempYC,tempZC,Exth,Eyth,Ezth);
  } else {
     grid->curlN2C_withghost(tempXC,tempYC,tempZC,Exth,Eyth,Ezth);
  }
 
  // update the magnetic field
  addscale(-c*dt,1,Bxc,tempXC,nxc,nyc);
  addscale(-c*dt,1,Byc,tempYC,nxc,nyc);
  addscale(-c*dt,1,Bzc,tempZC,nxc,nyc);
  // communicate ghost
  communicateCenter(nxc,nyc,Bxc,vct);
  communicateCenter(nxc,nyc,Byc,vct);
  communicateCenter(nxc,nyc,Bzc,vct);
 
  //BC are applied only on coarsest grid
  if (grid->getLevel() == 0) {
    // apply boundary conditions on B
    if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==0 ){ // perfect conductor
      BperfectConductorLeft(0,grid,vct);
    } else if (vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==1) {
      BmagneticMirrorLeft(0,grid,vct);
    }
    // boundary condition: Xright
    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==0 ){ // perfect conductor
      BperfectConductorRight(0,grid,vct);
    } else if (vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==1){
      BmagneticMirrorRight(0,grid,vct);
    } 
    // boundary condition: Yleft
    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==0){ // perfect conductor
      BperfectConductorLeft(1,grid,vct);
    } else if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==1) {
      BmagneticMirrorLeft(1,grid,vct);
    }
    // boundary condition: Yright
    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==0){ // perfect conductor
      BperfectConductorRight(1,grid,vct);
    } else if (vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==1){
      BmagneticMirrorRight(1,grid,vct);
    }
    }
  if (grid->getLevel() > 0) {
    if (nmessagerecuBC > 0){
      magneticfieldBC(Bxn,grid,Bxn_new,vct); 
      magneticfieldBC(Byn,grid,Byn_new,vct);
      magneticfieldBC(Bzn,grid,Bzn_new,vct);
    }
    //Recalculate magnetic field in the ghost centers
    i=0;
    if(vct->getXleft_neighbor()==MPI_PROC_NULL){ 
      for (j=0;j < nyc;j++){
	Bxc[i][j][0] =(Bxn[i][j][0] + Bxn[i][j+1][0] + Bxc[i+1][j][0])/3.;
	Byc[i][j][0] =(Byn[i][j][0] + Byn[i][j+1][0] + Byc[i+1][j][0])/3.;
	Bzc[i][j][0] =(Bzn[i][j][0] + Bzn[i][j+1][0] + Bzc[i+1][j][0])/3.;
      }
    }
    i=nxc-1;
    if(vct->getXright_neighbor()==MPI_PROC_NULL){ 
      for (j=0;j < nyc;j++){
	Bxc[i][j][0] = (Bxn[i+1][j][0] + Bxn[i+1][j+1][0] + Bxc[i-1][j][0])/3.;
	Byc[i][j][0] = (Byn[i+1][j][0] + Byn[i+1][j+1][0] + Byc[i-1][j][0])/3.;
	Bzc[i][j][0] = (Bzn[i+1][j][0] + Bzn[i+1][j+1][0] + Bzc[i-1][j][0])/3.;
      }
    }
    j=0;
    if(vct->getYleft_neighbor()==MPI_PROC_NULL){ 
      for (i=0;i < nxc;i++){
	Bxc[i][j][0] = (Bxn[i][j][0] + Bxn[i+1][j][0] + Bxc[i][j+1][0])/3.;
	Byc[i][j][0] = (Byn[i][j][0] + Byn[i+1][j][0] + Byc[i][j+1][0])/3.;
	Bzc[i][j][0] = (Bzn[i][j][0] + Bzn[i+1][j][0] + Bzc[i][j+1][0])/3.;
      }
    }
    j=nyc-1;
    if(vct->getYright_neighbor()==MPI_PROC_NULL){ 
      for (i=0;i < nxc;i++){
	Bxc[i][j][0] = (Bxn[i][j+1][0] + Bxn[i+1][j+1][0] + Bxc[i][j-1][0])/3.;
	Byc[i][j][0] = (Byn[i][j+1][0] + Byn[i+1][j+1][0] + Byc[i][j-1][0])/3.;
	Bzc[i][j][0] = (Bzn[i][j+1][0] + Bzn[i+1][j+1][0] + Bzc[i][j-1][0])/3.;
      }
    }
    //And now for the corners
    if(vct->getXleft_neighbor()==MPI_PROC_NULL && vct->getYleft_neighbor()==MPI_PROC_NULL){ 
      Bxc[0][0][0] = (Bxn[0][0][0] + Bxn[1][0][0] + Bxn[0][1][0])/3.;
      Byc[0][0][0] = (Byn[0][0][0] + Byn[1][0][0] + Byn[0][1][0])/3.;
      Bzc[0][0][0] = (Bzn[0][0][0] + Bzn[1][0][0] + Bzn[0][1][0])/3.;
    }
    if(vct->getXleft_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL){ 
      Bxc[0][nyc-1][0] = (Bxn[0][nyc][0] + Bxn[1][nyc][0] + Bxn[0][nyc-1][0])/3.;
      Byc[0][nyc-1][0] = (Byn[0][nyc][0] + Byn[1][nyc][0] + Byn[0][nyc-1][0])/3.;
      Bzc[0][nyc-1][0] = (Bzn[0][nyc][0] + Bzn[1][nyc][0] + Bzn[0][nyc-1][0])/3.;
    }
    if(vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYleft_neighbor()==MPI_PROC_NULL){ 
      Bxc[nxc-1][0][0] = (Bxn[nxc][0][0] + Bxn[nxc-1][0][0] + Bxn[nxc][1][0])/3.;
      Byc[nxc-1][0][0] = (Byn[nxc][0][0] + Byn[nxc-1][0][0] + Byn[nxc][1][0])/3.;
      Bzc[nxc-1][0][0] = (Bzn[nxc][0][0] + Bzn[nxc-1][0][0] + Bzn[nxc][1][0])/3.;
    }
    if(vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL){ 
      Bxc[nxc-1][nyc-1][0] = (Bxn[nxc][nyc][0] + Bxn[nxc-1][nyc][0] + Bxn[nxc][nyc-1][0])/3.;
      Byc[nxc-1][nyc-1][0] = (Byn[nxc][nyc][0] + Byn[nxc-1][nyc][0] + Byn[nxc][nyc-1][0])/3.;
      Bzc[nxc-1][nyc-1][0] = (Bzn[nxc][nyc][0] + Bzn[nxc-1][nyc][0] + Bzn[nxc][nyc-1][0])/3.;
    }

  }
  
    // interpolate C2N to update the active nodes. Ghost nodes have just been updated by the BCs.
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
		MPI_Barrier(vct->getCART_COMM());
		if (vct->getCartesian_rank()==0)
			cout << "LOADING EM FIELD FROM RESTART FILE in " + RestartDirName + "/restart.hdf" << endl;
		string name_file;
        stringstream ss;
	    ss << vct->getCartesian_rank();
		name_file = RestartDirName + "/restart" + ss.str() + ".hdf";
		cout << "nome file che apre: " << name_file << endl;
		// hdf stuff
		hid_t    file_id, dataspace;
		hid_t    datatype, dataset_id;
		herr_t   status;
		size_t   size;
		hsize_t     dims_out[2];           /* dataset dimensions */
		int status_n;
		MPI_Barrier(vct->getCART_COMM());
		// open the hdf file
		cout << "apri file HDF5" << endl;
		file_id = H5Fopen(name_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_id < 0){
			cout << "couldn't open file: " << name_file << endl;
			cout << "RESTART NOT POSSIBLE" << endl;
		}

		cout << "apri data set Bx" << endl;
		//dataset_id = H5Dopen(file_id,"/fields/Bx/cycle_0");
		dataset_id = H5Dopen1(file_id,"/fields/Bx/cycle_0");   //--- HDF5 1.8
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
		cout << "letto set Bx" << endl;

		// Byn
		cout << "apri data set By" << endl;
		//dataset_id = H5Dopen(file_id,"/fields/By/cycle_0");
		dataset_id = H5Dopen1(file_id,"/fields/By/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Byn[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Byn,vct);
		status = H5Dclose(dataset_id);


		// Bzn
		cout << "apri data set Bz" << endl;
		//dataset_id = H5Dopen(file_id,"/fields/Bz/cycle_0");
		dataset_id = H5Dopen1(file_id,"/fields/Bz/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Bzn[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Bzn,vct);
		status = H5Dclose(dataset_id);

		// Ex
		cout << "apri data set Ex" << endl;
		//dataset_id = H5Dopen(file_id,"/fields/Ex/cycle_0");
		dataset_id = H5Dopen1(file_id,"/fields/Ex/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Ex[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Ex,vct);
		status = H5Dclose(dataset_id);


		// Ey
		cout << "apri data set Ey" << endl;
		//dataset_id = H5Dopen(file_id,"/fields/Ey/cycle_0");
		dataset_id = H5Dopen1(file_id,"/fields/Ey/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Ey[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Ey,vct);
		status = H5Dclose(dataset_id);

		// Ez
		cout << "apri data set Ez" << endl;
		//dataset_id = H5Dopen(file_id,"/fields/Ez/cycle_0");
		dataset_id = H5Dopen1(file_id,"/fields/Ez/cycle_0");
		status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,
						 H5S_ALL,H5P_DEFAULT,temp_storage);
		k=0;
		for (int i=1; i < nxn-1; i++)
			for (int j=1; j <nyn-1; j++)
				Ez[i][j][0] = temp_storage[k++];
		communicateNode(nxn, nyn, Ez,vct);
		status = H5Dclose(dataset_id);

		// open the charge density for species
		cout << "apri data set charge density for species" << endl;

		stringstream *species_name = new stringstream[ns];
		for (int is=0; is < ns;is++){
			species_name[is] << is;
			string name_dataset = "/moments/species_" + species_name[is].str() + "/rho/cycle_0";
			//dataset_id = H5Dopen(file_id,name_dataset.c_str());
			dataset_id = H5Dopen1(file_id,name_dataset.c_str());
			status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL,H5S_ALL,H5P_DEFAULT,temp_storage);
			k=0;
			for (int i=1; i < nxn-1; i++)
				for (int j=1; j <nyn-1; j++)
					rhons[is][i][j][0] = temp_storage[k++];
			communicateNode(nxn, nyn, rhons,is,vct);
			status = H5Dclose(dataset_id);

		}
		// deallocate
		delete[] species_name;

		cout << "close set charge density for species" << endl;
		cout << "interpolate B" << endl;
		MPI_Barrier(vct->getCART_COMM());
		// initialize B on centers
		grid->interpN2C(Bxc,Bxn);
		grid->interpN2C(Byc,Byn);
		grid->interpN2C(Bzc,Bzn);
		cout << "interpolate rho" << endl;
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons); // calculate density on the center
											 // close the hdf file
		status = H5Fclose(file_id);
		MPI_Barrier(vct->getCART_COMM());
        cout << "close file" << endl;
	}
	MPI_Barrier(vct->getCART_COMM());
	cout << " initialization restart fields" << endl;
	cout << " comm restart fields" << endl;
	// communicate ghost
	//communicateNode(nxn,nyn,Bxn,vct);
	//communicateNode(nxn,nyn,Byn,vct);
	//communicateNode(nxn,nyn,Bzn,vct);
	//communicateNode(nxn,nyn,Ex,vct);
	//communicateNode(nxn,nyn,Ey,vct);
	//communicateNode(nxn,nyn,Ez,vct);
	// communicate ghost
	//communicateCenter(nxc,nyc,Bxc,vct);
	//communicateCenter(nxc,nyc,Byc,vct);
	//communicateCenter(nxc,nyc,Bzc,vct);


	cout << "fine initialization restart fields" << endl;
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

inline void EMfields::initLightwave(VirtualTopology *vct, Grid *grid){
    int i;
    double intensity,k,kx,ky,globalx,globaly,startingT,coarsedx,coarsedy;

    coarsedx = grid->getDX()*pow(grid->getRatio(),grid->getLevel()); 
    coarsedy = grid->getDY()*pow(grid->getRatio(),grid->getLevel()); 

    intensity =0.1;
    kx = PI/2./(Lx+2*coarsedx);
    ky = PI/2./(Ly+2*coarsedy);
    k = sqrt(kx*kx+ky*ky);
    Ox = 0;
    Oy = 0;
    for (int i=1; i < grid->getLevel()+1; i++){
        Ox += grid->getOx(i);   
        Oy += grid->getOy(i);   
    }
    startingT = 0.25/k;
    if (vct->getCartesian_rank_COMMTOTAL()==0){
        cout << "Init TE waves" << endl;
        cout << "kx = "<<kx<<" ky= "<<ky<<" k= "<<k<<" Starting T = "<<startingT<<endl;
    }

 
    for (i=0; i <nxn; i++)
        for (int j=0; j <nyn; j++){
            // initialize the density for species
            for (int is=0; is < ns; is++){
            	rhons[is][i][j][0] = rhoINIT[is]/FourPI;
            }
            globalx = grid->getXN(i,j,0)+ Ox + coarsedx;
            globaly = grid->getYN(i,j,0)+ Oy + coarsedy;
	    Ex[i][j][0] =  intensity*ky*cos(kx*globalx)*cos(ky*globaly)*sin(k*startingT)/kx;
	    Ey[i][j][0] =  intensity*sin(kx*globalx)*sin(ky*globaly)*sin(k*startingT);
	    Ez[i][j][0] =  0.0;
	    Bxn[i][j][0] =  0.0;
	    Byn[i][j][0] =  0.0;
	    Bzn[i][j][0] =  intensity*k*cos(kx*globalx)*sin(ky*globaly)*cos(k*startingT)/kx;
        }
		// communicate ghost
		communicateNode(nxn,nyn,Bxn,vct);
                communicateNode(nxn,nyn,Byn,vct);
                communicateNode(nxn,nyn,Bzn,vct);

        for (i=0; i <nxc; i++)
        	for (int j=0; j <nyc; j++){
                    globalx = grid->getXC(i,j,0)+ Ox + coarsedx;
                    globaly = grid->getYC(i,j,0)+ Oy + coarsedy;
	            Bxc[i][j][0] =  0.0;
	            Byc[i][j][0] =  0.0;
	            Bzc[i][j][0] =  intensity*k*cos(kx*globalx)*sin(ky*globaly)/kx;
                }
		communicateCenter(nxc, nyc,Bxc,vct);
		communicateCenter(nxc, nyc,Byc,vct);
		communicateCenter(nxc, nyc,Bzc,vct);
		communicateNode(nxn,nyn,Ex,vct);
		communicateNode(nxn,nyn,Ey,vct);
		communicateNode(nxn,nyn,Ez,vct);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C(rhocs,is,rhons);
              //Apply boundary conditions on coarse grid for electric fields. It is done for magnetic fields in Maxwell source.

       if (grid->getLevel()==0){
	    // boundary condition: Xleft
           if(vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==0 ){
	    	perfectConductorLeft(Ex,Ey,Ez,0,grid);
            } else if (vct->getXleft_neighbor()==MPI_PROC_NULL && bcEMfaceXleft ==1 ){
	    	magneticMirrorLeft(Ex,Ey,Ez,0,grid);
            }
	    // boundary condition: Xright
	    if(vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright==0 ){
	    	perfectConductorRight(Ex,Ey,Ez,0,grid);
            } else if (vct->getXright_neighbor()==MPI_PROC_NULL && bcEMfaceXright ==1 ){
	    	magneticMirrorRight(Ex,Ey,Ez,0,grid);
            }

	    // boundary condition: Yleft
	    if(vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==0 ){
	    	perfectConductorLeft(Ex,Ey,Ez,1,grid);
            } else if (vct->getYleft_neighbor()==MPI_PROC_NULL && bcEMfaceYleft ==1 ){
	    	magneticMirrorLeft(Ex,Ey,Ez,1,grid);
            }

	    // boundary condition: Yright
	    if(vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright==0 ){
	    	perfectConductorRight(Ex,Ey,Ez,1,grid);
            } else if (vct->getYright_neighbor()==MPI_PROC_NULL && bcEMfaceYright ==1 ){
	    	magneticMirrorRight(Ex,Ey,Ez,1,grid);
            }

       }
}

/**
*
 *
 * Initialize the EM for the GEM Challenge. The equilibrium chosen for the reconnection challenge problem is a Harris equilibrium with
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
				Bzn[i][j][0] = B0z;
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
/** Double Harris Equilibrium **/
inline void EMfields::initDoubleHarris(VirtualTopology *vct, Grid *grid){
        int i;
        double globalx, globaly,coarsedx,coarsedy;
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
        Ox = 0;
        Oy = 0;
        for (int i=1; i < grid->getLevel()+1; i++){
            Ox += grid->getOx(i);   
            Oy += grid->getOy(i);   
        }
        coarsedx = grid->getDX()*pow(grid->getRatio(),grid->getLevel()); 
        coarsedy = grid->getDY()*pow(grid->getRatio(),grid->getLevel()); 
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
		    globalx = grid->getXN(i,j,0)+ Ox +coarsedx ;
		    globaly = grid->getYN(i,j,0)+ Oy +coarsedy ;
		    // initialize the density for species
		    for (int is=0; is < ns; is++){
		      //if (DriftSpecies[is])
		      if (is == 0 || is == 1) // at the beginning no drift is added at the beginning
			{
			  rhons[is][i][j][0] = ((rhoINIT[is]/(cosh((globaly-Ly/4)/delta)*cosh((globaly-Ly/4)/delta))))/FourPI;
			  rhons[is][i][j][0] += ((rhoINIT[is]/(cosh((globaly-3*Ly/4)/delta)*cosh((globaly-3*Ly/4)/delta))))/FourPI;
			}
		      else
			rhons[is][i][j][0] = rhoINIT[is]/FourPI;
		    }
		    // electric field
		    Ex[i][j][0] =  0.0;
		    Ey[i][j][0] =  0.0;
		    Ez[i][j][0] =  0.0;
		    // Magnetic field
		    Bxn[i][j][0] = B0x*(-1.0+tanh((globaly - Ly/4)/delta)+tanh(-(globaly - 3*Ly/4)/delta));
		    // add the initial GEM perturbation
		    Bxn[i][j][0] +=(B0x*pertGEM)*(M_PI/Ly)*cos(2*M_PI*globalx/Lx)*sin(M_PI*(globaly- Ly/2)/Ly  );
		    Byn[i][j][0] = B0y -(B0x*pertGEM)*(2*M_PI/Lx)*sin(2*M_PI*globalx/Lx)*cos(M_PI*(globaly- Ly/2)/Ly);
		    // add the initial X perturbation
		    xpert = globalx- Lx/4;
		    ypert = globaly- Ly/4;
		    exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));
		    
		    Bxn[i][j][0] +=(B0x*pertX)*exp_pert*(
							 -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
							 -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
							 );
		    
		    Byn[i][j][0] +=(B0x*pertX)*exp_pert*(
							 cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
							 +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
							 );
		    // add the second initial X perturbation
		    xpert = globalx- 3*Lx/4;
		    ypert = globaly- 3*Ly/4;
		    exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));
		    
		    Bxn[i][j][0] +=(-B0x*pertX)*exp_pert*(
							  -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
							  -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
							  );
		    
		    Byn[i][j][0] +=(-B0x*pertX)*exp_pert*(
							  cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
							  +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
							  );
		    // guide field
		    Bzn[i][j][0] = B0z;
		  }
		// communicate ghost
		communicateNode(nxn,nyn,Bxn,vct);
		communicateNode(nxn,nyn,Byn,vct);
		communicateNode(nxn,nyn,Bzn,vct);

		//ME
		//for (int is=0; is<ns; is++)
		//  communicateNode(nxn,nyn,rhons,is,vct);
		//end ME
		for (int i=0; i <nxc; i++)
			for (int j=0; j <nyc; j++){
                                globalx = grid->getXC(i,j,0)+ Ox + coarsedx;
                                globaly = grid->getYC(i,j,0)+ Oy + coarsedy;
				Bxc[i][j][0] = B0x*(-1.0+tanh((globaly - Ly/4)/delta)+tanh(-(globaly - 3*Ly/4)/delta));
				// add the initial GEM perturbation
				Bxc[i][j][0] +=(B0x*pertGEM)*(M_PI/Ly)*cos(2*M_PI*globalx/Lx)*sin(M_PI*(globaly- Ly/2)/Ly  );
				Byc[i][j][0] = B0y -(B0x*pertGEM)*(2*M_PI/Lx)*sin(2*M_PI*globalx/Lx)*cos(M_PI*(globaly- Ly/2)/Ly);
				// add the initial X perturbation
				xpert = globalx- Lx/4;
				ypert = globaly- Ly/4;
				exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

				Bxc[i][j][0] +=(B0x*pertX)*exp_pert*(
													 -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
													 -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
													 );

				Byc[i][j][0] +=(B0x*pertX)*exp_pert*(
													 cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
													 +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
													 );
				// add the second initial X perturbation
								xpert = globalx- 3*Lx/4;
								ypert = globaly- 3*Ly/4;
								exp_pert = exp(-(xpert/delta)*(xpert/delta)-(ypert/delta)*(ypert/delta));

								Bxc[i][j][0] +=(-B0x*pertX)*exp_pert*(
																	 -cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*ypert/delta
																	 -cos(M_PI*xpert/10.0/delta)*sin(M_PI*ypert/10.0/delta)*M_PI/10.0
																	 );

								Byc[i][j][0] +=(-B0x*pertX)*exp_pert*(
																	 cos(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*2.0*xpert/delta
																	 +sin(M_PI*xpert/10.0/delta)*cos(M_PI*ypert/10.0/delta)*M_PI/10.0
																	 );
				// guide field
				Bzc[i][j][0] = B0z;
			}
				communicateCenter(nxc, nyc,Bxc,vct);
		communicateCenter(nxc, nyc,Byc,vct);
		communicateCenter(nxc, nyc,Bzc,vct);
		communicateNode(nxn,nyn,Ex,vct);
		communicateNode(nxn,nyn,Ey,vct);
		communicateNode(nxn,nyn,Ez,vct);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C_alsoGC(rhocs,is,rhons);
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
				grid->interpN2C_alsoGC(Bxc,Bxn);
		grid->interpN2C_alsoGC(Byc,Byn);
		grid->interpN2C_alsoGC(Bzc,Bzn);
		// communicate ghost
		communicateNode(nxn,nyn,Bxn,vct);
		communicateNode(nxn,nyn,Byn,vct);
		communicateNode(nxn,nyn,Bzn,vct);
		communicateNode(nxn,nyn,Ex,vct);
		communicateNode(nxn,nyn,Ey,vct);
		communicateNode(nxn,nyn,Ez,vct);
		for (int is=0 ; is<ns; is++)
			grid->interpN2C_alsoGC(rhocs,is,rhons);
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
/**
*
*  Init Force Free (JxB=0)
*/
inline void EMfields::initForceFree(VirtualTopology *vct, Grid *grid){
	if (restart1 ==0){

		// initialize
		if (vct->getCartesian_rank() ==0){
			cout << "-------------------------" << endl;
			cout << "Initialize Force Free " << endl;
			cout << "-------------------------" << endl;
			cout << "B0x                              = " << B0x << endl;
			cout << "B0y                              = " << B0y << endl;
			cout << "B0z                              = " << B0z << endl;
			cout << "Delta (current sheet thickness) = " << delta << endl;
			for (int i=0; i < ns; i++){
				cout << "rho species " << i <<" = " << rhoINIT[i];
			}
			cout << "-------------------------" << endl;
		}

		for (int i=0; i <nxn; i++)
			for (int j=0; j <nyn; j++){
				// initialize constant density
				for (int is=0; is < ns; is++){
						rhons[is][i][j][0] = rhoINIT[is]/FourPI;
				}
				// electric field
				Ex[i][j][0] =  0.0;
				Ey[i][j][0] =  0.0;
				Ez[i][j][0] =  0.0;
				// Magnetic field
				Bxn[i][j][0] = B0x*tanh((grid->getYN(i,j,0) - Ly/2)/delta);
				// add the initial GEM perturbation
				Bxn[i][j][0] +=(B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXN(i,j,0)/Lx)*sin(M_PI*(grid->getYN(i,j,0)- Ly/2)/Ly  );
				Byn[i][j][0] = B0y -(B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXN(i,j,0)/Lx)*cos(M_PI*(grid->getYN(i,j,0)- Ly/2)/Ly);
				// guide field
				Bzn[i][j][0] = B0z/cosh((grid->getYN(i,j,0) - Ly/2)/delta);
			}
				// communicate ghost
		//		communicateNode(nxn,nyn,Bxn,vct);
		//communicateNode(nxn,nyn,Byn,vct);
		//communicateNode(nxn,nyn,Bzn,vct);
		for (int i=0; i <nxc; i++)
			for (int j=0; j <nyc; j++){
				Bxc[i][j][0] = B0x*tanh((grid->getYC(i,j,0) - Ly/2)/delta);
				// add the initial GEM perturbation
				Bxc[i][j][0] +=(B0x/10.0)*(M_PI/Ly)*cos(2*M_PI*grid->getXC(i,j,0)/Lx)*sin(M_PI*(grid->getYC(i,j,0)- Ly/2)/Ly  );
				Byc[i][j][0] = B0y -(B0x/10.0)*(2*M_PI/Lx)*sin(2*M_PI*grid->getXC(i,j,0)/Lx)*cos(M_PI*(grid->getYC(i,j,0)- Ly/2)/Ly);
				// guide field
				Bzc[i][j][0] = B0z/cosh((grid->getYC(i,j,0) - Ly/2)/delta);
			}
		//communicateCenter(nxc, nyc,Bxc,vct);
		//communicateCenter(nxc, nyc,Byc,vct);
		//communicateCenter(nxc, nyc,Bzc,vct);
		//communicateNode(nxn,nyn,Ex,vct);
		//communicateNode(nxn,nyn,Ey,vct);
		//communicateNode(nxn,nyn,Ez,vct);
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
/** Initialize the weight for the BC projections between grids **/
/** Also initialize the number and size of messages that must be sent between grids every cycle for the BC **/
/** The first xnnl quadruplets are the weights for the y=Oy-finedy ghost nodes**/
/**The next xnnu quadruplets are the weights for the y=Oy+finely+finedy ghost nodes**/
/** The next ynnl quadruplets are the weights for the x=Ox-finedx ghost nodes**/
/** The last ynnu quadruplets are the weights for the x=Ox+finelx+finedx ghost nodes**/
/*inline int EMfields::initWeightBC(VirtualTopology *vct, Grid *grid, CollectiveIO* col){

int i,j,ix,iy,nproc;
double finedx, finedy,finelx,finely, xfirst, xlast,xfirstnext, Ox, Oy,xloc,yloc;
double coarsedx, coarsedy,coarselx,coarsely;
double finelxplusfinedx;
double xshift, yshift;

 nmessagerecuBC=0;
 nmessageBC=0;

// AMR, ME//
 targetBOTTOM=0;
 targetTOP=0;
 targetLEFT=0;
 targetRIGHT=0;
// end AMR, ME//

if ( grid->getLevel() < vct->getNgrids()-1) {                             //If this grid is considered as coarse by another grid
    finedx = grid->getDX()/col->getRatio();
    finelx = col->getLx()/pow(col->getRatio(),grid->getLevel()+1);
    finedy = grid->getDY()/col->getRatio();
    finely = col->getLy()/pow(col->getRatio(),grid->getLevel()+1);
    finelxplusfinedx = col->getLx()/(double)col->getNxc()*((double)col->getNxc()+1)/(double)pow((double)col->getRatio(),(double)grid->getLevel()+1.);
    Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid
    Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid
    j=0;
// Computation of the weights
// If the y=Oy-finedy ghost nodes of the fine grid intersect this process
    if (grid->getmodifiedYstart(vct) < Oy-finedy && grid->getmodifiedYend(vct) > Oy-finedy && grid->getmodifiedXstart(vct) < Ox+finelx+finedx && grid->getmodifiedXend(vct) > Ox-finedx) {
        xfirst = max(-finedx, ceil((grid->getmodifiedXstart(vct)-Ox)/finedx)*finedx); 
        xfirstnext = ceil((grid->getmodifiedXend(vct)-Ox)/finedx)*finedx; 
        xlast  = min(finelx+finedx, floor((grid->getmodifiedXend(vct)-xfirst-Ox)/finedx)*finedx+xfirst);

	// to remove                                 
	//if (vct->getCartesian_rank_COMMTOTAL()== 11 || vct->getCartesian_rank_COMMTOTAL()== 14)
	//  { 
	//    cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << " xlast1: " << xlast <<endl; 
	//    //cout << "grid->getXend() " <<grid->getXend() << " xlast " <<xlast << " Ox " <<Ox <<" fabs(grid->getXend()-xlast-Ox) " << fabs(grid->getXend()-xlast-Ox) << " DBL_EPSILON " << DBL_EPSILON << " fabs(grid->getXend()-xlast-Ox)*0.99999 "<<fabs(grid->getXend()-xlast-Ox)*0.99999  <<endl;                  
	//    }

        //if(grid->getXend()-xlast-Ox < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision
	if(fabs(grid->getXend()-xlast-Ox)*0.99999 < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision; 0.99999 because the obvious choices of ratios often bring to overlapping that this alone does not resolve          
            xlast = xlast - finedx;
        }

	if (vct->getCartesian_rank_COMMTOTAL()== 25 || vct->getCartesian_rank_COMMTOTAL()== 49) 
	  {
	    cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << " xlast bef correction: " << min(finelx+finedx, floor((grid->getmodifiedXend(vct)-xfirst-Ox)/finedx)*finedx+xfirst) <<endl;
	    cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << " xlast: " << xlast <<endl;       
	    cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << " xfirst: " << xfirst <<endl;       
	  }
        xnnl = floor((xlast-xfirst)/finedx+0.5)+1; //floor(x+0.5) used to round to closest integer in case the division is not working properly
        yloc = Oy - finedy;
        for (i =0;i<xnnl;i++) {
            xshift = xfirst + i * finedx ;
            xloc = xshift + Ox ;
            
            //Weights used for interpolation on nodes (for E and B)
            ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
            iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
            ixsent[i] = ix;
            iysent[i] = iy;
            weightBC[i][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
            weightBC[i][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
            weightBC[i][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
            weightBC[i][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
            nproc = col->getYLEN()*floor(xshift/((grid->getNXC()-2.)*finedx))+col->getXLEN()*col->getYLEN()*(grid->getLevel()+1); // rank of the proc on the fine grid receiving this point(in MPI_COMM_WORLD)
            nproc = max(nproc, col->getXLEN()*col->getYLEN()*(grid->getLevel()+1));
            nproc = min(nproc, col->getYLEN()*((grid->getLevel()+2)*col->getXLEN()-1));
            if (i==0){
                targetBC[0] = nproc;
                npointssent[0]=0;
                nmessageBC++;
		// AMR, ME
		BCSide[0]=0;
		targetBOTTOM++;
		// end AMR, ME
            }
            if(nproc != targetBC[j]){
                j++;
                nmessageBC++;
		// AMR, ME
		BCSide[j]=0;
		targetBOTTOM++;
		// end AMR, ME
                targetBC[j]=nproc; 
                npointssent[j]=0;
            }
                npointssent[j]++;
        }
        
    }
// If the y=Oy+finely+finedy ghost nodes of the fine grid intersect this process
    if (grid->getmodifiedYstart(vct) < Oy+finely+finedy && grid->getmodifiedYend(vct) > Oy+finely+finedy && grid->getmodifiedXstart(vct) < Ox+finelx+finedx && grid->getmodifiedXend(vct) > Ox-finedx) {
        xfirst = max(-finedx, ceil((grid->getmodifiedXstart(vct)-Ox)/finedx)*finedx); 
        xfirstnext = ceil((grid->getmodifiedXend(vct)-Ox)/finedx)*finedx; 
        xlast  = min(finelx+finedx, floor((grid->getmodifiedXend(vct)-xfirst-Ox)/finedx)*finedx+xfirst);
        //if(grid->getXend()-xlast-Ox < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision
	if(fabs(grid->getXend()-xlast-Ox)*0.99999  < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision ; 0.99999 because the obvious choices of ratios often\
 bring to overlapping that this alone does not resolve                                                              
            xlast = xlast - finedx;
        }
        xnnu = floor((xlast-xfirst)/finedx+0.5)+1; 
        yloc = Oy + finely + finedy;
        for (i =0;i<xnnu;i++) {
            xshift = xfirst + i * finedx ;
            xloc = xshift + Ox ;
            ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
            iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
            ixsent[xnnl+i] = ix;
            iysent[xnnl+i] = iy;
            weightBC[xnnl+i][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
            weightBC[xnnl+i][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
            weightBC[xnnl+i][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
            weightBC[xnnl+i][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
            nproc = col->getYLEN()*(floor(xshift/((grid->getNXC()-2.)*finedx))+1)-1+col->getXLEN()*col->getYLEN()*(grid->getLevel()+1); // rank of the proc on the fine grid receiving this point(in MPI_COMM_WORLD)
            nproc = max(nproc, col->getYLEN()*(col->getXLEN()*(grid->getLevel()+1)+1)-1);
            nproc = min(nproc, col->getYLEN()*col->getXLEN()*(grid->getLevel()+2)-1);
            if (i==0){
                if(nmessageBC>0){
                    j++;
                }
                targetBC[j] = nproc;
                npointssent[j]=0;
                nmessageBC++;
		// AMR, ME
		BCSide[j]=1;
		targetTOP++;
		// end AMR, ME
            }
            if(nproc != targetBC[j]){
                j++;
                nmessageBC++;
		// AMR, ME
		BCSide[j]=1;
		targetTOP++;
		// end AMR, ME
                targetBC[j]=nproc; 
                npointssent[j]=0;
            }
                npointssent[j]++;
       }
    }
// If the x=Ox-finedx ghost nodes of the fine grid intersect this process
    if (grid->getmodifiedYstart(vct) < Oy+finely+finedy && grid->getmodifiedYend(vct) > Oy-finedy && grid->getmodifiedXstart(vct) < Ox-finedx && grid->getmodifiedXend(vct) > Ox-finedx) {
        xfirst = max(0., ceil((grid->getYstart()-Oy)/finedy)*finedy); 
        xfirstnext = ceil((grid->getYend()-Oy)/finedy)*finedy; 
        xlast  = min(finely, floor((grid->getYend()-xfirst-Oy)/finedy)*finedy+xfirst);
        if(grid->getYend()-xlast-Oy < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision
            xlast = xlast - finedy;
        }
        ynnl = floor((xlast-xfirst)/finedy+0.5)+1; 
        xloc = Ox-finedx;
        for (i =0;i<ynnl;i++) {
            yshift = xfirst + i * finedy ;
            yloc = yshift + Oy ;
            ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
            iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
            ixsent[xnnl+xnnu+i] = ix;
            iysent[xnnl+xnnu+i] = iy;
            weightBC[xnnl+xnnu+i][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
            weightBC[xnnl+xnnu+i][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
            weightBC[xnnl+xnnu+i][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
            weightBC[xnnl+xnnu+i][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
            nproc =floor(yshift/((grid->getNYC()-2.)*finedy))+col->getXLEN()*col->getYLEN()*(grid->getLevel()+1); // rank of the proc on the fine grid receiving this point(in MPI_COMM_WORLD)
            nproc = max(nproc, col->getYLEN()*col->getXLEN()*(grid->getLevel()+1));
            nproc = min(nproc, col->getYLEN()*(1+col->getXLEN()*(grid->getLevel()+1))-1);
            if (i==0){
                if(nmessageBC>0){
                    j++;
                }
                targetBC[j] = nproc;
		// AMR, ME
		BCSide[j]=2;
		targetLEFT++;
		// end AMR, ME
                npointssent[j]=0;
                nmessageBC++;
            }
            if(nproc != targetBC[j]){
                j++;
                nmessageBC++;
                targetBC[j]=nproc;
		// AMR, ME
		BCSide[j]=2;
		targetLEFT++;
		// end AMR, ME
                npointssent[j]=0;
            }
                npointssent[j]++;
        }
    }
// If the x=Ox+finelx+finedx ghost nodes of the fine grid intersect this process
    if (grid->getmodifiedYstart(vct) < Oy+finely + finedy && grid->getmodifiedYend(vct) > Oy-finedy && grid->getmodifiedXstart(vct) < Ox+finelx+finedx && grid->getmodifiedXend(vct) > Ox+finelx+finedx) {
        xfirst = max(0., ceil((grid->getYstart()-Oy)/finedy)*finedy); 
        xfirstnext = ceil((grid->getYend()-Oy)/finedy)*finedy; 
        xlast  = min(finely, floor((grid->getYend()-xfirst-Oy)/finedy)*finedy+xfirst);
        if(grid->getYend()-xlast-Oy < DBL_EPSILON){ //If the fine subdivision overlap coarse subdivision
            xlast = xlast - finedy;
        }
        ynnu = floor((xlast-xfirst)/finedy+0.5)+1; 
        xloc = Ox+finelx+finedx;
        for (i =0;i<ynnu;i++) {
            yshift = xfirst + i * finedy ;
            yloc = yshift + Oy ;
            ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
            iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
            ixsent[xnnl+xnnu+ynnl+i] = ix;
            iysent[xnnl+xnnu+ynnl+i] = iy;
            weightBC[xnnl+xnnu+ynnl+i][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
            weightBC[xnnl+xnnu+ynnl+i][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
            weightBC[xnnl+xnnu+ynnl+i][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
            weightBC[xnnl+xnnu+ynnl+i][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
            nproc =floor(yshift/((grid->getNYC()-2.)*finedy))+col->getYLEN()*(col->getXLEN()*(grid->getLevel()+2)-1); // rank of the proc on the fine grid receiving this point(in MPI_COMM_WORLD)
            nproc = max(nproc, col->getYLEN()*(col->getXLEN()*(grid->getLevel()+2)-1));
            nproc = min(nproc, col->getYLEN()*col->getXLEN()*(grid->getLevel()+2)-1);
    	    if (i==0){
                if(nmessageBC>0){
                    j++;
                }
                targetBC[j] = nproc;
		// AMR, ME
		BCSide[j]=3;
		targetRIGHT++;
		// end AMR, ME
                npointssent[j]=0;
                nmessageBC++;
            }
            if(nproc != targetBC[j]){
                j++;
                nmessageBC++;
                targetBC[j]=nproc;
		// AMR, ME
		BCSide[j]=3;
		targetRIGHT++;
		// end AMR, ME
                npointssent[j]=0;
            }
                npointssent[j]++;
        }
    }
    for (i=0;i<nmessageBC;i++){
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << " sending BC "<<npointssent[i]<<" points to "<<targetBC[i] <<", side " << BCSide[i] <<endl;
    }
 }

// for debugging purposes                                                                                                                     
 int nmessagerecuBCLEFT=0;
 int nmessagerecuBCRIGHT=0;
 int nmessagerecuBCBOTTOM=0;
 int nmessagerecuBCTOP=0;


if ( grid->getLevel() > 0) {                             //If this grid is considered as fine by another grid
    coarsedx = grid->getDX()*col->getRatio();
    coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel()-1);
    coarsedy = grid->getDY()*col->getRatio();
    coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel()-1);
    Ox = grid->getOx(grid->getLevel()); //Origin x of the grid
    Oy = grid->getOy(grid->getLevel()); //Origin y of the grid
    j=0;
    if(vct->getCoordinates(1) == 0) {
        xfirst = grid->getmodifiedXstart(vct);
        if(vct->getCoordinates(0) == vct->getXLEN()-1) {
            xlast = grid->getXend()+grid->getDX(); 
        }else{
            xlast = grid->getXend()-grid->getDX();
        }
        
        yloc = max(Oy - grid->getDY(),0.);// Because when y < 0, it is on the same proc as if y=0  
        xnnl = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
        for (i=0; i< xnnl; i++) {
            xloc = max(Ox + xfirst +i*grid->getDX(),0.);// Because when x < 0, it is on the same proc as if x=0  
            nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
            if (i==0){
                fromBC[0] = nproc;
                nmessagerecuBC++;
		nmessagerecuBCBOTTOM++;
		BCSidecu[0]=0;
                npointsreceived[0]=0;
                xmrecv[0]=0;
                xprecv[0]=0;
                ymrecv[0]=0;
                yprecv[0]=0;
            }
            if(nproc != fromBC[j]){
                j++;
                nmessagerecuBC++;
		nmessagerecuBCBOTTOM++;
		BCSidecu[j]=0;
                fromBC[j]=nproc; 
                npointsreceived[j]=0;
                xmrecv[j]=0;
                xprecv[j]=0;
                ymrecv[j]=0;
                yprecv[j]=0;
            }
                npointsreceived[j]++;
                ixmrecvfirst[j]=floor((xloc-grid->getXstart()-Ox)/grid->getDX()+0.5)-xmrecv[j]+1; //+1 because 0 is the ghost cell at x=Xstart-dx.
                xmrecv[j]++;

        }
   } 
    if(vct->getCoordinates(1) == vct->getYLEN()-1) {
        xfirst = grid->getmodifiedXstart(vct);
        if(vct->getCoordinates(0) == vct->getXLEN()-1) {
            xlast = grid->getXend()+grid->getDX(); 
        }else{
            xlast = grid->getXend()-grid->getDX();
        }
        yloc = min(Oy + grid->getYend()+grid->getDY(),coarsely);  
        xnnu = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
        for (i=0; i< xnnu; i++) {
            xloc = max(Ox + xfirst +i*grid->getDX(),0.);
            nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
            if (i==0){
                if(nmessagerecuBC>0){
                    j++;
                }
                fromBC[j] = nproc;
                nmessagerecuBC++;
		nmessagerecuBCTOP++;
		BCSidecu[j]=1;
                npointsreceived[j]=0;
                xmrecv[j]=0;
                xprecv[j]=0;
                ymrecv[j]=0;
                yprecv[j]=0;
            }
            if(nproc != fromBC[j]){
                j++;
                nmessagerecuBC++;
		nmessagerecuBCTOP++;
		BCSidecu[j]=1;
                fromBC[j]=nproc; 
                npointsreceived[j]=0;
                xmrecv[j]=0;
                xprecv[j]=0;
                ymrecv[j]=0;
                yprecv[j]=0;
            }
                npointsreceived[j]++;
                ixprecvfirst[j]=floor((xloc-grid->getXstart()-Ox)/grid->getDX()+0.5)-xprecv[j]+1;
                xprecv[j]++;

        }
   } 
    if(vct->getCoordinates(0) == 0) {
        xfirst = grid->getYstart();
        if(vct->getCoordinates(1) == vct->getYLEN()-1) {
            xlast = grid->getYend(); 
        }else{
            xlast = grid->getYend()-grid->getDY();
        }
        xloc = max(Ox -grid->getDX(),0.);  
        ynnl = floor((xlast-xfirst)/grid->getDY()+0.5)+1; 
        for (i=0; i< ynnl; i++) {
            yloc = Oy + xfirst +i*grid->getDY();
            nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
            if (i==0){
                if(nmessagerecuBC>0){
                    j++;
                }
                fromBC[j] = nproc;
                nmessagerecuBC++;
		nmessagerecuBCLEFT++;
		BCSidecu[j]=2;
                npointsreceived[j]=0;
                xmrecv[j]=0;
                xprecv[j]=0;
                ymrecv[j]=0;
                yprecv[j]=0;
            }
            if(nproc != fromBC[j]){
                j++;
                nmessagerecuBC++;
		nmessagerecuBCLEFT++;
		BCSidecu[j]=2;
                fromBC[j]=nproc; 
                npointsreceived[j]=0;
                xmrecv[j]=0;
                xprecv[j]=0;
                ymrecv[j]=0;
                yprecv[j]=0;
            }
                npointsreceived[j]++;
                iymrecvfirst[j]=floor((yloc-grid->getYstart()-Oy)/grid->getDY()+0.5)-ymrecv[j]+1;
                ymrecv[j]++;

        }
   } 
    if(vct->getCoordinates(0) == vct->getXLEN()-1) {
        xfirst = grid->getYstart();
        if(vct->getCoordinates(1) == vct->getYLEN()-1) {
            xlast = grid->getYend(); 
        }else{
            xlast = grid->getYend()-grid->getDY();
        }
        xloc = min(Ox +grid->getXend()+grid->getDX(),coarselx);  
        ynnu = floor((xlast-xfirst)/grid->getDY()+0.5)+1; 
        for (i=0; i< ynnu; i++) {
            yloc = Oy + xfirst +i*grid->getDY();
            nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
            if (i==0){
                if(nmessagerecuBC>0){
                    j++;
                }
                fromBC[j] = nproc;
                nmessagerecuBC++;
		nmessagerecuBCRIGHT++;
		BCSidecu[j]=3;
                npointsreceived[j]=0;
                xmrecv[j]=0;
                xprecv[j]=0;
                ymrecv[j]=0;
                yprecv[j]=0;
            }
            if(nproc != fromBC[j]){
                j++;
                nmessagerecuBC++;
		nmessagerecuBCRIGHT++;
		BCSidecu[j]=3;
                fromBC[j]=nproc; 
                npointsreceived[j]=0;
                xmrecv[j]=0;
                xprecv[j]=0;
                ymrecv[j]=0;
                yprecv[j]=0;
            }
                npointsreceived[j]++;
                iyprecvfirst[j]=floor((yloc-grid->getYstart()-Oy)/grid->getDY()+0.5)-yprecv[j]+1;
                yprecv[j]++;

        }
   } 

    for (i=0;i<nmessagerecuBC;i++){
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() <<  ": receiving BC "<<npointsreceived[i]<<" points from "<<fromBC[i] <<", side " << BCSidecu[i]<<endl;
    }

}

//some checks, from ME
// check that the number of sends match the number of receives                                                                                                            

 int npointssentBOTTOM=0;
 int npointssentTOP=0;
 int npointssentLEFT=0;
 int npointssentRIGHT=0;

 for (int i=0; i< nmessageBC; i++)
   {
     if (BCSide[i]==0)
       npointssentBOTTOM+= npointssent[i];

     if (BCSide[i]==1)
       npointssentTOP+= npointssent[i];

     if (BCSide[i]==2)
       npointssentLEFT+= npointssent[i];

     if (BCSide[i]==3)
       npointssentRIGHT+= npointssent[i];
   }

 int npointsreceivedBOTTOM=0;
 int npointsreceivedTOP=0;
 int npointsreceivedLEFT=0;
 int npointsreceivedRIGHT=0;
 
 for (int i=0; i< nmessagerecuBC; i++)
   {
     if (BCSidecu[i]==0)
       npointsreceivedBOTTOM+=npointsreceived[i];
 
     if (BCSidecu[i]==1)
       npointsreceivedTOP+=npointsreceived[i];

     if (BCSidecu[i]==2)
       npointsreceivedLEFT+=npointsreceived[i];

     if (BCSidecu[i]==3)
       npointsreceivedRIGHT+=npointsreceived[i];
   }

 int ALL_npointssentBOTTOM;
 int ALL_npointssentTOP;
 int ALL_npointssentLEFT;
 int ALL_npointssentRIGHT;

 int ALL_npointsreceivedBOTTOM;
 int ALL_npointsreceivedTOP;
 int ALL_npointsreceivedLEFT;
 int ALL_npointsreceivedRIGHT;
 
 MPI_Allreduce ( &npointssentBOTTOM, &ALL_npointssentBOTTOM, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
 MPI_Allreduce ( &npointssentTOP, &ALL_npointssentTOP, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
 MPI_Allreduce ( &npointssentLEFT, &ALL_npointssentLEFT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
 MPI_Allreduce ( &npointssentRIGHT, &ALL_npointssentRIGHT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

 MPI_Allreduce ( &npointsreceivedBOTTOM, &ALL_npointsreceivedBOTTOM, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
 MPI_Allreduce ( &npointsreceivedTOP, &ALL_npointsreceivedTOP, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
 MPI_Allreduce ( &npointsreceivedLEFT, &ALL_npointsreceivedLEFT, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
 MPI_Allreduce ( &npointsreceivedRIGHT, &ALL_npointsreceivedRIGHT, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
 
 int ALL_TARGETS;
 int ALL_RECEIVERS;

 //appropriate initializations done                                                                                                                                         
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

 // to remove

 for (int i=0; i< nmessageBC; i++)
   {
     if (BCSide[i]==0)
       {
         cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << " FIELDS: sends BC bottom to R" << targetBC[i] <<" ns "<<ns <<endl;
       }
   }

 for (int i=0; i<nmessagerecuBC; i++)
   {
     if (BCSidecu[i]==0)
       {
         cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << " FIELDS: receives BC bottom from R" << fromBC[i] <<" ns "<<ns<<endl;
       }
   }
 MPI_Barrier(vct->getCART_COMM_TOTAL());
 // end to remove

 if (ALL_TARGETS!=ALL_RECEIVERS)
   {
     if (vct->getCartesian_rank_COMMTOTAL()==0)
       {
	 cout <<"FIELDS: ALL_TARGETS= " <<ALL_TARGETS <<"!= ALL_RECEIVERS= "<<ALL_RECEIVERS <<endl;
	 cout <<"FIELDS: ALLTARGETS_LEFT= " << ALLTARGETS_LEFT << " ALLRECEIVERS_LEFT= " << ALLRECEIVERS_LEFT <<endl;
	 cout <<"FIELDS: ALLTARGETS_RIGHT= " << ALLTARGETS_RIGHT << " ALLRECEIVERS_RIGHT= "<< ALLRECEIVERS_RIGHT <<endl;
	 cout <<"FIELDS: ALLTARGETS_BOTTOM= " << ALLTARGETS_BOTTOM << " ALLRECEIVERS_BOTTOM= "<< ALLRECEIVERS_BOTTOM <<endl;
	 cout <<"FIELDS: ALLTARGETS_TOP= " << ALLTARGETS_TOP << " ALLRECEIVERS_TOP= "<< ALLRECEIVERS_TOP <<endl;
	 return -1;
       }
   }
 
 if (ALL_npointssentBOTTOM!= ALL_npointsreceivedBOTTOM || ALL_npointssentTOP!= ALL_npointsreceivedTOP || ALL_npointssentLEFT!= ALL_npointsreceivedLEFT || ALL_npointssentRIGHT!= ALL_npointsreceivedRIGHT )
   {
     if (vct->getCartesian_rank_COMMTOTAL()==0)
       {
	 cout << "Mismatch in points sent and received for interpolation!!!"<<endl;
	 cout << "ALL_npointssentBOTTOM: " << ALL_npointssentBOTTOM << ", ALL_npointsreceivedBOTTOM: " << ALL_npointsreceivedBOTTOM <<endl;
	 cout << "ALL_npointssentTOP: " << ALL_npointssentTOP<< ", ALL_npointsreceivedTOP: " << ALL_npointsreceivedTOP<<endl;
	 cout << "ALL_npointssentLEFT: " << ALL_npointssentLEFT<< ", ALL_npointsreceivedLEFT: " << ALL_npointsreceivedLEFT<<endl;
	 cout << "ALL_npointssentRIGHT: " << ALL_npointssentRIGHT<< ", ALL_npointsreceivedRIGHT: " << ALL_npointsreceivedRIGHT<<endl;
	 return -1;
       }
   }

 // end checks, by ME
 return 1;
}*/
/*//Initialize weight for projection of refined fields between grids
inline void EMfields::initWeightProj(VirtualTopology *vct, Grid *grid, CollectiveIO* col){

MPI_Status status;
int X,Y,Xlocal,Ylocal,i,j,k,ix,iy,nprocfirst,nproclast,nybloc,ierr;
double coarsedx,coarsedy,coarselx,coarsely,finelx,finely,Ox,Oy,finedx,finedy,correctionx,correctiony;
double xloc, yloc,xfirst,xlast,yfirst,ylast,xlocalfirst,xlocallast,ylocalfirst,ylocallast,xstop,ystop;
int start, step;


if ( grid->getLevel() > 0) {                             //If this grid is considered as fine by another grid
    coarsedx = grid->getDX()*col->getRatio();
    coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel()-1);
    coarsedy = grid->getDY()*col->getRatio();
    coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel()-1);
    Ox = grid->getOx(grid->getLevel()); //Origin x of finer grid
    Oy = grid->getOy(grid->getLevel()); //Origin y of finer grid

    xfirst = floor((grid->getXstart()+Ox)/coarsedx)*coarsedx;//Coordinate of the first point in x direction of the coarse grid influenced by the projection of this proc expressed in the coarse grid frame
    X= floor((xfirst)/coarsedx)/(grid->getNXC()-2); // X coordinate of the proc of the coarse grid intersecting the origin of this proc   
    yfirst = floor((grid->getYstart()+Oy)/coarsedy)*coarsedy;
    Y= floor((yfirst)/coarsedy)/(grid->getNYC()-2); // X coordinate of the proc of the coarse grid intersecting the origin of this proc   
    targetProj[0] = X*vct->getYLEN()+Y+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
    targetProj[1] = (X+1)*vct->getYLEN()+Y+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
    targetProj[2] = X*vct->getYLEN()+Y+1+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
    targetProj[3] = (X+1)*vct->getYLEN()+Y+1+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
	if (vct->getXright_neighbor()==MPI_PROC_NULL){
        xlast = ceil((grid->getXend()+Ox)/coarsedx)*coarsedx;
        lastindicex = grid->getNXN()-2;
    } else {
        xlast = ceil((grid->getXend()+Ox)/coarsedx-1./col->getRatio())*coarsedx;
        lastindicex = grid->getNXN()-3;
    }
	if (vct->getYright_neighbor()==MPI_PROC_NULL){
        ylast = ceil((grid->getYend()+Oy)/coarsedy)*coarsedy;
        lastindicey = grid->getNYN()-2;
    } else {
        ylast = ceil((grid->getYend()+Oy)/coarsedy-1./col->getRatio())*coarsedy;
        lastindicey = grid->getNYN()-3;
    }
    xstop = (X+1)*(grid->getNXC()-2)*coarsedx; //End of domain of processor whith coordinate X
    ystop = (Y+1)*(grid->getNYC()-2)*coarsedy; //End of domain of processor with coordinate Y
    //nxmsend is the number of points sent in the x direction to the left part of x=xstop, and nxp to the right part of x=xstop
    if (xstop < xlast){
        nmessageProj=nmessageProj*2;
        nxpsend = floor((xlast-xstop)/coarsedx+0.5) + 1;
        nxmsend = floor((xlast-xfirst)/coarsedx+0.5)+1-nxpsend;
    } else {
        nxmsend = floor((xlast-xfirst)/coarsedx+0.5)+1;
        nxpsend = 0;
    }
    //nym is the number of points sent in the y direction to the lower part of y=ystop, and nxp to the upper part of y=ystop
    if (ystop < ylast){
        nmessageProj=nmessageProj*2;
        nypsend = floor((ylast-ystop)/coarsedy+0.5) + 1;
        nymsend = floor((ylast-yfirst)/coarsedy+0.5)+1-nypsend;
    } else {
        nymsend = floor((ylast-yfirst)/coarsedy+0.5) + 1;
        nypsend = 0 ;
    }

    for (ix =0; ix<nxmsend+nxpsend; ix++){
        for (iy =0; iy<nymsend+nypsend; iy++){
            normalizeProj[ix][iy][0]=0.;
        }
    }

    for (i =0;i<lastindicex;i++) {
        xloc = grid->getXstart() + i * grid->getDX() + Ox -xfirst;
        ixsentProj[i] = int(floor(xloc/coarsedx));
        for (j =0;j<lastindicey;j++) {
            yloc = grid->getYstart() + j * grid->getDY() + Oy -yfirst;
            iysentProj[j] = int(floor(yloc/coarsedy));

            weightProj[i][j][3][0] = ((xloc - ixsentProj[i]*coarsedx)/coarsedx)*((yloc - iysentProj[j]*coarsedy)/coarsedy); // weight +:+
            weightProj[i][j][2][0] = ((xloc - ixsentProj[i]*coarsedx)/coarsedx)*(((iysentProj[j]+1)*coarsedy - yloc)/coarsedy); //weight +:-
            weightProj[i][j][1][0] = (((ixsentProj[i]+1)*coarsedx - xloc)/coarsedx)*((yloc - iysentProj[j]*coarsedy)/coarsedy); // weight -:+
            weightProj[i][j][0][0] = (((ixsentProj[i]+1)*coarsedx - xloc)/coarsedx)*(((iysentProj[j]+1)*coarsedy - yloc)/coarsedy); // weight -:-

            normalizeProj[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0];
            normalizeProj[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0];
            normalizeProj[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0];
            normalizeProj[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0];
        }
    }

    cout<< vct->getCartesian_rank_COMMTOTAL() << " sends to target0 = " << targetProj[0]<<" nmessageProj= "<<nmessageProj<<" nxmsend= "<<nxmsend<<" nxpsend= "<<nxpsend<<" nymsend= "<<nymsend<<" nypsend= "<<nypsend<< "xfirst = "<<xfirst<<"xlast = "<<xlast<<" yfirst= "<<yfirst<<" ylast = "<<ylast<<" xstop= "<<xstop<<endl;
    start = 0;
    for (i=0;i<nxmsend*nymsend;i++){
        bufferProjsend[start] = normalizeProj[i/nymsend][i%nymsend][0];
        start++;
     //   cout << bufferProjsend[i]<<" ";
    }
       //cout << endl;
       cout << vct->getCartesian_rank_COMMTOTAL() << "sending to "<< targetProj[0] << endl;
    ierr = MPI_Send(bufferProjsend,nxmsend*nymsend,MPI_DOUBLE,targetProj[0],1,MPI_COMM_WORLD);

    if(nxpsend > 0){
        step=start;
        for (i=0;i<nxpsend*nymsend;i++){
            bufferProjsend[start] = normalizeProj[nxmsend+i/nymsend][i%nymsend][0];
            start++;
        }
        ierr = MPI_Send(bufferProjsend+step,nxpsend*nymsend,MPI_DOUBLE,targetProj[1],1,MPI_COMM_WORLD);
    }

    if(nypsend > 0){
        step=start;
        for (i=0;i<nxmsend*nypsend;i++){
            bufferProjsend[start] = normalizeProj[i/nypsend][nymsend+i%nypsend][0];
            start++;
        }
        ierr = MPI_Send(bufferProjsend+step,nxmsend*nypsend,MPI_DOUBLE,targetProj[2],1,MPI_COMM_WORLD);
    }
    if(nypsend > 0 && nxpsend > 0){
        step=start;
        for (i=0;i<nxpsend*nypsend;i++){
            bufferProjsend[start] = normalizeProj[nxmsend+i/nypsend][nymsend+i%nypsend][0];
            start++;
        }
        ierr = MPI_Send(bufferProjsend+step,nxpsend*nypsend,MPI_DOUBLE,targetProj[3],1,MPI_COMM_WORLD);
    }


}
if (grid->getLevel() < vct->getNgrids()-1) {
    Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid
    Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid
    coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel());
    coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel());
    finelx = col->getLx()/pow(col->getRatio(),grid->getLevel()+1);
    finely = col->getLy()/pow(col->getRatio(),grid->getLevel()+1);
    finedx = grid->getDX()/col->getRatio();
    finedy = grid->getDY()/col->getRatio();
    for (ix =0; ix<nxn; ix++){
        for (iy =0; iy<nyn; iy++){
            normalizerecvProj[ix][iy][0]=0.; 
        }
    }
    if(grid->getXstart()<Ox+finelx && grid->getXend()>Ox && grid->getYstart()<Oy+finely && grid->getYend()>Oy){
        xfirst = max(floor(Ox/grid->getDX())*grid->getDX(),grid->getXstart());//x of the first point of this proc being updated by the finer domain (in the coarse domain frame) 
        yfirst = max(floor(Oy/grid->getDY())*grid->getDY(),grid->getYstart());//y of the first point of this proc being updated by the finer domain (in the coarse domain frame)
        ixrecvfirstProjglobal=floor((xfirst-grid->getXstart())/grid->getDX()+0.5)+1;//ix of the first point of this proc being updated by the finer domain (in the coarse domain frame) 
        iyrecvfirstProjglobal=floor((yfirst-grid->getYstart())/grid->getDY()+0.5)+1;//iy of the first point of this proc being updated by the finer domain (in the coarse domain frame) 

        xlast = min(ceil((Ox+finelx)/grid->getDX())*grid->getDX(),grid->getXend());//x of the last point of this proc intersecting the finer domain (in the coarse domain frame)
        ylast = min(ceil((Oy+finely)/grid->getDY())*grid->getDY(),grid->getYend());//y of the last point of this proc intersecting the finer domain (in the coarse domain frame)
        ixrecvlastProjglobal=floor((xlast-grid->getXstart())/grid->getDX()+0.5)+1;//ix of the first point of this proc being updated by the finer domain (in the coarse domain frame) 
        iyrecvlastProjglobal=floor((ylast-grid->getYstart())/grid->getDY()+0.5)+1;//iy of the first point of this proc being updated by the finer domain (in the coarse domain frame) 
        //Correction not to update coarse grid close to the fine/coarse boundary
        if (xfirst <= Ox)
	  ixrecvfirstProjglobal += 5; // skipped cells 
        if (yfirst <= Oy)
            iyrecvfirstProjglobal += 5; 
        if (xlast >= Ox+finelx)
            ixrecvlastProjglobal -= 5; 
        if (ylast >= Oy+finely)
            iyrecvlastProjglobal -= 5; 
        nxmrecv=floor((xlast-xfirst)/grid->getDX()+0.5)+1;
        nymrecv=floor((ylast-yfirst)/grid->getDY()+0.5)+1;


    correctionx = 0.;
    correctiony = 0.;
    if (xfirst < Ox)
        correctionx = grid->getDX();
    if (yfirst < Oy)
        correctiony = grid->getDY();
      
     nprocfirst = min(floor((yfirst+correctiony-Oy)/((grid->getNYC()-2.)*finedy)),col->getYLEN()-1.)+col->getYLEN()*min(floor((xfirst+correctionx-Ox)/((grid->getNXC()-2.)*finedx)),col->getXLEN()-1.)+col->getXLEN()*col->getYLEN()*(grid->getLevel()+1); // rank of the proc on the fine grid sending this point(in MPI_COMM_WORLD)


    correctionx = 0.;
    correctiony = 0.;
    if (xlast > Ox+finelx)
        correctionx = grid->getDX();
    if (ylast > Oy+finely)
        correctiony = grid->getDY();

    nproclast = min(floor((ylast-correctiony-Oy)/((grid->getNYC()-2.)*finedy)),col->getYLEN()-1.)+col->getYLEN()*min(floor((xlast-correctionx-Ox)/((grid->getNXC()-2.)*finedx)),col->getXLEN()-1.)+col->getXLEN()*col->getYLEN()*( grid->getLevel()+1); // rank of the proc on the fine grid sending this point(in MPI_COMM_WORLD)
    X = (nproclast-col->getYLEN()*col->getXLEN()*(grid->getLevel()+1))/col->getYLEN();// X position of nproclast
    Y = nproclast%col->getYLEN();// Y position of nproclast
    cout << vct->getCartesian_rank_COMMTOTAL() << " before correction x first = "<<xfirst<<" y first = "<<yfirst<<" nproc first = "<<nprocfirst<<" nproc last = "<<nproclast<<" X= "<<X<<" Y= "<<Y << endl;

    if(xlast == grid->getXend() && grid->getXend() < Ox+finelx){
        nproclast -= col->getYLEN(); 
    }
    if(ylast == grid->getYend() && grid->getYend() < Oy+finely){
        nproclast -= 1; 
    }
    cout << vct->getCartesian_rank_COMMTOTAL() << " nproc first = "<<nprocfirst<<" nproc last = "<<nproclast << endl;

   X=nproclast/col->getYLEN()-nprocfirst/col->getYLEN()+1;//Number of proc sending to this one in the X direction
   Y=nproclast%col->getYLEN()-nprocfirst%col->getYLEN()+1;//Number of procs sending to this one in the Y direction
   nmessagerecvProj = X*Y;
   for (i=0;i<nmessagerecvProj;i++){
       fromProj[i]= nprocfirst+i%Y+(i/Y)*col->getYLEN();

       // This correction account for the fact that fine procs do not project their last points unless they are on the last row(y) or last colum(x).
       correctionx = 0.;
       correctiony = 0.;
       if ((fromProj[i]-col->getXLEN()*col->getYLEN()*(grid->getLevel()+1))/col->getYLEN() < col->getXLEN()-1)
           correctionx = -1./col->getRatio();
       if (fromProj[i]%col->getYLEN() < col->getYLEN()-1)
           correctiony = -1./col->getRatio();

          Xlocal=(fromProj[i]-col->getXLEN()*col->getYLEN()*( grid->getLevel()+1))/col->getYLEN();
          Ylocal=(fromProj[i]-col->getXLEN()*col->getYLEN()*( grid->getLevel()+1))%col->getYLEN();

       xlocalfirst = max(          grid->getXstart(), floor(   (Ox+Xlocal*finelx/col->getXLEN())/grid->getDX()   ) *  grid->getDX()                            );//First updated point by proc fromProj[i]
       ylocalfirst = max(          grid->getYstart(), floor(   (Oy+Ylocal*finely/col->getYLEN())/grid->getDY()   ) *  grid->getDY()                            );//First updated point by proc fromProj[i]
       xlocallast = min(          grid->getXend(), ceil(   (Ox+(Xlocal+1)*finelx/col->getXLEN())/grid->getDX()+correctionx   ) *  grid->getDX()                            );//First updated point by proc fromProj[i]
       ylocallast = min(          grid->getYend(), ceil(   (Oy+(Ylocal+1)*finely/col->getYLEN())/grid->getDY()+correctiony   ) *  grid->getDY()                            );//First updated point by proc fromProj[i]
       nxrecvProj[i]=floor((xlocallast-xlocalfirst)/grid->getDX()+1.5);
       nyrecvProj[i]=floor((ylocallast-ylocalfirst)/grid->getDY()+1.5);
       npointsreceivedProj[i]=nxrecvProj[i]*nyrecvProj[i];
       ixrecvfirstProj[i]=floor((xlocalfirst-grid->getXstart())/grid->getDX()+0.5)+1;
       iyrecvfirstProj[i]=floor((ylocalfirst-grid->getYstart())/grid->getDY()+0.5)+1;

       cout << vct->getCartesian_rank_COMMTOTAL() << "receiving from "<< fromProj[i] << " number of points = "<<npointsreceivedProj[i]<< "xlocalfirst= "<< xlocalfirst<<" xlocallast= "<<xlocallast<<" correctionx = "<<correctionx<<" ylocalfirst= "<< ylocalfirst<<" ylocallast = "<<ylocallast<<endl;
       ierr = MPI_Recv(bufferProj,npointsreceivedProj[i],MPI_DOUBLE,fromProj[i],1,MPI_COMM_WORLD, &status);
       for (j=0;j<npointsreceivedProj[i];j++){
           ix = ixrecvfirstProj[i] + j/nyrecvProj[i];
           iy = iyrecvfirstProj[i] + j%nyrecvProj[i];
           normalizerecvProj[ix][iy][0] += bufferProj[j];
       }
       
   }
   }
   // Completing normalizerecvProj with the normalization gathered by the other coarse procs
    if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //Send a message to proc  coorinateX-1
        for (i=1;i<nyn-1;i++){
                bufferProjsend[i-1] = normalizerecvProj[1][i][0];
            }
        if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
        //also receive a message
            cout << vct->getCartesian_rank_COMMTOTAL() << "SENd RECEIVE" <<endl;
            ierr = MPI_Sendrecv(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
            cout << vct->getCartesian_rank_COMMTOTAL() << "SENd only" <<endl;
            ierr = MPI_Send(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM());
        }
    } else {
        //Do not send ...
          if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
          //... but receive only
            cout << vct->getCartesian_rank_COMMTOTAL() << "RECEIVE only" <<endl;
              ierr = MPI_Recv(bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(), &status);
          }
    }
    if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
        //If a message from proc coorinateX+1 has been received
        for (i=1;i<nyn-1;i++){
                normalizerecvProj[nxn-2][i][0] += bufferProjrecv[i-1];
        }
    //Send a message to proc coordinateX+1
        for(i=1;i<nyn-1;i++){
            bufferProjsend[i-1] = normalizerecvProj[nxn-2][i][0];
        }
        if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //also receive a message
        ierr = MPI_Sendrecv(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
        cout << vct->getCartesian_rank_COMMTOTAL() <<" send to "<< vct->getXright_neighbor() << endl;
        ierr = MPI_Send(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM());
        }
    
    } else {
          //Do not send ...
           if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
           //...receive only
           cout << vct->getCartesian_rank_COMMTOTAL() <<" receive from "<< vct->getXleft_neighbor() << endl;
           ierr = MPI_Recv(bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(), &status);
           }
    }

    if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //A message from proc coordinateX-1 has been received
        for (i=1;i<nyn-1;i++){
            normalizerecvProj[1][i][0] =  bufferProjrecv[i-1];
        }
    }

    // In the Y direction now
    if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //Send a message to proc  coorinateY-1
        for (i=1;i<nxn-1;i++){
                bufferProjsend[i-1] = normalizerecvProj[i][1][0];
            }
        if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
        //also receive a message
            ierr = MPI_Sendrecv(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
            ierr = MPI_Send(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM());
        }
    } else {
        //Do not send ...
          if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
          //... but receive only
              ierr = MPI_Recv(bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(), &status);
          }
    }
    if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
        //If a message from proc coorinateY+1 has been received
        for (i=1;i<nxn-1;i++){
                normalizerecvProj[i][nyn-2][0] += bufferProjrecv[i-1];
        }
    //Send a message to proc coordinateY+1
        for(i=1;i<nxn-1;i++){
            bufferProjsend[i-1] = normalizerecvProj[i][nyn-2][0];
        }
        if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //also receive a message
        ierr = MPI_Sendrecv(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
        cout << vct->getCartesian_rank_COMMTOTAL() <<" send to "<< vct->getYright_neighbor() << endl;
        ierr = MPI_Send(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM());
        }
    
    } else {
          //Do not send ...
           if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
           //...receive only
           cout << vct->getCartesian_rank_COMMTOTAL() <<" receive from "<< vct->getYleft_neighbor() << endl;
           ierr = MPI_Recv(bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(), &status);
           }
    }
    if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //A message from proc coordinateX-1 has been received
        for (i=1;i<nxn-1;i++){
            normalizerecvProj[i][1][0] =  bufferProjrecv[i-1];
        }
    }
    

         //cout << "Proc " << vct->getCartesian_rank_COMMTOTAL() <<"receives from nproc first= " << nprocfirst<<" nproc last = "<<nproclast << "xfirst = "<<xfirst << "xlast = "<< xlast<< "yfirst = "<<yfirst << "ylast = "<< ylast << "nproc last = " << nproclast<< endl;
   for (i=0;i<nmessagerecvProj;i++){
       cout << " fromProj " << fromProj[i] << " nxrecvProj = "<<nxrecvProj[i]<<" nyrecvProj= "<<nyrecvProj[i]<< "ixrecvfirstProj= " <<ixrecvfirstProj[i]<<"iyrecvfirstProj= "<<iyrecvfirstProj[i]<<endl;
   }

   if (vct->getCartesian_rank_COMMTOTAL()==6){
       for(j=0;j<nyn;j++){
           for(i=0;i<nxn;i++){
               cout << normalizerecvProj[i][j][0] << " ";
           }
           cout << endl;
       }
   }


   // }

}
}*/
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

inline void EMfields::communicateGhostP2GOS(int ns,int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct){

  communicateInterpOS(nxn,nyn,ns,rhons,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,Jxs,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,Jys,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,Jzs,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,pXXsn,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,pXYsn,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,pXZsn,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,pYYsn,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,pYZsn,0,0,0,0,dx,dy,vct);
  communicateInterpOS(nxn,nyn,ns,pZZsn,0,0,0,0,dx,dy,vct);

  communicateNodeOS(nxn, nyn, rhons, ns,vct);
  communicateNodeOS(nxn, nyn, Jxs, ns,vct);
  communicateNodeOS(nxn, nyn, Jys, ns,vct);
  communicateNodeOS(nxn, nyn, Jzs, ns,vct);
  communicateNodeOS(nxn, nyn, pXXsn, ns,vct);
  communicateNodeOS(nxn, nyn, pXYsn, ns,vct);
  communicateNodeOS(nxn, nyn, pXZsn, ns,vct);
  communicateNodeOS(nxn, nyn, pYYsn, ns,vct);
  communicateNodeOS(nxn, nyn, pYZsn, ns,vct);
  communicateNodeOS(nxn, nyn, pZZsn, ns,vct);

}

/** adjust densities on boundaries that are not periodic */
inline  void EMfields::adjustNonPeriodicDensities(VirtualTopology *vct, int cycle){
 
  int cartesian_rank_total;
  MPI_Comm_rank   (vct->getCART_COMM_TOTAL(), &cartesian_rank_total);
  bool RefinedLevel;
  if (cartesian_rank_total > vct->getXLEN() * vct->getYLEN() -1)
    {
      RefinedLevel= true;
    }
  else
    {
      RefinedLevel= false;
    }
        
  //if (cycle ==0 || vct->getRefLevelAdj()==0 || (vct->getRefLevelAdj()==1 && !RefinedLevel))
  if (RefinedLevel)
    {
      for (int is=0; is < ns; is++)
	{
	  if (vct->getXleft_neighbor()==MPI_PROC_NULL)
	    {
	      for (int j=1; j < nyn-1; j++)   
		{
		  rhons[is][0][j][0]*= 2  ;
		  Jxs[is][0][j][0]  *= 2  ;
		  Jys[is][0][j][0]  *= 2  ;
		  Jzs[is][0][j][0]  *= 2  ;
		  pXXsn[is][0][j][0]  *=2 ;
		  pXYsn[is][0][j][0]  *=2 ;
		  pXZsn[is][0][j][0]  *=2 ;
		  pYYsn[is][0][j][0]  *=2 ;
		  pYZsn[is][0][j][0]  *=2 ;
		  pZZsn[is][0][j][0]  *=2 ;
		}
	    }
	  if (vct->getXright_neighbor()==MPI_PROC_NULL)
	    {
	      for (int j=1; j < nyn-1; j++)
		{
		  rhons[is][nxn-1][j][0]*= 2  ;
		  Jxs[is][nxn-1][j][0]  *= 2  ;
		  Jys[is][nxn-1][j][0]  *= 2  ;
		  Jzs[is][nxn-1][j][0]  *= 2  ;
		  pXXsn[is][nxn-1][j][0]  *=2 ;
		  pXYsn[is][nxn-1][j][0]  *=2 ;
		  pXZsn[is][nxn-1][j][0]  *=2 ;
		  pYYsn[is][nxn-1][j][0]  *=2 ;
		  pYZsn[is][nxn-1][j][0]  *=2 ;
		  pZZsn[is][nxn-1][j][0]  *=2 ;
		}
	    }
	  if (vct->getYleft_neighbor()==MPI_PROC_NULL)
	    {
	      for (int i=1; i < nxn-1; i++)
		{
		  rhons[is][i][0][0]*= 2  ;
		  Jxs[is][i][0][0]  *= 2  ;
		  Jys[is][i][0][0]  *= 2  ;
		  Jzs[is][i][0][0]  *= 2  ;
		  pXXsn[is][i][0][0]  *=2 ;
		  pXYsn[is][i][0][0]  *=2 ;
		  pXZsn[is][i][0][0]  *=2 ;
		  pYYsn[is][i][0][0]  *=2 ;
		  pYZsn[is][i][0][0]  *=2 ;
		  pZZsn[is][i][0][0]  *=2 ;
		}
	    }
	  if (vct->getYright_neighbor()==MPI_PROC_NULL)
	    {
	      for (int i=1; i < nxn-1; i++)
		{
		  rhons[is][i][nyn-1][0]*= 2  ;
		  Jxs[is][i][nyn-1][0]  *= 2  ;
		  Jys[is][i][nyn-1][0]  *= 2  ;
		  Jzs[is][i][nyn-1][0]  *= 2  ;
		  pXXsn[is][i][nyn-1][0]  *=2 ;
		  pXYsn[is][i][nyn-1][0]  *=2 ;
		  pXZsn[is][i][nyn-1][0]  *=2 ;
		  pYYsn[is][i][nyn-1][0]  *=2 ;
		  pYZsn[is][i][nyn-1][0]  *=2 ;
		  pZZsn[is][i][nyn-1][0]  *=2 ;
		}
	    } 
	  // now take care of the corners
	  
	  // bottom left
	  if (vct->getXleft_neighbor()==MPI_PROC_NULL && vct->getYleft_neighbor()==MPI_PROC_NULL)
	    {//bottom left corner proc 
	      rhons[is][0][0][0]*= 4  ;
	      Jxs[is][0][0][0]  *= 4  ;
	      Jys[is][0][0][0]  *= 4  ;
	      Jzs[is][0][0][0]  *= 4  ;
	      pXXsn[is][0][0][0]  *=4 ;
	      pXYsn[is][0][0][0]  *=4 ;
	      pXZsn[is][0][0][0]  *=4 ;
	      pYYsn[is][0][0][0]  *=4 ;
	      pYZsn[is][0][0][0]  *=4 ;
	      pZZsn[is][0][0][0]  *=4 ;
	    }
	  else if (vct->getXleft_neighbor()==MPI_PROC_NULL || vct->getYleft_neighbor()==MPI_PROC_NULL)
	    {//not bottom left corner proc
	      rhons[is][0][0][0]*= 2  ;
	      Jxs[is][0][0][0]  *= 2  ;
	      Jys[is][0][0][0]  *= 2  ;
	      Jzs[is][0][0][0]  *= 2  ;
	      pXXsn[is][0][0][0]  *=2 ;
	      pXYsn[is][0][0][0]  *=2 ;
	      pXZsn[is][0][0][0]  *=2 ;
	      pYYsn[is][0][0][0]  *=2 ;
	      pYZsn[is][0][0][0]  *=2 ;
	      pZZsn[is][0][0][0]  *=2 ;
	    }
	  if (vct->getXleft_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL)
	    {//upper left corner proc                                                                                   
	      rhons[is][0][nyn-1][0]*= 4  ;
	      Jxs[is][0][nyn-1][0]  *= 4  ;
	      Jys[is][0][nyn-1][0]  *= 4  ;
	      Jzs[is][0][nyn-1][0]  *= 4  ;
	      pXXsn[is][0][nyn-1][0]  *=4 ;
	      pXYsn[is][0][nyn-1][0]  *=4 ;
	      pXZsn[is][0][nyn-1][0]  *=4 ;
	      pYYsn[is][0][nyn-1][0]  *=4 ;
	      pYZsn[is][0][nyn-1][0]  *=4 ;
	      pZZsn[is][0][nyn-1][0]  *=4 ;
	    }
	  else if (vct->getXleft_neighbor()==MPI_PROC_NULL || vct->getYright_neighbor()==MPI_PROC_NULL)
	    {//not upper left corner proc                                                                               
	      rhons[is][0][nyn-1][0]*= 2  ;
	      Jxs[is][0][nyn-1][0]  *= 2  ;
	      Jys[is][0][nyn-1][0]  *= 2  ;
	      Jzs[is][0][nyn-1][0]  *= 2  ;
	      pXXsn[is][0][nyn-1][0]  *=2 ;
	      pXYsn[is][0][nyn-1][0]  *=2 ;
	      pXZsn[is][0][nyn-1][0]  *=2 ;
	      pYYsn[is][0][nyn-1][0]  *=2 ;
	      pYZsn[is][0][nyn-1][0]  *=2 ;
	      pZZsn[is][0][nyn-1][0]  *=2 ;
	    }
	  
	  if (vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYleft_neighbor()==MPI_PROC_NULL)
	    {//bottom right corner proc                                                                                 
	      rhons[is][nxn-1][0][0]*= 4  ;
	      Jxs[is][nxn-1][0][0]  *= 4  ;
	      Jys[is][nxn-1][0][0]  *= 4  ;
	      Jzs[is][nxn-1][0][0]  *= 4  ;
	      pXXsn[is][nxn-1][0][0]  *=4 ;
	      pXYsn[is][nxn-1][0][0]  *=4 ;
	      pXZsn[is][nxn-1][0][0]  *=4 ;
	      pYYsn[is][nxn-1][0][0]  *=4 ;
	      pYZsn[is][nxn-1][0][0]  *=4 ;
	      pZZsn[is][nxn-1][0][0]  *=4 ;
	    }
	  else if (vct->getXright_neighbor()==MPI_PROC_NULL || vct->getYleft_neighbor()==MPI_PROC_NULL)
	    {//not bottom right corner proc                                                                             
	      rhons[is][nxn-1][0][0]*= 2  ;
	      Jxs[is][nxn-1][0][0]  *= 2  ;
	      Jys[is][nxn-1][0][0]  *= 2  ;
	      Jzs[is][nxn-1][0][0]  *= 2  ;
	      pXXsn[is][nxn-1][0][0]  *=2 ;
	      pXYsn[is][nxn-1][0][0]  *=2 ;
	      pXZsn[is][nxn-1][0][0]  *=2 ;
	      pYYsn[is][nxn-1][0][0]  *=2 ;
	      pYZsn[is][nxn-1][0][0]  *=2 ;
	      pZZsn[is][nxn-1][0][0]  *=2 ;
	    }
	  if (vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL)
	    {//upper right corner proc                                                                                  
	      rhons[is][nxn-1][nyn-1][0]*= 4  ;
	      Jxs[is][nxn-1][nyn-1][0]  *= 4  ;
	      Jys[is][nxn-1][nyn-1][0]  *= 4  ;
	      Jzs[is][nxn-1][nyn-1][0]  *= 4  ;
	      pXXsn[is][nxn-1][nyn-1][0]  *=4 ;
	      pXYsn[is][nxn-1][nyn-1][0]  *=4 ;
	      pXZsn[is][nxn-1][nyn-1][0]  *=4 ;
	      pYYsn[is][nxn-1][nyn-1][0]  *=4 ;
	      pYZsn[is][nxn-1][nyn-1][0]  *=4 ;
	      pZZsn[is][nxn-1][nyn-1][0]  *=4 ;
	    }
	  else if (vct->getXright_neighbor()==MPI_PROC_NULL || vct->getYright_neighbor()==MPI_PROC_NULL)
	    {//not upper right corner proc                                                                            
	      rhons[is][nxn-1][nyn-1][0]*= 2  ;
	      Jxs[is][nxn-1][nyn-1][0]  *= 2  ;
	      Jys[is][nxn-1][nyn-1][0]  *= 2  ;
	      Jzs[is][nxn-1][nyn-1][0]  *= 2  ;
	      pXXsn[is][nxn-1][nyn-1][0]  *=2 ;
	      pXYsn[is][nxn-1][nyn-1][0]  *=2 ;
	      pXZsn[is][nxn-1][nyn-1][0]  *=2 ;
	      pYYsn[is][nxn-1][nyn-1][0]  *=2 ;
	      pYZsn[is][nxn-1][nyn-1][0]  *=2 ;
	      pZZsn[is][nxn-1][nyn-1][0]  *=2 ;
	    }
      	}// end species
    }
  

}

/** add an amount of charge density to charge density field at node X,Y */
inline void EMfields::addRho(double ***weight, int X, int Y, int Z, int is){
  	rhons[is][X-1][Y-1][0] += weight[0][0][0]*invVOL;
  	rhons[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
	rhons[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
	rhons[is][X][Y][0]    += weight[1][1][0]*invVOL;
}
/** add an amount of charge density to charge density field at node X,Y */
inline void EMfields::addRho_OS(double ***weight, int X, int Y, int Z, int is, Grid * g){
  if (X-1>=0 && Y-1 >=0 )
    rhons[is][X-1][Y-1][0] += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    rhons[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    rhons[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
  rhons[is][X][Y][0]    += weight[1][1][0]*invVOL;
}
/** add an amount of charge density to current density - direction X to current density field on the node*/
inline void EMfields::addJx(double ***weight, int X, int Y, int Z, int is){
	Jxs[is][X-1][Y-1][0] += weight[0][0][0]*invVOL;
	Jxs[is][X-1][Y][0]   += weight[0][1][0]*invVOL;
	Jxs[is][X][Y-1][0]   += weight[1][0][0]*invVOL;
	Jxs[is][X][Y][0]     += weight[1][1][0]*invVOL;}

/** add an amount of charge density to current density - direction X to current density field on the node*/
inline void EMfields::addJx_OS(double ***weight, int X, int Y, int Z, int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    Jxs[is][X-1][Y-1][0] += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    Jxs[is][X-1][Y][0]   += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    Jxs[is][X][Y-1][0]   += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
  Jxs[is][X][Y][0]     += weight[1][1][0]*invVOL;}
/** add an amount of current density - direction Y to current density field on the node */
inline  void EMfields::addJy(double ***weight, int X, int Y, int Z, int is){
	Jys[is][X-1][Y-1][0]  += weight[0][0][0]*invVOL;
	Jys[is][X-1][Y][0]    += weight[0][1][0]*invVOL;
	Jys[is][X][Y-1][0]    += weight[1][0][0]*invVOL;
	Jys[is][X][Y][0]      += weight[1][1][0]*invVOL;
}
/** add an amount of current density - direction Y to current density field on the node */
inline  void EMfields::addJy_OS(double ***weight, int X, int Y, int Z, int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    Jys[is][X-1][Y-1][0]  += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    Jys[is][X-1][Y][0]    += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    Jys[is][X][Y-1][0]    += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    Jys[is][X][Y][0]      += weight[1][1][0]*invVOL;
}
/** add an amount of current density - direction Z to current density field on the node */
inline  void EMfields::addJz(double ***weight, int X, int Y, int Z, int is){
	Jzs[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	Jzs[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	Jzs[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	Jzs[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of current density - direction Z to current density field on the node */
inline  void EMfields::addJz_OS(double ***weight, int X, int Y, int Z, int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    Jzs[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    Jzs[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    Jzs[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    Jzs[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XX to current density field on the node */
inline  void EMfields::addPxx(double ***weight, int X, int Y, int Z,int is){
	pXXsn[is][X-1][Y-1][0]   += weight[0][0][0]*invVOL;
	pXXsn[is][X-1][Y][0]     += weight[0][1][0]*invVOL;
	pXXsn[is][X][Y-1][0]     += weight[1][0][0]*invVOL;
	pXXsn[is][X][Y][0]       += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XX to current density field on the node */
inline  void EMfields::addPxx_OS(double ***weight, int X, int Y, int Z,int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    pXXsn[is][X-1][Y-1][0]   += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    pXXsn[is][X-1][Y][0]     += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    pXXsn[is][X][Y-1][0]     += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    pXXsn[is][X][Y][0]       += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XY to current density field on the node*/
inline  void EMfields::addPxy(double ***weight, int X, int Y, int Z,int is){
	pXYsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	pXYsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	pXYsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	pXYsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XY to current density field on the node*/
inline  void EMfields::addPxy_OS(double ***weight, int X, int Y, int Z,int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    pXYsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    pXYsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    pXYsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    pXYsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XZ to current density field on the node */
inline  void EMfields::addPxz(double ***weight, int X, int Y,int Z, int is){
	pXZsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	pXZsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	pXZsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	pXZsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction XZ to current density field on the node */
inline  void EMfields::addPxz_OS(double ***weight, int X, int Y,int Z, int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    pXZsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    pXZsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    pXZsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    pXZsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction YY to current density field on the node*/
inline  void EMfields::addPyy(double ***weight, int X, int Y, int Z, int is){
	pYYsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
	pYYsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
	pYYsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
	pYYsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction YY to current density field on the node*/
inline  void EMfields::addPyy_OS(double ***weight, int X, int Y, int Z, int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    pYYsn[is][X-1][Y-1][0]    += weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    pYYsn[is][X-1][Y][0]      += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    pYYsn[is][X][Y-1][0]      += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    pYYsn[is][X][Y][0]        += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction YZ to current density field on the node */
inline  void EMfields::addPyz(double ***weight, int X, int Y, int Z, int is){
	pYZsn[is][X-1][Y-1][0]+= weight[0][0][0]*invVOL;
	pYZsn[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
	pYZsn[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
	pYZsn[is][X][Y][0]    += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction YZ to current density field on the node */
inline  void EMfields::addPyz_OS(double ***weight, int X, int Y, int Z, int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    pYZsn[is][X-1][Y-1][0]+= weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    pYZsn[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    pYZsn[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    pYZsn[is][X][Y][0]    += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction ZZ to current density field on the node */
inline  void EMfields::addPzz(double ***weight, int X, int Y, int Z, int is){
	pZZsn[is][X-1][Y-1][0]+= weight[0][0][0]*invVOL;
	pZZsn[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
	pZZsn[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
	pZZsn[is][X][Y][0]    += weight[1][1][0]*invVOL;
}
/** add an amount of pressure density - direction ZZ to current density field on the node */
inline  void EMfields::addPzz_OS(double ***weight, int X, int Y, int Z, int is, Grid *g){
  if (X-1>=0 && Y-1 >=0 )
    pZZsn[is][X-1][Y-1][0]+= weight[0][0][0]*invVOL;
  if (X-1>=0 && Y<= g->getNYN()-1 )
    pZZsn[is][X-1][Y][0]  += weight[0][1][0]*invVOL;
  if (X <= g->getNXN()-1 && Y-1 >=0)
    pZZsn[is][X][Y-1][0]  += weight[1][0][0]*invVOL;
  if (X <= g->getNXN()-1 && Y<= g->getNYN()-1)
    pZZsn[is][X][Y][0]    += weight[1][1][0]*invVOL;
}


/**SPECIES: Sum the charge density of different species on NODES*/
inline void EMfields::sumOverSpecies(VirtualTopology *vct){
	// calculate the correct densities on the boundaries
	//adjustNonPeriodicDensities(vct);

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
  //AMR, ME: also ghostCell interpolated  
//grid->interpN2C(rhoc,rhon);
  grid->interpN2C_alsoGC(rhoc,rhon);
  communicateCenter(nxc,nyc,rhoc,vct); // this added by stef in KU Leuven
  //printRhoc(vct);
  
}

/** Calculate hat rho hat, Jx hat, Jy hat, Jz hat */
inline void EMfields::calculateHatFunctions(Grid *grid, VirtualTopology *vct){
  eqValue (0.0,tempXC,nxc,nyc);
  eqValue (0.0,tempYC,nxc,nyc);
  eqValue (0.0,tempZC,nxc,nyc);

  // crucial point: fixing rho in the ghost cells of the coarse grid using divE= 4 pi rho
  /* if (grid->getLevel()==1)
    {
      double compX, compY;
      double invdx= 1.0/grid->getDX();
      double invdy= 1.0/grid->getDY();
      for (int i=0; i< nxc; i++)
	{
	  if (vct->getYleft_neighbor()== MPI_PROC_NULL && bcEMfaceYleft==0)
	    {
	      int j=0;
	      //compX = .5*(Ex[i+1][j][0] - Ex[i][j][0])*invdx +  .5*(Ex[i+1][j+1][0] - Ex[i][j+1][0])*invdx;
	      //compY = .5*(Ey[i][j+1][0] - Ey[i][j][0])*invdy +  .5*(Ey[i+1][j+1][0] - Ey[i+1][j][0])*invdy;
	      //rhoc[i][j][0]=(compX + compY)/ FourPI;
	      rhoc[i][j][0]=0.;
	      rhon[i][j][0]=0.;
	    }
	  if (vct->getYright_neighbor()== MPI_PROC_NULL && bcEMfaceYright==0)
	    {
	      int j=nyc-1;
              //compX = .5*(Ex[i+1][j][0] - Ex[i][j][0])*invdx +  .5*(Ex[i+1][j+1][0] - Ex[i][j+1][0])*invdx;
              //compY = .5*(Ey[i][j+1][0] - Ey[i][j][0])*invdy +  .5*(Ey[i+1][j+1][0] - Ey[i+1][j][0])*invdy;
              //rhoc[i][j][0]=(compX + compY)/ FourPI;
              rhoc[i][j][0]=0.;
	      rhon[i][j][0]=0.;
	    } 
	}  
      for (int j=0; j<nyc; j++)
	{
	  if (vct->getXleft_neighbor()== MPI_PROC_NULL &&  bcEMfaceXleft==0)
	    {
	      int i=0;
	      //compX = .5*(Ex[i+1][j][0] - Ex[i][j][0])*invdx +  .5*(Ex[i+1][j+1][0] - Ex[i][j+1][0])*invdx;
              //compY = .5*(Ey[i][j+1][0] - Ey[i][j][0])*invdy +  .5*(Ey[i+1][j+1][0] - Ey[i+1][j][0])*invdy;
              //rhoc[i][j][0]=(compX + compY)/ FourPI;
              rhoc[i][j][0]=0.;
	      rhon[i][j][0]=0.;
	    }
          if (vct->getXright_neighbor()== MPI_PROC_NULL && bcEMfaceXright==0)
	    {
	      int i=nxc-1;
	      //compX = .5*(Ex[i+1][j][0] - Ex[i][j][0])*invdx +  .5*(Ex[i+1][j+1][0] - Ex[i][j+1][0])*invdx;
              //compY = .5*(Ey[i][j+1][0] - Ey[i][j][0])*invdy +  .5*(Ey[i+1][j+1][0] - Ey[i+1][j][0])*invdy;
              //rhoc[i][j][0]=(compX + compY)/ FourPI;
              rhoc[i][j][0]=0.;
	      rhon[i][j][0]=0.;
	    }
        }  
    }*/

  //smoothing
  smooth(Nvolte,Smooth,rhoc,0,grid,vct);
  
  // calculate j hat

  for (int is=0; is < ns; is++){
    grid->divSymmTensorN2C_alsoGC(tempXC,tempYC,tempZC,pXXsn,pXYsn,pXZsn,pYYsn,pYZsn,pZZsn,is);
    scale(tempXC,-dt/2.0,nxc,nyc);
    scale(tempYC,-dt/2.0,nxc,nyc);
    scale(tempZC,-dt/2.0,nxc,nyc);
    // communicate
    communicateCenter(nxc,nyc,tempXC,vct);
    communicateCenter(nxc,nyc,tempYC,vct);
    communicateCenter(nxc,nyc,tempZC,vct);
    grid->interpC2N_BC_alsoGN(tempXN,tempXC,vct);
    grid->interpC2N_BC_alsoGN(tempYN,tempYC,vct);
    grid->interpC2N_BC_alsoGN(tempZN,tempZC,vct);
    
    sum(tempXN,Jxs,nxn,nyn,is);
    sum(tempYN,Jys,nxn,nyn,is);
    sum(tempZN,Jzs,nxn,nyn,is);
    // PIDOT
    PIdot_alsoGN(Jxh,Jyh,Jzh,tempXN,tempYN,tempZN,is,grid);}
  
    //Here BC should be received for Jxyzh on ghost nodes of the fine grid so that in the following, ghost centers of rhoh
    // can be evaluated.
  
  communicateNode(nxn, nyn, Jxh, vct);
  communicateNode(nxn, nyn, Jyh, vct);
  communicateNode(nxn, nyn, Jzh, vct);

  smooth(Nvolte,Smooth,Jxh,1,grid,vct);
  smooth(Nvolte,Smooth,Jyh,1,grid,vct);
  smooth(Nvolte,Smooth,Jzh,1,grid,vct);
  
  // calculate rhoh
  grid->divN2C_plusghost(tempXC,Jxh,Jyh);
  scale(tempXC,-dt*th,nxc,nyc);
  sum(tempXC,rhoc,nxc,nyc);
  eq(rhoh,tempXC,nxc,nyc);
  // at this point is communicated

  if  (grid->getLevel()==1)                                                                             
    {
      int END;
      //if (grid->getLevel()==0)                                                                                 
      END=4;
      for (int i=0; i<nxc; i++)
        {
          if (vct->getYleft_neighbor()== MPI_PROC_NULL) // && bcEMfaceYleft ==0)// perfect conductor           
            {
              //int j=0;                                                                                       
              for (int j=0; j<END; j++)
                rhoh[i][j][0]=0.0;
            }
          if (vct->getYright_neighbor()== MPI_PROC_NULL)// && bcEMfaceYright ==0)                              
            {
              //int j=nyc-1;                                                                                   
              for (int j=nyc-1; j>nyc-1-END; j--)
                rhoh[i][j][0]=0.0;
            }
        }
      for (int j=0; j<nyc; j++)
        {
          if (vct->getXleft_neighbor()== MPI_PROC_NULL)//  && bcEMfaceXleft ==0)                               
            {
              //          int i=0;                                                                             
              for (int i=0; i<END; i++)
                rhoh[i][j][0]=0.0;
            }
          if (vct->getXright_neighbor()== MPI_PROC_NULL)//  && bcEMfaceXright ==0)                             
            {
              //          int i=nxc-1;                                                                         
              for (int i=nxc-1; i>nxc-1-END; i--)
                rhoh[i][j][0]=0.0;
            }
        }

	}


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

// calculate PI dot also on ghost nodes
inline void EMfields::PIdot_alsoGN(double ***PIdotX, double ***PIdotY, double ***PIdotZ, double ***vectX, double ***vectY, double ***vectZ, int ns, Grid *grid){
  double beta, edotb, omcx, omcy, omcz,denom;
  beta = .5*qom[ns]*dt/c;
  for(int i=0; i <nxn;i++)
    for(int j=0; j <nyn;j++){
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


/** Calculate MU dot (vectX, vectY, vectZ) including ghost cells*/
inline void EMfields::MUdot_plusghost(double ***MUdotX, double ***MUdotY, double ***MUdotZ, double ***vectX, double ***vectY, double ***vectZ, Grid *grid){
	double beta, edotb, omcx, omcy, omcz, denom;
	for(int i=0; i < nxn;i++)
		for(int j=0; j < nyn;j++){
			MUdotX[i][j][0] = 0.0;
			MUdotY[i][j][0] = 0.0;
			MUdotZ[i][j][0] = 0.0;
		}
			for (int is=0; is < ns; is++){
				beta = .5*qom[is]*dt/c;
				for(int i=0; i < nxn;i++)
					for(int j=0; j < nyn;j++){
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
			if (vct->getXleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int j=1; j<ny-1;j++)
					temp[1][j]=value*vector[1][j][0]+alpha*(vector[1][j-1][0]+vector[2][j][0]+vector[1][j+1][0]);}

			if (vct->getXright_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int j=1; j<ny-1;j++)
					temp[nx-2][j]=value*vector[nx-2][j][0]+alpha*(vector[nx-2][j-1][0]+vector[nx-3][j][0]+vector[nx-2][j+1][0]);}

			if (vct->getYleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int i=1; i<nx-1;i++)
					temp[i][1]=value*vector[i][1][0]+alpha*(vector[i-1][1][0]+vector[i+1][1][0]+vector[i][2][0]);}

			if (vct->getYright_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int i=1; i<nx-1;i++)
					temp[i][ny-2]=value*vector[i][ny-2][0]+alpha*(vector[i-1][ny-2][0]+vector[i][ny-3][0]+vector[i+1][ny-2][0]);}

			if (vct->getXleft_neighbor()==MPI_PROC_NULL && vct->getYleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/2;
				temp[1][1]=value*vector[1][1][0]+alpha*(vector[2][1][0]+vector[1][2][0]);
			}
			if (vct->getXleft_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/2;
				temp[1][ny-2]=value*vector[1][ny-2][0]+alpha*(vector[2][ny-1][0]+vector[1][ny-3][0]);
			}
			if (vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/2;
				temp[nx-2][1]=value*vector[nx-2][1][0]+alpha*(vector[nx-3][1][0]+vector[nx-2][2][0]);
			}
			if (vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL){
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
			if (vct->getXleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int j=1; j<ny-1;j++)
					temp[1][j]=value*vector[is][1][j][0]+alpha*(vector[is][1][j-1][0]+vector[is][2][j][0]+vector[is][1][j+1][0]);}

			if (vct->getXright_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int j=1; j<ny-1;j++)
					temp[nx-2][j]=value*vector[is][nx-2][j][0]+alpha*(vector[is][nx-2][j-1][0]+vector[is][nx-3][j][0]+vector[is][nx-2][j+1][0]);}

			if (vct->getYleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int i=1; i<nx-1;i++)
					temp[i][1]=value*vector[is][i][1][0]+alpha*(vector[is][i-1][1][0]+vector[is][i+1][1][0]+vector[is][i][2][0]);}

			if (vct->getYright_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/3;
				for (int i=1; i<nx-1;i++)
					temp[i][ny-2]=value*vector[is][i][ny-2][0]+alpha*(vector[is][i-1][ny-2][0]+vector[is][i][ny-3][0]+vector[is][i+1][ny-2][0]);}
			// corners
			if (vct->getXleftYleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/2;
				temp[1][1]=value*vector[is][1][1][0]+alpha*(vector[is][2][1][0]+vector[is][1][2][0]);
			}
			if (vct->getXleftYright_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/2;
				temp[1][ny-2]=value*vector[is][1][ny-2][0]+alpha*(vector[is][2][ny-1][0]+vector[is][1][ny-3][0]);
			}
			if (vct->getXrightYleft_neighbor()==MPI_PROC_NULL){
				alpha=(1.0-value)/2;
				temp[nx-2][1]=value*vector[is][nx-2][1][0]+alpha*(vector[is][nx-3][1][0]+vector[is][nx-2][2][0]);
			}
			if (vct->getXrightYright_neighbor()==MPI_PROC_NULL){
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
			for (int i=1; i <  nyn-1;i++){
				Bxn[0][i][0] = 0.0;
				Byn[0][i][0] = (Byc[0][i-1][0]+Byc[0][i][0])/2.-Jz[0][i][0]*dx/2.;
				Bzn[0][i][0] = (Bzc[0][i-1][0]+Bzc[0][i][0])/2.+Jy[0][i][0]*dx/2.;
			}
				Bxn[0][0][0] = 0.0;
				Byn[0][0][0] = Byc[0][0][0];//-Jz[0][0][0]*dx/2.;
				Bzn[0][0][0] = Bzc[0][0][0];//+Jy[0][0][0]*dx/2.;
				Bxn[0][nyn-1][0] = 0.0;
				Byn[0][nyn-1][0] = Byc[0][nyc-1][0];//-Jz[0][nyn-1][0]*dx/2.;
				Bzn[0][nyn-1][0] = Bzc[0][nyc-1][0];//+Jy[0][nyn-1][0]*dx/2.;
 
			break;
		case 1: // boundary condition on Y-DIRECTION left
                  	for (int i=1; i < nxn-1;i++){
				Bxn[i][0][0] = (Bxc[i][0][0]+Bxc[i-1][0][0])/2.+Jz[i][0][0]*dy/2.;
				Byn[i][0][0] = 0.0;         
				Bzn[i][0][0] = (Bzc[i][0][0]+Bzc[i-1][0][0])/2.-Jx[i][0][0]*dy/2.;
			}
				Byn[0][0][0] = 0.0;
				Bxn[0][0][0] = Bxc[0][0][0];
				Bzn[0][0][0] = Bzc[0][0][0];
				Byn[nxn-1][0][0] = 0.0;
				Bxn[nxn-1][0][0] = Bxc[nxc-1][0][0];
				Bzn[nxn-1][0][0] = Bzc[nxc-1][0][0];
                        // Stefano's complicated version for the GEM challenge
		/*	for (int i=0; i < nxc;i++){
				Bxc[i][0][0] = B0x*tanh((grid->getYC(i,0,0) - Ly/2)/delta);
				//Bxc[i][1][0] = B0x*tanh((grid->getYC(i,0,0) - Ly/2)/delta);
				//Bxc[i][2][0] = B0x*tanh((grid->getYC(i,0,0) - Ly/2)/delta);
				Byc[i][0][0] = B0y;
				Bzc[i][0][0] = B0z;///cosh((grid->getYC(i,0,0) - Ly/2)/delta);
				Bzc[i][1][0] = Bzc[i][0][0];
				Bzc[i][2][0] = Bzc[i][0][0];
			}*/
			break;
	}
}
/** Perfect conductor boundary conditions for magnetic field  */
inline  void EMfields::BperfectConductorRight(int dir,Grid *grid,VirtualTopology *vct){
	switch(dir){
		case 0: // boundary condition on X-DIRECTION right
			for (int i=1; i <  nyn-1;i++){                        
				Bxn[nxn-1][i][0] = 0.0;                     
				Byn[nxn-1][i][0] = (Byc[nxn-2][i][0]+Byc[nxn-2][i-1][0])/2.+Jz[nxn-1][i][0]*dx/2.;
				Bzn[nxn-1][i][0] = (Bzc[nxn-2][i][0]+Bzc[nxn-2][i-1][0])/2.-Jy[nxn-1][i][0]*dx/2.;
			}
				Bxn[nxn-1][0][0] = 0.0;                     
				Byn[nxn-1][0][0] = Byc[nxn-2][0][0];
				Bzn[nxn-1][0][0] = Bzc[nxn-2][0][0];
				Bxn[nxn-1][nyn-1][0] = 0.0;                     
				Byn[nxn-1][nyn-1][0] = Byc[nxn-2][nyn-2][0];
				Bzn[nxn-1][nyn-1][0] = Bzc[nxn-2][nyn-2][0];

                        break;
		case 1: // boundary condition on Y-DIRECTION right
                       for (int i=1; i < nxn-1;i++){                          
				Bxn[i][nyn-1][0] = (Bxc[i][nyn-2][0]+Bxc[i-1][nyn-2][0])/2.-Jz[i][nyn-1][0]*dy/2.;
				Byn[i][nyn-1][0] = 0.0;             
				Bzn[i][nyn-1][0] = (Bzc[i][nyn-2][0]+Bzc[i-1][nyn-2][0])/2.+Jx[i][nyn-1][0]*dy/2.;
			}
				Bxn[0][nyn-1][0] = Bxc[0][nyn-2][0];
				Byn[0][nyn-1][0] = 0.0;             
				Bzn[0][nyn-1][0] = Bzc[0][nyn-2][0];
				Bxn[nxn-1][nyn-1][0] = Bxc[nxn-2][nyn-2][0];
				Byn[nxn-1][nyn-1][0] = 0.0;            
				Bzn[nxn-1][nyn-1][0] = Bzc[nxn-2][nyn-2][0];
 
                        // Stefano's complicated version for the GEM challenge
			/*for (int i=0; i < nxc;i++){
				Bxc[i][nyc-1][0] = B0x*tanh((grid->getYC(i,nyc-1,0) - Ly/2)/delta);
				//Bxc[i][nyc-2][0] = B0x*tanh((grid->getYC(i,nyc-1,0) - Ly/2)/delta);
				//Bxc[i][nyc-3][0] = B0x*tanh((grid->getYC(i,nyc-1,0) - Ly/2)/delta);
				Byc[i][nyc-1][0] = B0y;
				Bzc[i][nyc-1][0] = B0z;///cosh((grid->getYC(i,0,0) - Ly/2)/delta);;
				Bzc[i][nyc-2][0] = Bzc[i][nyc-1][0];
				Bzc[i][nyc-3][0] = Bzc[i][nyc-1][0];
			}*/
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
                                   //Stefano's version below
				   vectorX[0][i][0] =  vectorX[1][i][0];//-rhon[0][i][0]*dx;// - rhon[1][i][0]*dy;
				   vectorY[0][i][0] =  0.0;
				   vectorZ[0][i][0] =  0.0;
			   }

			   break;
		   case 1: // boundary condition on Y-DIRECTION left
		           //stefano's version below
                   	   /*susxx = new double[nxn];
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
	                   delete[] susxx;
	                   delete[] susxy;
	                   delete[] susxz;
	                   delete[] susyx;
	                   delete[] susyy;
	                   delete[] susyz;
	                   delete[] suszx;
	                   delete[] suszy;
	                   delete[] suszz;*/
			   for (int i=0; i < nxn;i++){
			       vectorY[i][0][0] = vectorY[i][1][0];//-rhon[i][0][0]*dy;// - rhon[i][1][0]*dx;
			       vectorX[i][0][0] =  0.0;
			       vectorZ[i][0][0] =  0.0;
                           }

				   break;
	   }

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
                                   //Stefano's version below
				   vectorX[nxn-1][i][0] = vectorX[nxn-2][i][0];//+rhon[nxn-1][i][0]*dx;
				   vectorY[nxn-1][i][0] = 0.0;
				   vectorZ[nxn-1][i][0] = 0.0;
			   }
			   break;
		   case 1: // boundary condition on Y-DIRECTION RIGHT
                           //stefano's version below
			   /*susxx = new double[nxn];
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
	                   delete[] susxx;
	                   delete[] susxy;
	                   delete[] susxz;
	                   delete[] susyx;
	                   delete[] susyy;
	                   delete[] susyz;
	                   delete[] suszx;
	                   delete[] suszy;
	                   delete[] suszz;*/
			   for (int i=0; i < nxn;i++){
				   vectorX[i][nyn-1][0] = 0.0;
				   vectorZ[i][nyn-1][0] = 0.0;
				   vectorY[i][nyn-1][0] = vectorY[i][nyn-2][0];// +rhon[i][nyn-1][0]*dy;

			   }

				   break;
	   }
}

/** Magnetic mirror boundary conditions: LEFT WALL */
inline  void EMfields::magneticMirrorLeft(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){
	   switch(dir){
		   case 0: // boundary condition on X-DIRECTION left
			   for (int i=0; i <  nyn;i++){
                                   vectorX[0][i][0] = 0.0;
				   vectorY[0][i][0] = vectorY[1][i][0];
				   vectorZ[0][i][0] = vectorY[1][i][0];
			   }
			   break;
		   case 1: // boundary condition on Y-DIRECTION left
                            for (int i=0; i < nxn;i++){
                               vectorY[i][0][0] = 0.0;
			       vectorX[i][0][0] = vectorX[i][1][0];
			       vectorZ[i][0][0] = vectorZ[i][1][0];
                           }
                           break;
	   }

}
/** Magnetic mirror boundary conditions: RIGHT WALL */
inline  void EMfields::magneticMirrorRight(double ***vectorX, double ***vectorY, double ***vectorZ,int dir,Grid *grid){
	   switch(dir){
		   case 0: // boundary condition on X-DIRECTION left
			   for (int i=0; i <  nyn;i++){
                                   vectorX[nxn-1][i][0] = 0.0;
				   vectorY[nxn-1][i][0] = vectorY[nxn-2][i][0];
				   vectorZ[nxn-1][i][0] = vectorZ[nxn-2][i][0];
			   }
			   break;
		   case 1: // boundary condition on Y-DIRECTION left
                            for (int i=0; i < nxn;i++){
                               vectorY[i][nyn-1][0] = 0.0;
			       vectorX[i][nyn-1][0] = vectorX[i][nyn-2][0];
			       vectorZ[i][nyn-1][0] = vectorZ[i][nyn-2][0];
                           }
			   break;
	   }

}
inline  void EMfields::BmagneticMirrorLeft(int dir, Grid *grid, VirtualTopology *vct){

	switch(dir){
		case 0: // boundary condition on X-DIRECTION left
			for (int i=1; i <  nyn-1;i++){
				Bxn[0][i][0] = (Bxc[0][i][0]+Bxc[0][i-1][0])/2.;
				Byn[0][i][0] = 0.0;
				Bzn[0][i][0] = 0.0;
			}
			Bxn[0][0][0] = Bxc[0][0][0];
			Byn[0][0][0] = 0.0;
			Bzn[0][0][0] = 0.0;
			Bxn[0][nyn-1][0] = Bxc[0][nyc-1][0];
			Byn[0][nyn-1][0] = 0.0;
			Bzn[0][nyn-1][0] = 0.0;
                         
			break;
                        
		case 1: // boundary condition on Y-DIRECTION left
                  	for (int i=1; i < nxn-1;i++){
				Bxn[i][0][0] = 0.0;
				Byn[i][0][0] = (Byc[i][0][0]+Byc[i-1][0][0])/2.;
				Bzn[i][0][0] = 0.0;
			}
                        break;
				Bxn[0][0][0] = 0.0;
				Byn[0][0][0] = Byc[0][0][0];
				Bzn[0][0][0] = 0.0;
				Bxn[nxn-1][0][0] = 0.0;
				Byn[nxn-1][0][0] = Byc[nxc-1][0][0];
				Bzn[nxn-1][0][0] = 0.0;
        }
}
inline  void EMfields::BmagneticMirrorRight(int dir, Grid *grid, VirtualTopology *vct){

	switch(dir){
		case 0: // boundary condition on X-DIRECTION right
			for (int i=1; i <  nyn-1;i++){
				Bxn[nxn-1][i][0] = (Bxc[nxn-2][i][0]+Bxc[nxn-2][i-1][0])/2.;
				Byn[nxn-1][i][0] = 0.0;
				Bzn[nxn-1][i][0] = 0.0;
			}
			Bxn[nxn-1][0][0] = Bxc[nxn-2][0][0];
			Byn[nxn-1][0][0] = 0.0;
			Bzn[nxn-1][0][0] = 0.0;
			Bxn[nxn-1][nyn-1][0] = Bxc[nxn-2][nyn-2][0];
			Byn[nxn-1][nyn-1][0] = 0.0;
			Bzn[nxn-1][nyn-1][0] = 0.0;
			break;
		case 1: // boundary condition on Y-DIRECTION right
                  	for (int i=1; i < nxn-1;i++){
				Bxn[i][nyn-1][0] = 0.0;
				Byn[i][nyn-1][0] = (Byc[i][nyn-2][0]+Byc[i-1][nyn-2][0])/2.;
				Bzn[i][nyn-1][0] = 0.0;
			}
			Bxn[0][nyn-1][0] = 0.0;
			Byn[0][nyn-1][0] = Byc[0][nyn-2][0];
			Bzn[0][nyn-1][0] = 0.0;
			Bxn[nxn-1][nyn-1][0] = 0.0;
			Byn[nxn-1][nyn-1][0] = Byc[nxn-2][nyn-2][0];
			Bzn[nxn-1][nyn-1][0] = 0.0;
                        break;
        }
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
/** Continuous BC for PHI  left side*/
inline void EMfields::bcPHI_Left_continuous(double ***im,int dir){
	switch(dir){
		case 0:
			for (int i=0; i <  nyc;i++){
				im[0][i][0] = im[1][i][0];

			}
			break;
		case 1:
			for (int i=1; i < nxc;i++){
				im[i][0][0] = im[i][1][0];

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
/** Continuous BC for PHI  right side*/
inline void EMfields::bcPHI_Right_continuous(double ***im,int dir){
	switch(dir){
		case 0:
			for (int i=0; i <  nyc;i++){
				im[nxc-1][i][0] = im[nxc-2][i][0];

			}
			break;
		case 1:
			for (int i=0; i < nxc;i++){
				im[i][nyc-1][0] = im[i][nyc-2][0];
			}
			break;
	}
}

/** Send the refined fields to coarser level*/
inline void EMfields::sendProjection(Grid *grid, VirtualTopology *vct){
int i,j,ierr;
MPI_Request request;
int start,step;

//cout << "SENDING projection"<<endl;

//buffer - / -
if (nmessageProj > 0) {
	eqValue (0.0, reducedEx,  nxnproj,nynproj);
	eqValue (0.0, reducedEy,  nxnproj,nynproj);
	eqValue (0.0, reducedEz,  nxnproj,nynproj);
	eqValue (0.0, reducedBxn, nxnproj,nynproj);
	eqValue (0.0, reducedByn, nxnproj,nynproj);
	eqValue (0.0, reducedBzn, nxnproj,nynproj);
    for (i=0;i<lastindicex;i++){
        for (j=0;j<lastindicey;j++){
            reducedEx[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0]  *Ex[i+1][j+1][0];
            reducedEy[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0]  *Ey[i+1][j+1][0];
            reducedEz[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0]  *Ez[i+1][j+1][0];
            reducedBxn[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0]*Bxn[i+1][j+1][0];
            reducedByn[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0]*Byn[i+1][j+1][0];
            reducedBzn[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0]*Bzn[i+1][j+1][0];
            reducedEx[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0]    *Ex[i+1][j+1][0];
            reducedEy[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0]    *Ey[i+1][j+1][0];
            reducedEz[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0]    *Ez[i+1][j+1][0];
            reducedBxn[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0]  *Bxn[i+1][j+1][0];
            reducedByn[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0]  *Byn[i+1][j+1][0];
            reducedBzn[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0]  *Bzn[i+1][j+1][0];
            reducedEx[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0]    *Ex[i+1][j+1][0];
            reducedEy[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0]    *Ey[i+1][j+1][0];
            reducedEz[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0]    *Ez[i+1][j+1][0];
            reducedBxn[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0]  *Bxn[i+1][j+1][0];
            reducedByn[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0]  *Byn[i+1][j+1][0];
            reducedBzn[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0]  *Bzn[i+1][j+1][0];
            reducedEx[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0]      *Ex[i+1][j+1][0];
            reducedEy[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0]      *Ey[i+1][j+1][0];
            reducedEz[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0]      *Ez[i+1][j+1][0];
            reducedBxn[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0]    *Bxn[i+1][j+1][0];
            reducedByn[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0]    *Byn[i+1][j+1][0];
            reducedBzn[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0]    *Bzn[i+1][j+1][0];
        }
    }
    /*if(vct->getCartesian_rank_COMMTOTAL()==25){
        cout << "Reduced Ex "<<endl;
        for (i=0;i<nxnproj;i++){
            for (j=0;j<nynproj;j++){
                cout << reducedEx[i][j][0]<< " ";
            }
            cout << endl;
        }
    }*/

    //Load reduced quantities in the buffer
    start=0;
       /*if(vct->getCartesian_rank_COMMTOTAL()==25){
               cout << "Proc 25 bufferProj send "<< nxmsend*nymsend<< " points ";
           }*/
    for (i=0;i<nxmsend*nymsend;i++){
        bufferProj[start] = reducedEx[i/nymsend][i%nymsend][0];
        bufferProj[start+1] = reducedEy[i/nymsend][i%nymsend][0];
        bufferProj[start+2] = reducedEz[i/nymsend][i%nymsend][0];
        bufferProj[start+3] = reducedBxn[i/nymsend][i%nymsend][0];
        bufferProj[start+4] = reducedByn[i/nymsend][i%nymsend][0];
        bufferProj[start+5] = reducedBzn[i/nymsend][i%nymsend][0];
        start += 6;
       }
       /*if(vct->getCartesian_rank_COMMTOTAL()==25){
               cout << endl;
           }*/
    
    
    
    //ierr = MPI_Isend(bufferProj,nxmsend*nymsend*6,MPI_DOUBLE,targetProj[0],1,MPI_COMM_WORLD,&request);
    ierr = MPI_Ssend(bufferProj,nxmsend*nymsend*6,MPI_DOUBLE,targetProj[0],1,MPI_COMM_WORLD);
    //cout << "Proc "<< vct->getCartesian_rank_COMMTOTAL() << "sends "<< nxmsend*nymsend<< " points to " << targetProj[0] <<endl;
}
//buffer + / -
if (nmessageProj > 1 && nxpsend > 0){
    step = start;
    //Load reduced quantities in the buffer
    for (i=0;i<nxpsend*nymsend;i++){
        bufferProj[start] = reducedEx[nxmsend+i/nymsend][i%nymsend][0];
        bufferProj[start+1] = reducedEy[nxmsend+i/nymsend][i%nymsend][0];
        bufferProj[start+2] = reducedEz[nxmsend+i/nymsend][i%nymsend][0];
        bufferProj[start+3] = reducedBxn[nxmsend+i/nymsend][i%nymsend][0];
        bufferProj[start+4] = reducedByn[nxmsend+i/nymsend][i%nymsend][0];
        bufferProj[start+5] = reducedBzn[nxmsend+i/nymsend][i%nymsend][0];
        start += 6;
    }
    
    //ierr = MPI_Isend(bufferProj+start,nxpsend*nymsend*6,MPI_DOUBLE,targetProj[1],1,MPI_COMM_WORLD,&request);
    ierr = MPI_Ssend(bufferProj+step,nxpsend*nymsend*6,MPI_DOUBLE,targetProj[1],1,MPI_COMM_WORLD);
    //cout << "Proc "<< vct->getCartesian_rank_COMMTOTAL() << "sends "<< nxpsend*nymsend<< " points to " << targetProj[1] <<endl;
}
//buffer - / +
if (nmessageProj > 1 && nypsend > 0){
    step = start;
    //Load reduced quantities in the buffer
    for (i=0;i<nxmsend*nypsend;i++){
        bufferProj[start] = reducedEx[i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+1] = reducedEy[i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+2] = reducedEz[i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+3] = reducedBxn[i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+4] = reducedByn[i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+5] = reducedBzn[i/nypsend][nymsend+i%nypsend][0];
        start += 6;
    }
 
    //ierr = MPI_Isend(bufferProj+start,nxmsend*nypsend*6,MPI_DOUBLE,targetProj[2],1,MPI_COMM_WORLD,&request);
    ierr = MPI_Ssend(bufferProj+step,nxmsend*nypsend*6,MPI_DOUBLE,targetProj[2],1,MPI_COMM_WORLD);
    //cout << "Proc "<< vct->getCartesian_rank_COMMTOTAL() << "sends "<< nxmsend*nypsend<< " points to " << targetProj[2] <<endl;
}
//buffer + / +
if (nmessageProj > 2){
    step = start;
    //Load reduced quantities in the buffer
    for (i=0;i<nxpsend*nypsend;i++){
        bufferProj[start] = reducedEx[nxmsend+i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+1] = reducedEy[nxmsend+i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+2] = reducedEz[nxmsend+i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+3] = reducedBxn[nxmsend+i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+4] = reducedByn[nxmsend+i/nypsend][nymsend+i%nypsend][0];
        bufferProj[start+5] = reducedBzn[nxmsend+i/nypsend][nymsend+i%nypsend][0];
        start += 6;
    }
    
    //ierr = MPI_Isend(bufferProj+start,nxpsend*nypsend*6,MPI_DOUBLE,targetProj[3],1,MPI_COMM_WORLD,&request);
    ierr = MPI_Ssend(bufferProj+step,nxpsend*nypsend*6,MPI_DOUBLE,targetProj[3],1,MPI_COMM_WORLD);
    //cout << "Proc "<< vct->getCartesian_rank_COMMTOTAL() << "sends "<< nxpsend*nypsend<< " points to " << targetProj[3] <<endl;
}
}

/** Send the boundary conditions to finer level*/
inline void EMfields::sendBC(Grid *grid, VirtualTopology *vct){

int i,j,ierr,start,startc;
MPI_Status status;
MPI_Request request;
// Send data on nodes
start = 0;
for (i=0;i<nmessageBC;i++) {
    for (j=0;j<npointssent[i];j++){
        bufferBC[6*(j+start)  ] = weightBC[j+start][0][0]* Ex[ixsent[j+start]-1][iysent[j+start]-1][0] + weightBC[j+start][1][0]* Ex[ixsent[j+start]-1][iysent[j+start]][0] + weightBC[j+start][2][0]* Ex[ixsent[j+start]][iysent[j+start]-1][0] + weightBC[j+start][3][0]*Ex[ixsent[j+start]][iysent[j+start]][0];
        bufferBC[6*(j+start)+1] = weightBC[j+start][0][0]* Ey[ixsent[j+start]-1][iysent[j+start]-1][0] + weightBC[j+start][1][0]* Ey[ixsent[j+start]-1][iysent[j+start]][0] + weightBC[j+start][2][0]* Ey[ixsent[j+start]][iysent[j+start]-1][0] + weightBC[j+start][3][0]*Ey[ixsent[j+start]][iysent[j+start]][0];
        bufferBC[6*(j+start)+2] = weightBC[j+start][0][0]* Ez[ixsent[j+start]-1][iysent[j+start]-1][0] + weightBC[j+start][1][0]* Ez[ixsent[j+start]-1][iysent[j+start]][0] + weightBC[j+start][2][0]* Ez[ixsent[j+start]][iysent[j+start]-1][0] + weightBC[j+start][3][0]*Ez[ixsent[j+start]][iysent[j+start]][0];
        bufferBC[6*(j+start)+3] = weightBC[j+start][0][0]*Bxn[ixsent[j+start]-1][iysent[j+start]-1][0] + weightBC[j+start][1][0]*Bxn[ixsent[j+start]-1][iysent[j+start]][0] + weightBC[j+start][2][0]*Bxn[ixsent[j+start]][iysent[j+start]-1][0] + weightBC[j+start][3][0]*Bxn[ixsent[j+start]][iysent[j+start]][0];
        bufferBC[6*(j+start)+4] = weightBC[j+start][0][0]*Byn[ixsent[j+start]-1][iysent[j+start]-1][0] + weightBC[j+start][1][0]*Byn[ixsent[j+start]-1][iysent[j+start]][0] + weightBC[j+start][2][0]*Byn[ixsent[j+start]][iysent[j+start]-1][0] + weightBC[j+start][3][0]*Byn[ixsent[j+start]][iysent[j+start]][0];
        bufferBC[6*(j+start)+5] = weightBC[j+start][0][0]*Bzn[ixsent[j+start]-1][iysent[j+start]-1][0] + weightBC[j+start][1][0]*Bzn[ixsent[j+start]-1][iysent[j+start]][0] + weightBC[j+start][2][0]*Bzn[ixsent[j+start]][iysent[j+start]-1][0] + weightBC[j+start][3][0]*Bzn[ixsent[j+start]][iysent[j+start]][0];
    }        
    // ierr = MPI_Isend(bufferBC+6*start,npointssent[i]*6,MPI_DOUBLE,targetBC[i],1,MPI_COMM_WORLD,&request);
    ierr = MPI_Ssend(bufferBC+6*start,npointssent[i]*6,MPI_DOUBLE,targetBC[i],1,MPI_COMM_WORLD);
    start += npointssent[i];
}
}
/** Receive the refined fields from finer level*/
inline void EMfields::receiveProjection(CollectiveIO *col, Grid *grid, VirtualTopology *vct){
int ierr;
MPI_Status status;
MPI_Request request;
int i,j,k,ix,iy,step;
double Ox,Oy,finelx,finely;

//int ID= vct->getCartesian_rank_COMMTOTAL();   // just for DDT

Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid
Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid
finelx = col->getLx()/pow(col->getRatio(),grid->getLevel()+1);
finely = col->getLy()/pow(col->getRatio(),grid->getLevel()+1);
eqValue (0.0, Ex_recvbufferproj,nxn,nyn);
eqValue (0.0, Ey_recvbufferproj,  nxn,nyn);
eqValue (0.0, Ez_recvbufferproj,  nxn,nyn);
eqValue (0.0, Bxn_recvbufferproj, nxn,nyn);
eqValue (0.0, Byn_recvbufferproj, nxn,nyn);
eqValue (0.0, Bzn_recvbufferproj, nxn,nyn);

    for (i=0;i<nmessagerecvProj;i++) {
    //cout << "Proc "<< vct->getCartesian_rank_COMMTOTAL() << "receives "<< npointsreceivedProj[i]<< " points from " << fromProj[i] <<"init position = "<< ixrecvfirstProj[i]<<" "<<iyrecvfirstProj[i]<< endl;
    ierr = MPI_Recv(bufferProj,npointsreceivedProj[i]*6,MPI_DOUBLE,fromProj[i],1,MPI_COMM_WORLD, &status);
       for (j=0;j<npointsreceivedProj[i];j++){
           ix = ixrecvfirstProj[i] + j/nyrecvProj[i];
           iy = iyrecvfirstProj[i] + j%nyrecvProj[i];
           Ex_recvbufferproj[ix][iy][0] += bufferProj[6*j];
           Ey_recvbufferproj[ix][iy][0] += bufferProj[6*j+1];
           Ez_recvbufferproj[ix][iy][0] += bufferProj[6*j+2];
           Bxn_recvbufferproj[ix][iy][0] += bufferProj[6*j+3];
           Byn_recvbufferproj[ix][iy][0] += bufferProj[6*j+4];
           Bzn_recvbufferproj[ix][iy][0] += bufferProj[6*j+5];
       }
       /*if(vct->getCartesian_rank_COMMTOTAL()==10){
           cout << "proc 10 receives message from "<<fromProj[i]<<endl;
           for (j=0;j<npointsreceivedProj[i];j++){
               ix = ixrecvfirstProj[i] + j/nyrecvProj[i];
               iy = iyrecvfirstProj[i] + j%nyrecvProj[i];
               cout << Ex_recvbufferproj[ix][iy][0]<< " ";
           }
       }*/
    }
    
   // Completing recvprojbuffers with the values gathered by the other coarse procs
    if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //Send a message to proc  coorinateX-1
        step = 0;
        for (i=1;i<nyn-1;i++){
                bufferProjsend[step  ] =  Ex_recvbufferproj[1][i][0];
                bufferProjsend[step+1] =  Ey_recvbufferproj[1][i][0];
                bufferProjsend[step+2] =  Ez_recvbufferproj[1][i][0];
                bufferProjsend[step+3] = Bxn_recvbufferproj[1][i][0];
                bufferProjsend[step+4] = Byn_recvbufferproj[1][i][0];
                bufferProjsend[step+5] = Bzn_recvbufferproj[1][i][0];
                step +=6;
            }
        if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
        //also receive a message
            ierr = MPI_Sendrecv(bufferProjsend,6*(nyn-2),MPI_DOUBLE,vct->getXleft_neighbor(),1,bufferProjrecv,6*(nyn-2),MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
            ierr = MPI_Send(bufferProjsend,6*(nyn-2),MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM());
        }
    } else {
        //Do not send ...
          if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
          //... but receive only
              ierr = MPI_Recv(bufferProjrecv,6*(nyn-2),MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(), &status);
          }
    }
    if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
        //If a message from proc coorinateX+1 has been received
        step = 0;
        for (i=1;i<nyn-1;i++){
                Ex_recvbufferproj[nxn-2][i][0] += bufferProjrecv[step];
                Ey_recvbufferproj[nxn-2][i][0] += bufferProjrecv[step+1];
                Ez_recvbufferproj[nxn-2][i][0] += bufferProjrecv[step+2];
                Bxn_recvbufferproj[nxn-2][i][0] += bufferProjrecv[step+3];
                Byn_recvbufferproj[nxn-2][i][0] += bufferProjrecv[step+4];
                Bzn_recvbufferproj[nxn-2][i][0] += bufferProjrecv[step+5];
                step +=6;
        }
    //Send a message to proc coordinateX+1
        step =0;
        for(i=1;i<nyn-1;i++){
            bufferProjsend[step] = Ex_recvbufferproj[nxn-2][i][0];
            bufferProjsend[step+1] = Ey_recvbufferproj[nxn-2][i][0];
            bufferProjsend[step+2] = Ez_recvbufferproj[nxn-2][i][0];
            bufferProjsend[step+3] = Bxn_recvbufferproj[nxn-2][i][0];
            bufferProjsend[step+4] = Byn_recvbufferproj[nxn-2][i][0];
            bufferProjsend[step+5] = Bzn_recvbufferproj[nxn-2][i][0];
            step +=6;
        }
        if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //also receive a message
        ierr = MPI_Sendrecv(bufferProjsend,6*(nyn-2),MPI_DOUBLE,vct->getXright_neighbor(),1,bufferProjrecv,6*(nyn-2),MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
        //cout << vct->getCartesian_rank_COMMTOTAL() <<" send to "<< vct->getXright_neighbor() << endl;
        ierr = MPI_Send(bufferProjsend,6*(nyn-2),MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM());
        }
    
    } else {
          //Do not send ...
           if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
           //...receive only
           //cout << vct->getCartesian_rank_COMMTOTAL() <<" receive from "<< vct->getXleft_neighbor() << endl;
           ierr = MPI_Recv(bufferProjrecv,6*(nyn-2),MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(), &status);
           }
    }
    if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //A message from proc coordinateX-1 has been received
        step=0;
        for (i=1;i<nyn-1;i++){
            Ex_recvbufferproj[1][i][0] =  bufferProjrecv[step  ];
            Ey_recvbufferproj[1][i][0] =  bufferProjrecv[step+1];
            Ez_recvbufferproj[1][i][0] =  bufferProjrecv[step+2];
            Bxn_recvbufferproj[1][i][0] =  bufferProjrecv[step+3];
            Byn_recvbufferproj[1][i][0] =  bufferProjrecv[step+4];
            Bzn_recvbufferproj[1][i][0] =  bufferProjrecv[step+5];
            step +=6;
        }
    }

    // In the Y direction now
    if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //Send a message to proc  coorinateY-1
        step =0;
        for (i=1;i<nxn-1;i++){
                bufferProjsend[step  ] = Ex_recvbufferproj[i][1][0];
                bufferProjsend[step+1] = Ey_recvbufferproj[i][1][0];
                bufferProjsend[step+2] = Ez_recvbufferproj[i][1][0];
                bufferProjsend[step+3] = Bxn_recvbufferproj[i][1][0];
                bufferProjsend[step+4] = Byn_recvbufferproj[i][1][0];
                bufferProjsend[step+5] = Bzn_recvbufferproj[i][1][0];
                step +=6;
            }
        if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
        //also receive a message
            ierr = MPI_Sendrecv(bufferProjsend,6*(nxn-2),MPI_DOUBLE,vct->getYleft_neighbor(),1,bufferProjrecv,6*(nxn-2),MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
            ierr = MPI_Send(bufferProjsend,6*(nxn-2),MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM());
        }
    } else {
        //Do not send ...
          if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
          //... but receive only
              ierr = MPI_Recv(bufferProjrecv,6*(nxn-2),MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(), &status);
          }
    }
    if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
        //If a message from proc coorinateY+1 has been received
        step=0;
        for (i=1;i<nxn-1;i++){
                Ex_recvbufferproj[i][nyn-2][0] += bufferProjrecv[step  ];
                Ey_recvbufferproj[i][nyn-2][0] += bufferProjrecv[step+1];
                Ez_recvbufferproj[i][nyn-2][0] += bufferProjrecv[step+2];
                Bxn_recvbufferproj[i][nyn-2][0] += bufferProjrecv[step+3];
                Byn_recvbufferproj[i][nyn-2][0] += bufferProjrecv[step+4];
                Bzn_recvbufferproj[i][nyn-2][0] += bufferProjrecv[step+5];
                step +=6;
        }
    //Send a message to proc coordinateY+1
        step =0;
        for(i=1;i<nxn-1;i++){
            bufferProjsend[step  ] = Ex_recvbufferproj[i][nyn-2][0];
            bufferProjsend[step+1] = Ey_recvbufferproj[i][nyn-2][0];
            bufferProjsend[step+2] = Ez_recvbufferproj[i][nyn-2][0];
            bufferProjsend[step+3] = Bxn_recvbufferproj[i][nyn-2][0];
            bufferProjsend[step+4] = Byn_recvbufferproj[i][nyn-2][0];
            bufferProjsend[step+5] = Bzn_recvbufferproj[i][nyn-2][0];
            step +=6;
        }
        if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //also receive a message
        ierr = MPI_Sendrecv(bufferProjsend,6*(nxn-2),MPI_DOUBLE,vct->getYright_neighbor(),1,bufferProjrecv,6*(nxn-2),MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
        //send only
        //cout << vct->getCartesian_rank_COMMTOTAL() <<" send to "<< vct->getYright_neighbor() << endl;
        ierr = MPI_Send(bufferProjsend,6*(nxn-2),MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM());
        }
    
    } else {
          //Do not send ...
           if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
           //...receive only
           //cout << vct->getCartesian_rank_COMMTOTAL() <<" receive from "<< vct->getYleft_neighbor() << endl;
           ierr = MPI_Recv(bufferProjrecv,6*(nxn-2),MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(), &status);
           }
    }
    if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //A message from proc coordinateX-1 has been received
        step=0;
        for (i=1;i<nxn-1;i++){
            Ex_recvbufferproj[i][1][0] =  bufferProjrecv[step  ];
            Ey_recvbufferproj[i][1][0] =  bufferProjrecv[step+1];
            Ez_recvbufferproj[i][1][0] =  bufferProjrecv[step+2];
            Bxn_recvbufferproj[i][1][0] =  bufferProjrecv[step+3];
            Byn_recvbufferproj[i][1][0] =  bufferProjrecv[step+4];
            Bzn_recvbufferproj[i][1][0] =  bufferProjrecv[step+5];
            step +=6;
        }
    }

if (nmessagerecvProj>0) {
    //Average between the current fields and the projected ones stores in the recvbufferproj and normalization
            //cout << "receive ixrecvfirstProjglobal= "<<ixrecvfirstProjglobal<<" iyrecvfirstProjglobal= "<< iyrecvfirstProjglobal<<" ixrecvlastProjglobal= "<< ixrecvlastProjglobal<<" iyrecvlastProjglobal= "<<iyrecvlastProjglobal<<endl;
    for (i=ixrecvfirstProjglobal ;i<ixrecvlastProjglobal+1;i++){
        for (j=iyrecvfirstProjglobal;j<iyrecvlastProjglobal+1;j++){
	  Ex[i][j][0] = ( Ex[i][j][0] + Ex_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] )/2.;
	  Ey[i][j][0] = ( Ey[i][j][0] + Ey_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] )/2.;
	  Ez[i][j][0] = ( Ez[i][j][0] + Ez_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] )/2.;
	  //Bxn[i][j][0] = (Bxn[i][j][0]+ Bxn_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] )/2.;
	  //Byn[i][j][0] = (Byn[i][j][0]+ Byn_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] )/2.;
          //Bzn[i][j][0] = (Bzn[i][j][0]+ Bzn_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] )/2.;
	  // complete substitution
	  //Ex[i][j][0] =   Ex_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] ;   
	  //Ey[i][j][0] =   Ey_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] ;               
	  //Ez[i][j][0] =   Ez_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] ;               
	  ////	  Bxn[i][j][0] =  Bxn_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] ;               
	  //// Byn[i][j][0] =  Byn_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] ;               
	  //// Bzn[i][j][0] =  Bzn_recvbufferproj[i][j][0]/normalizerecvProj[i][j][0] ;
	   
	  //cout << "R" <<vct->getCartesian_rank_COMMTOTAL()<<" i " << i <<" j " <<j <<" normalizerecvProj[i][j][0]" << normalizerecvProj[i][j][0] << " coords x and y " << grid->getXN(i, j, 0) << " and " << grid->getYN(i, j, 0)  <<endl;
       }
    }
 } // edn if (nmessagerecvProj>0)
 
// this part is useful to fix the first active node of the coarse grid in case it falls in the PM / MP / PP area of the corresponding refined proc
// it just copies it from the proc to the left
// It is done on the entire coarse grid to simplify the if conditions and alson in the hope it helps with mischievous BC

// X dir
 
      if (vct->getXright_neighbor()!= MPI_PROC_NULL)
       {
	 //cout << "R" << vct->getCartesian_rank_COMMTOTAL() <<" wants to send to R" << vct->getXright_neighbor() << endl; 
	 for (int j=0; j<nyn; j++)
	   {
	     bufferProjsend[3*j]=   Ex[nxn-2][j][0];
	     bufferProjsend[3*j+1]= Ey[nxn-2][j][0];
	     bufferProjsend[3*j+2]= Ez[nxn-2][j][0];
	   }
	 
	 ierr= MPI_Isend(bufferProjsend,3*nyn,MPI_DOUBLE,vct->getXright_neighbor(),666,MPI_COMM_WORLD, &request); 
	 MPI_Wait(&request, &status); 
	 //cout <<"R" << vct->getCartesian_rank_COMMTOTAL()<< " ended send to R" << vct->getXright_neighbor() <<endl;
       }
     
     if (vct->getXleft_neighbor()!= MPI_PROC_NULL) 
       {
	 //cout <<"R" << vct->getCartesian_rank_COMMTOTAL() <<" wants to receive from R" << vct->getXleft_neighbor() <<endl;
	 ierr = MPI_Recv(bufferProjrecv,3*nyn,MPI_DOUBLE,vct->getXleft_neighbor(),666,MPI_COMM_WORLD, &status);      
	 //cout <<"R" << vct->getCartesian_rank_COMMTOTAL() <<" ended receive from R" << vct->getXleft_neighbor() <<endl;
	 for (int j=0; j<nyn; j++)
	   {
	     Ex[1][j][0] = bufferProjrecv[3*j];
	     Ey[1][j][0] = bufferProjrecv[3*j+1];
	     Ez[1][j][0] = bufferProjrecv[3*j+2];
	   }
	 
       }
     // Y dir
      if (vct->getYright_neighbor()!= MPI_PROC_NULL)
       {
	 //cout <<"R" << vct->getCartesian_rank_COMMTOTAL() <<" wants to send to R" << vct->getYright_neighbor() <<endl;
	 for (int i=0; i<nxn; i++)
	   {
	     bufferProjsend[3*i]=   Ex[i][nyn-2][0];
	     bufferProjsend[3*i+1]= Ey[i][nyn-2][0];
	     bufferProjsend[3*i+2]= Ez[i][nyn-2][0];
	   }
	 ierr = MPI_Isend(bufferProjsend,3*nxn,MPI_DOUBLE,vct->getYright_neighbor(),667,MPI_COMM_WORLD, &request);
	 MPI_Wait(&request, &status);
	 //cout <<"R" << vct->getCartesian_rank_COMMTOTAL() <<" ended send to R" << vct->getYright_neighbor() <<endl;
       }
 
     if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
       {
	 //cout <<"R" << vct->getCartesian_rank_COMMTOTAL()<< " wants to receive from R" << vct->getYleft_neighbor() <<endl;
	 ierr = MPI_Recv(bufferProjrecv,3*nxn,MPI_DOUBLE,vct->getYleft_neighbor(),667,MPI_COMM_WORLD, &status);
	 //cout <<"R" << vct->getCartesian_rank_COMMTOTAL() <<" ended receive from R" << vct->getYleft_neighbor() <<endl;
	 for (int i=0; i<nxn; i++)
	   {
	     Ex[i][1][0] = bufferProjrecv[3*i];
	     Ey[i][1][0] = bufferProjrecv[3*i+1];
	     Ez[i][1][0] = bufferProjrecv[3*i+2];
	   }
	 
       }
 
 // up to now, the ghost nodes have not been modified
 communicateNode(nxn,nyn,Ex,vct);
 communicateNode(nxn,nyn,Ey,vct);
 communicateNode(nxn,nyn,Ez,vct);
 
 calculateB_afterProj(grid, vct);
 
 //cout <<"END RECEIVE PROJ" << endl;
}
/** Receive the boundary conditions from coarser level*/
inline void EMfields::receiveBC(Grid *grid, VirtualTopology *vct, CollectiveIO *col){

int i,j,ierr,count;
double isend[12],irecv[12];
MPI_Status status;
MPI_Request request;


for (i=0;i<nmessagerecuBC;i++) {
    ierr = MPI_Recv(bufferBC,npointsreceived[i]*6,MPI_DOUBLE,fromBC[i],1,MPI_COMM_WORLD, &status);
    //Updating the electric fields on ghost nodes
    for (j=0;j<xmrecv[i];j++) {
        vectX[ixmrecvfirst[i]+j][0][0] =  bufferBC[6*j];
        vectY[ixmrecvfirst[i]+j][0][0] =  bufferBC[6*j+1];
        vectZ[ixmrecvfirst[i]+j][0][0] =  bufferBC[6*j+2];
        Bxn_new[ixmrecvfirst[i]+j][0][0] = bufferBC[6*j+3];//B_new on finer grid receives B(n+1) and stores it before it is used as BC in calculateB
        Byn_new[ixmrecvfirst[i]+j][0][0] = bufferBC[6*j+4];
        Bzn_new[ixmrecvfirst[i]+j][0][0] = bufferBC[6*j+5];
    }
    count = xmrecv[i];
    for (j=0;j<xprecv[i];j++) {
        vectX[ixprecvfirst[i]+j][grid->getNYN()-1][0] =  bufferBC[6*(count+j)];
        vectY[ixprecvfirst[i]+j][grid->getNYN()-1][0] =  bufferBC[6*(count+j)+1];
        vectZ[ixprecvfirst[i]+j][grid->getNYN()-1][0] =  bufferBC[6*(count+j)+2];
        Bxn_new[ixprecvfirst[i]+j][grid->getNYN()-1][0] = bufferBC[6*(count+j)+3];
        Byn_new[ixprecvfirst[i]+j][grid->getNYN()-1][0] = bufferBC[6*(count+j)+4];
        Bzn_new[ixprecvfirst[i]+j][grid->getNYN()-1][0] = bufferBC[6*(count+j)+5];
    }
    count = count+xprecv[i];
    for (j=0;j<ymrecv[i];j++) {
        vectX[0][iymrecvfirst[i]+j][0] =  bufferBC[6*(count+j)];
        vectY[0][iymrecvfirst[i]+j][0] =  bufferBC[6*(count+j)+1];
        vectZ[0][iymrecvfirst[i]+j][0] =  bufferBC[6*(count+j)+2];
        Bxn_new[0][iymrecvfirst[i]+j][0] = bufferBC[6*(count+j)+3];
        Byn_new[0][iymrecvfirst[i]+j][0] = bufferBC[6*(count+j)+4];
        Bzn_new[0][iymrecvfirst[i]+j][0] = bufferBC[6*(count+j)+5];
    }
    count = count+ymrecv[i];
    for (j=0;j<yprecv[i];j++) {
        vectX[grid->getNXN()-1][iyprecvfirst[i]+j][0] =  bufferBC[6*(count+j)];
        vectY[grid->getNXN()-1][iyprecvfirst[i]+j][0] =  bufferBC[6*(count+j)+1];
        vectZ[grid->getNXN()-1][iyprecvfirst[i]+j][0] =  bufferBC[6*(count+j)+2];
        Bxn_new[grid->getNXN()-1][iyprecvfirst[i]+j][0] = bufferBC[6*(count+j)+3];
        Byn_new[grid->getNXN()-1][iyprecvfirst[i]+j][0] = bufferBC[6*(count+j)+4];
        Bzn_new[grid->getNXN()-1][iyprecvfirst[i]+j][0] = bufferBC[6*(count+j)+5];
    }
}
// Now finish the communication of BC by updating the few ghost nodes that were not updated directly from the coarser grid
//communicate from right to left: 2 nodes
if (vct->getCoordinates(0)== 0 || vct->getCoordinates(0)==vct->getXLEN()-1){
    if (vct->getCoordinates(1) > 0) {
        isend[0]  = vectX[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][1][0];
        isend[1]  = vectY[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][1][0];
        isend[2]  = vectZ[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][1][0];
        isend[3]  = vectX[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][2][0];
        isend[4]  = vectY[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][2][0];
        isend[5]  = vectZ[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][2][0];
        isend[6] =  Bxn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][1][0];
        isend[7] =  Byn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][1][0];
        isend[8] =  Bzn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][1][0];
        isend[9] =  Bxn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][2][0];
        isend[10] = Byn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][2][0];
        isend[11] = Bzn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][2][0];
        if (vct->getCoordinates(1) == vct->getYLEN()-1){
            ierr = MPI_Send(isend,12,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM());
        } else {
            ierr = MPI_Sendrecv(isend,12,MPI_DOUBLE,vct->getYleft_neighbor(),1,irecv,12,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(),&status);
        }
    }

    if (vct->getCoordinates(1) < vct->getYLEN()-1) {
        if (vct->getCoordinates(1) == 0){
            ierr = MPI_Recv(irecv,12,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(),&status);
        }
        vectX[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-2][0]  = irecv[0]; 
        vectY[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-2][0]  = irecv[1];
        vectZ[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-2][0]  = irecv[2];
        vectX[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-1][0]  = irecv[3];
        vectY[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-1][0]  = irecv[4];
        vectZ[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-1][0]  = irecv[5];
        Bxn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-2][0] = irecv[6]; 
        Byn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-2][0] = irecv[7];
        Bzn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-2][0] = irecv[8];
        Bxn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-1][0] = irecv[9];
        Byn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-1][0] = irecv[10];
        Bzn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-1][0] = irecv[11];
    }
}
if (vct->getCoordinates(1)== 0 || vct->getCoordinates(1)==vct->getYLEN()-1){
    if (vct->getCoordinates(0) > 0){
        isend[0]  = vectX[1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[1]  = vectY[1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[2]  = vectZ[1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[3]  = vectX[2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[4]  = vectY[2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[5]  = vectZ[2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[6]  = Bxn_new[1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[7]  = Byn_new[1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[8]  = Bzn_new[1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[9]  = Bxn_new[2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[10] = Byn_new[2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[11] = Bzn_new[2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        if (vct->getCoordinates(0) == vct->getXLEN()-1){
            ierr = MPI_Send(isend,12,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM());
        } else {
            ierr = MPI_Sendrecv(isend,12,MPI_DOUBLE,vct->getXleft_neighbor(),1,irecv,12,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(),&status);
        }
    }
    if (vct->getCoordinates(0) < vct->getXLEN()-1) {
        if (vct->getCoordinates(0) == 0) {
            ierr = MPI_Recv(irecv,12,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(),&status);
        }
        vectX[grid->getNXN()-2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0]  = irecv[0];
        vectY[grid->getNXN()-2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0]  = irecv[1];
        vectZ[grid->getNXN()-2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0]  = irecv[2];
        vectX[grid->getNXN()-1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0]  = irecv[3];
        vectY[grid->getNXN()-1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0]  = irecv[4];
        vectZ[grid->getNXN()-1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0]  = irecv[5];
        Bxn_new[grid->getNXN()-2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[6];
        Byn_new[grid->getNXN()-2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[7];
        Bzn_new[grid->getNXN()-2][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[8];
        Bxn_new[grid->getNXN()-1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[9];
        Byn_new[grid->getNXN()-1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[10];
        Bzn_new[grid->getNXN()-1][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[11];
    }
}
//Communicate from left to right: 1 node
if (vct->getCoordinates(0)== 0 || vct->getCoordinates(0)==vct->getXLEN()-1){
    if(vct->getCoordinates(1) < vct->getYLEN()-1){
        isend[0] = vectX[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-3][0];
        isend[1] = vectY[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-3][0];
        isend[2] = vectZ[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-3][0];
        isend[3] = Bxn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-3][0];
        isend[4] = Byn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-3][0];
        isend[5] = Bzn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][grid->getNYN()-3][0];
        if (vct->getCoordinates(1) == 0){
            ierr = MPI_Send(isend,6,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM());
        } else {
            ierr = MPI_Sendrecv(isend,6,MPI_DOUBLE,vct->getYright_neighbor(),1,irecv,6,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(),&status);
        }
    }
    if (vct->getCoordinates(1) > 0) {
        if(vct->getCoordinates(1) == vct->getYLEN()-1){
            ierr = MPI_Recv(irecv,6,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(),&status);
        }
        vectX[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][0][0] = irecv[0]; 
        vectY[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][0][0] = irecv[1];
        vectZ[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][0][0] = irecv[2];
        Bxn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][0][0] = irecv[3]; 
        Byn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][0][0] = irecv[4];
        Bzn_new[vct->getCoordinates(0)/(vct->getXLEN()-1)*(grid->getNXN()-1)][0][0] = irecv[5];
    }
}
if (vct->getCoordinates(1)== 0 || vct->getCoordinates(1)==vct->getYLEN()-1){
    if(vct->getCoordinates(0) < vct->getXLEN()-1){
        isend[0] = vectX[grid->getNXN()-3][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[1] = vectY[grid->getNXN()-3][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[2] = vectZ[grid->getNXN()-3][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[3] =  Bxn_new[grid->getNXN()-3][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[4] = Byn_new[grid->getNXN()-3][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        isend[5] = Bzn_new[grid->getNXN()-3][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0];
        if (vct->getCoordinates(0) == 0){
	  //cout << "send to " << vct->getXright_neighbor() << endl ;
            ierr = MPI_Send(isend,6,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM());
        } else {
	  //cout << "send to " << vct->getXright_neighbor() <<" receive from " << vct->getXleft_neighbor() << endl ;
            ierr = MPI_Sendrecv(isend,6,MPI_DOUBLE,vct->getXright_neighbor(),1,irecv,6,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(),&status);
        }
    }
    if(vct->getCoordinates(0) > 0) {
        if(vct->getCoordinates(0) == vct->getXLEN()-1){
            //cout << "receive from  " << vct->getXleft_neighbor() << endl ;
            ierr = MPI_Recv(irecv,6,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(),&status);
        }
        vectX[0][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[0]; 
        vectY[0][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[1];
        vectZ[0][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[2];
        Bxn_new[0][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[3]; 
        Byn_new[0][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[4];
        Bzn_new[0][vct->getCoordinates(1)/(vct->getYLEN()-1)*(grid->getNYN()-1)][0] = irecv[5];
    }
}
if (nmessagerecuBC >0){
// Store boundary conditions in the appropritate buffers
    if (vct->getCoordinates(1)==0) { 
        for (j=0;j < nxn;j++){
            bufferBCExxm[j]  = vectX[j][0][0]*th+Ex[j][0][0]*(1-th);
            bufferBCEyxm[j]  = vectY[j][0][0]*th+Ey[j][0][0]*(1-th);
            bufferBCEzxm[j]  = vectZ[j][0][0]*th+Ez[j][0][0]*(1-th);
            
        }
    }
    if (vct->getCoordinates(1)==vct->getYLEN()-1){
        for (j=0;j < nxn;j++){
            bufferBCExxp[j]  = vectX[j][nyn-1][0]*th+Ex[j][nyn-1][0]*(1-th);
            bufferBCEyxp[j]  = vectY[j][nyn-1][0]*th+Ey[j][nyn-1][0]*(1-th);
            bufferBCEzxp[j]  = vectZ[j][nyn-1][0]*th+Ez[j][nyn-1][0]*(1-th);
        }
    }
    if (vct->getCoordinates(0)==0){
        for (j=0;j < nyn;j++){
            bufferBCExym[j]  = vectX[0][j][0]*th+Ex[0][j][0]*(1-th);
            bufferBCEyym[j]  = vectY[0][j][0]*th+Ey[0][j][0]*(1-th);
            bufferBCEzym[j]  = vectZ[0][j][0]*th+Ez[0][j][0]*(1-th);
        }
    }
    if (vct->getCoordinates(0)==vct->getXLEN()-1){
        for (j=0;j < nyn;j++){
            bufferBCExyp[j]  = vectX[nxn-1][j][0]*th+Ex[nxn-1][j][0]*(1-th);
            bufferBCEyyp[j]  = vectY[nxn-1][j][0]*th+Ey[nxn-1][j][0]*(1-th);
            bufferBCEzyp[j]  = vectZ[nxn-1][j][0]*th+Ez[nxn-1][j][0]*(1-th);
        }
    }
}
}
/** Output ghost values **/
inline void EMfields::outputghost(VirtualTopology *vct, CollectiveIO *col, int cycle){
int j;
ofstream outputfile;
stringstream strm;
string num;
const char * name;

strm << col->getSaveDirName();
strm << "/ghost";
strm << vct->getCartesian_rank_COMMTOTAL();

strm >> num;
name = num.c_str();

if (cycle > 0){
outputfile.open(name,ios::app);
} else {
outputfile.open(name);
}
for (j=0;j < nxn;j++){
    outputfile<<Bxn[j][0][0]<<"\n";
}
for (j=0;j < nxn;j++){
    outputfile<<Bxn[j][nyn-1][0]<<"\n";
}
for (j=1;j < nyn-1;j++){
    outputfile<<Bxn[0][j][0]<<"\n";
}
for (j=1;j < nyn-1;j++){
    outputfile<<Bxn[nxn-1][j][0]<<"\n";
    }

 /*if (cycle > 0){
   outputfile.open(name,ios::app);
 } else {
   outputfile.open(name);
 }
 for (j=0;j < nxn;j++){
   outputfile<<Ez[j][0][0]<<"\n";
 }
 for (j=0;j < nxn;j++){
   outputfile<<Ez[j][nyn-1][0]<<"\n";
 }
 for (j=1;j < nyn-1;j++){
   outputfile<<Ez[0][j][0]<<"\n";
 }
 for (j=1;j < nyn-1;j++){
   outputfile<<Ez[nxn-1][j][0]<<"\n";
 }

 /*if (cycle > 0){
   outputfile.open(name,ios::app);
 } else {
   outputfile.open(name);
 }
 
 for (j=0;j < nxn;j++){
   outputfile<<rhons[0][j][0][0]<<"\n";
 }
 for (j=0;j < nxn;j++){
   outputfile<<rhons[0][j][nyn-1][0]<<"\n";
 }
 for (j=1;j < nyn-1;j++){
   outputfile<<rhons[0][0][j][0]<<"\n";
 }
 for (j=1;j < nyn-1;j++){
   outputfile<<rhons[0][nxn-1][j][0]<<"\n";
   }*/
outputfile.close();

}

/** Apply direct electrif field boundary conditions*/
inline void EMfields::electricfieldBC(double ***componentX,double ***componentY, double ***componentZ, Grid *grid, VirtualTopology *vct){
int j;
    if (vct->getCoordinates(1)==0) { 
        for (j=0;j < nxn;j++){
           componentX[j][0][0]=     bufferBCExxm[j]  ; 
           componentY[j][0][0]=     bufferBCEyxm[j]  ; 
           componentZ[j][0][0]=     bufferBCEzxm[j]  ; 
        } 
    }
    if (vct->getCoordinates(1)==vct->getYLEN()-1) { 
        for (j=0;j < nxn;j++){
           componentX[j][nyn-1][0]= bufferBCExxp[j]  ; 
           componentY[j][nyn-1][0]= bufferBCEyxp[j]  ; 
           componentZ[j][nyn-1][0]= bufferBCEzxp[j]  ; 
        } 
    }
    if (vct->getCoordinates(0)==0) { 
        for (j=0;j < nyn;j++){
           componentX[0][j][0]    = bufferBCExym[j]; 
           componentY[0][j][0]    = bufferBCEyym[j]; 
           componentZ[0][j][0]    = bufferBCEzym[j]; 
        } 
    }
    if (vct->getCoordinates(0)==vct->getXLEN()-1) { 
        for (j=0;j < nyn;j++){
           componentX[nxn-1][j][0]= bufferBCExyp[j]; 
           componentY[nxn-1][j][0]= bufferBCEyyp[j]; 
           componentZ[nxn-1][j][0]= bufferBCEzyp[j]; 
        } 
    }
}
/** Apply magnetic field boundary conditions*/
inline void EMfields::magneticfieldBC(double ***component, Grid *grid, double ***from, VirtualTopology *vct){
int j;
    if (vct->getCoordinates(1)==0) { 
        for (j=0;j < nxn;j++){
           component[j][0][0]    = from[j][0][0]; 
        } 
    }
    if (vct->getCoordinates(1)==vct->getYLEN()-1) { 
        for (j=0;j < nxn;j++){
           component[j][nyn-1][0]= from[j][nyn-1][0]; 
        } 
    }
    if (vct->getCoordinates(0)==0) {
        for (j=0;j < nyn;j++){
           component[0][j][0]    = from[0][j][0] ; 
        } 
    }
    if (vct->getCoordinates(0)==vct->getXLEN()-1) {
        for (j=0;j < nyn;j++){
           component[nxn-1][j][0]= from[nxn-1][j][0]; 
        } 
    }
}

/** get Potential array ***/
inline double*** EMfields::getPHI() {return(PHI);}
/** get Ex(X,Y,Z)  */
inline double &EMfields::getEx(int indexX, int indexY, int indexZ) const{
	return(Ex[indexX][indexY][0]);}
/** get Electric field  component X array*/
inline double*** EMfields::getEx() {

return(Ex);

}
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
	return(Jxh[indexX][indexY][0]);}
/** get current array X component **/
inline double*** EMfields::getJx() {return(Jx);}
/** get current -Direction Y */
inline double &EMfields::getJy(int indexX, int indexY, int indexZ) const{
	return(Jyh[indexX][indexY][0]);}
/** get current array Y component **/
inline double*** EMfields::getJy() {return(Jy);}
/** get current -Direction Z */
inline double &EMfields::getJz(int indexX, int indexY, int indexZ) const{
	return(Jzh[indexX][indexY][0]);}
/** get current array Z component **/
inline double*** EMfields::getJz() {return(Jz);}
/**SPECIES: get current array X component */
inline double**** EMfields::getJxs() {return(Jxs);}
/** get Jxs(X,Y,Z,is) : density for species*/
inline double &EMfields::getJxs(int indexX, int indexY, int indexZ, int is) const{
	return(Jxs[is][indexX][indexY][0]);}
/**SPECIES: get current array Y component */
inline double**** EMfields::getJys() {return(Jys);}
/** get Jys(X,Y,Z,is) : density for species*/
inline double &EMfields::getJys(int indexX, int indexY, int indexZ, int is) const{
	return(Jys[is][indexX][indexY][0]);}
/**SPECIES: get current array Z component */
inline double**** EMfields::getJzs() {return(Jzs);}
/** get Jzs(X,Y,Z,is) : density for species*/
inline double &EMfields::getJzs(int indexX, int indexY, int indexZ, int is) const{
	return(Jzs[is][indexX][indexY][0]);}
//ME
inline double &EMfields::getpXXns(int indexX, int indexY, int indexZ, int is) const{
  return(pXXsn[is][indexX][indexY][0]);}
inline double &EMfields::getpXYns(int indexX, int indexY, int indexZ, int is) const{
  return(pXYsn[is][indexX][indexY][0]);}
inline double &EMfields::getpXZns(int indexX, int indexY, int indexZ, int is) const{
  return(pXZsn[is][indexX][indexY][0]);}
inline double &EMfields::getpYYns(int indexX, int indexY, int indexZ, int is) const{
  return(pYYsn[is][indexX][indexY][0]);}
inline double &EMfields::getpYZns(int indexX, int indexY, int indexZ, int is) const{
  return(pYZsn[is][indexX][indexY][0]);}
inline double &EMfields::getpZZns(int indexX, int indexY, int indexZ, int is) const{
  return(pZZsn[is][indexX][indexY][0]);}

//sets
inline void EMfields::setRHOns(double value, int indexX, int indexY,int indexZ,int is)
{rhons[is][indexX][indexY][0]=value;}
inline void EMfields::setJxs(double value, int indexX, int indexY,int indexZ,int is)
{Jxs[is][indexX][indexY][0]=value;}
inline void EMfields::setJys(double value, int indexX, int indexY,int indexZ,int is)
{Jys[is][indexX][indexY][0]=value;}
inline void EMfields::setJzs(double value, int indexX, int indexY,int indexZ,int is)
{Jzs[is][indexX][indexY][0]=value;}
inline void EMfields::setpXXns(double value, int indexX, int indexY,int indexZ,int is)
{pXXsn[is][indexX][indexY][0]=value;}
inline void EMfields::setpXYns(double value, int indexX, int indexY,int indexZ,int is)
{pXYsn[is][indexX][indexY][0]=value;}
inline void EMfields::setpXZns(double value, int indexX, int indexY,int indexZ,int is)
{pXZsn[is][indexX][indexY][0]=value;}
inline void EMfields::setpYYns(double value, int indexX, int indexY,int indexZ,int is)
{pYYsn[is][indexX][indexY][0]=value;}
inline void EMfields::setpYZns(double value, int indexX, int indexY,int indexZ,int is)
{pYZsn[is][indexX][indexY][0]=value;}
inline void EMfields::setpZZns(double value, int indexX, int indexY,int indexZ,int is)
{pZZsn[is][indexX][indexY][0]=value;}

/** Print info about electromagnetic field */
inline void EMfields::print(void) const{


}
/** constructor */
inline EMfields::EMfields(CollectiveIO *col,Grid *grid, VCtopology *vct){
        ratio = col->getRatio();
	nxc = grid->getNXC();
	nxn = grid->getNXN();
	nyc = grid->getNYC();
	nyn = grid->getNYN();
	// Arnaud
        //nxnproj = nxn /col->getRatio() +2;
        //nynproj = nyn /col->getRatio() +2;
	//cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << ": WithINT: nxnproj " << nxnproj <<" nynproj "<<nynproj << " nxn " << nxn << " nyn " << nyn  << " ratio " << col->getRatio()<<endl;
	// ME
	nxnproj = ceil((double)nxn /(double)col->getRatio()) +2;
        nynproj = ceil((double)nyn /(double)col->getRatio()) +2;
	cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << ": WithDOUBLE: nxnproj " << nxnproj <<" nynproj "<<nynproj << " nxn " << nxn << " nyn " << nyn  << " ratio " << col->getRatio()<<endl; 
	// end ME
	dx = grid->getDX();
	dy = grid ->getDY();
	invVOL = grid->getInvVOL();
	xStart = grid->getXstart();
	xEnd   = grid->getXend();
	yStart = grid->getYstart();
	yEnd   = grid->getYend();
	Lx = col->getLx();
	Ly = col->getLy();
        Ox = grid->getOx(grid->getLevel());
        Oy = grid->getOy(grid->getLevel());
	ns  = col->getNs();
	c = col->getC();
	dt = col->getDt();
	Smooth = col->getSmooth();
	Nvolte = col ->getNvolte();
	th = col->getTh();
	delt = c*th*dt;
        xnnl = 0;
        xnnu = 0;
        ynnl = 0;
        ynnu = 0;
        nmessageBC = 0;
        nmessagerecuBC = 0;
        nmessageProj = 1;
        nmessagerecvProj = 0;
        omega=1.;
        Vinj = col->getVinj();
        if (vct->getCoordinates(1)== col->getYLEN()-1){
            iyend = grid->getNYN()-1;
        } else {
            iyend = grid->getNYN()-2;
        }
        if (vct->getCoordinates(0)== col->getXLEN()-1){
            ixend = grid->getNXN()-1;
        } else {
            ixend = grid->getNXN()-2;
        }
	// FLAG ON POISSON CORRECTION
	// make always the divergence cleaning


	qom = new double[ns];
	targetBC  = new int[col->getXLEN()*col->getYLEN()];
	// AMR, ME
	BCSide  = new int[col->getXLEN()*col->getYLEN()];
	BCSidecu= new int[col->getXLEN()*col->getYLEN()];
	xfirst_COARSE= new double[4*col->getXLEN()*col->getYLEN()];// 4 are the sides
	yfirst_COARSE= new double[4*col->getXLEN()*col->getYLEN()];// 4 are the sides 
	// end AMR, ME
        targetProj = new int[4];
	npointssent  = new int[col->getXLEN()*col->getYLEN()];
	fromBC  = new int[col->getXLEN()*col->getYLEN()];//int[4];
        fromProj = new int[col->getXLEN()*col->getYLEN()];
	xmrecv  = new int[4];
	xprecv  = new int[4];
	ymrecv  = new int[4];
	yprecv  = new int[4];
	nxrecvProj  = new int[col->getXLEN()*col->getYLEN()];// Number of points received during Proj phase in the x direction by each of the finer grid process
	nyrecvProj  = new int[col->getXLEN()*col->getYLEN()];// Number of points received during Proj phase in the y direction by each of the finer grid process
	ixmrecvfirst  = new int[4];
	ixprecvfirst  = new int[4];
	iymrecvfirst  = new int[4];
	iyprecvfirst  = new int[4];
	npointsreceived  = new int[4];
        npointsreceivedProj = new int[col->getXLEN()*col->getYLEN()];
	ixrecvfirstProj  = new int[col->getXLEN()*col->getYLEN()];
	iyrecvfirstProj  = new int[col->getXLEN()*col->getYLEN()];
        ixsent = new int[2*(nxn+nyn-2)*(int)ceil(ratio)];
        iysent = new int[2*(nxn+nyn-2)*(int)ceil(ratio)];
        ixsentProj = new int[nxn];
        iysentProj = new int[nyn];
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
	Ex_old = newArr3(double,nxn,nyn,1);
	Ey_old = newArr3(double,nxn,nyn,1);
	Ez_old = newArr3(double,nxn,nyn,1);
	Ex_recvbufferproj = newArr3(double,nxn,nyn,1);
	Ey_recvbufferproj = newArr3(double,nxn,nyn,1);
	Ez_recvbufferproj = newArr3(double,nxn,nyn,1);
	Exth = newArr3(double,nxn,nyn,1);
	Eyth = newArr3(double,nxn,nyn,1);
	Ezth = newArr3(double,nxn,nyn,1);
	Bxc_old   = newArr3(double,nxc,nyc,1);
	Byc_old   = newArr3(double,nxc,nyc,1);
	Bzc_old   = newArr3(double,nxc,nyc,1);
	Bxn   = newArr3(double,nxn,nyn,1);
	Byn   = newArr3(double,nxn,nyn,1);
	Bzn   = newArr3(double,nxn,nyn,1);
	Bxn_new   = newArr3(double,nxn,nyn,1);
	Byn_new   = newArr3(double,nxn,nyn,1);
	Bzn_new   = newArr3(double,nxn,nyn,1);
	Bxn_recvbufferproj   = newArr3(double,nxn,nyn,1);
	Byn_recvbufferproj   = newArr3(double,nxn,nyn,1);
	Bzn_recvbufferproj   = newArr3(double,nxn,nyn,1);
	rhon  = newArr3(double,nxn,nyn,1);
	Jx    = newArr3(double,nxn,nyn,1);
	Jy    = newArr3(double,nxn,nyn,1);
	Jz    = newArr3(double,nxn,nyn,1);
	Jxh   = newArr3(double,nxn,nyn,1);
	Jyh   = newArr3(double,nxn,nyn,1);
	Jzh   = newArr3(double,nxn,nyn,1);
        // 2*(nxn+nyn-2) ghost nodes. 4 weights for each of them
	weightBC   = newArr3(double,2*(nxn+nyn-2)*(int)ceil(ratio),4,1);
        // maximum number of points from the coarse grid in the fine grid. 4 weights for each of them
	weightProj   = newArr4(double,nxn,nyn,4,1);
    normalizeProj = newArr3(double,nxnproj,nynproj,1);
    reducedEx = newArr3(double,nxnproj,nynproj,1);
    reducedEy = newArr3(double,nxnproj,nynproj,1);
    reducedEz = newArr3(double,nxnproj,nynproj,1);
    reducedBxn= newArr3(double,nxnproj,nynproj,1);
    reducedByn= newArr3(double,nxnproj,nynproj,1);
    reducedBzn= newArr3(double,nxnproj,nynproj,1);
    normalizerecvProj = newArr3(double,nxn,nyn,1);
        // 2*(nxn+nyn-2) ghost nodes. 6 components of electric field, new magnetic field
        bufferBC = new double[2*(nxn+nyn-2)*6*(int)ceil(ratio)]; //buffer sent by MPI
        //Buffers in which BC are stored
        bufferBCExxm = new double[nxn];
        bufferBCExxp = new double[nxn];
        bufferBCExym = new double[nyn];
        bufferBCExyp = new double[nyn];
        bufferBCEyxm = new double[nxn];
        bufferBCEyxp = new double[nxn];
        bufferBCEyym = new double[nyn];
        bufferBCEyyp = new double[nyn];
        bufferBCEzxm = new double[nxn];
        bufferBCEzxp = new double[nxn];
        bufferBCEzym = new double[nyn];
        bufferBCEzyp = new double[nyn];
                //nxnproj*nynproj projected nodes. 6 components of electric and magnetic field.
        bufferProj = new double[nxnproj*nynproj*6];
        bufferProjsend = new double[2*max(nxn,nyn)*6];
        bufferProjrecv = new double[2*max(nxn,nyn)*6];
	// ME
	// the 6 initial double are to store the map info
	INFObufferProjsend = new double[10+2*max(nxn,nyn)*6];
        INFObufferProj = new double[10+nxnproj*nynproj*6];// declared as bufferProj + the extra 10 double in the header
	// end ME
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
	delArr3(Ex_old,nxn,nyn);
	delArr3(Ey_old,nxn,nyn);
	delArr3(Ez_old,nxn,nyn);
	delArr3(Ex_recvbufferproj,nxn,nyn);
	delArr3(Ey_recvbufferproj,nxn,nyn);
	delArr3(Ez_recvbufferproj,nxn,nyn);
	delArr3(Exth,nxn,nyn);
	delArr3(Eyth,nxn,nyn);
	delArr3(Ezth,nxn,nyn);
	delArr3(Bxn,nxn,nyn);
	delArr3(Byn,nxn,nyn);
	delArr3(Bzn,nxn,nyn);
	delArr3(Bxc_old,nxc,nyc);
	delArr3(Byc_old,nxc,nyc);
	delArr3(Bzc_old,nxc,nyc);
	delArr3(Bxn_new,nxn,nyn);
	delArr3(Byn_new,nxn,nyn);
	delArr3(Bzn_new,nxn,nyn);
	delArr3(Bxn_recvbufferproj,nxn,nyn);
	delArr3(Byn_recvbufferproj,nxn,nyn);
	delArr3(Bzn_recvbufferproj,nxn,nyn);
	delArr3(rhon,nxn,nyn);
	delArr3(Jx,nxn,nyn);
	delArr3(Jy,nxn,nyn);
	delArr3(Jz,nxn,nyn);
	delArr3(Jxh,nxn,nyn);
	delArr3(Jyh,nxn,nyn);
	delArr3(Jzh,nxn,nyn);
	delArr3(weightBC,2*(nxn+nyn-2)*ratio,4);
	delArr4(weightProj,nxn,nyn,4);
	delArr3(normalizeProj,nxnproj,nynproj);
	delArr3(reducedEx,nxnproj,nynproj);
	delArr3(reducedEy,nxnproj,nynproj);
	delArr3(reducedEz,nxnproj,nynproj);
	delArr3(reducedBxn,nxnproj,nynproj);
	delArr3(reducedByn,nxnproj,nynproj);
	delArr3(reducedBzn,nxnproj,nynproj);
	delArr3(normalizerecvProj,nxn,nyn);
        delete[] bufferBC;
        delete[] bufferProj;
        delete[] bufferBCExxm ;
        delete[] bufferBCExxp ;
        delete[] bufferBCExym ;
        delete[] bufferBCExyp ;
        delete[] bufferBCEyxm ;
        delete[] bufferBCEyxp ;
        delete[] bufferBCEyym ;
        delete[] bufferBCEyyp ;
        delete[] bufferBCEzxm ;
        delete[] bufferBCEzxp ;
        delete[] bufferBCEzym ;
        delete[] bufferBCEzyp ;
       
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

	// AMR, ME
	delete[] BCSide;
	// end AMR, ME

}

// AMR, ME
// used by PRASend
inline int EMfields::getNmessageBC(){return nmessageBC;}
inline int *EMfields::getTargetBC(){return targetBC;}
inline int *EMfields::getBCSide(){return BCSide;}
inline int EMfields::getTargetBOTTOM(){return targetBOTTOM;}
inline int EMfields::getTargetTOP(){return targetTOP;}
inline int EMfields::getTargetLEFT(){return targetLEFT;}
inline int EMfields::getTargetRIGHT(){return targetRIGHT;}
// used by PRAReceive
inline int EMfields::getNmessagerecuBC(){return nmessagerecuBC;}
// end AMR,ME

inline void EMfields::printRhons(int is, VirtualTopology *vct)
{
  cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << ", local rho species "<< is << endl;
  for (int i=0; i< nxn; i++)
  {
    cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<"\t";
    for (int j=0; j<nyn; j++)
    {
      cout << "(" << i << ", " <<j <<"): " << rhons[is][i][j][0] << "\t";
    }
    cout << "\n"; 
  }
}

inline void EMfields::printRhoc(VirtualTopology *vct)
{
  cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << ", local rho center " << endl;
  for (int i=0; i< nxc; i++)
  {
    cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<"\t";
    for (int j=0; j<nyc; j++)
    {
      cout << "(" << i << ", " <<j <<"): " << rhoc[i][j][0] << "\t";
    }
    cout << "\n";
  }
}

inline void EMfields::setRhoc(double val)
{
  for (int i=0; i< nxc; i++)
  {
    for (int j=0; j<nyc; j++)
    {
      rhoc[i][j][0]= val;
    }
  } 
}

inline void EMfields::printVecN(double ***vec, VirtualTopology *vct)
{
  //cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<" vec on N" << endl;
  for (int i=0; i< nxn; i++)
  {
    cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<"\t";
    for (int j=0; j<nyn; j++)
    {
      cout << "(" << i << ", " <<j <<"): " << vec[i][j][0] << "\t";
    }
    cout << "\n";
  }
}

inline void EMfields::printVecC(double ***vec, VirtualTopology *vct)
{
  //cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << " vec on C" << endl;
  for (int i=0; i< nxc; i++)
  {
    cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<"\t";
    for (int j=0; j<nyc; j++)
    {
      cout << "(" << i << ", " <<j <<"): " << vec[i][j][0] << "\t";
    }
    cout << "\n";
  }
}

inline void EMfields::printVecN(double ****vec, VirtualTopology *vct, int is)
{
  //cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<" vec on N" << endl;
  for (int i=0; i< nxn; i++)
  {
    cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<"\t";
    for (int j=0; j<nyn; j++)
    {
      cout << "(" << i << ", " <<j <<"): " << vec[is][i][j][0] << "\t";
    }
    cout << "\n";
  }
}

inline void EMfields::printVecC(double ****vec, VirtualTopology *vct, int is)
{
  //cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << " vec on C" << endl;
  for (int i=0; i< nxc; i++)
    {
      cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<"\t";
      for (int j=0; j<nyc; j++)
	{
	  cout << "(" << i << ", " <<j <<"): " << vec[is][i][j][0] << "\t";
	}
      cout << "\n";
    }
}

inline void EMfields::printJxs(int ns, VirtualTopology *vct)
{
  cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << ", local Jx species "<< ns << endl;
  printVecN(Jxs, vct, ns);
}

inline void EMfields::printPxxsn(int ns, VirtualTopology *vct)
{
  cout << "R" <<vct->getCartesian_rank_COMMTOTAL() << ", local pXXns species "<< ns << endl;
  printVecN(pXXsn, vct, ns);
}

inline void EMfields::calculateB_afterProj(Grid *grid, VirtualTopology *vct){
  int i,j;


  //Recalculate the new Eth
  for (i=0;i<nxn;i++){
      for (j=0;j<nyn;j++){
	
	Exth[i][j][0] = (1-th)*Ex_old[i][j][0]+th*Ex[i][j][0];
	Eyth[i][j][0] = (1-th)*Ey_old[i][j][0]+th*Ey[i][j][0];
	Ezth[i][j][0] = (1-th)*Ez_old[i][j][0]+th*Ez[i][j][0];

      }
  }        


  //Go back to the old B at centers
  for (i=0;i<nxc;i++){
     for (j=0;j<nyc;j++){
	 Bxc[i][j][0] = Bxc_old[i][j][0];
	 Byc[i][j][0] = Byc_old[i][j][0];
	 Bzc[i][j][0] = Bzc_old[i][j][0];
     }
  }        
  
  //Recalculate B
  eqValue (0.0,tempXC,nxc,nyc);
  eqValue (0.0,tempYC,nxc,nyc);
  eqValue (0.0,tempZC,nxc,nyc);

  // calculate the curl of Eth                                                                                     
  grid->curlN2C(tempXC,tempYC,tempZC,Exth,Eyth,Ezth); // I don't care about BC; just about the internal domain

  // update the magnetic field                                                                                     
  addscale(-c*dt,1,Bxc,tempXC,nxc,nyc);
  addscale(-c*dt,1,Byc,tempYC,nxc,nyc);
  addscale(-c*dt,1,Bzc,tempZC,nxc,nyc);

  // communicate ghost                                                                           
  communicateCenter(nxc,nyc,Bxc,vct);
  communicateCenter(nxc,nyc,Byc,vct);
  communicateCenter(nxc,nyc,Bzc,vct);
  

  // update nodes, which are not received by projection
  // communicateNode is already in the intepolation
  grid->interpC2N(Bxn, Bxc, vct);
  grid->interpC2N(Byn, Byc, vct);
  grid->interpC2N(Bzn, Bzc, vct);

}

/** Initialize the weight for the BC projections between grids **/
/** Also initialize the number and size of messages that must be sent between grids every cycle for the BC **/
/** The first xnnl quadruplets are the weights for the y=Oy-finedy ghost nodes**/
/**The next xnnu quadruplets are the weights for the y=Oy+finely+finedy ghost nodes**/
/** The next ynnl quadruplets are the weights for the x=Ox-finedx ghost nodes**/
/** The last ynnu quadruplets are the weights for the x=Ox+finelx+finedx ghost nodes**/

/** modified by ME to avoid problems when nodes sent and receive do not match**/
inline int EMfields::initWeightBC(VirtualTopology *vct, Grid *grid, CollectiveIO* col){
  int i,j,ix,iy,nproc;
  double finedx, finedy,finelx,finely, xfirst, xlast,xfirstnext, Ox, Oy,xloc,yloc;
  double coarsedx, coarsedy,coarselx,coarsely;
  double finelxplusfinedx;
  double xshift, yshift;

  if (vct->getNgrids()<2)
    return 1;

  nmessagerecuBC=0;
  nmessageBC=0;

  // refined level part
  nmessagerecuBCLEFT=0;
  nmessagerecuBCRIGHT=0;
  nmessagerecuBCBOTTOM=0;
  nmessagerecuBCTOP=0;

  npointsreceivedBOTTOM=0;
  npointsreceivedTOP=0;
  npointsreceivedLEFT=0;
  npointsreceivedRIGHT=0;

  if ( grid->getLevel() > 0) {                             //If this grid is considered as fine by another grid
    coarsedx = grid->getDX()*col->getRatio();
    coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel()-1);
    coarsedy = grid->getDY()*col->getRatio();
    coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel()-1);
    Ox = grid->getOx(grid->getLevel()); //Origin x of the grid
    Oy = grid->getOy(grid->getLevel()); //Origin y of the grid
    j=0;
    // i spans the points, j the procs
    if(vct->getCoordinates(1) == 0) {
      xfirst = grid->getmodifiedXstart(vct);
      if(vct->getCoordinates(0) == vct->getXLEN()-1) {
	xlast = grid->getXend()+grid->getDX(); 
      }else{
	xlast = grid->getXend()-grid->getDX();
      }
        
      yloc = max(Oy - grid->getDY(),0.);// Because when y < 0, it is on the same proc as if y=0  
      xnnl = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
      for (i=0; i< xnnl; i++) {
	xloc = max(Ox + xfirst +i*grid->getDX(),0.);// Because when x < 0, it is on the same proc as if x=0  
	nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	if (i==0){  // first point
	  fromBC[0] = nproc;
	  nmessagerecuBC++;
	  nmessagerecuBCBOTTOM++;
	  BCSidecu[0]=0;
	  npointsreceived[0]=0;
	  xmrecv[0]=0;
	  xprecv[0]=0;
	  ymrecv[0]=0;
	  yprecv[0]=0;
	  // ME
	  xfirst_COARSE[0]=xloc;
	  yfirst_COARSE[0]=yloc;
	  // end ME
	}
	if(nproc != fromBC[j]){ //first point in a new proc
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCBOTTOM++;
	  BCSidecu[j]=0;
	  fromBC[j]=nproc; 
	  npointsreceived[j]=0;
	  xmrecv[j]=0;
	  xprecv[j]=0;
	  ymrecv[j]=0;
	  yprecv[j]=0;
	  // ME   
	  xfirst_COARSE[j]=xloc;
          yfirst_COARSE[j]=yloc;
          // end ME 

	}
	npointsreceived[j]++;
	ixmrecvfirst[j]=floor((xloc-grid->getXstart()-Ox)/grid->getDX()+0.5)-xmrecv[j]+1; //+1 because 0 is the ghost cell at x=Xstart-dx.
	xmrecv[j]++;

      }
    } 
    if(vct->getCoordinates(1) == vct->getYLEN()-1) {
      xfirst = grid->getmodifiedXstart(vct);
      if(vct->getCoordinates(0) == vct->getXLEN()-1) {
	xlast = grid->getXend()+grid->getDX(); 
      }else{
	xlast = grid->getXend()-grid->getDX();
      }
      yloc = min(Oy + grid->getYend()+grid->getDY(),coarsely);  
      xnnu = floor((xlast-xfirst)/grid->getDX()+0.5)+1; 
      for (i=0; i< xnnu; i++) {
	xloc = max(Ox + xfirst +i*grid->getDX(),0.);
	nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	if (i==0){  // first point
	  if(nmessagerecuBC>0){
	    j++;
	  }
	  fromBC[j] = nproc;
	  nmessagerecuBC++;
	  nmessagerecuBCTOP++;
	  BCSidecu[j]=1;
	  npointsreceived[j]=0;
	  xmrecv[j]=0;
	  xprecv[j]=0;
	  ymrecv[j]=0;
	  yprecv[j]=0;
	  // ME  
	  xfirst_COARSE[j]=xloc;
          yfirst_COARSE[j]=yloc;
          // end ME

	}
	if(nproc != fromBC[j]){ // first point in a new proc
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCTOP++;
	  BCSidecu[j]=1;
	  fromBC[j]=nproc; 
	  npointsreceived[j]=0;
	  xmrecv[j]=0;
	  xprecv[j]=0;
	  ymrecv[j]=0;
	  yprecv[j]=0;
	  // ME                                                                                                                                                                 
          xfirst_COARSE[j]=xloc;
          yfirst_COARSE[j]=yloc;
          // end ME  
	}
	npointsreceived[j]++;
	ixprecvfirst[j]=floor((xloc-grid->getXstart()-Ox)/grid->getDX()+0.5)-xprecv[j]+1;
	xprecv[j]++;

      }
    } 
    if(vct->getCoordinates(0) == 0) {
      xfirst = grid->getYstart();
      if(vct->getCoordinates(1) == vct->getYLEN()-1) {
	xlast = grid->getYend(); 
      }else{
	xlast = grid->getYend()-grid->getDY();
      }
      xloc = max(Ox -grid->getDX(),0.);  
      ynnl = floor((xlast-xfirst)/grid->getDY()+0.5)+1; 
      for (i=0; i< ynnl; i++) {
	yloc = Oy + xfirst +i*grid->getDY();
	nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	if (i==0){
	  if(nmessagerecuBC>0){
	    j++;
	  }
	  fromBC[j] = nproc;
	  nmessagerecuBC++;
	  nmessagerecuBCLEFT++;
	  BCSidecu[j]=2;
	  npointsreceived[j]=0;
	  xmrecv[j]=0;
	  xprecv[j]=0;
	  ymrecv[j]=0;
	  yprecv[j]=0;
	  // ME  
	  xfirst_COARSE[j]=xloc;
          yfirst_COARSE[j]=yloc;
          // end ME  
	}
	if(nproc != fromBC[j]){
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCLEFT++;
	  BCSidecu[j]=2;
	  fromBC[j]=nproc; 
	  npointsreceived[j]=0;
	  xmrecv[j]=0;
	  xprecv[j]=0;
	  ymrecv[j]=0;
	  yprecv[j]=0;
	  // ME                
	  xfirst_COARSE[j]=xloc;
          yfirst_COARSE[j]=yloc;
          // end ME  
	}
	npointsreceived[j]++;
	iymrecvfirst[j]=floor((yloc-grid->getYstart()-Oy)/grid->getDY()+0.5)-ymrecv[j]+1;
	ymrecv[j]++;

      }
    } 

    if(vct->getCoordinates(0) == vct->getXLEN()-1) {
      xfirst = grid->getYstart();
      if(vct->getCoordinates(1) == vct->getYLEN()-1) {
	xlast = grid->getYend(); 
      }else{
	xlast = grid->getYend()-grid->getDY();
      }
      xloc = min(Ox +grid->getXend()+grid->getDX(),coarselx);  
      ynnu = floor((xlast-xfirst)/grid->getDY()+0.5)+1; 
      for (i=0; i< ynnu; i++) {
	yloc = Oy + xfirst +i*grid->getDY();
	nproc =floor(yloc/((grid->getNYC()-2)*coarsedy))+floor(xloc/((grid->getNXC()-2)*coarsedx))*col->getYLEN()+col->getYLEN()*col->getXLEN()*(grid->getLevel()-1); // rank of the proc on the coarse grid sending this point(in MPI_COMM_WORLD)
	if (i==0){
	  if(nmessagerecuBC>0){
	    j++;
	  }
	  fromBC[j] = nproc;
	  nmessagerecuBC++;
	  nmessagerecuBCRIGHT++;
	  BCSidecu[j]=3;
	  npointsreceived[j]=0;
	  xmrecv[j]=0;
	  xprecv[j]=0;
	  ymrecv[j]=0;
	  yprecv[j]=0;
	  // ME                
	  xfirst_COARSE[j]=xloc;
          yfirst_COARSE[j]=yloc;
          // end ME  
	}
	if(nproc != fromBC[j]){
	  j++;
	  nmessagerecuBC++;
	  nmessagerecuBCRIGHT++;
	  BCSidecu[j]=3;
	  fromBC[j]=nproc; 
	  npointsreceived[j]=0;
	  xmrecv[j]=0;
	  xprecv[j]=0;
	  ymrecv[j]=0;
	  yprecv[j]=0;
	  // ME                                                                                                                                                                     
          xfirst_COARSE[j]=xloc;
          yfirst_COARSE[j]=yloc;
          // end ME  
	}
	npointsreceived[j]++;
	iyprecvfirst[j]=floor((yloc-grid->getYstart()-Oy)/grid->getDY()+0.5)-yprecv[j]+1;
	yprecv[j]++;

      }
    } 

    for (i=0;i<nmessagerecuBC;i++){
      cout << "R" << vct->getCartesian_rank_COMMTOTAL() <<  ": receiving BC "<<npointsreceived[i]<<" points from "<<fromBC[i] <<", side " << BCSidecu[i]<<endl;
    }

  } 

  // end refined level part

  // send/ receive

  // each refined proc send to all coarse grid procs its list of senders
  MPI_Request requestISend;
  MPI_Status status;
  int TagMap=0;

  // exchange the map
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

  int *buffer=new int[col->getXLEN()*col->getYLEN()*4];
  if (grid->getLevel() == 0 )
    {
      nmessageBC=0;


      for (int r=0; r< vct->getXLEN()*vct->getYLEN(); r++)
	{
	  MPI_Recv(buffer, col->getXLEN()*col->getYLEN()*4,MPI_INT, MPI_ANY_SOURCE, TagMap, vct->getCART_COMM_TOTAL(), &status);
	  for (int i=0; i< col->getXLEN()*col->getYLEN()*4; i++)
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

  MPI_Barrier(vct->getCART_COMM());
  if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())))
    {
      cout << "Level " << grid->getLevel() << " finished exchanging the map for BC\n";
    }

  // debug, can be commented
  int ALL_TARGETS;
  int ALL_RECEIVERS;

  MPI_Allreduce( &nmessageBC, &ALL_TARGETS, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
  MPI_Allreduce ( &nmessagerecuBC, &ALL_RECEIVERS, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

  if (ALL_TARGETS!=ALL_RECEIVERS)
    {
      if (vct->getCartesian_rank_COMMTOTAL()==0)
	{
	  cout <<"FIELDS: ALL_TARGETS= " <<ALL_TARGETS <<"!= ALL_RECEIVERS= "<<ALL_RECEIVERS <<endl;
	  return -1;
	}
    }
  else
    {
      if (vct->getCartesian_rank_COMMTOTAL()==0)
        {
	  cout << "FIELDS: everything OK with ALL_TARGETS, ALL_RECEIVERS" <<endl;
	}
    }
  // end debug, can be commented
  // end exchange the map


  // exchange the info
  
  // refined grid: prepare the info & send

  int TagInfo=1;
  int NInfo=5;   //how many doubles to send per proc                                                                                                                             
  double *info=new double[NInfo];

  if (grid->getLevel() > 0 )
    {
      for (int i=0; i< nmessagerecuBC; i++)
	{
	  info[0]= (double)vct->getCartesian_rank_COMMTOTAL();
	  info[1]= (double)npointsreceived[i];
	  info[2]= (double)BCSidecu[i];
	  info[3]= xfirst_COARSE[i];
	  info[4]= yfirst_COARSE[i];

	  MPI_Isend (info, NInfo, MPI_DOUBLE, fromBC[i], TagInfo, vct->getCART_COMM_TOTAL(), &requestISend);
	  MPI_Wait(&requestISend, &status);
	}
    }

  // coarse grid: receive the info
  npointssentBOTTOM=0;    // all the points sent to ANT proc as BOTTOM BC
  npointssentTOP=0;
  npointssentLEFT=0;
  npointssentRIGHT=0;

  if (grid->getLevel() == 0 )
    {
      for (int i=0; i<nmessageBC; i++)
	{
	  MPI_Recv(info, NInfo,MPI_DOUBLE, MPI_ANY_SOURCE, TagInfo, vct->getCART_COMM_TOTAL(), &status);

	  targetBC[i]= (int)info[0];
	  npointssent[i]= (int)info[1];
	  BCSide[i]= (int)info[2];
	  xfirst_COARSE[i]= info[3];
	  yfirst_COARSE[i]= info[4];
	  // this is actually debug
	  switch(BCSide[i])
	    {
	    case 0:
	      {
		npointssentBOTTOM+= npointssent[i];
		targetBOTTOM++;
		break;
	      }

	    case 1:
	      npointssentTOP+= npointssent[i];
	      targetTOP++;
	      break;
	      
	    case 2:
	      npointssentLEFT+= npointssent[i];
	      targetLEFT++;
	      break;

	    case 3:
	      npointssentRIGHT+= npointssent[i];
	      targetRIGHT++;
	      break;

	    default:
	      {
		cout << "Problems in initWeightBC, exiting..." << endl;
		return -1;
	      }
	    }// end switch  
	  // end this is actually debug
	}
    }
  // end coarse grid receiving the info

  MPI_Barrier(vct->getCART_COMM());
  if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())))
    {
      cout << "Level " << grid->getLevel() << " finished exchanging the info for BC\n";
    }

  // debug, can be commented

  if (grid->getLevel() > 0)
    {
      npointsreceivedBOTTOM=0;
      npointsreceivedTOP=0;
      npointsreceivedLEFT=0;
      npointsreceivedRIGHT=0;
      
      for (int i=0; i< nmessagerecuBC; i++)
	{
	  switch(BCSidecu[i])
	    {
	    case 0:
	      npointsreceivedBOTTOM+=npointsreceived[i];
	      break;
	      
	    case 1:
	      npointsreceivedTOP+=npointsreceived[i];
	      break;
	      
	    case 2:
	      npointsreceivedLEFT+=npointsreceived[i];
	      break;
	      
	    case 3:
	      npointsreceivedRIGHT+=npointsreceived[i];
	      break;
	      
	    default:
	      {
		cout << "Problems in initWeightBC, exiting..." << endl;
		return -1;
	      }
	    }// end switch
	}
    } // end if level

  int ALL_npointssentBOTTOM;
  int ALL_npointssentTOP;
  int ALL_npointssentLEFT;
  int ALL_npointssentRIGHT;

  int ALL_npointsreceivedBOTTOM;
  int ALL_npointsreceivedTOP;
  int ALL_npointsreceivedLEFT;
  int ALL_npointsreceivedRIGHT;
 
  MPI_Allreduce ( &npointssentBOTTOM, &ALL_npointssentBOTTOM, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
  MPI_Allreduce ( &npointssentTOP, &ALL_npointssentTOP, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
  MPI_Allreduce ( &npointssentLEFT, &ALL_npointssentLEFT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
  MPI_Allreduce ( &npointssentRIGHT, &ALL_npointssentRIGHT, 1,  MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

  MPI_Allreduce ( &npointsreceivedBOTTOM, &ALL_npointsreceivedBOTTOM, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
  MPI_Allreduce ( &npointsreceivedTOP, &ALL_npointsreceivedTOP, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
  MPI_Allreduce ( &npointsreceivedLEFT, &ALL_npointsreceivedLEFT, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());
  MPI_Allreduce ( &npointsreceivedRIGHT, &ALL_npointsreceivedRIGHT, 1,MPI_INT, MPI_SUM, vct->getCART_COMM_TOTAL());

  if (ALL_npointssentBOTTOM!= ALL_npointsreceivedBOTTOM || ALL_npointssentTOP!= ALL_npointsreceivedTOP || ALL_npointssentLEFT!= ALL_npointsreceivedLEFT || ALL_npointssentRIGHT!= ALL_npointsreceivedRIGHT )
    {
      if (vct->getCartesian_rank_COMMTOTAL()==0)
	{
	  cout << "Mismatch in points sent and received for interpolation!!!"<<endl;
	  cout << "ALL_npointssentBOTTOM: " << ALL_npointssentBOTTOM << ", ALL_npointsreceivedBOTTOM: " << ALL_npointsreceivedBOTTOM <<endl;
	  cout << "ALL_npointssentTOP: " << ALL_npointssentTOP<< ", ALL_npointsreceivedTOP: " << ALL_npointsreceivedTOP<<endl;
	  cout << "ALL_npointssentLEFT: " << ALL_npointssentLEFT<< ", ALL_npointsreceivedLEFT: " << ALL_npointsreceivedLEFT<<endl;
	  cout << "ALL_npointssentRIGHT: " << ALL_npointssentRIGHT<< ", ALL_npointsreceivedRIGHT: " << ALL_npointsreceivedRIGHT<<endl;
	  return -1;
	}
    }
  else
    {
      if (vct->getCartesian_rank_COMMTOTAL()==0)
        {
	  cout << "Points sent and received for BC checked! Congrats!!!" << endl;
	}
    }

  // end debug, can be commented

  // end exchange the info
  
  // now coarse grid has to built its weights
   
  if (grid->getLevel() < vct->getNgrids()-1)
    {

      //If this grid is considered as coarse by another grid
      finedx = grid->getDX()/col->getRatio();
      finelx = col->getLx()/pow(col->getRatio(),grid->getLevel()+1);
      finedy = grid->getDY()/col->getRatio();
      finely = col->getLy()/pow(col->getRatio(),grid->getLevel()+1);
      finelxplusfinedx = col->getLx()/(double)col->getNxc()*((double)col->getNxc()+1)/(double)pow((double)col->getRatio(),(double)grid->getLevel()+1.);
      Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid
      Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid
      j=0;

      int already_done=0; 
      // Computation of the weights
      for (int r=0; r< nmessageBC; r++)
	{
	  xloc= xfirst_COARSE[r];
	  yloc= yfirst_COARSE[r];

	  if (xloc < grid->getXstart() || xloc > grid->getXend() || yloc < grid->getYstart() || yloc > grid->getYend())
            {
	      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << "Problems in initWeightBC, exiting..." <<endl;
	      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << "xloc: " <<xloc <<", yloc: " <<yloc <<endl;
	      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<"grid->getXstart() " <<grid->getXstart() <<"grid->getXend() " <<grid->getXend() <<endl;
	      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<"grid->getYstart() " <<grid->getYstart() <<"grid->getYend() " <<grid->getYend() <<endl;
	      return -1;
            }
	  switch (BCSide[r])
	    {
	    case 0:
	      {
		// BOTTOM
		xfirst= xloc ;
		for (int i=0; i< npointssent[r]; i++)
		  {
		    xshift = xfirst + i * finedx ;
		    //xloc = xshift + Ox ;
		    xloc = xshift;  
		    if (xloc < grid->getXstart() || xloc > grid->getXend() )
		      {
			cout <<"R" << vct->getCartesian_rank_COMMTOTAL()<< " xloc " << xloc<< " Problems in initWeightBC, exiting..." <<endl;
			cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<"grid->getXstart() " <<grid->getXstart() <<"grid->getXend() " <<grid->getXend() <<endl;
			return -1;
		      }

		    //Weights used for interpolation on nodes (for E and B)
		    ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
		    iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
		    ixsent[already_done] = ix;
		    iysent[already_done] = iy;
		    weightBC[already_done][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
		    weightBC[already_done][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
		    weightBC[already_done][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
		    weightBC[already_done][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
		    already_done++;
		  }
		break;
	      }// end case 0
	    case 1:
	      {
		// TOP
		xfirst= xloc;
		for (int i=0; i< npointssent[r]; i++)
                  {
		    xshift = xfirst + i * finedx ;
		    //xloc = xshift + Ox ;
		    xloc = xshift ;
		    if (xloc < grid->getXstart() || xloc > grid->getXend() )
                      {
                        cout <<"R" << vct->getCartesian_rank_COMMTOTAL()<< " xloc " << xloc<< " Problems in initWeightBC, exiting..." <<endl;
			cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<"grid->getXstart() " <<grid->getXstart() <<"grid->getXend() " <<grid->getXend() <<endl;
                        return -1;
                      }

		    ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
		    iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
		    ixsent[already_done] = ix;
		    iysent[already_done] = iy;
		    weightBC[already_done][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
		    weightBC[already_done][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
		    weightBC[already_done][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
		    weightBC[already_done][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
		    already_done++;
		  }
		break;
	      } // end case 1
	    case 2:
	      {
		// LEFT
		xfirst= yloc;
		for (int i=0; i< npointssent[r]; i++)
                  {
		    yshift = xfirst + i * finedy ;
		    //yloc = yshift + Oy ;
		    yloc = yshift ; 
		    if (yloc < grid->getYstart() || yloc > grid->getYend())
		      {
			cout <<"R" << vct->getCartesian_rank_COMMTOTAL()<< " yloc " << yloc<< " Problems in initWeightBC, exiting..." <<endl;
                        cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<"grid->getYstart() " <<grid->getYstart() <<"grid->getYend() " <<grid->getYend() <<endl;
			return -1;
		      }
		    
		    ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
		    iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
		    ixsent[already_done] = ix;
		    iysent[already_done] = iy;
		    weightBC[already_done][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
		    weightBC[already_done][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
		    weightBC[already_done][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
		    weightBC[already_done][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
		    already_done++;
		  }
		break;
	      } // end case 2
	    case 3:
	      {
		// RIGHT
		xfirst = yloc;
		for (int i=0; i< npointssent[r]; i++)
                  {
		    yshift = xfirst + i * finedy ;
		    //yloc = yshift + Oy ;
		    yloc = yshift  ;
		    if (yloc < grid->getYstart() || yloc > grid->getYend())
		    {
		      cout <<"R" << vct->getCartesian_rank_COMMTOTAL()<< " yloc " << yloc<< " Problems in initWeightBC, exiting..." <<endl;
		      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<"grid->getYstart() " <<grid->getYstart() <<"grid->getYend() " <<grid->getYend() <<endl;
		      return -1;
		    }
		    ix = 2 +  int(floor((xloc-grid->getXstart())/grid->getDX()));
		    iy = 2 +  int(floor((yloc-grid->getYstart())/grid->getDY()));
		    ixsent[already_done] = ix;
		    iysent[already_done] = iy;
		    weightBC[already_done][3][0] = ((xloc - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((yloc - grid->getYN(ix-1,iy-1,0))/grid->getDY()); // weight +:+
		    weightBC[already_done][2][0] = ((xloc - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - yloc)/grid->getDY()); //weight +:-
		    weightBC[already_done][1][0] = ((grid->getXN(ix,iy-1,0) - xloc)/grid->getDX())*((yloc - grid->getYN(ix,iy-1,0))/grid->getDY()); // weight -:+
		    weightBC[already_done][0][0] = ((grid->getXN(ix,iy,0) - xloc)/grid->getDX())*((grid->getYN(ix,iy,0) - yloc)/grid->getDY()); // weight -:-
		    already_done++;
		  }
		break;
	      } // end case 3
	      
	    } // end switch
	  
	}  // end for on targets

      // end coarse grid builds its weights
      
      // debug, possible to remove
      int total_sent=0;
      for (int i=0; i< nmessageBC; i++)
	{
	  total_sent+= npointssent[i];
	}
      if (already_done != total_sent)
	{
	  cout << "R" << vct->getCartesian_rank_COMMTOTAL() << ": mess in initWeightBC: exiting..." <<endl;
	  return -1;
	}
      // end debug, possible to remvoe
      
    }// end coarse grod
  
  // end send/receive


  MPI_Barrier(vct->getCART_COMM());
  if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN())))
    {
      cout << "Level " << grid->getLevel() << " finished initWeightBC\n";
    }  

  MPI_Barrier(vct->getCART_COMM_TOTAL());
  if (vct->getCartesian_rank_COMMTOTAL()==0)
    {
      cout << "everybody finished initWeightBC" << endl;
    }

    
  delete []buffer;
  delete []info;
  return 1;
}

//Initialize weight for projection of refined fields between grids
inline int EMfields::initWeightProj(VirtualTopology *vct, Grid *grid, CollectiveIO* col){

  if (vct->getNgrids()<2)
    return 1;


  MPI_Status status;
  int X,Y,Xlocal,Ylocal,i,j,k,ix,iy,nprocfirst,nproclast,nybloc,ierr;
  double coarsedx,coarsedy,coarselx,coarsely,finelx,finely,Ox,Oy,finedx,finedy,correctionx,correctiony;
  double xloc, yloc,xfirst,xlast,yfirst,ylast,xlocalfirst,xlocallast,ylocalfirst,ylocallast,xstop,ystop;
  int start, step;

  int ID= vct->getCartesian_rank_COMMTOTAL();

  int *ProjCoarseGrid= new int[vct->getXLEN()*vct->getYLEN()];
  int *TOTALProjCoarseGrid= new int[vct->getXLEN()*vct->getYLEN()];
  int coordXprojMM, coordYprojMM; 
  int coordXprojMP, coordYprojMP;
  int coordXprojPM, coordYprojPM;
  int coordXprojPP, coordYprojPP;
  // to see if eliminate cells in the projection
  int xleft=0.;
  int xright=0;
  int yleft=0;
  int yright=0;
  /*MPI_Barrier(vct->getCART_COMM_TOTAL());
  if (vct->getCartesian_rank_COMMTOTAL()==0)
  cout << "AFter the Barrier at the beginnign in initWeightProj" <<endl;*/

  for (int i=0; i< vct->getXLEN()*vct->getYLEN(); i++)
    {
      ProjCoarseGrid[i]=0;
    }

  /*MPI_Barrier(vct->getCART_COMM_TOTAL());
  if (vct->getCartesian_rank_COMMTOTAL()==0)
  cout << "AFter the Barrier at the beginnign in initWeightProj" <<endl;*/

  if ( grid->getLevel() > 0) {                             //If this grid is considered as fine by another grid
    coarsedx = grid->getDX()*col->getRatio();
    coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel()-1);
    coarsedy = grid->getDY()*col->getRatio();
    coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel()-1);
    Ox = grid->getOx(grid->getLevel()); //Origin x of finer grid
    Oy = grid->getOy(grid->getLevel()); //Origin y of finer grid

    xfirst = floor((grid->getXstart()+Ox)/coarsedx)*coarsedx;//Coordinate of the first point in x direction of the coarse grid influenced by the projection of this proc expressed in the coarse grid frame
    X= floor((xfirst)/coarsedx)/(grid->getNXC()-2); // X coordinate of the proc of the coarse grid intersecting the origin of this proc   
    yfirst = floor((grid->getYstart()+Oy)/coarsedy)*coarsedy;
    Y= floor((yfirst)/coarsedy)/(grid->getNYC()-2); // X coordinate of the proc of the coarse grid intersecting the origin of this proc   
    targetProj[0] = X*vct->getYLEN()+Y+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
    targetProj[1] = (X+1)*vct->getYLEN()+Y+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
    targetProj[2] = X*vct->getYLEN()+Y+1+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
    targetProj[3] = (X+1)*vct->getYLEN()+Y+1+vct->getXLEN()*vct->getYLEN()*(grid->getLevel()-1);
    if (vct->getXright_neighbor()==MPI_PROC_NULL){
      xlast = ceil((grid->getXend()+Ox)/coarsedx)*coarsedx;
      lastindicex = grid->getNXN()-2; 
    } else {
      xlast = ceil((grid->getXend()+Ox)/coarsedx-1./col->getRatio())*coarsedx;
      lastindicex = grid->getNXN()-3; 
    }
    if (vct->getYright_neighbor()==MPI_PROC_NULL){
      ylast = ceil((grid->getYend()+Oy)/coarsedy)*coarsedy;
      lastindicey = grid->getNYN()-2; 
    } else {
      ylast = ceil((grid->getYend()+Oy)/coarsedy-1./col->getRatio())*coarsedy;
      lastindicey = grid->getNYN()-3; 
    }
    xstop = (X+1)*(grid->getNXC()-2)*coarsedx; //End of domain of processor whith coordinate X
    ystop = (Y+1)*(grid->getNYC()-2)*coarsedy; //End of domain of processor with coordinate Y
    //nxmsend is the number of points sent in the x direction to the left part of x=xstop, and nxp to the right part of x=xstop

    if (xstop < xlast){
      nmessageProj=nmessageProj*2;
      nxpsend = floor((xlast-xstop)/coarsedx+0.5) ; 
      nxmsend = floor((xlast-xfirst)/coarsedx+0.5)+1-nxpsend; // the common node is assigned to the first core; the node on the right core is fixed at the end of receiveProj
    } else {
      nxmsend = floor((xlast-xfirst)/coarsedx+0.5)+1;
      nxpsend = 0;
    }
    //nym is the number of points sent in the y direction to the lower part of y=ystop, and nxp to the upper part of y=ystop
    if (ystop < ylast){
      nmessageProj=nmessageProj*2;
      nypsend = floor((ylast-ystop)/coarsedy+0.5) ; 
      nymsend = floor((ylast-yfirst)/coarsedy+0.5)+1-nypsend; // the common node is assigned to the first core; the node on the right core is fixed at the end of receiveProj
    } else {
      nymsend = floor((ylast-yfirst)/coarsedy+0.5) + 1;
      nypsend = 0 ;
    }

    for (ix =0; ix<nxmsend+nxpsend; ix++){
      for (iy =0; iy<nymsend+nypsend; iy++){
	normalizeProj[ix][iy][0]=0.;
      }
    }
    // xfirst, x last are physical coords, I need the indexes
    double CoarseXLen= coarselx/ vct->getXLEN();
    double CoarseYLen= coarsely/ vct->getYLEN();

    coordXprojMM= floor((xfirst-X*CoarseXLen)/coarsedx+0.5)+1; // +1 to take into account the ghost node; without the +1 it's a disaster, keep it here     
    coordYprojMM= floor((yfirst-Y*CoarseYLen)/coarsedy+0.5)+1;
    coordXprojMP= floor((xfirst-X*CoarseXLen)/coarsedx+0.5)+1;
    coordYprojMP=2; 
    coordXprojPM=2; 
    coordYprojMP=2; 
    coordXprojPM=2; 
    coordYprojPM= floor((yfirst-Y*CoarseYLen)/coarsedy+0.5)+1;
    coordXprojPP=2; 
    coordYprojPP=2; 

    //cout << "R" << vct->getCartesian_rank_COMMTOTAL() << ", coordXprojMM "<< coordXprojMM << ", coordYprojMM "<< coordYprojMM << ", coordXprojMP "<< coordXprojMP << ", coordYprojMP "<< coordYprojMP << ", coordXprojPM "<< coordXprojPM << ", coordYprojPM "<< coordYprojPM << ", coordXprojPP "<< coordXprojPP << ", coordYprojPP "<< coordYprojPP <<endl;
    
    for (i =0;i<lastindicex;i++) {
      xloc = grid->getXstart() + i * grid->getDX() + Ox -xfirst;
      ixsentProj[i] = int(floor(xloc/coarsedx));
      for (j =0;j<lastindicey;j++) {
	yloc = grid->getYstart() + j * grid->getDY() + Oy -yfirst;
	iysentProj[j] = int(floor(yloc/coarsedy));

	weightProj[i][j][3][0] = ((xloc - ixsentProj[i]*coarsedx)/coarsedx)*((yloc - iysentProj[j]*coarsedy)/coarsedy); // weight +:+
	weightProj[i][j][2][0] = ((xloc - ixsentProj[i]*coarsedx)/coarsedx)*(((iysentProj[j]+1)*coarsedy - yloc)/coarsedy); //weight +:-
	weightProj[i][j][1][0] = (((ixsentProj[i]+1)*coarsedx - xloc)/coarsedx)*((yloc - iysentProj[j]*coarsedy)/coarsedy); // weight -:+
	weightProj[i][j][0][0] = (((ixsentProj[i]+1)*coarsedx - xloc)/coarsedx)*(((iysentProj[j]+1)*coarsedy - yloc)/coarsedy); // weight -:-

	normalizeProj[ixsentProj[i]+1][iysentProj[j]+1][0] += weightProj[i][j][3][0];
	normalizeProj[ixsentProj[i]+1][iysentProj[j]][0] += weightProj[i][j][2][0];
	normalizeProj[ixsentProj[i]][iysentProj[j]+1][0] += weightProj[i][j][1][0];
	normalizeProj[ixsentProj[i]][iysentProj[j]][0] += weightProj[i][j][0][0];
      }
    }
    /*if (vct->getCartesian_rank_COMMTOTAL()== 24){
        for (i=0;i<nxmsend;i++){
            for (j=0;j<nymsend;j++){
                cout<< normalizeProj[i][j][0]<< " ";
            }
            cout << endl;
        }
	}*/
    //cout<< vct->getCartesian_rank_COMMTOTAL() << " sends to target0 = " << targetProj[0]<<" nmessageProj= "<<nmessageProj<<" nxmsend= "<<nxmsend<<" nxpsend= "<<nxpsend<<" nymsend= "<<nymsend<<" nypsend= "<<nypsend<< "xfirst = "<<xfirst<<"xlast = "<<xlast<<" yfirst= "<<yfirst<<" ylast = "<<ylast<<" xstop= "<<xstop <<" lastindicey "<< lastindicey << "lastindicex  "<<lastindicex <<endl;
    
    ProjCoarseGrid[targetProj[0]]++;   

    if(nxpsend > 0){
      ProjCoarseGrid[targetProj[1]]++;
    }

    if(nypsend > 0){
      ProjCoarseGrid[targetProj[2]]++;
    }
    if(nypsend > 0 && nxpsend > 0){
      ProjCoarseGrid[targetProj[3]]++;
    }


  }
  //MPI_Barrier(vct->getCART_COMM_TOTAL());
  //if (vct->getCartesian_rank_COMMTOTAL()==0)
  //  cout << "AFter the Barrier in initWeightProj" <<endl;

  if (grid->getLevel()==0) // init the normalizerecvProj, while the refined grid does something else
    {
      for (ix =0; ix<nxn; ix++){
        for (iy =0; iy<nyn; iy++){
	  normalizerecvProj[ix][iy][0]=0.; 
        }
      }
    }

  MPI_Allreduce(ProjCoarseGrid,TOTALProjCoarseGrid,vct->getXLEN()*vct->getYLEN(),MPI_INT,MPI_SUM,vct->getCART_COMM_TOTAL());
  //if (vct->getCartesian_rank_COMMTOTAL()==0)
  //  {
  //    cout << "initWeightProj: after MPI_Reduce\n";
  //    for (int i=0; i<vct->getXLEN()*vct->getYLEN(); i++)
  //    {
  //	  cout << "Messages to R" << i << ": " << TOTALProjCoarseGrid[i] <<endl;
  // 	}
  //  }

  // refined grid builds the weights and  sends the messages

  int TagProj=5;
  //double prova=0.0;
  if (grid->getLevel()==1)
    {

      // send to the coarse grid to eliminate cells in the projection
      if (vct->getXleft_neighbor()== MPI_PROC_NULL)
	xleft=1;
      if (vct->getXright_neighbor()== MPI_PROC_NULL)
	xright=1;
      if (vct->getYleft_neighbor()== MPI_PROC_NULL)
	yleft=1;
      if (vct->getYright_neighbor()== MPI_PROC_NULL)
	yright=1;

      //by now, messages are just one int
      start = 0;
      for (i=0;i<nxmsend*nymsend;i++){
        bufferProjsend[start] = normalizeProj[i/nymsend][i%nymsend][0];
        start++;
	//   cout << bufferProjsend[i]<<" ";
      }
      // build the message - minus
      INFObufferProjsend[0]= (double)vct->getCartesian_rank_COMMTOTAL(); // rank of the refined grid sender
      INFObufferProjsend[1]= (double)nxmsend; //number of points in the x dir
      INFObufferProjsend[2]= (double)nymsend; //number of points in the y dir
      INFObufferProjsend[3]= (double)(nxmsend*nymsend); // total number of points
      INFObufferProjsend[4]= (double)coordXprojMM; // x index of the first point
      INFObufferProjsend[5]= (double)coordYprojMM; // y index of the first point 
      //cout << "R" << vct->getCartesian_rank_COMMTOTAL() << " targetProj[0]  " << targetProj[0] << " coordXprojMM " << coordXprojMM <<" coordYprojMM " << coordYprojMM << endl;
      INFObufferProjsend[6]= (double)xleft;
      INFObufferProjsend[7]= (double)xright;
      INFObufferProjsend[8]= (double)yleft;
      INFObufferProjsend[9]= (double)yright;
      memcpy (INFObufferProjsend+10, bufferProjsend, nxmsend*nymsend*sizeof(double) );
      // end build the message - minus
      //ierr = MPI_Send(bufferProjsend,nxmsend*nymsend,MPI_DOUBLE,targetProj[0],1,MPI_COMM_WORLD); //original, Arnaud
      //ierr = MPI_Send(&prova,1,MPI_DOUBLE,targetProj[0],TagProj,vct->getCART_COMM_TOTAL()); // to test comm
      ierr = MPI_Send(INFObufferProjsend,nxmsend*nymsend+10,MPI_DOUBLE,targetProj[0],TagProj,vct->getCART_COMM_TOTAL());

      if(nxpsend > 0){
        step=start;
        for (i=0;i<nxpsend*nymsend;i++){
	  bufferProjsend[start] = normalizeProj[nxmsend+i/nymsend][i%nymsend][0];
	  start++;
        }
	// build the message - x plus y minus
	INFObufferProjsend[0]= (double)vct->getCartesian_rank_COMMTOTAL(); // rank of the refined grid sender
	INFObufferProjsend[1]= (double)nxpsend; //number of points in the x dir
	INFObufferProjsend[2]= (double)nymsend; //number of points in the y dir 
	INFObufferProjsend[3]= (double)(nxpsend*nymsend); // total number of points 
	INFObufferProjsend[4]= (double)coordXprojPM; // x index of the first point 
	INFObufferProjsend[5]= (double)coordYprojPM; // y index of the first point   
	//cout << "R" << vct->getCartesian_rank_COMMTOTAL() << " targetProj[1]  " << targetProj[1] << " coordXprojPM " << coordXprojPM <<" coordYprojPM " << coordYprojPM << endl;
	INFObufferProjsend[6]= (double)xleft;
	INFObufferProjsend[7]= (double)xright;
	INFObufferProjsend[8]= (double)yleft;
	INFObufferProjsend[9]= (double)yright;
	memcpy (INFObufferProjsend+10, bufferProjsend+step, nxpsend*nymsend*sizeof(double) );
	// end build the message - x plus y minus
        //ierr = MPI_Send(bufferProjsend+step,nxpsend*nymsend,MPI_DOUBLE,targetProj[1],1,MPI_COMM_WORLD); //original Arnaud
	//ierr = MPI_Send(&prova,1,MPI_DOUBLE,targetProj[1],TagProj,vct->getCART_COMM_TOTAL()); // to test comm
	ierr = MPI_Send(INFObufferProjsend,nxpsend*nymsend+10,MPI_DOUBLE,targetProj[1],TagProj,vct->getCART_COMM_TOTAL()); 
	
      }

      if(nypsend > 0){
        step=start;
        for (i=0;i<nxmsend*nypsend;i++){
	  bufferProjsend[start] = normalizeProj[i/nypsend][nymsend+i%nypsend][0];
	  start++;
        }
	// build the message - x minus y plus
	INFObufferProjsend[0]= (double)vct->getCartesian_rank_COMMTOTAL(); // rank of the refined grid sender
        INFObufferProjsend[1]= (double)nxmsend; //number of points in the x dir  
        INFObufferProjsend[2]= (double)nypsend; //number of points in the y dir   
        INFObufferProjsend[3]= (double)(nxmsend*nypsend); // total number of points     
        INFObufferProjsend[4]= (double)coordXprojMP; // x index of the first point     
        INFObufferProjsend[5]= (double)coordYprojMP; // y index of the first point        
	//cout << "R" << vct->getCartesian_rank_COMMTOTAL() << " targetProj[2]  " << targetProj[2] << " coordXprojMP " << coordXprojMP <<" coordYprojMP " << coordYprojMP << endl;
	INFObufferProjsend[6]= (double)xleft;
	INFObufferProjsend[7]= (double)xright;
	INFObufferProjsend[8]= (double)yleft;
	INFObufferProjsend[9]= (double)yright;
        memcpy (INFObufferProjsend+10, bufferProjsend+step, nxmsend*nypsend*sizeof(double) );
	// end build the message - x minus y plus 
        //ierr = MPI_Send(bufferProjsend+step,nxmsend*nypsend,MPI_DOUBLE,targetProj[2],1,MPI_COMM_WORLD); // original from Arnaud
	//ierr = MPI_Send(&prova,1,MPI_DOUBLE,targetProj[2],TagProj,vct->getCART_COMM_TOTAL()); // to test comm
	ierr = MPI_Send(INFObufferProjsend,nxmsend*nypsend+10,MPI_DOUBLE,targetProj[2],TagProj,vct->getCART_COMM_TOTAL());
      }
      if(nypsend > 0 && nxpsend > 0){
        step=start;
        for (i=0;i<nxpsend*nypsend;i++){
	  bufferProjsend[start] = normalizeProj[nxmsend+i/nypsend][nymsend+i%nypsend][0];
	  start++;
        }
	// build the message - xpplus y plus
	INFObufferProjsend[0]= (double)vct->getCartesian_rank_COMMTOTAL(); // rank of the refined grid sender
        INFObufferProjsend[1]= (double)nxpsend; //number of points in the x dir  
        INFObufferProjsend[2]= (double)nypsend; //number of points in the y dir         
        INFObufferProjsend[3]= (double)(nxpsend*nypsend); // total number of points        
        INFObufferProjsend[4]= (double)coordXprojPP; // x index of the first point     
        INFObufferProjsend[5]= (double)coordYprojPP; // y index of the first point          
	//cout << "R" << vct->getCartesian_rank_COMMTOTAL() << " targetProj[3]  " << targetProj[3] << " coordXprojPP " << coordXprojPP <<" coordYprojPP " << coordYprojPP << endl;
	INFObufferProjsend[6]= (double)xleft;
	INFObufferProjsend[7]= (double)xright;
	INFObufferProjsend[8]= (double)yleft;
	INFObufferProjsend[9]= (double)yright;
        memcpy (INFObufferProjsend+10, bufferProjsend+step, nxpsend*nypsend*sizeof(double) );
	// end build - x plus y plus
        //ierr = MPI_Send(bufferProjsend+step,nxpsend*nypsend,MPI_DOUBLE,targetProj[3],1,MPI_COMM_WORLD); // original from Arnaud
	//ierr = MPI_Send(&prova,1,MPI_DOUBLE,targetProj[3],TagProj,vct->getCART_COMM_TOTAL());
	ierr = MPI_Send(INFObufferProjsend,nxpsend*nypsend+10,MPI_DOUBLE,targetProj[3],TagProj,vct->getCART_COMM_TOTAL());
      }
    }

  if (grid->getLevel()==0)
    {
      nmessagerecvProj=TOTALProjCoarseGrid[vct->getCartesian_rank_COMMTOTAL()];
      Ox = grid->getOx(grid->getLevel()+1); //Origin x of finer grid
      Oy = grid->getOy(grid->getLevel()+1); //Origin y of finer grid
      coarselx = col->getLx()/pow(col->getRatio(),grid->getLevel());
      coarsely = col->getLy()/pow(col->getRatio(),grid->getLevel());
      finelx = col->getLx()/pow(col->getRatio(),grid->getLevel()+1);
      finely = col->getLy()/pow(col->getRatio(),grid->getLevel()+1);
      finedx = grid->getDX()/col->getRatio();
      finedy = grid->getDY()/col->getRatio();
      
      bool skip_startX= false;
      bool skip_endX=false;
      bool skip_startY=false;
      bool skip_endY=false;
      int min_startX= nxn+2;
      int min_startY= nyn+2;
      int max_endX= -2;
      int max_endY= -2;

      for (i=0; i< nmessagerecvProj; i++)
	{
	  //ierr = MPI_Recv(&prova,1,MPI_DOUBLE,MPI_ANY_SOURCE,TagProj,vct->getCART_COMM_TOTAL(), &status); // to test comm
	  //ierr = MPI_Recv(bufferProj,npointsreceivedProj[i],MPI_DOUBLE,fromProj[i],1,MPI_COMM_WORLD, &status); // original from Arnaud
	  //ierr = MPI_Recv(INFObufferProj,6+nxnproj*nynproj*6,MPI_DOUBLE,MPI_ANY_SOURCE,TagProj,vct->getCART_COMM_TOTAL(), &status); //wrong
	  ierr = MPI_Recv(INFObufferProj,10+nxnproj*nynproj*6,MPI_DOUBLE,MPI_ANY_SOURCE,TagProj,vct->getCART_COMM_TOTAL(), &status);

	  fromProj[i]=(int)INFObufferProj[0];   //sending proc
	  nxrecvProj[i]=(int)INFObufferProj[1]; // points in the x dir
	  nyrecvProj[i]=(int)INFObufferProj[2]; // points in the y dir
	  npointsreceivedProj[i]= (int)INFObufferProj[3];// total number of points
	  ixrecvfirstProj[i] = (int)INFObufferProj[4];// x coord first point
	  iyrecvfirstProj[i] = (int)INFObufferProj[5];// x coord first point
	  // for global
	  min_startX= min(min_startX, ixrecvfirstProj[i]);
	  max_endX= max(max_endX, ixrecvfirstProj[i]+ nxrecvProj[i]-1);
	  min_startY= min(min_startY, iyrecvfirstProj[i]);
	  max_endY= max(max_endY, iyrecvfirstProj[i]+ nyrecvProj[i]-1);
	  //cout << "R" << vct->getCartesian_rank_COMMTOTAL() << " num message " << i <<" of " << nmessagerecvProj << " ixrecvfirstProj[i]  " <<  ixrecvfirstProj[i] << " iyrecvfirstProj[i] "<< iyrecvfirstProj[i] << " nxrecvProj[i] " << nxrecvProj[i] << " nyrecvProj[i] " << nyrecvProj[i] << " npointsreceivedProj[i] " << npointsreceivedProj[i] << endl;

	  if (fabs(INFObufferProj[6] -1.)<0.1  ) // to avoid mess with the conversion int / double
	    {
	      skip_startX= true;
	    }
	  if (fabs(INFObufferProj[7] -1.)<0.1  )// to avoid mess with the conversion int / double 
            {
	      skip_endX= true;
	    }
	  if (fabs(INFObufferProj[8] -1.)<0.1  )// to avoid mess with the conversion int / double    
	    {
	      skip_startY= true;
	    }
          if (fabs(INFObufferProj[9] -1.)<0.1  )// to avoid mess with the conversion int / double 
	    {
	      skip_endY= true;
	    }
	  memcpy (bufferProj, INFObufferProj+10, npointsreceivedProj[i]*sizeof(double)); // bufferProj
	  // unpackage the message
	  
	  // filling normalizerecvProj
	  for (j=0;j<npointsreceivedProj[i];j++){
	    ix = ixrecvfirstProj[i] + j/nyrecvProj[i];
	    iy = iyrecvfirstProj[i] + j%nyrecvProj[i];
	    normalizerecvProj[ix][iy][0] += bufferProj[j];
	    if (vct->getCartesian_rank_COMMTOTAL() == 27)
	      {
		cout << "test27 " << "ix " << ix << " iy " << iy << " normalizerecvProj[ix][iy][0] " << normalizerecvProj[ix][iy][0] <<endl;
	      }
	    if (vct->getCartesian_rank_COMMTOTAL() == 35)
              {
                cout << "test35 " << "ix " << ix << " iy " << iy << " normalizerecvProj[ix][iy][0] " << normalizerecvProj[ix][iy][0] <<endl;
              }
	    }

	  

	}// end cycle on messages to receive

      // for global
      int Skip=1;
      if (nmessagerecvProj>0)
	{
	  ixrecvfirstProjglobal= min_startX;
	  if (skip_startX)
	    ixrecvfirstProjglobal+= Skip;

	  iyrecvfirstProjglobal= min_startY;
          if (skip_startY)
            iyrecvfirstProjglobal+= Skip;

	  ixrecvlastProjglobal= max_endX;
          if (skip_endX)
            ixrecvlastProjglobal-= Skip;

	  iyrecvlastProjglobal= max_endY;
          if (skip_endY)
            iyrecvlastProjglobal-= Skip;

	  //cout <<"R" << vct->getCartesian_rank_COMMTOTAL()<< " ixrecvfirstProjglobal " <<ixrecvfirstProjglobal << " iyrecvfirstProjglobal " <<iyrecvfirstProjglobal << " ixrecvlastProjglobal " <<ixrecvlastProjglobal << " iyrecvlastProjglobal " <<iyrecvlastProjglobal <<endl;
	}


      // still the coarse grid internal exchanges: the first/ last active node is shared among two processors:
      // sum the contributions from both
      // Completing normalizerecvProj with the normalization gathered by the other coarse procs
      if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //Send a message to proc  coorinateX-1
        for (i=1;i<nyn-1;i++){
	  bufferProjsend[i-1] = normalizerecvProj[1][i][0];
	}
        if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
	  //also receive a message
	  //cout << vct->getCartesian_rank_COMMTOTAL() << "SENd RECEIVE" <<endl;
	  ierr = MPI_Sendrecv(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
	  //send only
	  //cout << vct->getCartesian_rank_COMMTOTAL() << "SENd only" <<endl;
	  ierr = MPI_Send(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM());
        }
      } else {
        //Do not send ...
	if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
          //... but receive only
	  //cout  << "RECEIVE only" <<endl;
	  ierr = MPI_Recv(bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM(), &status);
	}
      }
      if (Ox <=grid->getXend() && Ox+finelx >= grid->getXend()) {
        //If a message from proc coorinateX+1 has been received
        for (i=1;i<nyn-1;i++){
	  normalizerecvProj[nxn-2][i][0] += bufferProjrecv[i-1];
        }
	//Send a message to proc coordinateX+1
        for(i=1;i<nyn-1;i++){
	  bufferProjsend[i-1] = normalizerecvProj[nxn-2][i][0];
        }
        if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
	  //also receive a message
	  ierr = MPI_Sendrecv(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
	  //send only
	  //cout << vct->getCartesian_rank_COMMTOTAL() <<" send to "<< vct->getXright_neighbor() << endl;
	  ierr = MPI_Send(bufferProjsend,nyn-2,MPI_DOUBLE,vct->getXright_neighbor(),1,vct->getCART_COMM());
        }
    
      } else {
	//Do not send ...
	if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
	  //...receive only
	  //cout << vct->getCartesian_rank_COMMTOTAL() <<" receive from "<< vct->getXleft_neighbor() << endl;
	  ierr = MPI_Recv(bufferProjrecv,nyn-2,MPI_DOUBLE,vct->getXleft_neighbor(),1,vct->getCART_COMM(), &status);
	}
      }
      if (Ox <=grid->getXstart() && Ox+finelx >= grid->getXstart()) {
        //A message from proc coordinateX-1 has been received
        for (i=1;i<nyn-1;i++){
	  normalizerecvProj[1][i][0] =  bufferProjrecv[i-1];
        }
	/*if (vct->getCartesian_rank_COMMTOTAL() == 35)
          {
	    for (i=1;i<nyn-1;i++){
	      cout << "R35, i and j " << 1 <<"  and " << i <<" normalizerecvProj[1][j][0] " << normalizerecvProj[1][i][0] <<endl;
	    }
	    }*/
      }

      // In the Y direction now
      if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //Send a message to proc  coorinateY-1
        for (i=1;i<nxn-1;i++){
	  bufferProjsend[i-1] = normalizerecvProj[i][1][0];
	}
        if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
	  //also receive a message
	  ierr = MPI_Sendrecv(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
	  //send only
	  ierr = MPI_Send(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM());
        }
      } else {
        //Do not send ...
	if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
          //... but receive only
	  ierr = MPI_Recv(bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM(), &status);
	}
      }
      if (Oy <=grid->getYend() && Oy+finely >= grid->getYend()) {
        //If a message from proc coorinateY+1 has been received
        for (i=1;i<nxn-1;i++){
	  normalizerecvProj[i][nyn-2][0] += bufferProjrecv[i-1];
        }
	//Send a message to proc coordinateY+1
        for(i=1;i<nxn-1;i++){
	  bufferProjsend[i-1] =  normalizerecvProj[i][nyn-2][0];
        }
        if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
	  //also receive a message
	  ierr = MPI_Sendrecv(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(), &status);
        } else {
	  //send only
	  cout << vct->getCartesian_rank_COMMTOTAL() <<" send to "<< vct->getYright_neighbor() << endl;
	  ierr = MPI_Send(bufferProjsend,nxn-2,MPI_DOUBLE,vct->getYright_neighbor(),1,vct->getCART_COMM());
        }
	 
      } else {
	//Do not send ...
	if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
	  //...receive only
	  cout << vct->getCartesian_rank_COMMTOTAL() <<" receive from "<< vct->getYleft_neighbor() << endl;
	  ierr = MPI_Recv(bufferProjrecv,nxn-2,MPI_DOUBLE,vct->getYleft_neighbor(),1,vct->getCART_COMM(), &status);
	}
      }
      if (Oy <=grid->getYstart() && Oy+finely >= grid->getYstart()) {
        //A message from proc coordinateX-1 has been received
        for (i=1;i<nxn-1;i++){
	  normalizerecvProj[i][1][0] =  bufferProjrecv[i-1];
        }
      }


      /*cout <<"R" << vct->getCartesian_rank_COMMTOTAL() << ", normalizerecvProj: " << endl;
      for (int i=0; i< nxn; i++)
	{
	  cout << "R" << vct->getCartesian_rank_COMMTOTAL() ;
	  for (int j=0; j< nyn; j++)
	    {
	      cout << " i, j : " << i << " " << j << ": " << normalizerecvProj[i][j][0] <<endl;
	    }
	  cout <<endl;
	  }*/

    }// end coarse grid
    

  MPI_Barrier(vct->getCART_COMM_TOTAL());
  if (vct->getCartesian_rank_COMMTOTAL()==0)
    cout << "AFter the final Barrier in initWeightProj" <<endl;

  delete []ProjCoarseGrid;
  delete []TOTALProjCoarseGrid;
  
  return 1;
}

#endif

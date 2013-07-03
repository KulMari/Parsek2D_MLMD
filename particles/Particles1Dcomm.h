/*******************************************************************************************
Particles1Dcommcomm.h  -  Class for particles of the same species, in a 2D space and 3component velocity with communications methods
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef Part1DCOMM_H
#define Part1DCOMM_H

#include "Particles.h"
/**
* 
* Abstract class for particles of the same species, in a 2D space and 3component velocity with communications methods
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Particles1Dcomm : public Particles {
  public:
     /** constructor */
     Particles1Dcomm();
     /** destructor */
     ~Particles1Dcomm();
     /** allocate particles */
     void allocate(int species, CollectiveIO* col,VirtualTopology* vct, Grid* grid);
    
     /** interpolation method GRID->PARTICLE order 1: CIC */
     void interpP2G(Field* EMf, Grid *grid, VirtualTopology* vct);
     /** method for communicating exiting particles to X-RIGHT, X-LEFT, Y-RIGHT, Y-LEFT, Z-RIGHT, Z-LEFT processes */
     void communicate(VirtualTopology* ptVCT);
     /** put a particle exiting to X-LEFT in the bufferXLEFT for communication and check if you're sending the particle to the right subdomain*/
     void bufferXleft(double *b_, int np, VirtualTopology* vct);
     /** put a particle exiting to X-RIGHT in the bufferXRIGHT for communication and check if you're sending the particle to the right subdomain*/
     void bufferXright(double *b_, int np, VirtualTopology* vct);
     /** put a particle exiting to Y-LEFT in the bufferYLEFT for communication and check if you're sending the particle to the right subdomain*/
     void bufferYleft(double *b_, int np, VirtualTopology* vct);
     /** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication and check if you're sending the particle to the right subdomain*/
     void bufferYright(double *b_, int np, VirtualTopology* vct);
     /** Delete the a particle from a list(array) and pack the list(array) */
     void del_pack(int np, int *nplast);
     /** method to debuild the buffer received */
     void unbuffer(double *b_);
     /** resize the receiving buffer */
     void resize_buffers(int new_buffer_size);
     /** a method to compute how many particles are not in the right domain */
     int isMessagingDone(VirtualTopology* ptVCT);
     /** calculate the maximum number exiting from this domain */
     int maxNpExiting();
     /** calculate the weights given the position of particles */
  //   void calculateWeights(double*** weight, double xp, double yp, double zp,int ix, int iy, int iz, Grid* grid);
     /** get X-position array for all the particles */
     double* getXall() const;
     /** get Y-position array for all the particles */
     double* getYall() const;
     /** get Z-position array for all the particles */
     double* getZall() const;
     /** get u (X-velocity) array for all the particles */
     double* getUall() const;
     /** get v (Y-velocity) array for all the particles */
     double* getVall() const;
     /** get w (Z-velocity) array for all the particles */
     double* getWall() const;
     /** get the ID array   */
     unsigned long* getParticleIDall() const;
     /** get X-position of particle with label indexPart */
     double getX(int indexPart)const;
     /** get Y-position of particle with label indexPart */
     double getY(int indexPart)const;
     /** get Z-position of particle with label indexPart */
     double getZ(int indexPart)const;
     /** get u (X-velocity) of particle with label indexPart */
     double getU(int indexPart) const;
     /** get v (Y-velocity) of particle with label indexPart */
     double getV(int indexPart) const;
     /** get w (Z-velocity) of particle with label indexPart */
     double getW(int indexPart) const;
     /** get ID of particle with label indexPart */
     unsigned long getParticleID(int indexPart) const;
     /**get charge of particle with label indexPart */
     double getQ(int indexPart) const;
     /** get charge of array for ID particles */
     double* getQall() const;
     /** get the number of particles of this subdomain */
     int getNOP() const;
     /** Print particles info: positions, velocities */
     void Print(VirtualTopology* ptVCT) const;
     /** Print the number of particles of this subdomain */
     void PrintNp(VirtualTopology* ptVCT) const;
     
  protected:
    /** number of species */
    int ns;
    /** maximum number of particles of this species on this domain. used for memory allocation */
    int npmax;
    /** nxmax        */
    int nxmax;
    /** number of particles of this species on this domain */
    int nop;
    /** number of particles per cell */
    int npcel;
    /** number of particles per cell - X direction */
    int npcelx;

    /** charge to mass ratio */
    double qom;
    /** thermal velocity  - Direction X*/
    double uth;
	/** u0 Drift velocity - Direction X */
    double u0;
    /** Positions arra - X component */
    double* x;
    /** Velocities array - X component */
    double*  u;
    /** TrackParticleID */
    bool TrackParticleID;
    /** ParticleID */
    unsigned long* ParticleID;
    /** rank of processor in which particle is created (for ID) */
    int  BirthRank[2];
    /** number of variables to be stored in buffer for communication for each particle  */
    int nVar;
    /** Charge array */
    double* q;
    /** Simulation domain lengths */
    double xstart, xend,invVOL;
    /** time step */
    double dt;
    /** Lx = simulation box length - x direction   */
    double Lx;
    /** number of nodes */
    int nxn;
    /** buffers for communication */
    /** size of sending buffers for exiting particles, DEFINED IN METHOD "COMMUNICATE" */
    int buffer_size;
    /** buffer with particles going to the right processor - Direction X */
    double *b_XDX;
    /** pointer to the buffer for resizing */
    double *b_XDX_ptr;
    /** buffer with particles going to the left processor - Direction X */
    double *b_XSN;
    /** pointer to the buffer for resizing */
    double *b_XSN_ptr;
    /** number of particles exiting to X-RIGHT per cycle*/
    int npExitXright;
    /** number of particles exiting to X-LEFT per cycle*/
    int npExitXleft;
	/** total number of particles exiting per cycle */
    int npExit;
    /** number of particles not in the right domain   */
    int rightDomain;
    /** bool for communication verbose */
    bool cVERBOSE;
    /** Boundary condition on particles:
    <ul>
    	<li>0 = exit</li>
    	<li>1 = perfect mirror</li>
    	<li>2 = riemission</li>
    	<li>3 = periodic condition </li>
    </ul>
    */
    /** Boundary Condition Particles: FaceXright */
    int bcPfaceXright;
    /** Boundary Condition Particles: FaceXleft */
    int bcPfaceXleft;
    /** Boundary Condition Particles: FaceYright */
    int bcPfaceYright;
    /** Boundary Condition Particles: FaceYleft */
    int bcPfaceYleft;
    /** restart variable for loading particles from restart file */
    int restart;
	
	/** BEAM parameter */
	/** species weight */
	double w;
	/** beam kinetic energy */
	double ekin;
	/** beam current */
	double ibeam;
	/** charge of ion */
	int zion;
	/** a of ion */
	int aion;
	
	/** constants */
	/** c = light speed */
    double c;
    /** atomic mass unit [kg] */
    double amu ;
    /** proton charge */
    double echarge;
    /** electron mass */
    double emass;
    /** Permittivity of free space*/
    double eps0;
    /** Permeability of free space */
    double mu0;
    /** Conversion factor, Joules per eV */
    double jperev;
    /** Boltzmann's constant */
    double boltzmann;
    /** pi */
    double pi;
};


#endif






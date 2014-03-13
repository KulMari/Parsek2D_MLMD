/*******************************************************************************************
Particles2Dcommcomm.h  -  Class for particles of the same species, in a 2D space and 3component velocity with communications methods
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef Part2DCOMM_H
#define Part2DCOMM_H

#include "Particles.h"
// AMR, ME
#include <stdio.h>
#include <string.h>
// end AMR, ME

/**
* 
* Abstract class for particles of the same species, in a 2D space and 3component velocity with communications methods
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Particles2Dcomm : public Particles {
  public:

  double x_sum_CC;
  double u_sum_CC;
  double roundPrec(double x, int prec);
     /** constructor */
     Particles2Dcomm();
     /** destructor */
     ~Particles2Dcomm();
     /** allocate particles */
     void allocate(int species, CollectiveIO* col,VirtualTopology* vct, Grid* grid);
    
     /** interpolation method GRID->PARTICLE order 1: CIC */
     void interpP2G(Field* EMf, Grid *grid, VirtualTopology* vct);
     /** method for communicating exiting particles to X-RIGHT, X-LEFT, Y-RIGHT, Y-LEFT, Z-RIGHT, Z-LEFT processes */
     // AMR, ME: I need the Grid object for the level
     int communicate(VirtualTopology* ptVCT);
     int communicate(VirtualTopology* ptVCT, Grid* grid, int BC_partCommunicate); 
     // end AMR, ME
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
     /** method to debuild the buffer received, not to be used anymore */
     int unbuffer(double *b_);
     /** method to debuild the buffer received */
     int unbuffer(double *b_, VirtualTopology *ptVCT);
     /** resize the receiving buffer */
     bool resize_buffers(int new_buffer_size);
     /** resize the buffer for SP send/ receives*/
     bool resize_buffersSP(int new_np);
     /** resize the buffer for OS send/ receives*/
     bool resize_buffersOS(int new_np);
     /** a method to compute how many particles are not in the right domain */
     int isMessagingDone(VirtualTopology* ptVCT);
     /** a method to compute if refined particle boundary communication has to be repeated */
     int isMessagingDoneSP(VirtualTopology* ptVCT);
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

     /**AMR methods, ME*/
     /**init operations connected with the repopulation of particles...*/
     int initPRAVariables(int species, CollectiveIO* col,VirtualTopology* vct, Grid* grid, Field* EMf);
     /**prints some info about the PRA area; only in verbose mode*/
     void checkAfterInitPRAVariables(int species, CollectiveIO* col,VirtualTopology* vct, Grid* grid); 
     /*apply BC conditions on particles in communicate*/
     bool applyParticleBC(int BC_partCommunicate, int np_current, int *nplast, int *npDeletedBoundary, Grid * grid, VirtualTopology * ptVCT);
     /** method for communicating exiting particles to X-RIGHT, X-LEFT, Y-RIGHT, Y-LEFT, Z-RIGHT, Z-LEFT processes,
	 /** other variation of communicate, to be used when communicated splitted particles between finer grid procs*/
     int communicateSP(VirtualTopology* ptVCT, Grid* grid); 
     // to communicate particles for OS operations (patching the "other side" of the moments on the refined grid) to the right processor
     int communicateOS(VirtualTopology* ptVCT, Grid* grid);
     // init the vector for communication of PRA particles from coarser to finer grids
     void initPRABuffers(Grid* grid, VirtualTopology* vct);
     // decide if a coarser grid particle falls in the RPA area and therefore needs to be communicated
     // in case, put in the comm vectors
     // ParticleNumber is the current number of the particle to deal with in the species object
     int PRARepopulationAdd(int ParticleNumber);
     // resize PRA send buffers; TO DO; resize the receive ones as well
     void resizePRAbuffers(int new_max_nop);
     // send PRA particles to finer grids
     int PRASend(Grid* grid, VirtualTopology* vct);
     // receive PRA particles from coarser grids
     int PRAReceive(Grid* grid, VirtualTopology* vct, Field* EMf);
     // when subcycling, to inject repopulated particles in the refined grid
     int SubCyclingParticles(Grid* grid, bool CoarseOp, VirtualTopology* ptVCT);
     // used for internal ops for particle repopulation, the size of the vector must be MAX_NP_REPOP_SIZE*nVar
     void setToMINVAL(double *vec);
     // used for internal ops for particle repopulation                             
     void setToMINVAL_SP(double *vec, int newNP);
     // used for internal ops for patching the refined grid boundary values
     void setToMINVAL_OS(double *vec, int newNP);
     // used for internal ops in communicate; the first buffer_size elements of the vector to MINVAL
     void setToMINVAL_comm(double *vec);
     // split a coarse grid particle into a finer grid particle
     // ret: <0, if the size of SplittedP needs to be increase, 0 otherwise
     // arg: n_rec_p: coarse particle number in REPOP_receive_b
     //      n_ref_p_: finer particle number in SplittedParticles_Comm (array of the splitted particles to communicate)
     // modifies nop
     // BE CAREFUL TO NVAR FOR ARRAYS!!!
     int SplitCoarseParticle(VirtualTopology* vct, Grid* grid, int n_rec_p); 
     // add splitted particles to the appropriate comm vector; dim check first and index increment
     int pack_SP(double *vec, int part_index, double x_F, double y_F, double u_F, double v_F, double w_F, double q_F, unsigned long ParticleID_F);
     // count how many particles sit in the PRA area -- very expensive, just for test --
     int CountPRAParticles(VirtualTopology* vct);
     int printPRAparticles(VirtualTopology* vct, Grid* grid);
     void setPRACollectionMethod(int safe);
     int getPRACollectionMethod();
     // safe way of collective repopulation particles in the coarser levels
     int CollectivePRARepopulationAdd(VirtualTopology* ptVCT, Grid* grid);
     // Randomize position for particles in the PRA, to see if it affects the noise in the refiend grid with high RFs                                          
     int RandomizePositionPRAParticles(int is, VirtualTopology* vct,  Grid* grid);
     // OS operations (different from the normal buffer / del_pack/ unbuffer because tilde values are not passed)
     /** put a OS particle exiting to X-LEFT in the bufferXLEFT for communication*/
     void bufferXleftOS(double *b_, int np, VirtualTopology* vct);
     /** put a OS particle exiting to X-RIGHT in the bufferXRIGHT for communication*/
     void bufferXrightOS(double *b_, int np, VirtualTopology* vct);
     /** put a particle exiting to Y-LEFT in the bufferYLEFT for communication*/
     void bufferYleftOS(double *b_, int np, VirtualTopology* vct);
     /** put a particle exiting to Y-RIGHT in the bufferYRIGHT for communication*/
     void bufferYrightOS(double *b_, int np, VirtualTopology* vct);
     /** Delete the OS particle from a list(array) and pack the list(array) */
     void del_packOS(int np_current, int *nplast);
     /*unbuffer OS particles*/
     int unbufferOS(double *b_, VirtualTopology *ptVCT);
     /** interpolation method GRID->PARTICLE order 1: CIC, refined grid, OS particles */
     int interpP2G_OS(Field* EMf, Grid *grid, VirtualTopology* vct);
     /** to store the values of corner nodes to the original value before interpNode, in case OS operations are done**/
     void storeCornerOsValues(Field * EMf, VirtualTopology* vct);
     /** to restore the values of corner nodes to the original value before interpNode, in case OS operations are done (otherwise InterpNode operations done twice there)**/
     void REstoreCornerOsValues(Field * EMf, VirtualTopology* vct);
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
    /** number of particles per cell - Y direction */
    int npcely;
    /** charge to mass ratio */
    double qom;
    /** thermal velocity  - Direction X*/
    double uth;
    /** thermal velocity  - Direction Y*/
    double vth;
    /** thermal velocity  - Direction Z*/
    double wth;
    /** u0 Drift velocity - Direction X */
    double u0;
    /** v0 Drift velocity - Direction Y */
    double v0;
    /** w0 Drift velocity - Direction Z */
    double w0;
    /** Positions arra - X component */
    double* x;
    /** Positions array - Y component */
    double*  y;
    /** Implicit positions array - X component */
    double*  xptilde;
    /** Implicit positions array - Y component */
    double*  yptilde;
    /** Velocities array - X component */
    double*  u;
    /** Implicit velocities array - X component */
    double*   uptilde;
    /** Velocities array - Y component */
    double*   v;
    /** Implicit velocities array - Y component */
    double*   vptilde;
    /** Velocities array - Z component */
    double*  w;
    /** Implicit velocities array - Z component */
    double*  wptilde;
    /** TrackParticleID */
    bool TrackParticleID;
    /** for MLMD operations, only for coarse grid with PRACollectionMethod=0; if 0, coarse particle not yet collected for repopulation ops, otherwise already collected **/
    bool* AlreadyAccumulated;
    /** ParticleID */
    unsigned long* ParticleID;
    /** rank of processor in which particle is created (for ID) */
    int  BirthRank[2];
    /** number of variables to be stored in buffer for communication for each particle  */
    int nVar;
    /** Charge array */
    double* q;
    /** Simulation domain lengths */
    double xstart, xend, ystart, yend;
	/** cell size */
	double dx, dy, invVOL;
    /** time step */
    double dt;
    /** am I subcycling? */
    int SubCycling;
    /** Lx = simulation box length - x direction   */
    double Lx;
    /** Ly = simulation box length - y direction   */
    double Ly;
    /** number of nodes */
    int nxn;
    /** number of nodes */
    int nyn;
    /** buffers for communication */
    /** size of sending buffers for exiting particles, DEFINED IN METHOD "COMMUNICATE" */
    int buffer_size;
    /** max length of communication buffers*/
    int MAX_BUFFER_SIZE;
    // used for some communicate internal ops; size MAX_BUFFER_SIZE
    double *MIN_VAL_VEC_COMM;
    /** buffer with particles going to the right processor - Direction X */
    double *b_XDX;
    /** pointer to the buffer for resizing */
    double *b_XDX_ptr;
    /** buffer with particles going to the left processor - Direction X */
    double *b_XSN;
    /** pointer to the buffer for resizing */
    double *b_XSN_ptr;
    /** buffer with particles going to the right processor - Direction Y */
    double *b_YDX;
    /** pointer to the buffer for resizing */
    double *b_YDX_ptr;
    /** buffer with particles going to the left processor - Direction Y */
    double *b_YSN;
    /** pointer to the buffer for resizing */
    double *b_YSN_ptr;
    /** number of particles exiting to X-RIGHT per cycle*/
    int npExitXright;
    /** number of particles exiting to X-LEFT per cycle*/
    int npExitXleft;
    /** number of particles exiting to Y-RIGHT per cycle*/
    int npExitYright;
    /** number of particles exiting to Y-LEFT per cycle*/
    int npExitYleft;
    /** total number of particles exiting per cycle */
    int npExit;
	/** total number of particles eliminated at the boundaries*/
    int npDeletedBoundary;
	 /** total number of particles eliminated at the boundaries*/
    int npDeletedDipole;
    /** number of particles not in the right domain   */
    int rightDomain;
    /** number of particles not in the right domain, X or Y dir   */
    int rightDomainX, rightDomainY;
    /** flags for rightDomainX, rigthDomainY */
    bool firstrightDomainX, firstrightDomainY;
    /** Bool for communication verbose */
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
    /** speed of light in vacuum */
    double c;
    /** reconnection thickness */
    double delta;
    /** restart variable for loading particles from restart file */
    int restart;
    /** number of iteration for PC mover */
    int NiterMover;
    /** velocity of the injection of the particles */
    double Vinj;
    /** AMR variables, ME */
    /**number of cells, ghost cell INCLUDED, for particle repopulation; x left*/
    int PRA_Xleft;
    /**number of cells, ghost cell INCLUDED, for particle repopulation; x right*/
    int PRA_Xright;
    /**number of cells, ghost cell INCLUDED, for particle repopulation; y left*/
    int PRA_Yleft;
    /**number of cells, ghost cell INCLUDED, for particle repopulation; y right*/
    int PRA_Yright;
    /**coordinates on the local grid for particle repopulation*/
    /**x start left*/
    double PRA_oxStartLeft;
    /**x end left*/
    double PRA_oxEndLeft;
    /**x start right*/
    double PRA_oxStartRight;
    /**x end right*/
    double PRA_oxEndRight;
    /**y start left*/
    double PRA_oyStartLeft;
    /**y end left*/
    double PRA_oyEndLeft;
    /**y start right*/
    double PRA_oyStartRight;
    /**y end right*/
    double PRA_oyEndRight;

    /**if there is a child, coordinates on the local grid for the child's PRA*/
    /*the dx/ dy for children particle repopop is already taken into account*/
    /**x start left*/
    double PRA_CoxStartLeft;
    /**x end left*/
    double PRA_CoxEndLeft;
    /**x start right*/
    double PRA_CoxStartRight;
    /**x end right*/
    double PRA_CoxEndRight;
    /**y start left*/
    double PRA_CoyStartLeft;
    /**y end left*/
    double PRA_CoyEndLeft;
    /**y start right*/
    double PRA_CoyStartRight;
    /**y end right*/
    double PRA_CoyEndRight;

    //for coarse grid, to see if a refined grid is intersecting you for PRA purposes
    bool PRAIntersection;

    /* needed 
       -- for level 0 to mark wether the communicate is the last one, with definitive position, or one in the inner mover
          (in the other levels difference with communicateFirst and Second)
       -- in the unbuffer of all levels to decide wether PRA ops must be initiated
          (it's the unbuffer after definitive positions or an intermediate one) 
       set at the appropriate time in the particle mover*/
    int LastCommunicate; // moved to public for debug issues
    // needed  to decide wether to pack PRA particles for finer grids (also checked in the coarser grid)
    // practically, tells if there is a finer grids without having to pass the Grid class
    // init in allocate
    int FinerLevels_PRAOps;

    // which particles need to be sent from this level to the finer ones for particle repop
    // bottom, top, left, right mark the directions
    // the buffers are populated one particle at the time during communicate ops
    // allocated in allocate with MAX_REPOP_SIZE
    /** buffers for communication */
    /** buffer with particles going to the finer procs in the bottom (0) direction */
    int MAX_NP_REPOP_SIZE;  // MAX number of particles for repopulation  (for vector sizes, be careful to nVar)
    
    double *MIN_VAL_VEC;  // useful for some intermediate ops
    double *MIN_VAL_VEC_SP;  // useful for some intermediate ops 

    double *REPOP_b_BOTTOM;
    /** pointer to the buffer for resizing */
    double *REPOP_b_BOTTOM_ptr;
    /** size, in number of particles*/
    int np_REPOP_b_BOTTOM;  // number of particles

    /** buffer with particles going to the finer procs in the top (1) direction */
    double *REPOP_b_TOP;
    /** pointer to the buffer for resizing */
    double *REPOP_b_TOP_ptr;
    /** size*/
    int np_REPOP_b_TOP;

    /** buffer with particles going to the finer procs in the left (2) direction */
    double *REPOP_b_LEFT;
    /** pointer to the buffer for resizing */
    double *REPOP_b_LEFT_ptr;
    /** size*/
    int np_REPOP_b_LEFT;
    
    /** buffer with particles going to the finer procs in the right (3) direction */
    double *REPOP_b_RIGHT;
    /** pointer to the buffer for resizing */
    double *REPOP_b_RIGHT_ptr;
    /** size*/
    int np_REPOP_b_RIGHT;

    // buffer for receiving PRA particles, used only if level>0
    double *REPOP_receive_b;
    // buffer for generated finer grid particles not belonging to the current processor (to communicate away)
    // for size, nVar * npart
    double *SplittedParticles_Comm_BOTTOM;
    double *SplittedParticles_Comm_BOTTOM_ptr;//alias for resize
    double *SplittedParticles_Comm_TOP;
    double *SplittedParticles_Comm_TOP_ptr;//alias for resize  
    double *SplittedParticles_Comm_LEFT;
    double *SplittedParticles_Comm_LEFT_ptr;//alias for resize  
    double *SplittedParticles_Comm_RIGHT;
    double *SplittedParticles_Comm_RIGHT_ptr;//alias for resize  
    // max number of particles in SplittedParticles_Comm (for vector sizes, be careful to nVar)
    int max_np_SplitPartComm;
    // max_np_SplitPartComm can be increased through resize_buffersSP
    // if, during a resize, max_np_SplitPartComm > MAX_NP_SPLIPARTCOMM, the simulation is stopped
    int MAX_NP_SPLIPARTCOMM;

    // buffer for OS particles
    // for size, nVarOS * npart
    double *OSParticles_Comm_BOTTOM;
    double *OSParticles_Comm_BOTTOM_ptr; //alias for resizing
    double *OSParticles_Comm_TOP;
    double *OSParticles_Comm_TOP_ptr; //alias for resizing
    double *OSParticles_Comm_LEFT;
    double *OSParticles_Comm_LEFT_ptr; //alias for resizing  
    double *OSParticles_Comm_RIGHT;
    double *OSParticles_Comm_RIGHT_ptr; //alias for resizing  
    double *MIN_VAL_VEC_OS;
    // max number of particles in OSParticles_Comm_ (for vector sizes, be careful to nVar) 
    int max_np_OsPartComm;
    // max_np_OsPartComm can be increased through resize_buffersOS 
    // if, during a resize, max_np_OSPartComm > MAX_NP_OSPARTCOMM, the simulation is stopped 
    int MAX_NP_OSPARTCOMM;
    int nVarOS;
    // fix corner for OS particles -- only if OS operations
    double OSfix_upperLeft_rhons;
    double OSfix_upperRight_rhons;
    double OSfix_lowerLeft_rhons;
    double OSfix_lowerRight_rhons;

    double OSfix_upperLeft_Jx;
    double OSfix_upperRight_Jx;
    double OSfix_lowerLeft_Jx;
    double OSfix_lowerRight_Jx;

    double OSfix_upperLeft_Jy;
    double OSfix_upperRight_Jy;
    double OSfix_lowerLeft_Jy;
    double OSfix_lowerRight_Jy;
    
    double OSfix_upperLeft_Jz;
    double OSfix_upperRight_Jz;
    double OSfix_lowerLeft_Jz;
    double OSfix_lowerRight_Jz;

    double OSfix_upperLeft_Pxx;
    double OSfix_upperRight_Pxx;
    double OSfix_lowerLeft_Pxx;
    double OSfix_lowerRight_Pxx;

    double OSfix_upperLeft_Pxy;
    double OSfix_upperRight_Pxy;
    double OSfix_lowerLeft_Pxy;
    double OSfix_lowerRight_Pxy;

    double OSfix_upperLeft_Pxz;
    double OSfix_upperRight_Pxz;
    double OSfix_lowerLeft_Pxz;
    double OSfix_lowerRight_Pxz;

    double OSfix_upperLeft_Pyy;
    double OSfix_upperRight_Pyy;
    double OSfix_lowerLeft_Pyy;
    double OSfix_lowerRight_Pyy;

    double OSfix_upperLeft_Pyz;
    double OSfix_upperRight_Pyz;
    double OSfix_lowerLeft_Pyz;
    double OSfix_lowerRight_Pyz;
    
    double OSfix_upperLeft_Pzz;
    double OSfix_upperRight_Pzz;
    double OSfix_lowerLeft_Pzz;
    double OSfix_lowerRight_Pzz;
    // end fix corner
	
    // these variables have the same name as in EMfields, but they are NOT the same; calculated in initPRAVariables
    // used by PRASend
    int nmessageBC;
    int *targetBC;
    int *BCSide;
    int targetBOTTOM;
    int targetTOP;
    int targetLEFT;
    int targetRIGHT;
    // used by PRAReceive
    int nmessagerecuBC;
    int *fromBC;
    int *BCSidecu;
    // these variables have the same name as in EMfields, but they are NOT the same; calculated in initPRAVariables

    // variables used in splitting particles, which I don't want to recalculate every time
    int ratio;

    // variable for debugging PRASend/ PRAReceive
    int GeneratedRP;    // how many refined particles are generated by splitting
    //int LocallyAcceptedRP; // how many of the generated particles are accepted by the current processor
    //int ExternallyGeneratedAcceptedRP; // how many of the particles generated by other processors are accepted by the current processor
    int AcceptedRPAfterSplit; // how many of the generated particles are preliminarly accepted (fall in the PRA)     
    int AcceptedRPAfterComm; // how many of the generated particles are definitively accepted (after the communicateSP) 
    int RejectedRP;       // how many particles are rejected (do not fall in the PRA)
    int CommunicatedRP;   // particles generated by the local proc, accepetd and sent to other procs

    // if Subcycling and refined grid, the # of particles from the coarse grid, to reinject in the refined grid at each cycle when they are not directly received from the coarse grid
    int RP_nop;

    // variables for distribute splitted particles in the right processor
    int n_ref_p_comm_BOTTOM; 
    int n_ref_p_comm_TOP; 
    int n_ref_p_comm_LEFT; 
    int n_ref_p_comm_RIGHT;

    // used in communicateSP
    double Modified_xstart;
    double Modified_xend;
    double Modified_ystart;
    double Modified_yend;

    // 0: coarse grid particles to send to refined grid collected one by one during the mover
    // 1: coarse grid particles to send to refined grid collected collectively AFTER the mover - safer but slower
    int PRACollectionMethod;

    // to provide "other side" moments for the refined grid
    int nop_OS;
    double *OS_x;
    double *OS_y;      
    double *OS_u;
    double *OS_v;
    double *OS_w;
    double *OS_q;
    unsigned long* OS_ParticleID;

    int npmax_OS; // max number of particles in the OS buffers

    // to hos repopulated particles when subcyclig
    double *RP_x;
    double *RP_y;
    double *RP_u;
    double *RP_v;
    double *RP_w;
    double *RP_q;
    unsigned long *RP_ParticleID;

    int SizeRP_Sub;


    // for the mover
    double ***XN;
    double ***YN;
}
;


#endif






/***************************************************************************
                  VCtopology.h  -  Virtual cartesian topology
                  A virtual topology is a mechanism for naming the processes
                  in a communicator in a way that fits the communication
                  pattern better. Since our processes will communicate mainly
                  with the nearest neighbours after the fashion of a two-dimensional
                  grid, we create a virtual topology to reflect this fact
                             -------------------
    begin                : Fri Jun 4 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef VCtopology_H
#define VCtopology_H

#include "mpi.h"
#include "VirtualTopology.h"
#include <iostream>



using std::cout;
using std::endl;

/**
*  
* Virtual cartesian topology
* A virtual topology is a mechanism for naming the processes
* in a communicator in a way that fits the communication
* pattern better. Since our processes will communicate mainly
* with the nearest neighbours after the fashion of a two-dimensional
* grid, we create a virtual topology to reflect this fact
* @date Fri Jun 4 2004
* @par Copyright:
* (C) 2004 Los Alamos National Laboratory
* @author Stefano Markidis, Giovanni Lapenta
* @version 1.0
*/


class VCtopology : public VirtualTopology {
  public:
    /** constructor: Define topology parameters: dimension, domain decomposition,... */
    VCtopology();
    /** destructor */
    ~VCtopology();
    /** Find the neighbors in the new communicator  */
    void setup_vctopology(MPI_Comm comm_old);
    /** Print topology info */
    void Print();
    /** Print the mapping of topology */
    void PrintMapping();
    /** get the CART_COMM communicator, associated to one grid **/
    MPI_Comm getCART_COMM();
    /** get the CART_COMM communicator, associated to all grid **/
    MPI_Comm getCART_COMM_TOTAL();
    /** get and set XLEN */
    int getXLEN();
    void setXLEN(int XLENinput);
    /** get and set YLEN */
    int getYLEN();
    void setYLEN(int YLENinput);
    /** set divisions **/
    void setDivisions();
    /** get and set nprocs */
    int getNprocs();
    void setNprocs();
    /** get and set ngridss */
    int getNgrids();
    void setNgrids(int ngridsinput);
    /** get periodicity on boundaries - DIRECTION X*/
    bool getPERIODICX();
    /** get periodicity on boundaries - DIRECTION Y*/
    bool getPERIODICY();
    /** set periodicity */
    void setPeriodicity(int xleft, int xright, int yleft, int yright);
    /** get the cartesian rank of the process */
    int getCartesian_rank();
    // AMR ME
    // gets the rank relative to ALL levels, in MPI_COMM_TOTAL                                                      
    int getCartesian_rank_COMMTOTAL();
    // get the rank in MPI_COMM_WORLD                                                                                         
    int getRank_MPI_COMM_WORLD();
    // end AMR ME
    /** get the cartesian rank of XLEFT neighbor */
    int getXleft_neighbor();
    /** get the cartesian rank of XRIGHT neighbor */
    int getXright_neighbor();
    /** get the cartesian rank of YLEFT neighbor */
    int getYleft_neighbor();
    /** get the cartesian rank of YRIGHT neighbor */
    int getYright_neighbor();
    /** get the cartesian rank of XLEFT(-) YLEFT(-) SAME Z neighbor  */
    int getXleftYleft_neighbor();
    /** get the cartesian rank of XLEFT(-) YRIGHT(+) SAME Z neighbor */
    int getXleftYright_neighbor();
    /** get the cartesian rank of XRIGHT(+) YLEFT(-) SAME Z neighbor */
    int getXrightYleft_neighbor();
    /** get the cartesian rank of XRIGHT(+) YRIGHT(+) SAME Z neighbor */
    int getXrightYright_neighbor();
       
    /** get the coordinates in dir direction of process*/
    int getCoordinates(int dir);
    /** get Periodicity condition in dir direction */
    int getPeriods(int dir);
    /** if cVERBOSE == true, print to the screen all the comunication */
    bool getcVERBOSE();
    // set RefLevelAd
    void setRefLevelAdj(int val);
    // get RefLevelAd
    int getRefLevelAdj();
  private:
    /** New communicator with virtual cartesian topology */
    MPI_Comm CART_COMM;
    MPI_Comm CART_COMM_test;
    /** New communicator with virtual cartesian topology including grid levels */
    MPI_Comm CART_COMM_TOTAL;
    /** MPI status during sending and receiving communication */
    MPI_Status status;
    /** Direction X for shift MPI_Cart_Shift*/
    int XDIR;
    /** Direction Y for shift MPI_Cart_Shift*/
    int YDIR;
    /** RIGHT = +    upwards   shift */
    int RIGHT;
    /** LEFT  = -    downwards shift */
    int LEFT;
    /** dimension of virtual topology */
    int PROCDIM;
    /** number of subdomains - Direction X */
    int XLEN;
    /** number of subdomains - Direction Y */
    int YLEN;
    /** nprocs = number of processors */
    int nprocs;
    /** ngrids = number of grids */
    int ngrids;
    /** periodicity on boundaries - DIRECTION X*/
    bool PERIODICX;
    /** periodicity on boundaries - DIRECTION Y*/
    bool PERIODICY;
    /** rank may be reordered     */
    int reorder;
    /** arrays for Create_Cart_create  */
    int *divisions;
    int *periodic_divisions;
    int *periods;
    int *twoDperiods;
    int *global_coordinates;

    /** cartesian rank */
    int cartesian_rank;
    /** coordinates on processors grid */
    int* coordinates;
    /** cartesian rank of XLEFT neighbor */
    int xleft_neighbor;
    /** cartesian rank of XRIGHT neighbor */
    int xright_neighbor;
    /** cartesian rank of YLEFT neighbor */
    int yleft_neighbor;
    /** cartesian rank of YRIGHT neighbor */
    int yright_neighbor;

    /** cartesian rank of 12 neighbors: exchange a line during ghost transmission */
    /** cartesian rank of XLEFT(-) YLEFT(+) SAME Z  neighbor */
    int XleftYleft_neighbor;
    /** cartesian rank of XLEFT(-) YRIGHT(+) SAME Z neighbor */
    int XleftYright_neighbor;
    /** cartesian rank of XRIGHT(+) YLEFT(-) SAME Z neighbor */
    int XrightYleft_neighbor;
    /** cartesian rank of XRIGHT(+) YRIGHT(+) SAME Z neighbor */
    int XrightYright_neighbor;
    
    
   
    /** if cVERBOSE == true, print to the screen all the comunication */
    bool cVERBOSE;

    // ME
    int RefLevelAdj; /* options:                                                              
		       --0: same adjust for coarse and refined level (multiply)
		       --1: interp of OS particles for the refined level */
};

#endif


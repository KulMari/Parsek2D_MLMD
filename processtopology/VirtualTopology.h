/***************************************************************************
    VirtualTopology.h - Abstract Base class for virtual process topologies
                             -------------------
    begin                : Wed Jun 2 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef VirtualTopology_H
#define VirtualTopology_H

#include "mpi.h"
/**
*  Abstract base class for virtual process topologies
*
* @date Fri Jun 4 2004
* @par Copyright:
* (C) 2004 Los Alamos National Laboratory
* @author Stefano Markidis, Giovanni Lapenta
* @version 1.0
*/
class VirtualTopology {
  public:
    /** Find the neighbors in the new communicator  */
    virtual void setup_vctopology(MPI_Comm comm_old)=0;
    /** Get CART_COMM communicators **/
    virtual MPI_Comm getCART_COMM()=0;
    /** Get CART_COMM_TOTAL communicators **/
    virtual MPI_Comm getCART_COMM_TOTAL()=0;
    /** get boundary communicators **/
    virtual MPI_Comm getCOMM_B_LEFT()=0;
    virtual MPI_Comm getCOMM_B_RIGHT()=0;
    virtual MPI_Comm getCOMM_B_BOTTOM()=0;
    virtual MPI_Comm getCOMM_B_TOP()=0;
    virtual MPI_Comm getCOMM_B_ALL()=0;
    /** get ngrids */
    virtual int getNgrids()=0;
    /** Print topology info */
    virtual void Print()=0;
    /** Print the mapping of topology */
    virtual void PrintMapping()=0;
    /** get XLEN */
    virtual int getXLEN()=0;
    /** get YLEN */
    virtual int getYLEN()=0;
    /** get nprocs */
    virtual int getNprocs()=0;
    /** get periodicity on boundaries - DIRECTION X*/
    virtual bool getPERIODICX()=0;
    /** get periodicity on boundaries - DIRECTION Y*/
    virtual bool getPERIODICY()=0;
    /** get the cartesian rank of the process */
    virtual int getCartesian_rank()=0;
    /** get the cartesian rank of the process in the global communicator (on the multiple levels)*/
    virtual int getCartesian_rank_COMMTOTAL()=0;
    // get the rank in MPI_COMM_WORLD (not necessarely the same as in commtotal)
    virtual int getRank_MPI_COMM_WORLD()=0;
    /** get the rank of the process on the boundary communicators */
    virtual int getRank_COMM_B_LEFT()=0;
    virtual int getRank_COMM_B_RIGHT()=0;
    virtual int getRank_COMM_B_BOTTOM()=0;
    virtual int getRank_COMM_B_TOP()=0;
    /** get the cartesian rank of XLEFT neighbor */
    virtual int getXleft_neighbor()=0;
    /** get the cartesian rank of XRIGHT neighbor */
    virtual int getXright_neighbor()=0;
    /** get the cartesian rank of YLEFT neighbor */
    virtual int getYleft_neighbor()=0;
    /** get the cartesian rank of YRIGHT neighbor */
    virtual int getYright_neighbor()=0;
    /** get the cartesian rank of XLEFT(-) YLEFT(-) SAME Z neighbor  */
    virtual int getXleftYleft_neighbor()=0;
    /** get the cartesian rank of XLEFT(-) YRIGHT(+) SAME Z neighbor */
    virtual int getXleftYright_neighbor()=0;
    /** get the cartesian rank of XRIGHT(+) YLEFT(-) SAME Z neighbor */
    virtual int getXrightYleft_neighbor()=0;
    /** get the cartesian rank of XRIGHT(+) YRIGHT(+) SAME Z neighbor */
    virtual int getXrightYright_neighbor()=0;
    /** get the neighbors in the boundary communicators */
    virtual int getLeftNeighbor_COMM_B_LEFT()=0;
    virtual int getRightNeighbor_COMM_B_LEFT()=0;
    virtual int getLeftNeighbor_COMM_B_RIGHT()=0;
    virtual int getRightNeighbor_COMM_B_RIGHT()=0;
    virtual int getLeftNeighbor_COMM_B_BOTTOM()=0;
    virtual int getRightNeighbor_COMM_B_BOTTOM()=0;
    virtual int getLeftNeighbor_COMM_B_TOP()=0;
    virtual int getRightNeighbor_COMM_B_TOP()=0;
    /** get the coordinates in dir direction of process*/
    virtual int getCoordinates(int dir)=0;
    /** get Periodicity condition in dir direction */
    virtual int getPeriods(int dir)=0;
    /** if cVERBOSE == true, print to the screen all the comunication */
    virtual bool getcVERBOSE()=0;

    //ME 
    // set RefLevelAd                                                                                                                                   
    virtual void setRefLevelAdj(int val)=0;
    // get RefLevelAd                                                                                                                                   
    virtual int getRefLevelAdj()=0;
};

#endif

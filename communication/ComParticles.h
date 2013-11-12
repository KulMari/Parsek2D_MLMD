/***************************************************************************
    ComParticles.h  -  Library to manage communication of particles among processors
                       -------------------
    begin                : Fri Jun 4 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef ComParticles_H
#define ComParticles_H

#include <mpi.h>
#include "ComBasic.h"


/** comunicate particles and receive particles in X direction */
inline void communicateParticles(int buffer_size,double *b_Xleft, double *b_Xright,VirtualTopology *vct){
      communicateParticlesDIR(vct->getCART_COMM(), buffer_size,vct->getCartesian_rank(),vct->getXright_neighbor(),vct->getXleft_neighbor(),0,vct->getXLEN(),vct->getYLEN(),b_Xright,b_Xleft);
      
}
/** comunicate particles and receive particles X and Y directions direction */
inline void communicateParticles(int buffer_size,double *b_Xleft, double *b_Xright,double *b_Yleft, double *b_Yright,VirtualTopology *vct){
      
      communicateParticlesDIR(vct->getCART_COMM(), buffer_size,vct->getCartesian_rank(),vct->getXright_neighbor(),vct->getXleft_neighbor(),0,vct->getXLEN(),vct->getYLEN(),b_Xright,b_Xleft);
      communicateParticlesDIR(vct->getCART_COMM(), buffer_size,vct->getCartesian_rank(),vct->getYright_neighbor(),vct->getYleft_neighbor(),1,vct->getXLEN(),vct->getYLEN(),b_Yright,b_Yleft);
      
}
/** communicate the number of particles are not in the right domain*/
inline int reduceNumberParticles(MPI_Comm CART_COMM, int rightDomain){
  int result=0;
  MPI_Barrier(CART_COMM);
  MPI_Allreduce(&rightDomain, &result, 1, MPI_INT, MPI_SUM, CART_COMM);
  return(result);
}
/** communicate the maximum number of particles from a domain */
inline int reduceMaxNpExiting(MPI_Comm CART_COMM, int npExitingMax){
  int result=0;
  MPI_Barrier(CART_COMM);
  MPI_Allreduce(&npExitingMax, &result, 1, MPI_INT, MPI_MAX, CART_COMM);
  return(result);
}

/** communicate particles on the boundary communicators */
/** commID:  0:left, 1:right, 2:bottom, 3:top */

inline int communicateParticles_BoundaryComm(int comm_ID, int buffer_size,double *b_Xleft, double *b_Xright,VirtualTopology *vct){

  MPI_Comm Comm;
  int myrank; //on the appropriate boundary communicator                                  
  int leftNeighbor, rightNeighbor;//on the appropriate   
  int CommSize; 
  int dir=0;

  switch(comm_ID){
  case 0: //left
    myrank= vct->getRank_COMM_B_LEFT();
    Comm= vct->getCOMM_B_LEFT();
    leftNeighbor=vct->getLeftNeighbor_COMM_B_LEFT();
    rightNeighbor=vct->getRightNeighbor_COMM_B_LEFT();
    CommSize= vct->getYLEN();
    break;
  case 1: //right
    myrank= vct->getRank_COMM_B_RIGHT();
    Comm= vct->getCOMM_B_RIGHT();
    leftNeighbor=vct->getLeftNeighbor_COMM_B_RIGHT();
    rightNeighbor=vct->getRightNeighbor_COMM_B_RIGHT();
    CommSize= vct->getYLEN();
    break;
  case 2: //bottom
    myrank= vct->getRank_COMM_B_BOTTOM();
    Comm= vct->getCOMM_B_BOTTOM();
    leftNeighbor=vct->getLeftNeighbor_COMM_B_BOTTOM();
    rightNeighbor=vct->getRightNeighbor_COMM_B_BOTTOM();
    CommSize= vct->getXLEN();
    break;
  case 3: //top
    myrank= vct->getRank_COMM_B_TOP();
    Comm= vct->getCOMM_B_TOP();
    leftNeighbor=vct->getLeftNeighbor_COMM_B_TOP();
    rightNeighbor=vct->getRightNeighbor_COMM_B_TOP();
    CommSize= vct->getXLEN();
    break;
  default:
    cout <<"ERROR: wrong input in communicateParticles_BC, ...\n";
    cerr <<"ERROR: wrong input in communicateParticles_BC, ...\n";
    return -1;
  }
    
  
  communicateParticlesDIR(Comm, buffer_size,myrank,rightNeighbor,leftNeighbor,0,CommSize,1,b_Xright,b_Xleft);

  return 1;
}

#endif

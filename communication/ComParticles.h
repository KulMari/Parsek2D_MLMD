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

#endif

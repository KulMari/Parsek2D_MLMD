/***************************************************************************
    ComBasic.h  -  Library to handle Basic Communication
                       -------------------
    begin                : Fri Jun 4 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef ComBasic_H
#define ComBasic_H

#include <mpi.h>
#include <math.h>

#include "ComParser.h"
/** communicate particles along a direction **/
inline void communicateParticlesDIR(MPI_Comm CART_COMM, int buffer_size, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, double *b_right, double *b_left){
  MPI_Status status;
  double *LEN = new double[2];
  LEN[0] = XLEN, LEN[1] = YLEN;
  switch (DIR){
    case 0:
     myrank = (int) floor( (double) (myrank/(YLEN)));
     break;
    case 1:
      myrank = (int) floor( (double) (myrank) );
      break;
  }
  
  if (myrank%2==0){
            // On the boundaries send e receive only if you have periodic condition: send to X-RIGHT
            if (right_neighbor != MPI_PROC_NULL){
              if (LEN[DIR] > 1){
		MPI_Sendrecv_replace(&b_right[0],buffer_size,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, CART_COMM, &status);
		}
		else
                 swapBuffer(buffer_size,b_left,b_right);
            }
         } else {
            // On the boundaries send e receive only if you have periodic condition: send to X-LEFT
            if (left_neighbor != MPI_PROC_NULL){
              if (LEN[DIR] > 1){
                MPI_Sendrecv_replace(&b_left[0],buffer_size,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, CART_COMM, &status);
                }
		else
                swapBuffer(buffer_size,b_left,b_right);
            }
         }
   if (myrank%2==1){
            // On the boundaries send e receive only if you have periodic condition: send to X-RIGHT
            if (right_neighbor != MPI_PROC_NULL){
              if (LEN[DIR] > 1)
                MPI_Sendrecv_replace(&b_right[0],buffer_size,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, CART_COMM, &status);

            }
         } else  {
            // On the boundaries send e receive only if you have periodic condition: send to X-LEFT
            if (left_neighbor != MPI_PROC_NULL){
              if (LEN[DIR] > 1)
                MPI_Sendrecv_replace(&b_left[0],buffer_size,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, CART_COMM, &status);
              }
         }
  delete[] LEN;
}
/** communicate ghost along a direction (DIR=0 --> x-direction ; DIR=1 --> y direction) **/
inline void communicateGhostFace(MPI_Comm CART_COMM, int b_len, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, double *ghostRightFace, double *ghostLeftFace){

  MPI_Status status;
  double *LEN = new double[2];
  int rankF,ierr;
  LEN[0] = XLEN, LEN[1] = YLEN;
  switch (DIR){
    case 0:
        rankF = (int) floor( (double) (myrank/(YLEN)) ); 
        break;
    case 1:
        rankF = myrank;
        break;
  }
  if (rankF%2==0 && right_neighbor != MPI_PROC_NULL && LEN[DIR] > 1){ // SEND-RECEIVE RIGHT{
              ierr = MPI_Sendrecv_replace(&ghostRightFace[0],b_len,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, CART_COMM, &status);
                }
  else if (rankF%2== 1 && left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1){ // SEND-RECEIVE LEFT
              ierr = MPI_Sendrecv_replace(&ghostLeftFace[0],b_len,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, CART_COMM, &status);
                }

  if (rankF%2==1 && right_neighbor != MPI_PROC_NULL &&  LEN[DIR] > 1){  // SEND-RECEIVE RIGHT
              ierr = MPI_Sendrecv_replace(&ghostRightFace[0],b_len,MPI_DOUBLE,right_neighbor,1,right_neighbor,1, CART_COMM, &status);
                }
  else if (rankF%2==0 && left_neighbor != MPI_PROC_NULL  && LEN[DIR] > 1 ){   // SEND-RECEIVE LEFT
              ierr = MPI_Sendrecv_replace(&ghostLeftFace[0],b_len,MPI_DOUBLE,left_neighbor,1,left_neighbor,1, CART_COMM, &status);
  }
  // just swap the buffer if you have just 1 processor in 1 direction
  if (LEN[DIR] == 1 && right_neighbor != MPI_PROC_NULL && left_neighbor != MPI_PROC_NULL)
    swapBuffer(b_len,ghostLeftFace,ghostRightFace);


  delete[] LEN;
}


inline void communicateGhostCorner(MPI_Comm CART_COMM, int myrank, int right_neighborC, int left_neighborC, int XLEN, int YLEN, double *ghostRightCorner, double *ghostLeftCorner){
  MPI_Status status;
  int rankC, ierr;
  rankC = (int) floor( (double) (myrank/(YLEN)) );
  bool comNotDone = true;

  if (right_neighborC == myrank && left_neighborC == myrank){
     swapBuffer(ghostLeftCorner,ghostRightCorner);
     comNotDone = false;
  }
  // if it's possible communicate corners

  if (rankC%2==0 && right_neighborC != MPI_PROC_NULL && comNotDone){
              ierr=MPI_Sendrecv_replace(ghostRightCorner,1,MPI_DOUBLE,right_neighborC,1,right_neighborC,1, CART_COMM, &status);
  }
  else if (rankC%2== 1 && left_neighborC != MPI_PROC_NULL && comNotDone){           
              ierr=MPI_Sendrecv_replace(ghostLeftCorner,1,MPI_DOUBLE,left_neighborC,1,left_neighborC,1, CART_COMM, &status);
  }

  if (rankC%2==1 && right_neighborC != MPI_PROC_NULL && comNotDone){             
              ierr=MPI_Sendrecv_replace(ghostRightCorner,1,MPI_DOUBLE,right_neighborC,1,right_neighborC,1, CART_COMM, &status);
  }
  else if (rankC%2==0 && left_neighborC != MPI_PROC_NULL && comNotDone){   
              ierr=MPI_Sendrecv_replace(ghostLeftCorner,1,MPI_DOUBLE,left_neighborC,1,left_neighborC,1, CART_COMM, &status);
  }
  
}

inline void communicateGhostCorner_forcolumn(MPI_Comm CART_COMM, int myrank, int right_neighborC, int left_neighborC, int XLEN, int YLEN, double *ghostRightCorner, double *ghostLeftCorner){
  MPI_Status status;
  int rankC, ierr;
  rankC = (int) floor( (double) (myrank) );
  bool comNotDone = true;

  if (right_neighborC == myrank && left_neighborC == myrank){
     swapBuffer(ghostLeftCorner,ghostRightCorner);
     comNotDone = false;
  }
  // if it's possible communicate corners

  if (rankC%2==0 && right_neighborC != MPI_PROC_NULL && comNotDone){
              ierr=MPI_Sendrecv_replace(ghostRightCorner,1,MPI_DOUBLE,right_neighborC,1,right_neighborC,1, CART_COMM, &status);
  }
  else if (rankC%2== 1 && left_neighborC != MPI_PROC_NULL && comNotDone){           
              ierr=MPI_Sendrecv_replace(ghostLeftCorner,1,MPI_DOUBLE,left_neighborC,1,left_neighborC,1, CART_COMM, &status);
  }

  if (rankC%2==1 && right_neighborC != MPI_PROC_NULL && comNotDone){             
              ierr=MPI_Sendrecv_replace(ghostRightCorner,1,MPI_DOUBLE,right_neighborC,1,right_neighborC,1, CART_COMM, &status);
  }
  else if (rankC%2==0 && left_neighborC != MPI_PROC_NULL && comNotDone){   
              ierr=MPI_Sendrecv_replace(ghostLeftCorner,1,MPI_DOUBLE,left_neighborC,1,left_neighborC,1, CART_COMM, &status);
  }
  
}














#endif

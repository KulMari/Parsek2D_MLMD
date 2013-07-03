/***********************************************************************************
BcFields.h  -  series of functions to set the boundary condition for Fields
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta,  David Burgess
***************************************************************************/


#ifndef BcFields_H
#define BcFields_H

/**
* 
* series of functions to set the boundary condition for Fields 
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/

/** set the boundary condition on faces in 1D*/
inline void BCface(int nx,double ***vector,int bcFaceXright, int bcFaceXleft,VirtualTopology *vct){
   // XLEFT
   if (vct->getXleft_neighbor()==MPI_PROC_NULL){
      switch(bcFaceXleft){
        case 0:   // Dirichilet = 0 Second Order  attention --> BC are applied here also on ghost corners, and modified below 
	                    //  for corners of the whole system (?!?)
            vector[0][0][0] = -vector[1][0][0];
        break;

      case 1:  // Dirichilet = 0 First Order
        
            vector[0][0][0] = 0;
        break;

      case 2:   // Neumann = 0 First Order
      
            vector[0][0][0] = vector[1][0][0];

        break;


      }

   }
   // XRIGHT
   if (vct->getXright_neighbor()==MPI_PROC_NULL){
     switch(bcFaceXright){

      case 0:    // Dirichilet = 0 Second Order
       
           vector[nx-1][0][0] = -vector[nx-2][0][0];
        break;

      case 1:   // Dirichilet = 0 First Order
        
           vector[nx-1][0][0] = 0;
        break;

      case 2:    // Neumann = 0 First Order
        
           vector[nx-1][0][0] = vector[nx-2][0][0];
        break;

      }

   }
  
   
   
	   
}

/** set the boundary condition on faces in 2D*/
inline void BCface(int nx, int ny, double ***vector,int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct){
   // XLEFT
   if (vct->getXleft_neighbor()==MPI_PROC_NULL){
      switch(bcFaceXleft){

      case 0:   // Dirichilet = 0 Second Order  attention --> BC are applied here also on ghost corners, and modified below 
	for(int i=0; i < ny; i++)                         //  for corners of the whole system (?!?)
            vector[0][i][0] = -vector[1][i][0];
        break;

      case 1:  // Dirichilet = 0 First Order
        for(int i=0; i < ny; i++)
            vector[0][i][0] = 0;  
        break;

      case 2:   // Neumann = 0 First Order
        for(int i=0; i < ny; i++)
            vector[0][i][0] = vector[1][i][0];

        break;


      }

   }
   // XRIGHT
   if (vct->getXright_neighbor()==MPI_PROC_NULL){
     switch(bcFaceXright){

      case 0:    // Dirichilet = 0 Second Order
        for(int i=0; i < ny; i++)
           vector[nx-1][i][0] = -vector[nx-2][i][0];
        break;

      case 1:   // Dirichilet = 0 First Order
         for(int i=0; i < ny; i++)
           vector[nx-1][i][0] = 0; 
        break;

      case 2:    // Neumann = 0 First Order
        for(int i=0; i < ny; i++)
           vector[nx-1][i][0] = vector[nx-2][i][0];
        break;

      }

   }
   // YLEFT
   if (vct->getYleft_neighbor()==MPI_PROC_NULL){
     switch(bcFaceYleft){

      case 0:  // Dirichilet = 0 Second Order
      for(int i=0; i < nx; i++)
          vector[i][0][0] = -vector[i][1][0];
      break;

      case 1: // Dirichilet = 0 First Order
      for(int i=0; i < nx; i++)
          vector[i][0][0] = 0;       
      break;

      case 2:  // Neumann = 0 First Order
      for(int i=0; i < nx; i++)
          vector[i][0][0] = vector[i][1][0];
      break;

     }
   }
   // YRIGHT
   if (vct->getYright_neighbor()==MPI_PROC_NULL){
      switch(bcFaceYright){

      case 0: // Dirichilet = 0 Second Order
      for(int i=0; i < nx; i++)
          vector[i][ny-1][0] = -vector[i][ny-2][0];
      break;

      case 1: // Dirichilet = 0 First Order
      for(int i=0; i < nx; i++)
          vector[i][ny-1][0] = 0;   
      break;

      case 2: // Neumann = 0 First Order
      for(int i=0; i < nx; i++)
          vector[i][ny-1][0] = vector[i][ny-2][0];
      break;

      }

   }
   
   // CORNER XRIGHT + YRIGHT
   if (vct->getXright_neighbor()==MPI_PROC_NULL & vct->getYright_neighbor()==MPI_PROC_NULL){
	   vector[nx-1][ny-1][0]=0.5*(vector[nx-1][ny-2][0]+vector[nx-2][ny-1][0]);}      
   
   // CORNER XRIGHT + YLEFT	
     if (vct->getXright_neighbor()==MPI_PROC_NULL & vct->getYleft_neighbor()==MPI_PROC_NULL){
	   vector[nx-1][0][0]=0.5*(vector[nx-1][1][0]+vector[nx-2][0][0]);}      
   
   // CORNER XLEFT + YLEFT	
   if (vct->getXleft_neighbor()==MPI_PROC_NULL & vct->getYleft_neighbor()==MPI_PROC_NULL){
	   vector[0][0][0]=0.5*(vector[0][1][0]+vector[1][0][0]);}      

   // CORNER XLEFT + YRIGHT	
   if (vct->getXleft_neighbor()==MPI_PROC_NULL & vct->getYright_neighbor()==MPI_PROC_NULL){
	   vector[0][ny-1][0]=0.5*(vector[0][ny-2][0]+vector[1][ny-1][0]);}      
	   
}


#endif

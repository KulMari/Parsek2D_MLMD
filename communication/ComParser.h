/***********************************************************************************
ComParser.h  -  series of functions to 1) swap buffers in 1D cases 2)prepare ghost cells in communication buffers 3)parse the communication to ghost cells
                            -------------------
 developers: Stefano Markidis, Giovanni Lapenta
***************************************************************************/

#ifndef ComParser_H
#define ComParser_H

#include "../utility/Alloc.h"
/**
* 
* Prepare the ghost cells with nodes values, communicate the ghost cells, and replace the values on the ghost cells
* @date Fri Jun 4 2008
* @author Stefano Markidis, Giovanni Lapenta
* @version 4.0
*
*/


/** swap the communication buffer instead of communicating, useful if you have 1 processor in one direction*/
inline void swapBuffer(int buffer_size, double *b_left,double *b_right){
   double *temp = new double[buffer_size];
  
   for (register int i=0;i < buffer_size;i++){
      temp[i] = b_left[i];
      b_left[i] = b_right[i];
      b_right[i] = temp[i];
   }
  
   delete[] temp;
}
/** swap the communication buffer instead of communicating, useful if you have 1 processor in one direction*/
inline void swapBuffer(double *b_left,double *b_right){
  double temp = *b_left;
  *b_left = *b_right;
  *b_right = temp;

}/** swap the ghost cell face instead of communicating, useful if you have 1 processor in one direction*/
inline void swapGhostFace(int n1, int n2, double **ghostFaceLeft, double **ghostFaceRight){
   double **temp = newArr(double,n1,n2);
   for (register int i=0;i < n1;i++){
     for (register int j=0;j < n2;j++){
        temp[i][j] = ghostFaceLeft[i][j];
        ghostFaceLeft[i][j] = ghostFaceRight[i][j];
        ghostFaceRight[i][j] = temp[i][j];
      }
   }
   delArr(temp,n1);
}

/** NODE SOLVER: prepare node ghost cells in a buffer for communication - X direction*/
inline void makeNodeFaceX(int nxn, int nyn, double ***vector, double *ghostXrightFace, double *ghostXleftFace, VirtualTopology *vct){
   if (vct->getXleft_neighbor()!=MPI_PROC_NULL){
       for (register int j=1; j < nyn-1; j++){
            ghostXleftFace[j-1]  = vector[2][j][0];
        }
    } else {
       for (register int j=1; j < nyn-1; j++){
            ghostXleftFace[j-1]  = vector[0][j][0];
        }
    }

   if (vct->getXright_neighbor()!=MPI_PROC_NULL){
       for (register int j=1; j < nyn-1; j++){
            ghostXrightFace[j-1] = vector[nxn-3][j][0];
        }
    } else {
       for (register int j=1; j < nyn-1; j++){
            ghostXrightFace[j-1] = vector[nxn-1][j][0];
        }
    }
}

/** prepare center  ghost cells in a buffer for communication - X direction*/
inline void makeCenterFaceX(int nxc, int nyc, double ***vector, double *ghostXrightFace, double *ghostXleftFace, VirtualTopology *vct){
   if (vct->getXleft_neighbor()!= MPI_PROC_NULL){
       for (register int j=1; j < nyc-1; j++){
            	ghostXleftFace[j-1]  = vector[1][j][0];
       }
   } else {
       for (register int j=1; j < nyc-1; j++){
            	ghostXleftFace[j-1]  = vector[0][j][0];
       }

   }
   if (vct->getXright_neighbor()!= MPI_PROC_NULL){
       for (register int j=1; j < nyc-1; j++){
            	ghostXrightFace[j-1] = vector[nxc-2][j][0];
       }
   } else {
       for (register int j=1; j < nyc-1; j++){
            	ghostXrightFace[j-1] = vector[nxc-1][j][0];
       }
   }
}

/** NODE SOLVER prepare node ghost cells in a buffer for communication - Y direction*/
inline void makeNodeFaceY(int nxn, int nyn, double ***vector, double *ghostYrightFace, double *ghostYleftFace,VirtualTopology *vct){
   // YLEFT
   if (vct->getYleft_neighbor()== MPI_PROC_NULL){
      for (register int j=1; j < nxn-1; j++)
        ghostYleftFace[j-1] = vector[j][0][0];
	} else {
	   for (register int j=1; j < nxn-1; j++)
        ghostYleftFace[j-1] = vector[j][2][0];
	}
	// YRIGHT
   if (vct->getYright_neighbor()== MPI_PROC_NULL){
     for (register int j=1; j < nxn-1; j++)
	   ghostYrightFace[j-1]  = vector[j][nyn-1][0];
   } else {
      for (register int j=1; j < nxn-1; j++)
	   ghostYrightFace[j-1]  = vector[j][nyn-3][0];
   }
}

/** prepare center ghost cells in a buffer for communication - Y direction*/
inline void makeCenterFaceY(int nxc, int nyc, double ***vector, double *ghostYrightFace, double *ghostYleftFace,VirtualTopology *vct){
   if(vct->getYleft_neighbor()== MPI_PROC_NULL){
       for (register int j=1; j < nxc-1; j++){
            	ghostYleftFace[j-1] = vector[j][0][0];
       }
   } else {
       for (register int j=1; j < nxc-1; j++){
            	ghostYleftFace[j-1] = vector[j][1][0];
       }
   }
   if(vct->getYright_neighbor()== MPI_PROC_NULL){
       for (register int j=1; j < nxc-1; j++){
            	ghostYrightFace[j-1]  = vector[j][nyc-1][0];
       }
   } else {
       for (register int j=1; j < nxc-1; j++){
            	ghostYrightFace[j-1]  = vector[j][nyc-2][0];
       }
   }
}

/** NODE SPECIES: prepare node ghost cells in a buffer for communication - X direction*/
inline void makeNodeFaceX(int nxn, int nyn, double ****vector,  int ns, double *ghostXrightFace, double *ghostXleftFace,VirtualTopology *vct){
   if (vct->getXleft_neighbor()!=MPI_PROC_NULL){
       for (register int j=1; j < nyn-1; j++){
            ghostXleftFace[j-1]  = vector[ns][2][j][0];
        }
    } else {
       for (register int j=1; j < nyn-1; j++){
            ghostXleftFace[j-1]  = vector[ns][0][j][0];
        }
    }

   if (vct->getXright_neighbor()!=MPI_PROC_NULL){
       for (register int j=1; j < nyn-1; j++){
            ghostXrightFace[j-1] = vector[ns][nxn-3][j][0];
        }
    } else {
       for (register int j=1; j < nyn-1; j++){
            ghostXrightFace[j-1] = vector[ns][nxn-1][j][0];
        }
    } 
}
/** SPECIES: prepare center  ghost cells in a buffer for communication - X direction*/
inline void makeCenterFaceX(int nxc, int nyc, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace,VirtualTopology *vct){
   if (vct->getXleft_neighbor()!= MPI_PROC_NULL){
       for (register int j=1; j < nyc-1; j++){
            	ghostXleftFace[j-1]  = vector[ns][1][j][0];
       }
   } else {
       for (register int j=1; j < nyc-1; j++){
            	ghostXleftFace[j-1]  = vector[ns][0][j][0];
       }

   }
   if (vct->getXright_neighbor()!= MPI_PROC_NULL){
       for (register int j=1; j < nyc-1; j++){
            	ghostXrightFace[j-1] = vector[ns][nxc-2][j][0];
       }
   } else {
       for (register int j=1; j < nyc-1; j++){
            	ghostXrightFace[j-1] = vector[ns][nxc-1][j][0];
       }
   }

}

/** SPECIES: prepare node ghost cells in a buffer for communication - Y direction*/
inline void makeNodeFaceY(int nxn, int nyn, double ****vector, int ns,double *ghostYrightFace, double *ghostYleftFace,VirtualTopology *vct){
   // YLEFT
   if (vct->getYleft_neighbor()== MPI_PROC_NULL){
      for (register int j=1; j < nxn-1; j++)
        ghostYleftFace[j-1] = vector[ns][j][0][0];
	} else {
	   for (register int j=1; j < nxn-1; j++)
        ghostYleftFace[j-1] = vector[ns][j][2][0];
	}
	// YRIGHT
   if (vct->getYright_neighbor()== MPI_PROC_NULL){
     for (register int j=1; j < nxn-1; j++)
	   ghostYrightFace[j-1]  = vector[ns][j][nyn-1][0];
   } else {
      for (register int j=1; j < nxn-1; j++)
	   ghostYrightFace[j-1]  = vector[ns][j][nyn-3][0];
   }

}

/** SPECIES: prepare center ghost cells in a buffer for communication - Y direction*/
inline void makeCenterFaceY(int nxc, int nyc, double ****vector, int ns, double *ghostYrightFace, double *ghostYleftFace,VirtualTopology *vct){
   if(vct->getYleft_neighbor()== MPI_PROC_NULL){
       for (register int j=1; j < nxc-1; j++){
            	ghostYleftFace[j-1] = vector[ns][j][0][0];
       }
   } else {
       for (register int j=1; j < nxc-1; j++){
            	ghostYleftFace[j-1] = vector[ns][j][1][0];
       }
   }
   if(vct->getYright_neighbor()== MPI_PROC_NULL){
       for (register int j=1; j < nxc-1; j++){
            	ghostYrightFace[j-1]  = vector[ns][j][nyc-1][0];
       }
   } else {
       for (register int j=1; j < nxc-1; j++){
            	ghostYrightFace[j-1]  = vector[ns][j][nyc-2][0];
       }
   }

}

/** SOLVER NODES: prepare ghost Node corners for communication*/
inline void makeNodeCorner(int nxn, int nyn, double ***vector, double *ghostXrightYrightCorner, double *ghostXleftYrightCorner, double *ghostXrightYleftCorner, double *ghostXleftYleftCorner,VirtualTopology *vct){

  //upper left
  if (vct->getXleftYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[2][nyn-3][0];}
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[0][nyn-3][0];}
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[2][nyn-1][0];}
  else
    {*ghostXleftYrightCorner = vector[0][nyn-1][0];}
  
  // upper right
  if (vct->getXrightYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[nxn-3][nyn-3][0];} 
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[nxn-1][nyn-3][0];}
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[nxn-3][nyn-1][0];}
  else
    {*ghostXrightYrightCorner = vector[nxn-1][nyn-1][0];}


  // lower left
  if (vct->getXleftYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[2][2][0];}
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[0][2][0];}	
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[2][0][0];}
  else
    {*ghostXleftYleftCorner = vector[0][0][0];}


  // lower right
  if (vct->getXrightYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[nxn-3][2][0];}
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[nxn-1][2][0];}
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)                                                                          
    {*ghostXrightYleftCorner = vector[nxn-3][0][0];}
  else
    {*ghostXrightYleftCorner = vector[nxn-1][0][0];}
}

/** SOLVER NODES: prepare ghost Node corners for communication for OS particles*/
inline void makeNodeCornerInterpOS(int nxn, int nyn, double ****vector, int ns, double* ghostXrightYrightCorner, double* ghostXleftYrightCorner, double* ghostXrightYleftCorner, double* ghostXleftYleftCorner,VirtualTopology *vct){
  //upper left
  if (! (vct->getXleft_neighbor()== MPI_PROC_NULL) != ! (vct->getYright_neighbor()== MPI_PROC_NULL))
    {
      //cout << "R" << vct->getCartesian_rank_COMMTOTAL() << "in makeNodeCornerInterpOS\n"; 
      if (vct->getXleft_neighbor()== MPI_PROC_NULL)
	{*ghostXleftYrightCorner= vector[ns][0][nyn-2][0];}
      else if (vct->getYright_neighbor()== MPI_PROC_NULL)
	{*ghostXleftYrightCorner= vector[ns][1][nyn-1][0];}
    }
  else{*ghostXleftYrightCorner=0.0;}
  //upper right                                                                                                                                             
  if (! (vct->getXright_neighbor()== MPI_PROC_NULL) != ! (vct->getYright_neighbor()== MPI_PROC_NULL ))                                                     
    {//cout << "R" << vct->getCartesian_rank_COMMTOTAL()<< "in makeNodeCornerInterpOS\n";
      if (vct->getXright_neighbor()== MPI_PROC_NULL)                                                                                                        
        {*ghostXrightYrightCorner= vector[ns][nxn-1][nyn-2][0];}                                                                                         
      else if (vct->getYright_neighbor()== MPI_PROC_NULL)                                                                                                   
        {*ghostXrightYrightCorner= vector[ns][nxn-2][nyn-1][0];}                                                                                          
    }                                                                                                                                                       
  else{*ghostXrightYrightCorner=0.0;}                                                                                                                       
  // lower left
  if (! (vct->getXleft_neighbor()== MPI_PROC_NULL) != ! (vct->getYleft_neighbor()== MPI_PROC_NULL ))
    {//cout << "R" << vct->getCartesian_rank_COMMTOTAL()<< "in makeNodeCornerInterpOS\n";
      if (vct->getXleft_neighbor()== MPI_PROC_NULL)
	{*ghostXleftYleftCorner = vector[ns][0][1][0];}
      else if (vct->getYleft_neighbor()== MPI_PROC_NULL)
	{*ghostXleftYleftCorner = vector[ns][1][0][0];}
    }else{ghostXleftYleftCorner[0]=0.0;}
    // lower right                                                                                                                                        
  if (! (vct->getXright_neighbor()== MPI_PROC_NULL) != !( vct->getYleft_neighbor()== MPI_PROC_NULL )) 
    {//cout << "R" << vct->getCartesian_rank_COMMTOTAL()<< "in makeNodeCornerInterpOS\n";
      if (vct->getXright_neighbor()== MPI_PROC_NULL)
	{*ghostXrightYleftCorner = vector[ns][nxn-1][1][0];}
      else if(vct->getYleft_neighbor()== MPI_PROC_NULL)
	{*ghostXrightYleftCorner = vector[ns][nxn-2][0][0];}
    }else{*ghostXrightYleftCorner =0.0;}
  
}

/** prepare ghost center corners for communication*/
inline void makeCenterCorner(int nxc, int nyc, double ***vector, double *ghostXrightYrightCorner, double *ghostXleftYrightCorner, double *ghostXrightYleftCorner, double *ghostXleftYleftCorner,VirtualTopology *vct){

  //upper left     
  if (vct->getXleftYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[1][nyc-2][0];}
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[0][nyc-2][0];}
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[1][nyc-1][0];} 
  else
    {*ghostXleftYrightCorner = vector[0][nyc-1][0];}

  // upper right                                                                                                                                 
  if (vct->getXrightYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[nxc-2][nyc-2][0];}
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[nxc-1][nyc-2][0];}
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[nxc-2][nyc-1][0];}
  else
    {*ghostXrightYrightCorner = vector[nxc-1][nyc-1][0];}

  // lower left                                                                                                                                 
  if (vct->getXleftYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[1][1][0];}
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[0][1][0];}
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[1][0][0];}
  else
    {*ghostXleftYleftCorner = vector[0][0][0];}

  // lower right                                                                                                                                
  if (vct->getXrightYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[nxc-2][1][0];}
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[nxc-1][1][0];}
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[nxc-2][0][0];}
  else
    {*ghostXrightYleftCorner = vector[nxc-1][0][0];}

}

/** SPECIES: prepare ghost Node corners for communication*/
inline void makeNodeCorner(int nxn, int nyn, double ****vector, int ns, double* ghostXrightYrightCorner, double* ghostXleftYrightCorner, double* ghostXrightYleftCorner, double* ghostXleftYleftCorner,VirtualTopology *vct){
  
  //upper left
  if (vct->getXleftYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[ns][2][nyn-3][0];}
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[ns][0][nyn-3][0];}
  // new
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[ns][2][nyn-1][0];}
  // end new
  else
    {*ghostXleftYrightCorner = vector[ns][0][nyn-1][0];}
  
  // upper right
  if (vct->getXrightYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[ns][nxn-3][nyn-3][0];} 
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[ns][nxn-1][nyn-3][0];}
  // new
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[ns][nxn-3][nyn-1][0];}
  // end new
  else
    {*ghostXrightYrightCorner = vector[ns][nxn-1][nyn-1][0];}


  // lower left
  if (vct->getXleftYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[ns][2][2][0];}
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[ns][0][2][0];}
  // new
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[ns][2][0][0];}
  // end new
  else
    {*ghostXleftYleftCorner = vector[ns][0][0][0];}

  // lower right
  if (vct->getXrightYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[ns][nxn-3][2][0];}
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[ns][nxn-1][2][0];}
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)                                                                          
    {*ghostXrightYleftCorner = vector[ns][nxn-3][0][0];}
  else
    {*ghostXrightYleftCorner = vector[ns][nxn-1][0][0];}

}

/** SPECIES: prepare ghost Center corners for communication*/
inline void makeCenterCorner(int nxc, int nyc, double ****vector, int ns, double *ghostXrightYrightCorner, double *ghostXleftYrightCorner, double *ghostXrightYleftCorner, double *ghostXleftYleftCorner,VirtualTopology *vct){
 
  //upper left                                                                             
  if (vct->getXleftYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[ns][1][nyc-2][0];}
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[ns][0][nyc-2][0];}
  // new                                                                                                                                        
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[ns][1][nyc-1][0];}
  // end new                                                                                                                                    
  else
    {*ghostXleftYrightCorner = vector[ns][0][nyc-1][0];}

  // upper right                                                                                                                                 
  if (vct->getXrightYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[ns][nxc-2][nyc-2][0];}
  else if (vct->getYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[ns][nxc-1][nyc-2][0];}
  // new                                                                                                                                        
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYrightCorner = vector[ns][nxc-2][nyc-1][0];}
  // end new                                                                                                                                    
  else
    {*ghostXrightYrightCorner = vector[ns][nxc-1][nyc-1][0];}

  // lower left                                                                                                                                 
  if (vct->getXleftYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[ns][1][1][0];}
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[ns][0][1][0];}
  // new                                                                                                                                        
  else if (vct->getXleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYleftCorner = vector[ns][1][0][0];}
  // end new                                                                                                                                    
  else
    {*ghostXleftYleftCorner = vector[ns][0][0][0];}

  // lower right                                                                                                                                
  if (vct->getXrightYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[ns][nxc-2][1][0];}
  /*else if (vct->getXright_neighbor()!= MPI_PROC_NULL)                                                                                         
    {*ghostXrightYleftCorner = vector[nxn-1][2][0];}*/
  else if (vct->getYleft_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[ns][nxc-1][1][0];}
  else if (vct->getXright_neighbor()!= MPI_PROC_NULL)
    {*ghostXrightYleftCorner = vector[ns][nxc-2][0][0];}
  else
    {*ghostXrightYleftCorner = vector[ns][nxc-1][0][0];}


}


/** SPECIES: prepare Node ghost cells for interpolation  when there is periodicity*/
inline void makeNodeFaceInterp(int nx, int ny, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace,VirtualTopology *vct){
//X DIRECTION
   //LEFT
   if (vct->getXleft_neighbor()!=MPI_PROC_NULL){
       for (register int j=1; j < ny-1; j++){
            ghostXleftFace[j-1]  = vector[ns][1][j][0];
        }
    } else {
       for (register int j=1; j < ny-1; j++){
            ghostXleftFace[j-1]  = vector[ns][0][j][0];
        }
    }
   //RIGHT
   if (vct->getXright_neighbor()!=MPI_PROC_NULL){
       for (register int j=1; j < ny-1; j++){
            ghostXrightFace[j-1] = vector[ns][nx-2][j][0];
        }
    } else {
       for (register int j=1; j < ny-1; j++){
            ghostXrightFace[j-1] = vector[ns][nx-1][j][0];
        }
    } 
// Y DIRECTION
   //LEFT
   if (vct->getYleft_neighbor()== MPI_PROC_NULL){
      for (register int j=1; j < nx-1; j++)
        ghostYleftFace[j-1] = vector[ns][j][0][0];
	} else {
	   for (register int j=1; j < nx-1; j++)
        ghostYleftFace[j-1] = vector[ns][j][1][0];
	}
   //RIGHT
   if (vct->getYright_neighbor()== MPI_PROC_NULL){
     for (register int j=1; j < nx-1; j++)
	   ghostYrightFace[j-1]  = vector[ns][j][ny-1][0];
   } else {
      for (register int j=1; j < nx-1; j++)
	   ghostYrightFace[j-1]  = vector[ns][j][ny-2][0];
   }
}

/** SPECIES: prepare Corner node ghost cells for interpolation  when there is periodicity*/
inline void makeNodeCornerInterp(int nxn, int nyn, double ****vector, int ns, double* ghostXrightYrightCorner, double* ghostXleftYrightCorner, double* ghostXrightYleftCorner, double* ghostXleftYleftCorner,VirtualTopology *vct){
  
  //upper left                                                                                                                                                      
  if (vct->getXleftYright_neighbor()!= MPI_PROC_NULL)
    {*ghostXleftYrightCorner = vector[ns][1][nyn-2][0];} 
  // upper right                                                                                                                                                    
  if (vct->getXrightYright_neighbor()!= MPI_PROC_NULL)                                                                         
    {*ghostXrightYrightCorner = vector[ns][nxn-2][nyn-2][0];}  
  // lower left                                                                                                                                                     
  if (vct->getXleftYleft_neighbor()!= MPI_PROC_NULL)                                                                           
    {*ghostXleftYleftCorner = vector[ns][1][1][0];}
  // lower right                                                                                                                                                    
  if (vct->getXrightYleft_neighbor()!= MPI_PROC_NULL)                                                                          
    {*ghostXrightYleftCorner = vector[ns][nxn-2][1][0];}
}


/** take the communication buffer and parse it to ghost cell - DIRECTION X*/
inline void parseFaceX(int nxn, int nyn, double ***vector, double *ghostXrightFace, double *ghostXleftFace){
  // parse ghost: ghostXrightFace(last X column) and ghostXleftFace (first X column)
    for (register int j=1; j < nyn-1; j++){
        vector[nxn-1][j][0] = ghostXrightFace[j-1];
        vector[0][j][0]    = ghostXleftFace[j-1];
    }

}
/** take the communication buffer and parse it to ghost cell - DIRECTION Y*/
inline void parseFaceY(int nxn, int nyn, double ***vector,  double *ghostYrightFace, double *ghostYleftFace){
  // parse ghost: ghostYrightFace(last Y column) and ghostYleftFace (first Y column)
    for (register int j=1; j < nxn-1; j++){
        vector[j][nyn-1][0] = ghostYrightFace[j-1];
        vector[j][0][0] = ghostYleftFace[j-1];
    }

}


/** SPECIES: take the communication buffer and parse it to ghost cell - DIRECTION X*/
inline void parseFaceX(int nxn, int nyn, double ****vector, int ns, double *ghostXrightFace, double *ghostXleftFace){
  // parse ghost: ghostXrightFace(last X column) and ghostXleftFace (first X column)
	for (register int j=1; j < nyn-1; j++){
		vector[ns][nxn-1][j][0] = ghostXrightFace[j-1];
		vector[ns][0][j][0]     = ghostXleftFace[j-1];
	}

}
/** SPECIES: take the communication buffer and parse it to ghost cell - DIRECTION Y*/
inline void parseFaceY(int nxn, int nyn, double ****vector,  int ns, double *ghostYrightFace, double *ghostYleftFace){
  // parse ghost: ghostYrightFace(last Y column) and ghostYleftFace (first Y column)
	for (register int j=1; j < nxn-1; j++){
		vector[ns][j][nyn-1][0] = ghostYrightFace[j-1];
		vector[ns][j][0][0] = ghostYleftFace[j-1];
	}

}


/** SPECIES: add ghost cells values for corner values (FOR INTERPOLATION)*/
inline void addGhostCorner(int nx, int ny, int ns, double ****vector, double *ghostXrightYrightCorner, double *ghostXleftYrightCorner, double *ghostXrightYleftCorner,double *ghostXleftYleftCorner, VirtualTopology *vct){
	if (vct->getXrightYright_neighbor() != MPI_PROC_NULL )	
	  vector[ns][nx-2][ny-2][0]               += *ghostXrightYrightCorner;
	if (vct->getXleftYright_neighbor() != MPI_PROC_NULL )	
	  vector[ns][1][ny-2][0]    		     += *ghostXleftYrightCorner;
	if (vct->getXrightYleft_neighbor() != MPI_PROC_NULL )	
	  vector[ns][nx-2][1][0]           	     += *ghostXrightYleftCorner;
	if (vct->getXleftYleft_neighbor() != MPI_PROC_NULL )	
	  vector[ns][1][1][0]                        += *ghostXleftYleftCorner;
}

/** SPECIES: add ghost cells values for corner values (FOR INTERPOLATION)*/
inline void addGhostCornerOS(int nx, int ny, int ns, double ****vector, double *ghostXrightYrightCorner, double *ghostXleftYrightCorner, double *ghostXrightYleftCorner,double *ghostXleftYleftCorner, VirtualTopology *vct){
  
  if (! (vct->getXleft_neighbor()== MPI_PROC_NULL) != ! (vct->getYright_neighbor()== MPI_PROC_NULL))
    {
      if (vct->getXleft_neighbor()== MPI_PROC_NULL)
	{vector[ns][0][ny-2][0]+= *ghostXleftYrightCorner;}
      else if (vct->getYright_neighbor()== MPI_PROC_NULL)
	{vector[ns][1][ny-1][0]+= *ghostXleftYrightCorner;}
    }

  //upper right                                                                                                                                             
  if (! (vct->getXright_neighbor()== MPI_PROC_NULL) != ! (vct->getYright_neighbor()== MPI_PROC_NULL ))           
    {
      if (vct->getXright_neighbor()== MPI_PROC_NULL)
	{vector[ns][nx-1][ny-2][0]+=*ghostXrightYrightCorner;}
      else if (vct->getYright_neighbor()== MPI_PROC_NULL)
	{vector[ns][nx-2][ny-1][0]+=*ghostXrightYrightCorner;}
    }
  // lower left
  if (! (vct->getXleft_neighbor()== MPI_PROC_NULL) != ! (vct->getYleft_neighbor()== MPI_PROC_NULL ))
    {
      if (vct->getXleft_neighbor()== MPI_PROC_NULL)
	{vector[ns][0][1][0]+=*ghostXleftYleftCorner;}
      else if (vct->getYleft_neighbor()== MPI_PROC_NULL)
	{vector[ns][1][0][0]+=*ghostXleftYleftCorner;}
    }
  // lower right                                                                                                                                        
  if (! (vct->getXright_neighbor()== MPI_PROC_NULL) != !( vct->getYleft_neighbor()== MPI_PROC_NULL )) 
    {
      if (vct->getXright_neighbor()== MPI_PROC_NULL)
	{vector[ns][nx-1][1][0]+=*ghostXrightYleftCorner;}
      else if(vct->getYleft_neighbor()== MPI_PROC_NULL)
	{vector[ns][nx-2][0][0]+=*ghostXrightYleftCorner;}
    }
}
/** take the communication buffer and parse it to corner cells*/
inline void parseCorner(int nxn, int nyn, double ***vector, double *ghostXrightYrightCorner, double *ghostXleftYrightCorner, double *ghostXrightYleftCorner, double *ghostXleftYleftCorner){
	vector[nxn-1][nyn-1][0]         = *ghostXrightYrightCorner;
	vector[0][nyn-1][0]             = *ghostXleftYrightCorner;
	vector[nxn-1][0][0]             = *ghostXrightYleftCorner;
	vector[0][0][0]                 = *ghostXleftYleftCorner;

}

/** SPECIES: take the communication buffer and parse it to corner cells*/
inline void parseCorner(int nxn, int nyn, double ****vector, int ns,double *ghostXrightYrightCorner, double *ghostXleftYrightCorner, double *ghostXrightYleftCorner, double *ghostXleftYleftCorner){
	vector[ns][nxn-1][nyn-1][0]         = *ghostXrightYrightCorner;
	vector[ns][0][nyn-1][0]             = *ghostXleftYrightCorner;
	vector[ns][nxn-1][0][0]             = *ghostXrightYleftCorner;
	vector[ns][0][0][0]                 = *ghostXleftYleftCorner;
}

/** SPECIES:  add ghost center cells values for corner values (FOR INTERPOLATION)*/
inline void addGhostFace(int nx, int ny,int ns, double ****vector, double *ghostXrightFace, double *ghostXleftFace, double *ghostYrightFace, double *ghostYleftFace,VirtualTopology *vct){
  // add ghost values: ghostXrightFace(last X column) and ghostXleftFace (first X column)
	if (vct->getXright_neighbor() != MPI_PROC_NULL ){
		for (register int j=1; j < ny-1; j++)
			vector[ns][nx-2][j][0]   += ghostXrightFace[j-1];

	}
  
	if (vct->getXleft_neighbor() != MPI_PROC_NULL){
		for (register int j=1; j < ny-1; j++)
			vector[ns][1][j][0]      += ghostXleftFace[j-1];
   
	}
  
  // add ghost value: ghostYrightFace(last Y column) and ghostYleftFace (first Y column)
	if (vct->getYright_neighbor() != MPI_PROC_NULL){
		for (register  int j=1; j <nx-1; j++)
			vector[ns][j][ny-2][0] += ghostYrightFace[j-1];

	}
	if (vct->getYleft_neighbor() != MPI_PROC_NULL ){
		for (register  int j=1; j < nx-1; j++)
			vector[ns][j][1][0] += ghostYleftFace[j-1];
	}

}

inline void addRefinedGhostFace(int nxn, int nyn, int startingPoint, int FinalPoint,int ns, double ****vector, int DIR, double *ghostrightFace, double *ghostleftFace,VirtualTopology *vct){

  if (DIR==0)
  {
    // X DIR
    if (vct->getXright_neighbor() != MPI_PROC_NULL)
    {
      int j=0;
      while (startingPoint+j<= FinalPoint)
	{vector[ns][nxn-2][startingPoint+j][0] += ghostrightFace[j]; j++;}
    }
    if (vct->getXleft_neighbor() != MPI_PROC_NULL)
    {
      int j=0;
      while (startingPoint+j<= FinalPoint)
	{vector[ns][1][startingPoint+j][0]+= ghostleftFace[j]; j++;}
    }
  }
  else
  {// Y DIR
    int i=0;
    if (vct->getYright_neighbor() != MPI_PROC_NULL)
    {
      int i=0;
      while (startingPoint+i<= FinalPoint)
	{vector[ns][startingPoint+i][nyn-2][0]+= ghostrightFace[i]; i++;}
    }
    if (vct->getYleft_neighbor() != MPI_PROC_NULL)
    {
      int i=0;
      while (startingPoint+i<= FinalPoint)
	{vector[ns][startingPoint+i][1][0]+= ghostleftFace[i]; i++;}
    }
  }

}

/** prepare ghost Node corners for communication in 1D*/
inline void makeNodeCorner(int nxn,double ***vector, double *ghostXright, double *ghostXleft){
        *ghostXright = vector[nxn-3][0][0];
	    *ghostXleft  = vector[2][0][0];
}
/** prepare ghost center corners for communication in 1D, only for interpolation*/
inline void makeCenterCornerInterp(int nxn,double ****vector, int ns, double *ghostXright, double *ghostXleft){
        *ghostXright = vector[ns][nxn-2][0][0];
	    *ghostXleft  = vector[ns][1][0][0];
}
/** prepare ghost center corners for communication in 1D*/
inline void makeCenterCorner(int nxc,double ***vector, double *ghostXright, double *ghostXleft){
        *ghostXright = vector[nxc-2][0][0];
	    *ghostXleft  = vector[1][0][0];
}
/** SPECIES: prepare ghost Node corners for communication in 1D*/
inline void makeNodeCorner(int nxn,double ****vector, int ns, double* ghostXright,double* ghostXleft){
	*ghostXright = vector[ns][nxn-2][0][0];
	*ghostXleft  = vector[ns][1][0][0];
}
/** SPECIES: prepare ghost Center corners for communication in 1D*/
inline void makeCenterCorner(int nxc, double ****vector, int ns, double *ghostXright, double *ghostXleft){
	*ghostXright = vector[ns][nxc-2][0][0];
	*ghostXleft  = vector[ns][1][0][0];
}
/** SPECIES: prepare Corner node ghost cells for interpolation  in 1D*/
inline void makeNodeCornerInterp(int nxn,double ****vector, int ns, double* ghostXright, double* ghostXleft){
    *ghostXright = vector[ns][nxn-2][0][0];
	*ghostXleft  = vector[ns][1][0][0];
}
/** SPECIES: add ghost cells values for corner values (FOR INTERPOLATION) in 1D*/
inline void addGhostCorner(int nx,int ns,double ****vector,double *ghostXright,double *ghostXleft,VirtualTopology *vct){
	if (vct->getXright_neighbor() != MPI_PROC_NULL)	
	  vector[ns][nx-2][0][0]+= *ghostXright;
	if (vct->getXleft_neighbor() != MPI_PROC_NULL)	
	  vector[ns][1][0][0]   += *ghostXleft;
        
}
/** SPECIES: add ghost cells values for corner values (FOR INTERPOLATION) in 1D*/
inline void addGhostCornerCenter(int nx,int ns,double ****vector,double *ghostXright,double *ghostXleft,VirtualTopology *vct){
	if (vct->getXright_neighbor() != MPI_PROC_NULL )	
	  vector[ns][nx-2][0][0]+= *ghostXright;
	if (vct->getXleft_neighbor() != MPI_PROC_NULL )	
	  vector[ns][1][0][0]   += *ghostXleft;
        
}
/** take the communication buffer and parse it to corner cells in 1D*/
inline void parseCorner(int nxn,double ***vector, double *ghostXright, double *ghostXleft){
	vector[nxn-1][0][0]         = *ghostXright;
	vector[0][0][0]             = *ghostXleft;
}
/** SPECIES: take the communication buffer and parse it to corner cells in 1D*/
inline void parseCorner(int nxn,double ****vector, int ns,double *ghostXright,double *ghostXleft){
	vector[ns][nxn-1][0][0]             = *ghostXright;
	vector[ns][0][0][0]                 = *ghostXleft;
}

inline void makeNodeRefinedFaceInterp(int nxn, int nyn, int finalNode, double ****vector, int ns, double *ghostrightFace, double *ghostleftFace, int DIRFace, int startingNode, VirtualTopology * vct){

  if (DIRFace==0) // DIR X
    {
      int j=0;
      while(j+startingNode<=finalNode)
	{
	  ghostleftFace[j]  = vector[ns][1][j+startingNode][0];
	  ghostrightFace[j]  = vector[ns][nxn-2][j+startingNode][0];
	  j++;
	}
    }else // DIR Y
    {
      int i=0;
      while(i+startingNode<=finalNode)
	{
	  ghostleftFace[i]  = vector[ns][i+startingNode][1][0];
	  ghostrightFace[i]  = vector[ns][i+startingNode][nyn-2][0];
	  i++;
	}
    }
}

#endif

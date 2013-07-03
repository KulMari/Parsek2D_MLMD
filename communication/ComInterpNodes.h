/***********************************************************************************
ComInterpNodes.h  -  series of functions to communicate interpolation data among processors
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta,  David Burgess
***************************************************************************/

#ifndef ComInterpNodes_H
#define ComInterpNodes_H

#include "ComBasic.h"

/**
* 
* series of functions to communicate interpolation data among processors 
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/

/**SPECIES: communicate ghost cells and sum the contribution with a index indicating the number of species */
inline void communicateInterp(int nx, int ny, int ns, double ****vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft,double dx, double dy, VirtualTopology *vct){
  
  int startDirX=0;
  int startDirY=0;
  int finalNodeDirX=nx-1;
  int finalNodeDirY=ny-1;

  
  double ghostXrightYrightCorner, ghostXleftYrightCorner, ghostXrightYleftCorner, ghostXleftYleftCorner;
  
  
  int lenY=finalNodeDirY-startDirY+1;
  int lenX=finalNodeDirX-startDirX+1; 
  double *ghostXrightFace = new double[lenY];
  double *ghostXleftFace  = new double[lenY];
  double *ghostYrightFace = new double[lenX];
  double *ghostYleftFace  = new double[lenX];
  
  // DIR X=0, DIR Y=1
  // final node included
  makeNodeRefinedFaceInterp(nx, ny, finalNodeDirY, vector, ns, ghostXrightFace, ghostXleftFace, 0, startDirY, vct);
  makeNodeRefinedFaceInterp(nx, ny, finalNodeDirX, vector, ns, ghostYrightFace, ghostYleftFace, 1, startDirX, vct);
  
  makeNodeCornerInterp(nx,ny, vector, ns, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner, vct);
  communicateGhostFace(vct->getCART_COMM(),lenY,vct->getCartesian_rank(),vct->getXright_neighbor(),vct->getXleft_neighbor(),0,vct->getXLEN(),vct->getYLEN(),ghostXrightFace, ghostXleftFace);
  communicateGhostFace(vct->getCART_COMM(),lenX,vct->getCartesian_rank(),vct->getYright_neighbor(),vct->getYleft_neighbor(),1,vct->getXLEN(),vct->getYLEN(),ghostYrightFace, ghostYleftFace);
  
    // I just need to pass the diagonal corner, the rest is OK
  int left;
  int right;
  if(vct->getXrightYright_neighbor()!=MPI_PROC_NULL)
    {right=vct->getXrightYright_neighbor();}
  else
    {right=MPI_PROC_NULL;}
  if (vct->getXleftYleft_neighbor()!=MPI_PROC_NULL)
    {left= vct->getXleftYleft_neighbor();}
  else
    {left=MPI_PROC_NULL;}

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);
  
  if(vct->getXrightYleft_neighbor()!=MPI_PROC_NULL)
    {right=vct->getXrightYleft_neighbor();}
  else
    {right=MPI_PROC_NULL;}
  if (vct->getXleftYright_neighbor()!=MPI_PROC_NULL)
    {left= vct->getXleftYright_neighbor();}
  else
    {left=MPI_PROC_NULL;}
  
  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
  
  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()  << "AFter Communicate GhostFace: \n";
  addRefinedGhostFace(nx, ny, startDirY, finalNodeDirY, ns, vector, 0, ghostXrightFace, ghostXleftFace, vct);
  addRefinedGhostFace(nx, ny, startDirX, finalNodeDirX, ns, vector, 1, ghostYrightFace, ghostYleftFace, vct);
  
  addGhostCorner(nx,ny, ns, vector, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner, vct);
  delete[] ghostXrightFace;
  delete[] ghostXleftFace;
  delete[] ghostYrightFace;
  delete[] ghostYleftFace;  
  
}

inline void communicateInterpOS(int nx, int ny, int ns, double ****vector, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft,double dx, double dy, VirtualTopology *vct){
  
  double  ghostXrightYrightCorner;
  double  ghostXleftYrightCorner;
  double  ghostXrightYleftCorner; 
  double  ghostXleftYleftCorner;

  makeNodeCornerInterpOS(nx,ny, vector, ns, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner, vct);

  int left;
  int right;
  if(vct->getXrightYright_neighbor()!=MPI_PROC_NULL)
    {
      if (vct->getXright_neighbor()!=MPI_PROC_NULL)
	{right= vct->getXright_neighbor();}
      else if (vct->getYright_neighbor()!=MPI_PROC_NULL)
	{right= vct->getYright_neighbor();}
      else{right= MPI_PROC_NULL;}
    }
  else
    {right=MPI_PROC_NULL;}
  
  if (vct->getXleftYleft_neighbor()!=MPI_PROC_NULL)
    {
      if (vct->getXleft_neighbor()!=MPI_PROC_NULL)
	{left=vct->getXleft_neighbor();}
      else if (vct->getYleft_neighbor()!=MPI_PROC_NULL)
	{left= vct->getYleft_neighbor();}
      else
	{left=MPI_PROC_NULL;}
    }
  else
    {left=MPI_PROC_NULL;}

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);

  if(vct->getXrightYleft_neighbor()!=MPI_PROC_NULL)
    {
      if (vct->getXright_neighbor()!=MPI_PROC_NULL)
	{right=vct->getXright_neighbor();}
      else if(vct->getYleft_neighbor()!=MPI_PROC_NULL)
	{right=vct->getYleft_neighbor();}
      else
	{right=MPI_PROC_NULL;}
    }
  else
    {right=MPI_PROC_NULL;}

  if (vct->getXleftYright_neighbor()!=MPI_PROC_NULL)
    {
      if (vct->getXleft_neighbor()!=MPI_PROC_NULL)
	{left=vct->getXleft_neighbor();}
      else if(vct->getYright_neighbor()!=MPI_PROC_NULL)
	{left=vct->getYright_neighbor();}
      else
	{left=MPI_PROC_NULL;}
    }
  else
    {left=MPI_PROC_NULL;}

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
  
  addGhostCornerOS(nx,ny, ns, vector, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner, vct);

}
#endif

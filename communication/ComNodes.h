/***********************************************************************************
ComNodes.h  -  series of functions to communicate the values on nodes on processors
                            -------------------
developers: Stefano Markidis, Giovanni Lapenta
***************************************************************************/

#ifndef ComNodes_H
#define ComNodes_H

#include "ComBasic.h"
#include "../bc/BcFields.h"
#include "mpi.h"


/**
* 
* Prepare the ghost cells with nodes values, communicate the ghost cells, and replace the values on the ghost cells
* @date Fri Jun 4 2008
* @author Stefano Markidis, Giovanni Lapenta
* @version 2.0
*
*/

/** communicate ghost node values among processors in 2D*/
inline void communicateNode(int nx, int ny, double ***vector, VirtualTopology *vct){
  // allocate 4 ghost cell Faces
  double *ghostXrightFace = new double[ny-2];
  double *ghostXleftFace  = new double[ny-2];
  double *ghostYrightFace = new double[nx-2];
  double *ghostYleftFace  = new double[nx-2];
  
  // allocate 4 ghost cell corner
  double ghostXrightYrightCorner,ghostXleftYrightCorner,ghostXrightYleftCorner,ghostXleftYleftCorner;
  // prepare the values an put in a buffer
  makeNodeFaceX(nx,ny,vector,ghostXrightFace,ghostXleftFace,vct);
  makeNodeFaceY(nx,ny,vector,ghostYrightFace,ghostYleftFace,vct);
  makeNodeCorner(nx,ny, vector, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner,vct);
  //cout <<vct->getCartesian_rank()<< " make node done"<<endl;
  // communication of edges
  communicateGhostFace(vct->getCART_COMM(),(ny-2),vct->getCartesian_rank(),vct->getXright_neighbor(),vct->getXleft_neighbor(),0,vct->getXLEN(),vct->getYLEN(),ghostXrightFace, ghostXleftFace);
  communicateGhostFace(vct->getCART_COMM(),(nx-2),vct->getCartesian_rank(),vct->getYright_neighbor(),vct->getYleft_neighbor(),1,vct->getXLEN(),vct->getYLEN(),ghostYrightFace, ghostYleftFace);
  //  cout << "communicate node done"<<endl;
  // take the values from the buffer and put it back to the edge ghost cells
  parseFaceX(nx,ny,vector,ghostXrightFace,ghostXleftFace);
  parseFaceY(nx,ny,vector,ghostYrightFace,ghostYleftFace);
  // cout << "parse node done"<<endl;


  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYright_neighbor(), vct->getXleftYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);
  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYleft_neighbor(), vct->getXleftYright_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
  //cout << "interior done"<<endl;
  if (vct->getYright_neighbor() == MPI_PROC_NULL ){ //Top line, non periodic along Y
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYrightCorner);
  }
  //cout << "top line done"<<endl;
  if (vct->getYleft_neighbor() == MPI_PROC_NULL ){ //Bottom line, non periodic along Y
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYleftCorner);
  }
  //cout << "bottom line done"<<endl;
  if (vct->getXleft_neighbor() == MPI_PROC_NULL ){ //Left column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXleftYrightCorner, &ghostXleftYleftCorner);
  }
  //cout << "left column done"<<endl;
  if (vct->getXright_neighbor() == MPI_PROC_NULL ){ //Right column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXrightYleftCorner);
  }
  //cout << "communicate ghost done" << endl;
/*
  int left;
  int right;
  if(vct->getXrightYright_neighbor()!=MPI_PROC_NULL)
    {right=vct->getXrightYright_neighbor();}
  else if (vct->getYright_neighbor()!=MPI_PROC_NULL)
    {right=vct->getYright_neighbor();}
  // new
  else if (vct->getXright_neighbor()!=MPI_PROC_NULL)
    {right=vct->getXright_neighbor();}
  // end new
  else
    {right=MPI_PROC_NULL;}

  if (vct->getXleftYleft_neighbor()!=MPI_PROC_NULL)
    {left= vct->getXleftYleft_neighbor();}
  else if (vct->getYleft_neighbor()!=MPI_PROC_NULL)
    {left=vct->getYleft_neighbor();}
  // new
  else if (vct->getXleft_neighbor()!=MPI_PROC_NULL)
    {left=vct->getXleft_neighbor();}
  // end new
  else
    {left=MPI_PROC_NULL;}

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);
  cout << "communicate ghost corner done" << endl;
  if(vct->getXrightYleft_neighbor()!=MPI_PROC_NULL)
    {right=vct->getXrightYleft_neighbor();}
  else if(vct->getYleft_neighbor()!=MPI_PROC_NULL)
    {right=vct->getYleft_neighbor();}
  // new
  else if(vct->getXright_neighbor()!=MPI_PROC_NULL)
    {right=vct->getXright_neighbor();}
  // end new
  else
    {right=MPI_PROC_NULL;}
  if (vct->getXleftYright_neighbor()!=MPI_PROC_NULL)
    {left= vct->getXleftYright_neighbor();}
  else if(vct->getYright_neighbor()!=MPI_PROC_NULL)
    {left=vct->getYright_neighbor();}
  // new
  else if(vct->getXleft_neighbor()!=MPI_PROC_NULL)
    {left=vct->getXleft_neighbor();}
  // end new
  else
    {left=MPI_PROC_NULL;}

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
*/


  // take the values from the buffer and put it back to the edge ghost cells
  parseCorner(nx,ny,vector, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner);
	  

  delete[] ghostXrightFace;
  delete[] ghostXleftFace;
  delete[] ghostYrightFace;
  delete[] ghostYleftFace;
    
}

/**SPECIES: communicate ghost node values among processors in 2D*/
inline void communicateNode(int nx, int ny, double ****vector, int ns, VirtualTopology *vct){
  // allocate 4 ghost cell Faces
	double *ghostXrightFace = new double[ny-2];
	double *ghostXleftFace  = new double[ny-2];
	double *ghostYrightFace = new double[nx-2];
	double *ghostYleftFace  = new double[nx-2];
  
  // allocate 4 ghost cell corner
	double ghostXrightYrightCorner,ghostXleftYrightCorner,ghostXrightYleftCorner,ghostXleftYleftCorner;
  
	makeNodeFaceX(nx,ny,vector, ns,ghostXrightFace,ghostXleftFace,vct);
	makeNodeFaceY(nx,ny,vector, ns,ghostYrightFace,ghostYleftFace,vct);
  
	makeNodeCorner(nx,ny, vector, ns, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner,&ghostXleftYleftCorner,vct);
  
	communicateGhostFace(vct->getCART_COMM(),(ny-2),vct->getCartesian_rank(),vct->getXright_neighbor(),vct->getXleft_neighbor(),0,vct->getXLEN(),vct->getYLEN(),ghostXrightFace, ghostXleftFace);
	communicateGhostFace(vct->getCART_COMM(),(nx-2),vct->getCartesian_rank(),vct->getYright_neighbor(),vct->getYleft_neighbor(),1,vct->getXLEN(),vct->getYLEN(),ghostYrightFace, ghostYleftFace);
  
	parseFaceX(nx,ny,vector,ns,ghostXrightFace,ghostXleftFace);
	parseFaceY(nx,ny,vector,ns,ghostYrightFace,ghostYleftFace);

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYright_neighbor(), vct->getXleftYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);
  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYleft_neighbor(), vct->getXleftYright_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
  if (vct->getYright_neighbor() == MPI_PROC_NULL ){ //Top line
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYrightCorner);
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL ){ //Bottom line
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYleftCorner);
  }
  if (vct->getXleft_neighbor() == MPI_PROC_NULL ){ //Left column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXleftYrightCorner, &ghostXleftYleftCorner);
  }
  if (vct->getXright_neighbor() == MPI_PROC_NULL ){ //Right column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXrightYleftCorner);
  }

	/*int left;
	int right;
	if(vct->getXrightYright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXrightYright_neighbor();}
	else if (vct->getYright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getYright_neighbor();}
	// new
	else if (vct->getXright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXright_neighbor();}
	// end new
	else
	  {right=MPI_PROC_NULL;}

	if (vct->getXleftYleft_neighbor()!=MPI_PROC_NULL)
	  {left= vct->getXleftYleft_neighbor();}
	else if (vct->getYleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getYleft_neighbor();}
	// new
	else if (vct->getXleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getXleft_neighbor();}
	// end new
	else
	  {left=MPI_PROC_NULL;}

	communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);

	if(vct->getXrightYleft_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXrightYleft_neighbor();}
	else if(vct->getYleft_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getYleft_neighbor();}
	// new
	else if(vct->getXright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXright_neighbor();}
	// end new
	else
	  {right=MPI_PROC_NULL;}
	if (vct->getXleftYright_neighbor()!=MPI_PROC_NULL)
	  {left= vct->getXleftYright_neighbor();}
	else if(vct->getYright_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getYright_neighbor();}
	// new
	else if(vct->getXleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getXleft_neighbor();}
	// end new
	else
	  {left=MPI_PROC_NULL;}

	communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);   
*/
	parseCorner(nx,ny,vector, ns, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner);
	  

	delete[] ghostXrightFace;
	delete[] ghostXleftFace;
	delete[] ghostYrightFace;
	delete[] ghostYleftFace;
    
}


///** communicate ghost center values among processors in 2D*/
inline void communicateCenter(int nx, int ny, double ***vector, VirtualTopology *vct){
   // allocate 4 ghost cell Faces
	double *ghostXrightFace = new double[ny-2];
	double *ghostXleftFace  = new double[ny-2];
	double *ghostYrightFace = new double[nx-2];
	double *ghostYleftFace  = new double[nx-2];
  
  // allocate 4 ghost cell corner
	double ghostXrightYrightCorner,ghostXleftYrightCorner,ghostXrightYleftCorner,ghostXleftYleftCorner;
  
	makeCenterFaceX(nx,ny,vector,ghostXrightFace,ghostXleftFace,vct);
	makeCenterFaceY(nx,ny,vector,ghostYrightFace,ghostYleftFace,vct);
	makeCenterCorner(nx,ny, vector, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner,vct);

  
	communicateGhostFace(vct->getCART_COMM(),(ny-2),vct->getCartesian_rank(),vct->getXright_neighbor(),vct->getXleft_neighbor(),0,vct->getXLEN(),vct->getYLEN(),ghostXrightFace, ghostXleftFace);
	communicateGhostFace(vct->getCART_COMM(),(nx-2),vct->getCartesian_rank(),vct->getYright_neighbor(),vct->getYleft_neighbor(),1,vct->getXLEN(),vct->getYLEN(),ghostYrightFace, ghostYleftFace);
	
	parseFaceX(nx,ny,vector,ghostXrightFace,ghostXleftFace);
	parseFaceY(nx,ny,vector,ghostYrightFace,ghostYleftFace);

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYright_neighbor(), vct->getXleftYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);
  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYleft_neighbor(), vct->getXleftYright_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
  if (vct->getYright_neighbor() == MPI_PROC_NULL ){ //Top line
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYrightCorner);
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL ){ //Bottom line
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYleftCorner);
  }
  if (vct->getXleft_neighbor() == MPI_PROC_NULL ){ //Left column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXleftYrightCorner, &ghostXleftYleftCorner);
  }
  if (vct->getXright_neighbor() == MPI_PROC_NULL ){ //Right column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXrightYleftCorner);
  }
	
/*	int left;
	int right;
	if(vct->getXrightYright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXrightYright_neighbor();}
	else if (vct->getYright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getYright_neighbor();}
	// new
	else if (vct->getXright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXright_neighbor();}
	// end new
	else
	  {right=MPI_PROC_NULL;}

	if (vct->getXleftYleft_neighbor()!=MPI_PROC_NULL)
	  {left= vct->getXleftYleft_neighbor();}
	else if (vct->getYleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getYleft_neighbor();}
	// new
	else if (vct->getXleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getXleft_neighbor();}
	// end new
	else
	  {left=MPI_PROC_NULL;}

	communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);

	if(vct->getXrightYleft_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXrightYleft_neighbor();}
	else if(vct->getYleft_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getYleft_neighbor();}
	// new
	else if(vct->getXright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXright_neighbor();}
	// end new
	else
	  {right=MPI_PROC_NULL;}
	if (vct->getXleftYright_neighbor()!=MPI_PROC_NULL)
	  {left= vct->getXleftYright_neighbor();}
	else if(vct->getYright_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getYright_neighbor();}
	// new
	else if(vct->getXleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getXleft_neighbor();}
	// end new
	else
	  {left=MPI_PROC_NULL;}

	communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
*/
	parseCorner(nx,ny,vector, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner);
	  

	delete[] ghostXrightFace;
	delete[] ghostXleftFace;
	delete[] ghostYrightFace;
	delete[] ghostYleftFace;
    
}

///** SPECIES: communicate ghost center values among processors in 2D*/
inline void communicateCenter(int nx, int ny, double ****vector, int ns,VirtualTopology *vct){
  // allocate 4 ghost cell Faces
	double *ghostXrightFace = new double[ny-2];
	double *ghostXleftFace  = new double[ny-2];
	double *ghostYrightFace = new double[nx-2];
	double *ghostYleftFace  = new double[nx-2];
  
  // allocate 4 ghost cell corner
	double ghostXrightYrightCorner,ghostXleftYrightCorner,ghostXrightYleftCorner,ghostXleftYleftCorner;
  
	makeCenterFaceX(nx,ny,vector,ns,ghostXrightFace,ghostXleftFace,vct);
	makeCenterFaceY(nx,ny,vector,ns,ghostYrightFace,ghostYleftFace,vct);
	makeCenterCorner(nx,ny, vector, ns,&ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner,vct);

  
	communicateGhostFace(vct->getCART_COMM(),(ny-2),vct->getCartesian_rank(),vct->getXright_neighbor(),vct->getXleft_neighbor(),0,vct->getXLEN(),vct->getYLEN(),ghostXrightFace, ghostXleftFace);
	communicateGhostFace(vct->getCART_COMM(),(nx-2),vct->getCartesian_rank(),vct->getYright_neighbor(),vct->getYleft_neighbor(),1,vct->getXLEN(),vct->getYLEN(),ghostYrightFace, ghostYleftFace);
	
	parseFaceX(nx,ny,vector, ns, ghostXrightFace,ghostXleftFace);
	parseFaceY(nx,ny,vector, ns, ghostYrightFace,ghostYleftFace);

  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYright_neighbor(), vct->getXleftYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);
  communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXrightYleft_neighbor(), vct->getXleftYright_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner);
  if (vct->getYright_neighbor() == MPI_PROC_NULL ){ //Top line
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYrightCorner);
  }
  if (vct->getYleft_neighbor() == MPI_PROC_NULL ){ //Bottom line
      communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getXright_neighbor(), vct->getXleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYleftCorner);
  }
  if (vct->getXleft_neighbor() == MPI_PROC_NULL ){ //Left column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXleftYrightCorner, &ghostXleftYleftCorner);
  }
  if (vct->getXright_neighbor() == MPI_PROC_NULL ){ //Right column
      communicateGhostCorner_forcolumn(vct->getCART_COMM(),vct->getCartesian_rank(), vct->getYright_neighbor(), vct->getYleft_neighbor(), vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXrightYleftCorner);
  }

/*	int left;
	int right;
	if(vct->getXrightYright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXrightYright_neighbor();}
	else if (vct->getYright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getYright_neighbor();}
	// new
	else if (vct->getXright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXright_neighbor();}
	// end new
	else
	  {right=MPI_PROC_NULL;}

	if (vct->getXleftYleft_neighbor()!=MPI_PROC_NULL)
	  {left= vct->getXleftYleft_neighbor();}
	else if (vct->getYleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getYleft_neighbor();}
	// new
	else if (vct->getXleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getXleft_neighbor();}
	// end new
	else
	  {left=MPI_PROC_NULL;}

	communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYrightCorner, &ghostXleftYleftCorner);

	if(vct->getXrightYleft_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXrightYleft_neighbor();}
	else if(vct->getYleft_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getYleft_neighbor();}
	// new
	else if(vct->getXright_neighbor()!=MPI_PROC_NULL)
	  {right=vct->getXright_neighbor();}
	// end new
	else
	  {right=MPI_PROC_NULL;}
	if (vct->getXleftYright_neighbor()!=MPI_PROC_NULL)
	  {left= vct->getXleftYright_neighbor();}
	else if(vct->getYright_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getYright_neighbor();}
	// new
	else if(vct->getXleft_neighbor()!=MPI_PROC_NULL)
	  {left=vct->getXleft_neighbor();}
	// end new
	else
	  {left=MPI_PROC_NULL;}

	communicateGhostCorner(vct->getCART_COMM(),vct->getCartesian_rank(), right, left, vct->getXLEN(),vct->getYLEN(), &ghostXrightYleftCorner, &ghostXleftYrightCorner); 
*/
	parseCorner(nx,ny,vector, ns, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner);
	  

	delete[] ghostXrightFace;
	delete[] ghostXleftFace;
	delete[] ghostYrightFace;
	delete[] ghostYleftFace;
    
}
/**SPECIES: communicate ghost node values among processors in 2D*/
inline void communicateNodeOS(int nx, int ny, double ****vector, int ns, VirtualTopology *vct){
  // allocate 4 ghost cell Faces
  
  // allocate 4 ghost cell corner
  double ghostXrightYrightCorner,ghostXleftYrightCorner,ghostXrightYleftCorner,ghostXleftYleftCorner;
  
  makeNodeCorner(nx,ny, vector, ns, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner,&ghostXleftYleftCorner,vct);
  
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
  parseCorner(nx,ny,vector, ns, &ghostXrightYrightCorner, &ghostXleftYrightCorner, &ghostXrightYleftCorner, &ghostXleftYleftCorner);
    
}

#endif

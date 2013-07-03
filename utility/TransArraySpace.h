/***************************************************************************
TransArraySpace.h  -  functions to convert 3field to 1D field
                         
developers: Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
   
 ***************************************************************************/

#ifndef TransArraySpace_H
#define TransArraySpace_H

/**
* 
* Utility functions for moving array to physical space to krylocv space for solvers
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/

/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//              PHYS TO SOLVER            //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/** method to convert a 1D field in a 1D field not considering guard cells*/
inline void phys2solver(double* vectSolver, double*** vectPhys, int nx){
    for (register int i=1; i < nx-1; i++)
      *vectSolver++ = vectPhys[i][0][0]; 
}
   
/** method to convert a 2D field in a 1D field not considering guard cells*/
inline void phys2solver(double* vectSolver, double*** vectPhys, int nx, int ny){
    for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        *vectSolver++ = vectPhys[i][j][0]; 
}
/** method to convert a 2D field in a 1D field CONSIDERING guard cells*/
inline void phys2solver_plusghost(double* vectSolver, double*** vectPhys, int nx, int ny)
{
  for (register int i=0; i < nx; i++)
    for (register int j=0; j < ny; j++)
      *vectSolver++ = vectPhys[i][j][0];
}
/** method to convert a 1D field in a 1D field not considering guard cells, added by Enrico*/
 inline void phys2solver(double* vectSolver, double*** vectPhys1,double*** vectPhys2,int nx){
	for (register int i=1; i < nx-1; i++){
	    *vectSolver++ =  vectPhys1[i][0][0];
	    *vectSolver++ =  vectPhys2[i][0][0];
	  }
}
/** method to convert a 2D field in a 1D field not considering guard cells, added by Enrico*/
 inline void phys2solver(double* vectSolver, double*** vectPhys1,double*** vectPhys2,int nx, int ny){
	for (register int i=1; i < nx-1; i++)
	  for (register int j=1; j < ny-1; j++){
	    *vectSolver++ =  vectPhys1[i][j][0];
	    *vectSolver++ =  vectPhys2[i][j][0];
	  }
   }
/** method to convert a 1D field in a 1D field not considering guard cells, added by Enrico*/
 inline void phys2solver(double* vectSolver, double*** vectPhys1,double*** vectPhys2,double*** vectPhys3, int nx){
	for (register int i=1; i < nx-1; i++){
	    *vectSolver++ =  vectPhys1[i][0][0];
	    *vectSolver++ =  vectPhys2[i][0][0];
	    *vectSolver++ =  vectPhys3[i][0][0];
	}
}
/** method to convert a 2D field in a 1D field not considering guard cells, added by Enrico*/
 inline void phys2solver(double* vectSolver, double*** vectPhys1,double*** vectPhys2,double*** vectPhys3, int nx, int ny){
	for (register int i=1; i < nx-1; i++)
	  for (register int j=1; j < ny-1; j++){
	    *vectSolver++ =  vectPhys1[i][j][0];
	    *vectSolver++ =  vectPhys2[i][j][0];
	    *vectSolver++ =  vectPhys3[i][j][0];
	  }
   }

/** method to convert a 3D field in a 1D field not considering guard cells*/
inline void phys2solver(double* vectSolver, double*** vectPhys, int nx, int ny, int nz){
     for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        for (register int k=1; k < nz-1; k++)
           *vectSolver++ = vectPhys[i][j][k];
           

}

/** method to convert a 3D field in a 1D field not considering guard cells*/
inline void phys2solver(double* vectSolver, double*** vectPhys1,double*** vectPhys2,double*** vectPhys3, int nx, int ny, int nz){
     for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        for (register int k=1; k < nz-1; k++){
           *vectSolver++ =  vectPhys1[i][j][k];
           *vectSolver++ =  vectPhys2[i][j][k];
           *vectSolver++ =  vectPhys3[i][j][k];
        }
}


/** method to convert a 3D field in a 1D field considering only active cells */
inline void phys2solver(double* vectSolver, double*** vectPhys1,double*** vectPhys2,double*** vectPhys3, int xStartActiveN, int xEndActiveN, int yStartActiveN, int yEndActiveN, int zStartActiveN, int zEndActiveN){
     for (register int i=xStartActiveN; i <= xEndActiveN;i++)
      for (register int j=yStartActiveN; j <= yEndActiveN;j++)
        for (register int k=zStartActiveN; k <= zEndActiveN;k++){
          *vectSolver++ =  vectPhys1[i][j][k];
          *vectSolver++ =  vectPhys2[i][j][k];
          *vectSolver++ =  vectPhys3[i][j][k];
        }
}
/////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//              SOLVER TO PHYS            //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/** method to convert a 1D field in a 1D field not considering guard cells, added by Enrico*/
inline void solver2phys(double*** vectPhys, double* vectSolver, int nx){
    for (register int i=1; i < nx-1; i++)
           vectPhys[i][0][0] = *vectSolver++;
          
}
/** method to convert a 1D field in a 2D field not considering guard cells, added by Enrico*/
inline void solver2phys(double*** vectPhys, double* vectSolver, int nx, int ny){
    for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
           vectPhys[i][j][0] = *vectSolver++;
          
}
/** method to convert a 1D field in a 2D field CONSIDERING guard cells*/
inline void solver2phys_plusghost(double*** vectPhys, double* vectSolver, int nx, int ny)
{
  for (register int i=0; i < nx; i++)
    for (register int j=0; j < ny; j++)
      vectPhys[i][j][0] = *vectSolver++;

}
/** method to convert a 1D field in a 2D field not considering guard cells, added by Enrico*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double* vectSolver, int nx, int ny){
    for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++){
          vectPhys1[i][j][0] = *vectSolver++;
          vectPhys2[i][j][0] = *vectSolver++;
        }
}
/** method to convert a 1D field in a 1D field not considering guard cells, added by Enrico*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double* vectSolver, int nx){
    for (register int i=1; i < nx-1; i++){
          vectPhys1[i][0][0] = *vectSolver++;
          vectPhys2[i][0][0] = *vectSolver++;
    }
}
/** method to convert a 1D field in a 3D field considering only active cells, added by Enrico*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double*** vectPhys3, double* vectSolver,int nx, int ny){
	for (register int i=1; i < nx-1;i++)
	 for (register int j=1; j < ny-1;j++){
	  vectPhys1[i][j][0] = *vectSolver++;
	  vectPhys2[i][j][0] = *vectSolver++;
	  vectPhys3[i][j][0] = *vectSolver++;
	 }
}
/** method to convert a 1D field in a 3D field considering only active cells, added by Enrico*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double*** vectPhys3, double* vectSolver,int nx){
	for (register int i=1; i < nx-1;i++){
	  vectPhys1[i][0][0] = *vectSolver++;
	  vectPhys2[i][0][0] = *vectSolver++;
	  vectPhys3[i][0][0] = *vectSolver++;
	 }
}

/** method to convert a 1D field in a 3D field not considering guard cells*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double*** vectPhys3, double* vectSolver, int nx, int ny, int nz){
    for (register int i=1; i < nx-1; i++)
      for (register int j=1; j < ny-1; j++)
        for (register int k=1; k < nz-1; k++){
          vectPhys1[i][j][k] = *vectSolver++;
          vectPhys2[i][j][k] = *vectSolver++;
          vectPhys3[i][j][k] = *vectSolver++;
        }
}




/** method to convert a 1D field in a 3D field considering only active cells*/
inline void solver2phys(double*** vectPhys1, double*** vectPhys2, double*** vectPhys3, double* vectSolver,int xStartActiveN, int xEndActiveN, int yStartActiveN, int yEndActiveN, int zStartActiveN, int zEndActiveN){
    for (register int i=xStartActiveN; i <= xEndActiveN;i++)
      for (register int j=yStartActiveN; j <= yEndActiveN;j++)
        for (register int k=zStartActiveN; k <= zEndActiveN;k++){
          vectPhys1[i][j][k] = *vectSolver++;
          vectPhys2[i][j][k] = *vectSolver++;
          vectPhys3[i][j][k] = *vectSolver++;
        }
}





#endif

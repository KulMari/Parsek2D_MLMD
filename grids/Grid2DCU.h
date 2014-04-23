/*******************************************************************************************************
 Grid2DCU.h  -  uniform cartesian 2D local grid for each process, including che guard cells
p                             -------------------
 developers: Stefano Markidis, Giovanni Lapenta
 *******************************************************************************************************/

#ifndef GRID2DCU_H
#define GRID2DCU_H

#include <iostream>

#include "Grid.h"
#include "../communication/ComInterpNodes.h"
#include "../communication/ComNodes.h"
#include "../utility/Alloc.h"
#include "../mathlib/Basic.h"

using std::cout;
using std::endl;

/**
* 2D Uniform cartesian local(each processor has its own grid) grid 
*
* Note that GRID2DCU is implementing the abstract class Grid.h, therefore
* it must implement all the virtual methods in Grid.h
*
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta
* @version 2.0
*
*/
class Grid2DCU : public Grid {
  public:
      /** constructor */
      Grid2DCU(CollectiveIO *col, VirtualTopology *vct, int GridLevel);
      /** destructor */
      ~Grid2DCU();
      /** allocate grid arrays for this domain */
      void allocate(CollectiveIO *ptC, VirtualTopology *ptVCT);
      /** deallocate grid arrays for this domain */
      void deallocate();
      /** print grid info */
      void print(VirtualTopology* ptVCT);
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradC2N(double*** gradXN,double*** gradYN,double***scFieldC,VirtualTopology *vct);
      /** calculate gradient on nodes, given a scalar field defined on central points with one sided derivative on the BC */
      void gradC2N_onesided_derBC(double*** gradXN,double*** gradYN,double***scFieldC,VirtualTopology *vct);
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradN2C(double ***gradXC, double ***gradYC, double ***scFieldN);
      /** calculate gradient on nodes, given a scalar field defined on central points including ghost centers  */
      void gradN2C_plusghost(double ***gradXC, double ***gradYC, double ***scFieldN);
      /** calculate divergence on central points, given a vector field defined on nodes  */
      void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN);
      /** calculate divergence on central points, given a vector field defined on nodes including ghost centers calculated using ghost nodes */
      void divN2C_plusghost(double ***divC, double ***vecFieldXN, double ***vecFieldYN);
      /** calculate divergence on nodes, given a vector field defined on central points  */
      void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC,VirtualTopology *vct);
	  /** calculate divergence on nodes, given a vector field defined on central points with one sided derivative on the BC */
      void divC2N_onesided_derBC(double ***divN, double ***vecFieldXC, double ***vecFieldYC,VirtualTopology *vct);
      /** calculate curl on nodes, given a vector field defined on central points  */
      void curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC,VirtualTopology *vct);
      /** calculate curl on central points, given a vector field defined on nodes  */
      void curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN);
      /** calculate curl on central points, given a vector field defined on nodes including ghost cells */
      void curlN2C_withghost(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN);

      /** calculate divergence on central points, given a Tensor field defined on nodes  */
      void divSymmTensorN2C(double ***divCX, double ***divCY, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ,int ns);
      void divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ,int ns);

      void divSymmTensorN2C_alsoGC(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ,int ns);
      
      /** calculate laplacian on nodes, given a scalar field defined on nodes */
      void lapN2N(double*** lapN,double ***scFieldN,VirtualTopology *vct);
      /** calculate laplacian on active nodes, given a scalar field defined on nodes including ghost nodes*/
      void lapN2N_plusghost(double*** lapN,double ***scFieldN,VirtualTopology *vct);
      /** calculate laplacian on central points, given a scalar field defined on central points */
      void lapC2C(double*** lapC,double ***scFieldC,VirtualTopology *vct);
      /** calculate laplacian on central points, given a scalar field defined on central points; patch for ghost centers */
      void lapC2C_plusghost(double*** lapC,double ***scFieldC,VirtualTopology *vct);
	  /** interpolate on nodes from central points */
      void interpC2N(double ***vecFieldN, double ***vecFieldC,VirtualTopology *vct);
	  /** interpolate on nodes from central points */
      void interpC2N_BC(double ***vecFieldN, double ***vecFieldC,VirtualTopology *vct);
      // interpolate on nodes from central points, also GN
      void interpC2N_BC_alsoGN(double ***vecFieldN, double ***vecFieldC, VirtualTopology *vct);
      /** interpolate on central points from nodes */
      void interpN2C(double ***vecFieldC, double ***vecFieldN);
      /** interpolate on central points from nodes per species*/
      void interpN2C(double ****vecFieldC, int ns, double ****vecFieldN);
      /** interpolate on central points from nodes per species including ghost cells*/
      void interpN2C_alsoGC(double ****vecFieldC, int ns, double ****vecFieldN);
      // interpN2C also on ghost cells, make sure ghost nodes are updated                                                                       
      void interpN2C_alsoGC(double ***vecFieldC, double ***vecFieldN);
      /** return grid level */
      int getLevel(); 
      /** return ratio */
      double getRatio(); 
      /** return origin in x **/
      double &getOx(int level);
      /** return origin in y **/
      double &getOy(int level);
      /** return number of cells in direction X*/
      int getNXC();
      /** return number of nodes in direction X*/
      int getNXN();
      /** return number of cells in direction Y*/
      int getNYC();
      /** return number of nodes in direction Y*/
      int getNYN();
      /** return number of cells in direction Z*/
      int getNZC();
      /** return number of nodes in direction Z*/
      int getNZN();
      /** return grid spacing dx */
      double getDX();
      /** return grid spacing dy */
      double getDY();
      /** return grid spacing dz */
      double getDZ();
      /** get x-coordinate of the node(X,Y,Z) */
      double &getXN(int indexX, int indexY, int indexZ);
      /** get y-coordinate of the node(X,Y,Z) */
      double &getYN(int indexX, int indexY, int indexZ);
      /** get x-coordinate of the node(X,Y,Z) */
      double getModifiedXN(int indexX, int indexY, int indexZ);
      /** get y-coordinate of the node(X,Y,Z) */
      double getModifiedYN(int indexX, int indexY, int indexZ);
      /** get z-coordinate of the node(X,Y,Z) */
      double &getZN(int indexX, int indexY, int indexZ);
      /** get x-coordinate of the cell(X,Y,Z) */
      double &getXC(int indexX, int indexY, int indexZ);
      /** get y-coordinate of the cell(X,Y,Z) */
      double &getYC(int indexX, int indexY, int indexZ);
      /** get z-coordinate of the cell(X,Y,Z) */
      double &getZC(int indexX, int indexY, int indexZ);


      /** get the whole vector xc with the x-coordinate of center cells*/
      double*** getXC();
      /** get the whole vector yc with the y-coordinate of center cells*/
      double*** getYC();
       /** get the whole vector yc with the z-coordinate of center cells*/
      double*** getZC();
      /** get x coordinate of the first node of Simulation box - X direction */
      double getXstart();
      /** get x coordinate of the first node of Simulation box - X direction or of node 0 if this proc is on the x=0 border */
      double getmodifiedXstart(VirtualTopology *vct);
      /** get x coordinate of the last node of Simulation box - X direction*/
      double getXend();
      /** get x coordinate of the last node of Simulation box - X direction or of the VERY last node if this proc is on the x=Lx border*/
      double getmodifiedXend(VirtualTopology *vct);
      /** get y coordinate of the first node of Simulation box - Y direction */
      double getYstart();
      /** get y coordinate of the first node of Simulation box - Y direction or of node 0 if this proc is on the y=0 border */
      double getmodifiedYstart(VirtualTopology *vct);
      /** get y coordinate of the last node of Simulation box - Y direction */
      double getYend();
      /** get y coordinate of the first node of Simulation box - Y direction or of the VERY last node if this proc is on the y=Ly border */
      double getmodifiedYend(VirtualTopology *vct);
      /** get z coordinate of the first node of Simulation box - Z direction */
      double getZstart();
      /** get z coordinate of the last node of Simulation box - Z direction */
      double getZend();

      /** get the inverse of volume */
      double getInvVOL();

  private:
     /** Grid Level **/
     int level;
     /** Ratio **/
     double ratio;
     /** Origins of the grids in x direction **/
     double* ox;
     /** Origins of the grids in y direction **/
     double* oy;
     /** number of cells - X direction, including + 2 (guard cells) */
     int nxc;
     /** number of nodes - X direction, including + 2 extra nodes for guard cells */
     int nxn;
     /** number of cell - Y direction, including + 2 (guard cells) */
     int nyc;
     /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
     int nyn;
     /** dx = grid spacing - X direction */
     double dx;
     /** dy = grid spacing - Y direction */
     double dy;
     /** invdx = 1/dx */
     double invdx;
     /** invdy = 1/dy */
     double invdy;
     /** invol = inverse of volume*/
     double invVOL;
     // put them in Grid.h, to make them publicly accessible
     /** node - X coordinate (indexX, indexY)   */
     //double ***xn;
     /** node - Y coordinate (indexX, indexY)   */
     //double ***yn;
     /** centre of cell - X coordinate */
     double ***xc;
     /** centre of cell - Y coordinate */
     double ***yc;
     /** local grid boundaries coordinates  */
     double xStart, xEnd, yStart, yEnd;
     
};
/** constructor */
inline Grid2DCU::Grid2DCU(CollectiveIO* col, VirtualTopology* vct, int GridLevel){
   int i;
   double lengthx, lengthy;
   ox = new double[col->getNgrids()];
   oy = new double[col->getNgrids()];
   level = GridLevel;
   ratio = col->getRatio();
   // add 2 for the guard cells
   nxc = (col->getNxc())/(vct->getXLEN()) + 2;
   nyc = (col->getNyc())/(vct->getYLEN()) + 2;
   nxn = nxc + 1;
   nyn = nyc + 1;
   dx  = col->getLx()/col->getNxc()/(double)pow(ratio,level);
   dy  = col->getLy()/col->getNyc()/(double)pow(ratio,level);
   invVOL = 1.0/(dx*dy);
   invdx = 1.0/dx;
   invdy = 1.0/dy;

   lengthx=col->getLx()/(double)pow(ratio,level); //Total length of the grid in x
   lengthy=col->getLy()/(double)pow(ratio,level); //Total length of the grid in y
       //ox and oy are only used when different grids interact
      //ox and oy are the coordinates of the the current grid orgine in the frame of the grid from level just above it (just coarser grid), NOT in the frame of the level 0 grid.
   // For a centered refined grid
   ox[0]=(double)0.;
   oy[0]=(double)0.;
   //if (vct->getCartesian_rank_COMMTOTAL()==0){
   //cout << "level " << level << " Ox = " << ox[0]<< " Oy = "<<oy[0]<< " dx = "<<dx<<" dy = "<<dy<<" nxn = "<<nxn<< "nyn = "<<nyn<<endl;}

   for (i=1;i<col->getNgrids();i++){
     // refined grid centered

       //Manually tuned
     //ox[i] = 15.10-col->getLx()/pow(ratio,i)/2.;
     //oy[i] = 15.10-col->getLy()/pow(ratio,i)/2.;

     // from inputfile
     ox[i] = col->getL1_CX()-col->getLx()/pow(ratio,i)/2.;
     oy[i] = col->getL1_CY()-col->getLy()/pow(ratio,i)/2.;

     if (vct->getCartesian_rank()==0 && level == i){
       cout << "Level " << i << " getL1_CX: " << col->getL1_CX() << ", L1_CY: " << col->getL1_CY()  << ", Ox: " << ox[i]<< ", Oy: "<<oy[i]<< ", dx: "<<dx<<", dy: "<<dy<<endl;
     }
   }
   // local grid dimensions and boundaries of active nodes
   xStart = vct->getCoordinates(0)*(lengthx/(double) vct->getXLEN());
   xEnd   = xStart + (lengthx/(double) vct->getXLEN());
   yStart = vct->getCoordinates(1)*(lengthy/(double) vct->getYLEN());
   yEnd   = yStart + (lengthy/(double) vct->getYLEN());
 
   // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
   //xn = newArr3(double,nxn,nyn,1);
   //yn = newArr3(double,nxn,nyn,1);
   allocArr3(&xn, nxn, nyn, 1);
   allocArr3(&yn, nxn, nyn, 1);
   for (int i=0; i < nxn; i++){
     for (int j=0; j < nyn; j++){
          xn[i][j][0] = xStart + (i-1)*dx;
          yn[i][j][0] = yStart + (j-1)*dy;
        }
     }
   
   // arrays allocation: cells ---> the first cell has index 1, the last has index ncn-2!
   //xc = newArr3(double,nxc,nyc,1);
   //yc = newArr3(double,nxc,nyc,1);
   allocArr3(&xc, nxc, nyc, 1);
   allocArr3(&yc, nxc, nyc, 1);
   for (int i=0; i < nxc; i++){
     for (int j=0; j < nyc; j++){
          xc[i][j][0] = .5*(xn[i][j][0] + xn[i+1][j][0]);
          yc[i][j][0] = .5*(yn[i][j][0] + yn[i][j+1][0]);
         
     }
   }
}

/** deallocate the local grid */
inline Grid2DCU::~Grid2DCU(){
   // deallocate nodes
   //delArr3(xn,nxn,nyn);
   //delArr3(yn,nxn,nyn);
  freeArr3(&xn);
  freeArr3(&yn);
   // centers cells
  //delArr3(xc,nxc,nyc);
  //delArr3(yc,nxc,nyc);
  freeArr3(&xc);
  freeArr3(&yc);
   
   
}
/** print the local grid info */
inline void Grid2DCU::print(VirtualTopology* ptVCT){
    cout << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ")"<< endl;
    cout <<  "Number of cell: -X=" << nxc-2 << " -Y=" << nyc-2 << endl;
    cout <<  "Xin = " << xn[1][1][0] << "; Xfin = " << xn[nxn-2][1][0] << endl;
    cout <<  "Yin = " << yn[1][1][0] << "; Yfin = " << yn[1][nyn-2][0] << endl;
    cout << endl;
    
}
/** calculate gradient on nodes, given a scalar field defined on central points  with one sided derivative on BC*/
inline void Grid2DCU::gradC2N(double***  gradXN, double***  gradYN, double*** scFieldC,VirtualTopology *vct){
    for (register int i=1;i < nxn-1;i++)
	    for (register int j=1;j < nyn-1; j++){
            gradXN[i][j][0] = .5*(scFieldC[i][j][0] - scFieldC[i-1][j][0])*invdx +  .5*(scFieldC[i][j-1][0] - scFieldC[i-1][j-1][0])*invdx;
            gradYN[i][j][0] = .5*(scFieldC[i][j][0] - scFieldC[i][j-1][0])*invdy +  .5*(scFieldC[i-1][j][0] - scFieldC[i-1][j-1][0])*invdy;
         }
  }

/** calculate gradient on nodes, given a scalar field defined on central points  with one sided derivative on BC*/
inline void Grid2DCU::gradC2N_onesided_derBC(double***  gradXN, double***  gradYN, double*** scFieldC,VirtualTopology *vct){
    for (register int i=1;i < nxn-1;i++)
	    for (register int j=1;j < nyn-1; j++){
            gradXN[i][j][0] = .5*(scFieldC[i][j][0] - scFieldC[i-1][j][0])*invdx +  .5*(scFieldC[i][j-1][0] - scFieldC[i-1][j-1][0])*invdx;
            gradYN[i][j][0] = .5*(scFieldC[i][j][0] - scFieldC[i][j-1][0])*invdy +  .5*(scFieldC[i-1][j][0] - scFieldC[i-1][j-1][0])*invdy;
         }
	// adjust the boundary
	// XLEFT
    if (vct->getXleft_neighbor()==MPI_PROC_NULL){
       for (register int j=2;j < nyn-2;j++){
         gradXN[1][j][0] = gradXN[2][j][0];
		 gradYN[1][j][0] = (scFieldC[1][j][0] - scFieldC[1][j-1][0])*invdy;
	   }
    }
    // XRIGHT
   if (vct->getXright_neighbor()==MPI_PROC_NULL){
      for (register int j=2;j < nyn-2;j++){
         gradXN[nxn-2][j][0] = gradXN[nxn-3][j][0];
		 gradYN[nxn-2][j][0] = (scFieldC[nxc-2][j][0] - scFieldC[nxc-2][j-1][0])*invdy;
	  }
   }
   // YLEFT
  if (vct->getYleft_neighbor()==MPI_PROC_NULL){
   for (register int i=2;i < nxn-2;i++){
         gradXN[i][1][0] = (scFieldC[i][1][0] - scFieldC[i-1][1][0])*invdx;
		 gradYN[i][1][0] = gradYN[i][2][0];
	}
  }
  // YRIGHT
  if (vct->getYright_neighbor()==MPI_PROC_NULL){
     for (register int i=2;i < nxn-2;i++){
         gradXN[i][nyn-2][0] = (scFieldC[i][nyc-2][0] - scFieldC[i-1][nyc-2][0])*invdx;
		 gradYN[i][nyn-2][0] = gradYN[i][nyn-3][0];
	 }
  }
  // XLEFTYLEFT
  if (vct->getXleftYleft_neighbor()==MPI_PROC_NULL){
     gradXN[1][1][0] = gradXN[2][2][0];
	 gradYN[1][1][0] = gradYN[2][2][0];
  }
  // XLEFTYRIGHT
  if (vct->getXleftYright_neighbor()==MPI_PROC_NULL){
     gradXN[1][nyn-2][0] = gradXN[2][nyn-3][0]; 
	 gradYN[1][nyn-2][0] = gradYN[2][nyn-3][0];
  }
  // XRIGHTYLEFT
  if (vct->getXrightYleft_neighbor()==MPI_PROC_NULL){
     gradXN[nxn-2][1][0] = gradXN[nxn-3][2][0];
	 gradYN[nxn-2][1][0] = gradYN[nxn-3][2][0]; 
  }
  // XRIGHTYRIGHT
  if (vct->getXrightYright_neighbor()==MPI_PROC_NULL){
     gradXN[nxn-2][nyn-2][0] = gradXN[nxn-3][nyn-3][0]; 
	 gradXN[nxn-2][nyn-2][0] = gradXN[nxn-3][nyn-3][0];
  }
}

/** calculate gradient on nodes, given a scalar field defined on central points  */
inline void Grid2DCU::gradN2C(double ***gradXC, double ***gradYC, double ***scFieldN){
   for (register int i=1;i < nxc-1; i++)
     for (register int j=1;j < nyc-1; j++) {
          gradXC[i][j][0] = .5*(scFieldN[i+1][j][0] - scFieldN[i][j][0])*invdx + .5*(scFieldN[i+1][j+1][0] - scFieldN[i][j+1][0])*invdx;
          gradYC[i][j][0] = .5*(scFieldN[i][j+1][0] - scFieldN[i][j][0])*invdy + .5*(scFieldN[i+1][j+1][0] - scFieldN[i+1][j][0])*invdy;
       }
}
/** calculate gradient on nodes, given a scalar field defined on central points. Includes ghost centers calculated using ghost nodes.  */
inline void Grid2DCU::gradN2C_plusghost(double ***gradXC, double ***gradYC, double ***scFieldN){
   for (register int i=0;i < nxc; i++)
     for (register int j=0;j < nyc; j++) {
          gradXC[i][j][0] = .5*(scFieldN[i+1][j][0] - scFieldN[i][j][0])*invdx + .5*(scFieldN[i+1][j+1][0] - scFieldN[i][j+1][0])*invdx;
          gradYC[i][j][0] = .5*(scFieldN[i][j+1][0] - scFieldN[i][j][0])*invdy + .5*(scFieldN[i+1][j+1][0] - scFieldN[i+1][j][0])*invdy;
       }
}

/** calculate divergence on central points, given a vector field defined on nodes  */
inline void Grid2DCU::divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN){
  double compX;
  double compY;
   for (register int i=1;i < nxc-1 ;i++)
     for (register int j=1;j < nyc-1 ;j++){
          compX = .5*(vecFieldXN[i+1][j][0] - vecFieldXN[i][j][0])*invdx +  .5*(vecFieldXN[i+1][j+1][0] - vecFieldXN[i][j+1][0])*invdx;
          compY = .5*(vecFieldYN[i][j+1][0] - vecFieldYN[i][j][0])*invdy +  .5*(vecFieldYN[i+1][j+1][0] - vecFieldYN[i+1][j][0])*invdy;
          divC[i][j][0] = compX + compY;
      }      
}
/** calculate divergence on central points, given a vector field defined on nodes. This version includes ghost centers calculated using ghost nodes  */
inline void Grid2DCU::divN2C_plusghost(double ***divC, double ***vecFieldXN, double ***vecFieldYN){
  double compX;
  double compY;
   for (register int i=0;i < nxc ;i++)
     for (register int j=0;j < nyc ;j++){
          compX = .5*(vecFieldXN[i+1][j][0] - vecFieldXN[i][j][0])*invdx +  .5*(vecFieldXN[i+1][j+1][0] - vecFieldXN[i][j+1][0])*invdx;
          compY = .5*(vecFieldYN[i][j+1][0] - vecFieldYN[i][j][0])*invdy +  .5*(vecFieldYN[i+1][j+1][0] - vecFieldYN[i+1][j][0])*invdy;
          divC[i][j][0] = compX + compY;
      }      
}

/** calculate divergence on central points, given a Tensor field defined on nodes */
inline void Grid2DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ****pXX, double**** pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
  double comp1X, comp2X, comp3X;
  double comp1Y, comp2Y, comp3Y;
   for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++){
          comp1X = .5*(pXX[ns][i+1][j][0] - pXX[ns][i][j][0])*invdx +  .5*(pXX[ns][i+1][j+1][0] - pXX[ns][i][j+1][0])*invdx;
          comp2X = .5*(pXY[ns][i+1][j][0] - pXY[ns][i][j][0])*invdx +  .5*(pXY[ns][i+1][j+1][0] - pXY[ns][i][j+1][0])*invdx;
          comp3X = .5*(pXZ[ns][i+1][j][0] - pXZ[ns][i][j][0])*invdx +  .5*(pXZ[ns][i+1][j+1][0] - pXZ[ns][i][j+1][0])*invdx;
          comp1Y = .5*(pXY[ns][i][j+1][0] - pXY[ns][i][j][0])*invdy +  .5*(pXY[ns][i+1][j+1][0] - pXY[ns][i+1][j][0])*invdy;
          comp2Y = .5*(pYY[ns][i][j+1][0] - pYY[ns][i][j][0])*invdy +  .5*(pYY[ns][i+1][j+1][0] - pYY[ns][i+1][j][0])*invdy;
          comp3Y = .5*(pYZ[ns][i][j+1][0] - pYZ[ns][i][j][0])*invdy +  .5*(pYZ[ns][i+1][j+1][0] - pYZ[ns][i+1][j][0])*invdy;
          
          divCX[i][j][0] = comp1X + comp2X + comp3X;
          divCY[i][j][0] = comp1Y + comp2Y + comp3Y;
          
      }


}
/** calculate divergence on central points, given a Tensor field defined on nodes */
inline void Grid2DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double**** pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
  double comp1X, comp2X, comp3X;
  double comp1Y, comp2Y, comp3Y;
   for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++){
          comp1X = .5*(pXX[ns][i+1][j][0] - pXX[ns][i][j][0])*invdx +  .5*(pXX[ns][i+1][j+1][0] - pXX[ns][i][j+1][0])*invdx;
          comp2X = .5*(pXY[ns][i+1][j][0] - pXY[ns][i][j][0])*invdx +  .5*(pXY[ns][i+1][j+1][0] - pXY[ns][i][j+1][0])*invdx;
          comp3X = .5*(pXZ[ns][i+1][j][0] - pXZ[ns][i][j][0])*invdx +  .5*(pXZ[ns][i+1][j+1][0] - pXZ[ns][i][j+1][0])*invdx;
          comp1Y = .5*(pXY[ns][i][j+1][0] - pXY[ns][i][j][0])*invdy +  .5*(pXY[ns][i+1][j+1][0] - pXY[ns][i+1][j][0])*invdy;
          comp2Y = .5*(pYY[ns][i][j+1][0] - pYY[ns][i][j][0])*invdy +  .5*(pYY[ns][i+1][j+1][0] - pYY[ns][i+1][j][0])*invdy;
          comp3Y = .5*(pYZ[ns][i][j+1][0] - pYZ[ns][i][j][0])*invdy +  .5*(pYZ[ns][i+1][j+1][0] - pYZ[ns][i+1][j][0])*invdy;
          
          divCX[i][j][0] = comp1X + comp1Y;
          divCY[i][j][0] = comp2X + comp2Y;
	  divCZ[i][j][0] = comp3X + comp3Y;
           
     }


}

/** calculate divergence on central points, given a Tensor field defined on nodes */
inline void Grid2DCU::divSymmTensorN2C_alsoGC(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double**** pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
  double comp1X, comp2X, comp3X;
  double comp1Y, comp2Y, comp3Y;
  for (register int i=0;i < nxc;i++)
    for (register int j=0;j < nyc;j++){
      comp1X = .5*(pXX[ns][i+1][j][0] - pXX[ns][i][j][0])*invdx +  .5*(pXX[ns][i+1][j+1][0] - pXX[ns][i][j+1][0])*invdx;
      comp2X = .5*(pXY[ns][i+1][j][0] - pXY[ns][i][j][0])*invdx +  .5*(pXY[ns][i+1][j+1][0] - pXY[ns][i][j+1][0])*invdx;
      comp3X = .5*(pXZ[ns][i+1][j][0] - pXZ[ns][i][j][0])*invdx +  .5*(pXZ[ns][i+1][j+1][0] - pXZ[ns][i][j+1][0])*invdx;
      comp1Y = .5*(pXY[ns][i][j+1][0] - pXY[ns][i][j][0])*invdy +  .5*(pXY[ns][i+1][j+1][0] - pXY[ns][i+1][j][0])*invdy;
      comp2Y = .5*(pYY[ns][i][j+1][0] - pYY[ns][i][j][0])*invdy +  .5*(pYY[ns][i+1][j+1][0] - pYY[ns][i+1][j][0])*invdy;
      comp3Y = .5*(pYZ[ns][i][j+1][0] - pYZ[ns][i][j][0])*invdy +  .5*(pYZ[ns][i+1][j+1][0] - pYZ[ns][i+1][j][0])*invdy;

      divCX[i][j][0] = comp1X + comp1Y;
      divCY[i][j][0] = comp2X + comp2Y;
      divCZ[i][j][0] = comp3X + comp3Y;

    }


}
/** calculate divergence on nodes, given a vector field defined on central points with one_sided derivative */
inline void Grid2DCU::divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC,VirtualTopology *vct){
  double compX;
  double compY;
  for (register int i=1; i < nxn-1 ;i++)
	  for (register int j=1; j < nyn-1 ;j++){
          compX = .5*(vecFieldXC[i][j][0] - vecFieldXC[i-1][j][0])*invdx + .5*(vecFieldXC[i][j-1][0] - vecFieldXC[i-1][j-1][0])*invdx;
          compY = .5*(vecFieldYC[i][j][0] - vecFieldYC[i][j-1][0])*invdy + .5*(vecFieldYC[i-1][j][0] - vecFieldYC[i-1][j-1][0])*invdy;
          divN[i][j][0] = compX + compY;
       }
}
/** calculate divergence on nodes, given a vector field defined on central points with one_sided derivative */
inline void Grid2DCU::divC2N_onesided_derBC(double ***divN, double ***vecFieldXC, double ***vecFieldYC,VirtualTopology *vct){
  double compX;
  double compY;
  for (register int i=1; i < nxn-1 ;i++)
	  for (register int j=1; j < nyn-1 ;j++){
          compX = .5*(vecFieldXC[i][j][0] - vecFieldXC[i-1][j][0])*invdx + .5*(vecFieldXC[i][j-1][0] - vecFieldXC[i-1][j-1][0])*invdx;
          compY = .5*(vecFieldYC[i][j][0] - vecFieldYC[i][j-1][0])*invdy + .5*(vecFieldYC[i-1][j][0] - vecFieldYC[i-1][j-1][0])*invdy;
          divN[i][j][0] = compX + compY;
       }
    // adjust the boundaries
	// XLEFT
    if (vct->getXleft_neighbor()==MPI_PROC_NULL){
	  for (register int j=2;j < nyn-2;j++){
         compX = .5*(vecFieldXC[2][j][0] - vecFieldXC[1][j][0])*invdx + .5*(vecFieldXC[2][j-1][0] - vecFieldXC[1][j-1][0])*invdx;
		 compY = (vecFieldYC[1][j][0] - vecFieldYC[1][j-1][0])*invdy;
		 divN[1][j][0] = compX + compY;
	   }
    }
    // XRIGHT
   if (vct->getXright_neighbor()==MPI_PROC_NULL){
      for (register int j=2;j < nyn-2;j++){
         compX = .5*(vecFieldXC[nxc-2][j][0] - vecFieldXC[nxc-3][j][0])*invdx + .5*(vecFieldXC[nxc-2][j-1][0] - vecFieldXC[nxc-3][j-1][0])*invdx;
		 compY = (vecFieldYC[nxc-2][j][0] - vecFieldYC[nxc-2][j-1][0])*invdy;
		 divN[nxn-2][j][0] = compX + compY;
	  }
   }
   // YLEFT
  if (vct->getYleft_neighbor()==MPI_PROC_NULL){
   for (register int i=2;i < nxn-2;i++){
		 compX = (vecFieldXC[i][1][0] - vecFieldXC[i-1][1][0])*invdx;
		 compY = .5*(vecFieldYC[i][2][0] - vecFieldYC[i][1][0])*invdy + .5*(vecFieldYC[i-1][2][0] - vecFieldYC[i-1][1][0])*invdy;
		 divN[i][1][0] = compX + compY;
	 }
  }
  // YRIGHT
  if (vct->getYright_neighbor()==MPI_PROC_NULL){
     for (register int i=2;i < nxn-2;i++){
		 compX = (vecFieldXC[i][nyc-2][0] - vecFieldXC[i-1][nyc-2][0])*invdx;
		 compY = .5*(vecFieldYC[i][nyc-2][0] - vecFieldYC[i][nyc-3][0])*invdy + .5*(vecFieldYC[i-1][nyc-2][0] - vecFieldYC[i-1][nyc-3][0])*invdy;
		 divN[i][nyn-2][0] = compX + compY;
	 } 
  }
  // XLEFTYLEFT
  if (vct->getXleftYleft_neighbor()==MPI_PROC_NULL){
     divN[1][1][0] = divN[2][2][0];
	 divN[1][1][0] = divN[2][2][0];
  }
  // XLEFTYRIGHT
  if (vct->getXleftYright_neighbor()==MPI_PROC_NULL){
     divN[1][nyn-2][0] = divN[2][nyn-3][0]; 
	 divN[1][nyn-2][0] = divN[2][nyn-3][0];
  }
  // XRIGHTYLEFT
  if (vct->getXrightYleft_neighbor()==MPI_PROC_NULL){
     divN[nxn-2][1][0] = divN[nxn-3][2][0];
	 divN[nxn-2][1][0] = divN[nxn-3][2][0]; 
  }
  // XRIGHTYRIGHT
  if (vct->getXrightYright_neighbor()==MPI_PROC_NULL){
     divN[nxn-2][nyn-2][0] = divN[nxn-3][nyn-3][0]; 
	 divN[nxn-2][nyn-2][0] = divN[nxn-3][nyn-3][0];
  }
 
 }
/** calculate curl on nodes, given a vector field defined on central points */
inline void Grid2DCU::curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC,VirtualTopology *vct){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
    for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++){
          // curl - X
          compZDY = .5*(vecFieldZC[i][j][0] - vecFieldZC[i][j-1][0])*invdy +  .5*(vecFieldZC[i-1][j][0] - vecFieldZC[i-1][j-1][0])*invdy;
          // curl - Y
          compZDX = .5*(vecFieldZC[i][j][0] - vecFieldZC[i-1][j][0])*invdx +  .5*(vecFieldZC[i][j-1][0] - vecFieldZC[i-1][j-1][0])*invdx;
          // curl - Z
          compYDX = .5*(vecFieldYC[i][j][0] - vecFieldYC[i-1][j][0])*invdx +  .5*(vecFieldYC[i][j-1][0] - vecFieldYC[i-1][j-1][0])*invdx;
	      compXDY = .5*(vecFieldXC[i][j][0] - vecFieldXC[i][j-1][0])*invdy +  .5*(vecFieldXC[i-1][j][0] - vecFieldXC[i-1][j-1][0])*invdy;

          curlXN[i][j][0] = compZDY;
          curlYN[i][j][0] = - compZDX;
          curlZN[i][j][0] = compYDX - compXDY;
          }
}
/** calculate curl on central points, given a vector field defined on nodes  (nuovo)*/
inline void Grid2DCU::curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
    for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++){
          // curl - X
          compZDY = .5*(vecFieldZN[i][j+1][0] - vecFieldZN[i][j][0])*invdy +  .5*(vecFieldZN[i+1][j+1][0] - vecFieldZN[i+1][j][0])*invdy;
          // curl - Y
          compZDX = .5*(vecFieldZN[i+1][j][0] - vecFieldZN[i][j][0])*invdx +  .5*(vecFieldZN[i+1][j+1][0] - vecFieldZN[i][j+1][0])*invdx;
          // curl - Z
          compYDX = .5*(vecFieldYN[i+1][j][0] - vecFieldYN[i][j][0])*invdx +  .5*(vecFieldYN[i+1][j+1][0] - vecFieldYN[i][j+1][0])*invdx;
          compXDY = .5*(vecFieldXN[i][j+1][0] - vecFieldXN[i][j][0])*invdy +  .5*(vecFieldXN[i+1][j+1][0] - vecFieldXN[i+1][j][0])*invdy;


          curlXC[i][j][0] = compZDY ;
          curlYC[i][j][0] = - compZDX;
          curlZC[i][j][0] = compYDX - compXDY;
      }
 

}
/** calculate curl on central points, given a vector field defined on nodes including ghost*/
inline void Grid2DCU::curlN2C_withghost(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
    for (register int i=0;i < nxn-1;i++)
     for (register int j=0;j < nyn-1;j++){
          // curl - X
          compZDY = .5*(vecFieldZN[i][j+1][0] - vecFieldZN[i][j][0])*invdy +  .5*(vecFieldZN[i+1][j+1][0] - vecFieldZN[i+1][j][0])*invdy;
          // curl - Y
          compZDX = .5*(vecFieldZN[i+1][j][0] - vecFieldZN[i][j][0])*invdx +  .5*(vecFieldZN[i+1][j+1][0] - vecFieldZN[i][j+1][0])*invdx;
          // curl - Z
          compYDX = .5*(vecFieldYN[i+1][j][0] - vecFieldYN[i][j][0])*invdx +  .5*(vecFieldYN[i+1][j+1][0] - vecFieldYN[i][j+1][0])*invdx;
          compXDY = .5*(vecFieldXN[i][j+1][0] - vecFieldXN[i][j][0])*invdy +  .5*(vecFieldXN[i+1][j+1][0] - vecFieldXN[i+1][j][0])*invdy;


          curlXC[i][j][0] = compZDY ;
          curlYC[i][j][0] = - compZDX;
          curlZC[i][j][0] = compYDX - compXDY;
      }
 

}

/** calculate laplacian on nodes, given a scalar field defined on nodes */
inline void Grid2DCU::lapN2N(double*** lapN,double ***scFieldN,VirtualTopology *vct){
   // calculate laplacian as divercence of gradient
   // allocate 3 gradients: defined on central points
   //double*** gradXC = newArr3(double,nxc,nyc,1);
   //double*** gradYC = newArr3(double,nxc,nyc,1);
  double*** gradXC;
  double*** gradYC;

  allocArr3(&gradXC, nxc, nyc, 1);
  allocArr3(&gradYC, nxc, nyc, 1);

   eqValue (0.0, gradXC,nxc,nyc);
   eqValue (0.0, gradYC,nxc,nyc);
   gradN2C(gradXC,gradYC,scFieldN);
   communicateCenter(nxc,nyc,gradXC,vct);
   communicateCenter(nxc,nyc,gradYC,vct);
   divC2N(lapN,gradXC,gradYC,vct);
   // here you need to fix the BC
   // XLEFT
    if (vct->getXleft_neighbor()==MPI_PROC_NULL){
	  for (register int j=1;j < nyn-1;j++)
	   lapN[1][j][0] = (scFieldN[0][j][0] -2*scFieldN[1][j][0] + scFieldN[2][j][0]  )*invdx*invdx + (scFieldN[1][j-1][0]-2*scFieldN[1][j][0] + scFieldN[1][j+1][0])*invdy*invdy;
	}
    // XRIGHT
   if (vct->getXright_neighbor()==MPI_PROC_NULL){
      for (register int j=1;j < nyn-1;j++)
	   lapN[nxn-2][j][0]  = (scFieldN[nxn-3][j][0]  -2*scFieldN[nxn-2][j][0] + scFieldN[nxn-1][j][0])*invdx*invdx + (scFieldN[nxn-2][j-1][0] -2*scFieldN[nxn-2][j][0] +scFieldN[nxn-2][j+1][0])*invdy*invdy;
	}
   // YLEFT
  if (vct->getYleft_neighbor()==MPI_PROC_NULL){
   for (register int i=1;i < nxn-1;i++)
     lapN[i][1][0]  = (scFieldN[i-1][1][0]  -2*scFieldN[i][1][0] + scFieldN[i+1][1][0])*invdx*invdx + (scFieldN[i][0][0]-2*scFieldN[i][1][0]+scFieldN[i][2][0])*invdy*invdy;
  }
  // YRIGHT
  if (vct->getYright_neighbor()==MPI_PROC_NULL){
     for (register int i=1;i < nxn-1;i++)
       lapN[i][nyn-2][0]  = (scFieldN[i-1][nyn-2][0] -2*scFieldN[i][nyn-2][0] +scFieldN[i+1][nyn-2][0] )*invdx*invdx + (scFieldN[i][nyn-3][0] -2*scFieldN[i][nyn-2][0] + scFieldN[i][nyn-1][0])*invdy*invdy;
  }
   // deallocate
   //delArr3(gradXC,nxc,nyc);
   //delArr3(gradYC,nxc,nyc);
  freeArr3(&gradXC);
  freeArr3(&gradYC);
}
/** calculate laplacian on nodes, given a scalar field defined on nodes */
inline void Grid2DCU::lapN2N_plusghost(double*** lapN,double ***scFieldN,VirtualTopology *vct){
   // calculate laplacian as divercence of gradient
   // allocate 3 gradients: defined on central points
   //double*** gradXC = newArr3(double,nxc,nyc,1);
   //double*** gradYC = newArr3(double,nxc,nyc,1);

  double*** gradXC;
  double*** gradYC;

  allocArr3(&gradXC, nxc, nyc, 1);
  allocArr3(&gradYC, nxc, nyc, 1);

   eqValue (0.0, gradXC,nxc,nyc);
   eqValue (0.0, gradYC,nxc,nyc);
   gradN2C_plusghost(gradXC,gradYC,scFieldN);
   communicateCenter(nxc,nyc,gradXC,vct);
   communicateCenter(nxc,nyc,gradYC,vct);
   divC2N(lapN,gradXC,gradYC,vct);
   // deallocate
   //delArr3(gradXC,nxc,nyc);
   //delArr3(gradYC,nxc,nyc);
   freeArr3(&gradXC);
   freeArr3(&gradYC);
}

/** calculate laplacian on central points, given a scalar field defined on central points */
inline void Grid2DCU::lapC2C(double*** lapC,double ***scFieldC,VirtualTopology *vct){
    //  calculate laplacian as divercence of gradient
    // allocate 3 gradients: defined on nodes
    //double*** gradXN = newArr3(double,nxn,nyn,1);
    //double*** gradYN = newArr3(double,nxn,nyn,1);

  double*** gradXN;
  double*** gradYN;

  allocArr3(&gradXN, nxn, nyn, 1);
  allocArr3(&gradYN, nxn, nyn, 1);

    gradC2N(gradXN,gradYN,scFieldC,vct);//if ghost centers OK, active nodes OK
    divN2C(lapC,gradXN,gradYN);//active centers OK, with ghost centers before OK
    // here you need to fix the BC
    // it fixes the first active centers, so it was assumed the ghost centers were not OK
   // XLEFT
    if (vct->getXleft_neighbor()==MPI_PROC_NULL){
	  for (register int j=1;j < nyc-1;j++)
	   lapC[1][j][0] = (scFieldC[0][j][0] -2*scFieldC[1][j][0] + scFieldC[2][j][0]  )*invdx*invdx + (scFieldC[1][j-1][0]-2*scFieldC[1][j][0] + scFieldC[1][j+1][0])*invdy*invdy;
	}
    // XRIGHT
   if (vct->getXright_neighbor()==MPI_PROC_NULL){
      for (register int j=1;j < nyc-1;j++)
	   lapC[nxc-2][j][0]  = (scFieldC[nxc-3][j][0]  -2*scFieldC[nxc-2][j][0] + scFieldC[nxc-1][j][0])*invdx*invdx + (scFieldC[nxc-2][j-1][0] -2*scFieldC[nxc-2][j][0] +scFieldC[nxc-2][j+1][0])*invdy*invdy;
	}
   // YLEFT
  if (vct->getYleft_neighbor()==MPI_PROC_NULL){
   for (register int i=1;i < nxc-1;i++)
     lapC[i][1][0]  = (scFieldC[i-1][1][0]  -2*scFieldC[i][1][0] + scFieldC[i+1][1][0])*invdx*invdx + (scFieldC[i][0][0]-2*scFieldC[i][1][0]+scFieldC[i][2][0])*invdy*invdy;
  }
  // YRIGHT
  if (vct->getYright_neighbor()==MPI_PROC_NULL){
     for (register int i=1;i < nxc-1;i++)
       lapC[i][nyc-2][0]  = (scFieldC[i-1][nyc-2][0] -2*scFieldC[i][nyc-2][0] + scFieldC[i+1][nyc-2][0] )*invdx*invdx + (scFieldC[i][nyc-3][0] -2*scFieldC[i][nyc-2][0] + scFieldC[i][nyc-1][0])*invdy*invdy;
  }
    // deallocate
    //delArr3(gradXN,nxn,nyn);
    //delArr3(gradYN,nxn,nyn);

  freeArr3(&gradXN);
  freeArr3(&gradYN);
}                              
/** calculate laplacian on central points, given a scalar field defined on central points */
/** (dubious) patch for ghost centers **/
inline void Grid2DCU::lapC2C_plusghost(double*** lapC,double ***scFieldC,VirtualTopology *vct){
  //  calculate laplacian as divercence of gradient                                                 
  // allocate 3 gradients: defined on nodes                                                         
  //double*** gradXN = newArr3(double,nxn,nyn,1);
  //double*** gradYN = newArr3(double,nxn,nyn,1);

  double*** gradXN;
  double*** gradYN;
  allocArr3(&gradXN, nxn, nyn, 1);
  allocArr3(&gradYN, nxn, nyn, 1);
  
  gradC2N(gradXN,gradYN,scFieldC,vct);//if ghost centers OK, active nodes OK                        
  //divN2C(lapC,gradXN,gradYN);//active centers OK, with ghost centers before OK
  //Sep12, ME
  // fixes the ghost nodes in the non-boundary procs
  communicateNode(nxn,nyn,gradXN,vct);
  communicateNode(nxn,nyn,gradYN,vct);
  // now ghost centers OK in the non boundary procs
  divN2C_plusghost(lapC,gradXN,gradYN);                      
  //end Sep12, ME
  // here you need to fix the BC of the boundary procs (the other already OK)                                                                    
  // it fixes the ghost centers             
  // XLEFT                     
  if (vct->getXleft_neighbor()==MPI_PROC_NULL){
    for (register int j=1;j < nyc-1;j++)
      //lapC[1][j][0] = (scFieldC[0][j][0] -2*scFieldC[1][j][0] + scFieldC[2][j][0]  )*invdx*invdx + (scFieldC[1][j-1][0]-2*scFieldC[1][j][0] + scFieldC[1][j+1][0])*invdy*invdy;
      lapC[0][j][0] = ( -scFieldC[0][j][0] + scFieldC[1][j][0]  )*invdx*invdx + (scFieldC[0][j-1][0]-2*scFieldC[0][j][0] + scFieldC[0][j+1][0])*invdy*invdy; 
  }
  // XRIGHT                                                                                         
  if (vct->getXright_neighbor()==MPI_PROC_NULL){
    for (register int j=1;j < nyc-1;j++)
      //lapC[nxc-2][j][0]  = (scFieldC[nxc-3][j][0]  -2*scFieldC[nxc-2][j][0] + scFieldC[nxc-1][j][0])*invdx*invdx + (scFieldC[nxc-2][j-1][0] -2*scFieldC[nxc-2][j][0] +scFieldC[nxc-2][j+1][0])*invdy*invdy;
      lapC[nxc-1][j][0]  = (scFieldC[nxc-2][j][0]  -scFieldC[nxc-1][j][0])*invdx*invdx + (scFieldC[nxc-1][j-1][0] -2*scFieldC[nxc-1][j][0] +scFieldC[nxc-1][j+1][0])*invdy*invdy;
  }
  // YLEFT                                                                                           
  if (vct->getYleft_neighbor()==MPI_PROC_NULL){
    for (register int i=1;i < nxc-1;i++)
      //lapC[i][1][0]  = (scFieldC[i-1][1][0]  -2*scFieldC[i][1][0] + scFieldC[i+1][1][0])*invdx*invdx + (scFieldC[i][0][0]-2*scFieldC[i][1][0]+scFieldC[i][2][0])*invdy*invdy;
      lapC[i][0][0]  = (scFieldC[i-1][0][0]  -2*scFieldC[i][0][0] + scFieldC[i+1][0][0])*invdx*invdx + (-scFieldC[i][0][0]+scFieldC[i][1][0])*invdy*invdy; 
  }
  // YRIGHT                                                                                           
  if (vct->getYright_neighbor()==MPI_PROC_NULL){
    for (register int i=1;i < nxc-1;i++)
      //lapC[i][nyc-2][0]  = (scFieldC[i-1][nyc-2][0] -2*scFieldC[i][nyc-2][0] + scFieldC[i+1][nyc-2][0] )*invdx*invdx + (scFieldC[i][nyc-3][0] -2*scFieldC[i][nyc-2][0] + scFieldC[i][nyc-1][0])*invdy*invdy;
      lapC[i][nyc-1][0]  = (scFieldC[i-1][nyc-1][0] -2*scFieldC[i][nyc-1][0] + scFieldC[i+1][nyc-1][0] )*invdx*invdx + (scFieldC[i][nyc-2][0] -scFieldC[i][nyc-1][0])*invdy*invdy; 
  }

  // for the ghost corner of all boundary procs, but the corner procs
  communicateCenter(nxc,nyc,lapC,vct);

  //corners
  //if (vct->getXleft_neighbor()==MPI_PROC_NULL || vct->getYleft_neighbor()==MPI_PROC_NULL)
    lapC[0][0][0] = ( -scFieldC[0][0][0] + scFieldC[1][0][0]  )*invdx*invdx + (-scFieldC[0][0][0] + scFieldC[0][1][0])*invdy*invdy;
    //if (vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL)
    lapC[nxc-1][nyc-1][0]  = (scFieldC[nxc-2][nyc-1][0]  -scFieldC[nxc-1][nyc-1][0])*invdx*invdx + (scFieldC[nxc-1][nyc-2][0] -scFieldC[nxc-1][nyc-1][0])*invdy*invdy;
    //if (vct->getXright_neighbor()==MPI_PROC_NULL && vct->getYleft_neighbor()==MPI_PROC_NULL)
    lapC[nxc-1][0][0]  = (scFieldC[nxc-2][0][0]  -scFieldC[nxc-1][0][0])*invdx*invdx + (-scFieldC[nxc-1][0][0]+scFieldC[nxc-1][1][0])*invdy*invdy;
    //if (vct->getXleft_neighbor()==MPI_PROC_NULL && vct->getYright_neighbor()==MPI_PROC_NULL)
    lapC[0][nyc-1][0]  = ( -scFieldC[0][nyc-1][0] + scFieldC[1][nyc-1][0] )*invdx*invdx + (scFieldC[0][nyc-2][0] -scFieldC[0][nyc-1][0])*invdy*invdy;
    
  // deallocate                                                                                     
  //delArr3(gradXN,nxn,nyn);
  //delArr3(gradYN,nxn,nyn);
    freeArr3(&gradXN);
    freeArr3(&gradYN);
}
/** interpolate on nodes from central points: do this for the magnetic field */
inline void Grid2DCU::interpC2N(double ***vecFieldN, double ***vecFieldC, VirtualTopology *vct){
  // internal point
  for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
         vecFieldN[i][j][0] = (vecFieldC[i][j][0] + vecFieldC[i-1][j][0] + vecFieldC[i][j-1][0] + vecFieldC[i-1][j-1][0])*.25;

  // communicate
  communicateNode(nxn,nyn,vecFieldN,vct);

}
/** interpolate on nodes from central points: do this for the magnetic field */
inline void Grid2DCU::interpC2N_BC(double ***vecFieldN, double ***vecFieldC, VirtualTopology *vct){
  // internal point
  for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
         vecFieldN[i][j][0] = (vecFieldC[i][j][0] + vecFieldC[i-1][j][0] + vecFieldC[i][j-1][0] + vecFieldC[i-1][j-1][0])*.25;
  // fix the boundaries
  // XLEFT
  if (vct->getXleft_neighbor()==MPI_PROC_NULL){
     for (register int j=2;j < nyn-2;j++)
         vecFieldN[1][j][0] = (vecFieldC[1][j][0] + vecFieldC[1][j-1][0])*.5;
  }
  // XRIGHT
  if (vct->getXright_neighbor()==MPI_PROC_NULL){
	for (register int j=2;j < nyn-2;j++)
         vecFieldN[nxn-2][j][0] = (vecFieldC[nxc-2][j][0] + vecFieldC[nxc-2][j-1][0])*.5;
  }
   // YLEFT
  if (vct->getYleft_neighbor()==MPI_PROC_NULL){
       for (register int i=2;i < nxn-2;i++)
         vecFieldN[i][1][0] = (vecFieldC[i][1][0] + vecFieldC[i-1][1][0])*.5;
  }
  // YRIGHT
  if (vct->getYright_neighbor()==MPI_PROC_NULL){
     for (register int i=2;i < nxn-2;i++)
         vecFieldN[i][nyn-2][0] = (vecFieldC[i][nyc-2][0] + vecFieldC[i-1][nyc-2][0])*.5;
  }
  // XLEFTYLEFT
  if (vct->getXleftYleft_neighbor()==MPI_PROC_NULL){
     vecFieldN[1][1][0] = .33333333*(vecFieldN[2][1][0] + vecFieldN[2][2][0] + vecFieldN[1][2][0]); // interpolate from nodes
  }
  // XLEFTYRIGHT
  if (vct->getXleftYright_neighbor()==MPI_PROC_NULL){
     vecFieldN[1][nyn-2][0] = .33333333*(vecFieldN[2][nyn-2][0] + vecFieldN[2][nyn-3][0] + vecFieldN[1][nyn-3][0]); // interpolate from nodes
  }
  // XRIGHTYLEFT
  if (vct->getXrightYleft_neighbor()==MPI_PROC_NULL){
     vecFieldN[nxn-2][1][0] = .33333333*(vecFieldN[nxn-2][2][0] + vecFieldN[nxn-3][2][0] + vecFieldN[nxn-3][1][0]); // interpolate from nodes
  }
  // XRIGHTYRIGHT
  if (vct->getXrightYright_neighbor()==MPI_PROC_NULL){
     vecFieldN[nxn-2][nyn-2][0] = .33333333*(vecFieldN[nxn-2][nyn-3][0] + vecFieldN[nxn-3][nyn-2][0] + vecFieldN[nxn-3][nyn-3][0]); // interpolate from nodes
  }
  // communicate
  communicateNode(nxn,nyn,vecFieldN,vct);

}


/** interpolate on nodes from central points: do this for the magnetic field */
inline void Grid2DCU::interpC2N_BC_alsoGN(double ***vecFieldN, double ***vecFieldC, VirtualTopology *vct){
  // modified to fix BC on the outer crown
  // internal point
  for (register int i=1;i < nxn-1;i++)
    for (register int j=1;j < nyn-1;j++)
      vecFieldN[i][j][0] = (vecFieldC[i][j][0] + vecFieldC[i-1][j][0] + vecFieldC[i][j-1][0] + vecFieldC[i-1][j-1][0])*.25;

  // fix the boundaries
  // XLEFT
  if (vct->getXleft_neighbor()==MPI_PROC_NULL){
    for (register int j=1;j < nyn-1;j++)
      vecFieldN[0][j][0] = (vecFieldC[0][j][0] + vecFieldC[0][j-1][0])*.5;
  }
  // XRIGHT
  if (vct->getXright_neighbor()==MPI_PROC_NULL){
    for (register int j=1;j < nyn-1;j++)
      vecFieldN[nxn-1][j][0] = (vecFieldC[nxc-1][j][0] + vecFieldC[nxc-1][j-1][0])*.5;
  }
  // YLEFT
  if (vct->getYleft_neighbor()==MPI_PROC_NULL){
    for (register int i=1;i < nxn-1;i++)
      vecFieldN[i][0][0] = (vecFieldC[i][0][0] + vecFieldC[i-1][0][0])*.5;
  }
  // YRIGHT
  if (vct->getYright_neighbor()==MPI_PROC_NULL){
    for (register int i=1;i < nxn-1;i++)
      vecFieldN[i][nyn-1][0] = (vecFieldC[i][nyc-1][0] + vecFieldC[i-1][nyc-1][0])*.5;
  }
  // XLEFTYLEFT
  if (vct->getXleftYleft_neighbor()==MPI_PROC_NULL){
    //has to be fixed in any case, not only for the corner processor 
  vecFieldN[0][0][0] = .33333333*(vecFieldN[1][0][0] + vecFieldN[1][1][0] + vecFieldN[0][1][0]); // interpolate from nodes
    }
  // XLEFTYRIGHT
  if (vct->getXleftYright_neighbor()==MPI_PROC_NULL){
  //has to be fixed in any case, not only for the corner processor
  vecFieldN[0][nyn-1][0] = .33333333*(vecFieldN[1][nyn-1][0] + vecFieldN[1][nyn-2][0] + vecFieldN[0][nyn-2][0]); // interpolate from nodes
    }
  // XRIGHTYLEFT
  if (vct->getXrightYleft_neighbor()==MPI_PROC_NULL){
  //has to be fixed in any case, not only for the corner processor
    vecFieldN[nxn-1][0][0] = .33333333*(vecFieldN[nxn-1][1][0] + vecFieldN[nxn-2][1][0] + vecFieldN[nxn-2][0][0]); // interpolate from nodes
    }
  // XRIGHTYRIGHT
  if (vct->getXrightYright_neighbor()==MPI_PROC_NULL){
    //has to be fixed in any case, not only for the corner processor
    vecFieldN[nxn-1][nyn-1][0] = .33333333*(vecFieldN[nxn-1][nyn-2][0] + vecFieldN[nxn-2][nyn-1][0] + vecFieldN[nxn-2][nyn-2][0]); // interpolate from nodes
  }
  // communicate
  communicateNode(nxn,nyn,vecFieldN,vct);

}


/** interpolate on central points from nodes*/
inline void Grid2DCU::interpN2C(double ***vecFieldC, double ***vecFieldN){
  
 for (register int i=1;i < nxc-1;i++)
	   for (register int j=1;j < nyc-1;j++){
	   vecFieldC[i][j][0] = .25*(vecFieldN[i+1][j][0] + vecFieldN[i][j][0] + vecFieldN[i+1][j+1][0] + vecFieldN[i][j+1][0]);
	   }
}

inline void Grid2DCU::interpN2C_alsoGC(double ***vecFieldC, double ***vecFieldN){

  for (register int i=0;i < nxc;i++)
    for (register int j=0;j < nyc;j++){
      vecFieldC[i][j][0] = .25*(vecFieldN[i+1][j][0] + vecFieldN[i][j][0] + vecFieldN[i+1][j+1][0] + vecFieldN[i][j+1][0]);
    }

}
/** interpolate on central points from nodes per species */
inline void Grid2DCU::interpN2C(double ****vecFieldC, int ns, double ****vecFieldN){
 
  for (register int i=1;i < nxc-1;i++)
     for (register int j=1;j < nyc-1;j++)
         vecFieldC[ns][i][j][0] = .25*(vecFieldN[ns][i+1][j][0] + vecFieldN[ns][i][j][0] + vecFieldN[ns][i+1][j+1][0] + vecFieldN[ns][i][j+1][0]);
 
}
// interpolate on central points from nodes including ghost cells*/
inline void Grid2DCU::interpN2C_alsoGC(double ****vecFieldC, int ns, double ****vecFieldN){

  for (register int i=0;i < nxc;i++)
    for (register int j=0;j < nyc;j++)
      vecFieldC[ns][i][j][0] = .25*(vecFieldN[ns][i+1][j][0] + vecFieldN[ns][i][j][0] + vecFieldN[ns][i+1][j+1][0] + vecFieldN[ns][i][j+1][0]);
}

/** get level - level of the current grid*/
inline int Grid2DCU::getLevel(){
  return(level);
}
/** get Ratio*/
inline double Grid2DCU::getRatio(){
  return(ratio);
}

/** get grid origin  in x - level of the current grid*/
inline double &Grid2DCU::getOx(int level){
  return(ox[level]);
}

/** get grid origin  in y - level of the current grid*/
inline double &Grid2DCU::getOy(int level){
  return(oy[level]);
}

/** get nxc - number of cells in the X-DIRECTION*/
inline int Grid2DCU::getNXC(){
  return(nxc);
}
/** get nxn - number of nodes in the X-DIRECTION*/
inline int Grid2DCU::getNXN(){
  return(nxn);
}
/** get nyc - number of cells in the Y-DIRECTION*/
inline int Grid2DCU::getNYC(){
  return(nyc);
}
/** get nyn - number of nodes in the Y-DIRECTION*/
inline int Grid2DCU::getNYN(){
  return(nyn);
}

/** get dx - grid spacing */
inline double Grid2DCU::getDX(){
  return(dx);
}
/** get dy - grid spacing*/
inline double Grid2DCU::getDY(){
  return(dy);
}

/** get x-coordinate of the node(X,Y,Z) */
inline double &Grid2DCU::getXN(int indexX, int indexY, int indexZ){
  return(xn[indexX][indexY][0]);
}
/** get y-coordinate of the node(X,Y,Z) */
inline double &Grid2DCU::getYN(int indexX, int indexY, int indexZ){
  return(yn[indexX][indexY][0]);
}
/** get x-coordinate of the node(X,Y,Z), treats also the node following the last ghost node (for interpP2G_OS) */
inline double Grid2DCU::getModifiedXN(int indexX, int indexY, int indexZ){
  //cout << "indexX " << indexX << "indexX " << indexX <<endl;
  if (indexX== nxn)
    return xEnd +2* dx; // xend is the end of the active part of the grid
  if (indexX== -1)
    return -2*dx;
  if (indexY == -1 || indexY== nyn)
    indexY=0; //some dummy value, the coord which matters is x
  return(xn[indexX][indexY][0]);
}
/** get y-coordinate of the node(X,Y,Z), treats also the node following the last ghost node (for interpP2G_OS) */
inline double Grid2DCU::getModifiedYN(int indexX, int indexY, int indexZ){
  //cout << "indexX " << indexX << "indexX " << indexX <<endl;
  if (indexY== nyn)
    return yEnd +2* dy;
  if (indexY== -1)
    return -2*dy;
  if (indexX == -1 || indexX== nxn)
    indexX=0; //some dummy value, the coord which matters is y 
  return(yn[indexX][indexY][0]);
}
/** get x-coordinate of the center of cell(X,Y,Z) */
inline double &Grid2DCU::getXC(int indexX, int indexY, int indexZ){
  return(xc[indexX][indexY][0]);
}
/** get y-coordinate of the center of cell(X,Y,Z) */
inline double &Grid2DCU::getYC(int indexX, int indexY, int indexZ){
  return(yc[indexX][indexY][0]);
}

/** get the whole vector xc with coordinate of center cells*/
inline double*** Grid2DCU::getXC(){
  return(xc);
}
/** get the whole vector yc with coordinate of center cells*/
inline double*** Grid2DCU::getYC(){
  return(yc);
}
/** get x coordinate of the first node of Simulation box - X direction */
inline double Grid2DCU::getXstart(){
 return(xStart);
}
/** get x coordinate of the first node of Simulation box - X direction or of node 0 if this proc is on the x=0 border */
inline double Grid2DCU::getmodifiedXstart(VirtualTopology *vct){
 double res;
 res = xStart;
 if (vct->getCoordinates(0)==0) {
 res = xStart - dx;
 }
 return(res);
}
/** get x coordinate of the first node of Simulation box - X direction or of the VERY last node if this proc is on the x=Lx border */
inline double Grid2DCU::getmodifiedXend(VirtualTopology *vct){
 double res;
 res = xEnd;
 if (vct->getCoordinates(0)==vct->getXLEN()-1) {
 res = xEnd + dx;
 }
 return(res);
}
/** get x coordinate of the last node of Simulation box - X direction*/
inline double Grid2DCU::getXend(){
 return(xEnd);
}
/** get y coordinate of the first node of Simulation box - Y direction */
inline double Grid2DCU::getYstart(){
 return(yStart);
}
/** get y coordinate of the last node of Simulation box - Y direction */
inline double Grid2DCU::getYend(){
 return(yEnd);
}
/** get y coordinate of the first node of Simulation box - Y direction or of node 0 if this proc is on the y=0 border */
inline double Grid2DCU::getmodifiedYstart(VirtualTopology *vct){
 double res;
 res = yStart;
 if (vct->getCoordinates(1)==0) {
 res = yStart - dy;
 }
 return(res);
}
/** get y coordinate of the first node of Simulation box - Y direction or of the VERY last node if this proc is on the x=Ly border */
inline double Grid2DCU::getmodifiedYend(VirtualTopology *vct){
 double res;
 res = yEnd;
 if (vct->getCoordinates(1)==vct->getYLEN()-1) {
 res = yEnd + dy;
 }
 return(res);
}

/** get the inverse of volume */
inline double Grid2DCU::getInvVOL(){
 return(invVOL);
}





/** get nyc - number of cells in the Z-DIRECTION*/
inline int Grid2DCU::getNZC(){
  cout << "No need to call getNZC(). You are using 2D grid XY" << endl;
  return(1);
}
/** get nzn - number of nodes in the Z-DIRECTION*/
inline int Grid2DCU::getNZN(){
   cout << "No need to call getNZN(). You are using 2D grid XY" << endl;
   return(1);
}
/** get dz - grid spacing*/
inline double Grid2DCU::getDZ(){
   cout << "No need to call getDZ(). You are using 2D grid XY" << endl;
   return(1.0);
}
/** get z-coordinate of the center of cell(X,Y,Z) */
inline double &Grid2DCU::getZC(int indexX, int indexY, int indexZ){
   cout << "No need to call getZC(int indexX, int indexY, int indexZ). You are using 2D grid XY" << endl;
   return(xc[0][0][0]);
}
/** get z-coordinate of the node(X,Y,Z) */
inline double &Grid2DCU::getZN(int indexX, int indexY, int indexZ){
    cout << "No need to call getZN(int indexX, int indexY, int indexZ). You are using 2D grid XY" << endl;
   return(xn[0][0][0]);
}
/** get the whole vector zc with coordinate of center cells*/
inline double*** Grid2DCU::getZC(){
  cout << "No need to call getZC(). You are using 2D grid XY" << endl;
  return(NULL);
}

/** get y coordinate of the first node of Simulation box - Z direction */
inline double Grid2DCU::getZstart(){
  cout << "No need to call getZstart(). You are using 2D grid XY" << endl;
  return(0.0);
}
/** get y coordinate of the last node of Simulation box - Z direction */
inline double Grid2DCU::getZend(){
  cout << "No need to call getZend(). You are using 2D grid XY" << endl;
  return(0.0);
}
#endif

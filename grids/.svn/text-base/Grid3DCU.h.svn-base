 /*******************************************************************************************************
 Grid3DCU.h  -  uniform cartesian 2D local grid for each process, including che guard cells
                             -------------------
 developers: Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
 *******************************************************************************************************/

#ifndef Grid3DCU_H
#define Grid3DCU_H

#include <iostream>

#include "Grid.h"
#include "../communication/ComInterpNodes.h"
#include "../communication/ComNodes.h"
#include "../utility/Alloc.h"

using std::cout;
using std::endl;

/**
* 3D Uniform cartesian local(each processor has its own grid) grid 
*
* Note that Grid3DCU is implementing the abstract class Grid.h, therefore
* it must implements all the virtual methods in Grid.h
*
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Grid3DCU : public Grid {
  public:
      /** constructor */
      Grid3DCU(CollectiveIO *col, VirtualTopology *vct);
      /** destructor */
      ~Grid3DCU();
      /** allocate grid arrays for this domain */
      void allocate(CollectiveIO *ptC, VirtualTopology *ptVCT);
      /** deallocate grid arrays for this domain */
      void deallocate();
      /** print grid info */
      void print(VirtualTopology* ptVCT);
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradC2N(double*** gradXN,double*** gradYN,double***scFieldC );
	  /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradC2N(double*** gradXN,double*** gradYN,double*** gradZN,double***scFieldC );
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradN2C(double ***gradXC, double ***gradYC, double ***scFieldN);
	   /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradN2C(double ***gradXC, double ***gradYC, double ***gradZC,double ***scFieldN);
      /** calculate divergence on central points, given a vector field defined on nodes  */
      void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN);
	  /** calculate divergence on central points, given a vector field defined on nodes  */
      void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN,  double ***vecFieldZN);
      /** calculate divergence on nodes, given a vector field defined on central points  */
      void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC);
	  /** calculate divergence on nodes, given a vector field defined on central points  */
      void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC);
      /** calculate curl on nodes, given a vector field defined on central points  */
      void curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC);
      /** calculate curl on central points, given a vector field defined on nodes  */
      void curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN);

      /** calculate divergence on central points, given a Tensor field defined on nodes  */
      void divSymmTensorN2C(double ***divCX, double ***divCY, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ,int ns);
      void divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ,int ns);
      
      /** calculate laplacian on nodes, given a scalar field defined on nodes */
      void lapN2N(double*** lapN,double ***scFieldN,VirtualTopology *vct);
      /** calculate laplacian on central points, given a scalar field defined on central points */
      void lapC2C(double*** lapC,double ***scFieldC,VirtualTopology *vct);

      /** calculate divergence on boundaries */
      void divBCleft(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dirDER);
      /** calculate divergence on boundaries */
      void divBCright(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ, int rightActiveNode, int dirDER);
      /** calculate derivative on boundaries */
      void derBC(double ***derBC, double ***vector, int leftActiveNode, int dirDER);
      
      
      /** interpolate on nodes from central points */
      void interpC2N(double ***vecFieldN, double ***vecFieldC);
      /** interpolate on central points from nodes */
      void interpN2C(double ***vecFieldC, double ***vecFieldN);
      /** interpolate on central points from nodes per species*/
      void interpN2C(double ****vecFieldC, int ns, double ****vecFieldN);
      
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
      /** get x coordinate of the last node of Simulation box - X direction*/
      double getXend();
      /** get y coordinate of the first node of Simulation box - Y direction */
      double getYstart();
      /** get y coordinate of the last node of Simulation box - Y direction */
      double getYend();
      /** get z coordinate of the first node of Simulation box - Z direction */
      double getZstart();
      /** get z coordinate of the last node of Simulation box - Z direction */
      double getZend();

      /** get the inverse of volume */
      double getInvVOL();

  private:
     /** number of cells - X direction, including + 2 (guard cells) */
     int nxc;
     /** number of nodes - X direction, including + 2 extra nodes for guard cells */
     int nxn;
     /** number of cell - Y direction, including + 2 (guard cells) */
     int nyc;
     /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
     int nyn;
	  /** number of cell - Z direction, including + 2 (guard cells) */
     int nzc;
     /** number of nodes - Z direction, including + 2 extra nodes for guard cells */
     int nzn;
     /** dx = grid spacing - X direction */
     double dx;
     /** dy = grid spacing - Y direction */
     double dy;
	  /** dz = grid spacing -Z direction */
     double dz;
     /** invdx = 1/dx */
     double invdx;
     /** invdy = 1/dy */
     double invdy;
	  /** invdz = 1/dz */
     double invdz;
     /** invol = inverse of volume*/
     double invVOL;
     /** node - X coordinate   */
     double ***xn;
     /** node - Y coordinate   */
     double ***yn;
	  /** node - Z coordinate   */
     double ***zn;
     /** centre of cell - X coordinate */
     double ***xc;
     /** centre of cell - Y coordinate */
     double ***yc;
	 /** centre of cell - Y coordinate */
     double ***yc;
     /** local grid boundaries coordinates  */
     double xStart, xEnd, yStart, yEnd, zStart, zEnd;
     
};
/** constructor */
inline Grid3DCU::Grid3DCU(CollectiveIO* col, VirtualTopology* vct){
   // add 2 for the guard cells
   nxc = (col->getNxc())/(vct->getXLEN()) + 2;
   nyc = (col->getNyc())/(vct->getYLEN()) + 2;
   nzc = col->getNzc() + 2; // we use only 2D topology
   nxn = nxc + 1;
   nyn = nyc + 1;
   nzn = nzc + 1;
   dx  = col->getLx()/col->getNxc();
   dy  = col->getLy()/col->getNyc();
   dz  = col->getLz()/col->getNzc();
   invVOL = 1.0/(dx*dy*dz);
   invdx = 1.0/dx;
   invdy = 1.0/dy;
   invdz = 1.0/dz;
   // local grid dimensions and boundaries of active nodes
   xStart = vct->getCoordinates(0)*(col->getLx()/(double) vct->getXLEN());
   xEnd   = xStart + (col->getLx()/(double) vct->getXLEN());
   yStart = vct->getCoordinates(1)*(col->getLy()/(double) vct->getYLEN());
   yEnd   = yStart + (col->getLy()/(double) vct->getYLEN());
   zStart = 0;
   zEnd   = Lz;
   // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
   xn = newArr3(double,nxn,nyn,nzn);
   yn = newArr3(double,nxn,nyn,nzn);
   zn = newArr3(double,nxn,nyn,nzn);
   for (int i=0; i < nxn; i++){
     for (int j=0; j < nyn; j++){
	    for (int k=0; k < nzn; k++){ 
            xn[i][j][k] = xStart + (i-1)*dx;
            yn[i][j][k] = yStart + (j-1)*dy;
			zn[i][j][k] = zStart + (k-1)*dz;
		  }
        }
     }
   
   // arrays allocation: cells ---> the first cell has index 1, the last has index ncn-2!
   xc = newArr3(double,nxc,nyc,nzc);
   yc = newArr3(double,nxc,nyc,nzc);
   zc = newArr3(double,nxc,nyc,nzc);
   for (int i=0; i < nxc; i++){
     for (int j=0; j < nyc; j++){
	   for (int k=0; k < nzc; k++){
          xc[i][j][k] = .5*(xn[i][j][k] + xn[i+1][j][k]);
          yc[i][j][k] = .5*(yn[i][j][k] + yn[i][j+1][k]);
		  zc[i][j][k] = .5*(zn[i][j][k] + zn[i][j][k+1]);
	   }
         
     }
   }
}

/** deallocate the local grid */
inline Grid3DCU::~Grid3DCU(){
   // deallocate nodes
   delArr3(xn,nxn,nyn);
   delArr3(yn,nxn,nyn);
   delArr3(zn,nxn,nyn);
   // centers cells
   delArr3(xc,nxc,nyc);
   delArr3(yc,nxc,nyc);
   delArr3(zc,nxc,nyc);
   
}
/** print the local grid info */
inline void Grid3DCU::print(VirtualTopology* ptVCT){
    cout << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ")"<< endl;
    cout <<  "Number of cell: -X=" << nxc-2 << " -Y=" << nyc-2 << " -Z =" << nzc -2 << endl;
    cout <<  "Xin = " << xn[1][1][1] << "; Xfin = " << xn[nxn-2][1][1] << endl;
    cout <<  "Yin = " << yn[1][1][1] << "; Yfin = " << yn[1][nyn-2][1] << endl;
	cout <<  "Zin = " << zn[1][1][1] << "; Zfin = " << zn[1][1][nzn-2] << endl;
    cout << endl;
    
}
/** calculate gradient on nodes in 3D, given a scalar field defined on central points*/
inline void Grid3DCU::gradC2N(double***  gradXN, double***  gradYN, double***  gradZN, double*** scFieldC){
    for (register int i=1;i < nxc;i++)
	    for (register int j=1;j < nyc; j++)
		   for (register int k=1;k < nzc; k++){
            gradXN[i][j][k] = .25*(scFieldC[i][j][k] - scFieldC[i-1][j][k])*invdx +  .25*(scFieldC[i][j-1][k] - scFieldC[i-1][j-1][k])*invdx +  .25*(scFieldC[i][j][k-1] - scFieldC[i-1][j][k-1])*invdx +  .25*(scFieldC[i][j-1][k-1] - scFieldC[i-1][j-1][k-1])*invdx;
            gradYN[i][j][k] = .25*(scFieldC[i][j][k] - scFieldC[i][j-1][k])*invdy +  .25*(scFieldC[i-1][j][k] - scFieldC[i-1][j-1][k])*invdy +  .25*(scFieldC[i][j][k-1] - scFieldC[i][j-1][k-1])*invdy +  .25*(scFieldC[i-1][j][k-1] - scFieldC[i-1][j-1][k-1])*invdy;
            gradZN[i][j][k] = .25*(scFieldC[i][j][k] - scFieldC[i][j][k-1])*invdz +  .25*(scFieldC[i-1][j][k] - scFieldC[i-1][j][k-1])*invdz +  .25*(scFieldC[i][j-1][k] - scFieldC[i][j-1][k-1])*invdz +  .25*(scFieldC[i-1][j-1][k] - scFieldC[i-1][j-1][k-1])*invdz
		 }
    

}

/** calculate gradient on nodes in 3D, given a scalar field defined on central points  */
inline void Grid3DCU::gradN2C(double ***gradXC, double ***gradYC, double ***gradZC,double ***scFieldN){
   for (register int i=0;i < nxn-1; i++)
     for (register int j=0;j < nyn-1; j++) 
	   for (register int k=0;k < nzn-1; k++){
          gradXC[i][j][k] = .25*(scFieldN[i+1][j][k] - scFieldN[i][j][k])*invdx + .25*(scFieldN[i+1][j+1][k] - scFieldN[i][j+1][k])*invdx + .25*(scFieldN[i+1][j][k+1] - scFieldN[i][j][k+1])*invdx + .25*(scFieldN[i+1][j+1][k+1] - scFieldN[i][j+1][k+1])*invdx;
          gradYC[i][j][k] = .25*(scFieldN[i][j+1][k] - scFieldN[i][j][k])*invdy + .25*(scFieldN[i+1][j+1][k] - scFieldN[i+1][j][k])*invdy + .25*(scFieldN[i][j+1][k+1] - scFieldN[i][j][k+1])*invdy + .25*(scFieldN[i+1][j+1][k+1] - scFieldN[i+1][j][k+1])*invdy;
          gradZC[i][j][k] = .25*(scFieldN[i][j][k+1] - scFieldN[i][j][k])*invdz + .25*(scFieldN[i+1][j][k+1] - scFieldN[i+1][j][k])*invdz + .25*(scFieldN[i][j+1][k+1] - scFieldN[i][j+1][k])*invdz + .25*(scFieldN[i+1][j+1][k+1] - scFieldN[i+1][j+1][k])*invdz;
	   }
}

/** calculate divergence on central points, given a vector field defined on nodes  */
inline void Grid3DCU::divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN){
  double compX;
  double compY;
  double compZ;
   for (register int i=0;i < nxn-1 ;i++)
     for (register int j=0;j < nyn-1 ;j++)
	   for (register int k=0;k < nzn-1 ;k++){
          compX = .25*(vecFieldXN[i+1][j][k] - vecFieldXN[i][j][k])*invdx +  .25*(vecFieldXN[i+1][j+1][k] - vecFieldXN[i][j+1][k])*invdx +  .25*(vecFieldXN[i+1][j][k+1] - vecFieldXN[i][j][k+1])*invdx +  .25*(vecFieldXN[i+1][j+1][k+1] - vecFieldXN[i][j+1][k+1])*invdx;
          compY = .25*(vecFieldYN[i][j+1][k] - vecFieldYN[i][j][k])*invdy +  .25*(vecFieldYN[i+1][j+1][k] - vecFieldYN[i+1][j][k])*invdy +  .25*(vecFieldYN[i][j+1][k+1] - vecFieldYN[i][j][k+1])*invdy +  .25*(vecFieldYN[i+1][j+1][k+1] - vecFieldYN[i+1][j][k+1])*invdy;
          compZ = .25*(vecFieldZN[i][j][k+1] - vecFieldZN[i][j][k])*invdz +  .25*(vecFieldZN[i+1][j][k+1] - vecFieldZN[i+1][j][k])*invdz +  .25*(vecFieldZN[i][j+1][k+1] - vecFieldZN[i][j+1][k])*invdz +  .25*(vecFieldZN[i+1][j+1][k+1] - vecFieldZN[i+1][j+1][k])*invdz;
		  divC[i][j][k] = compX + compY + compZ;
      }      
}
/** calculate divergence on central points, given a Tensor field defined on nodes */
inline void Grid3DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double**** pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
  double comp1X, comp2X, comp3X;
  double comp1Y, comp2Y, comp3Y;
  double comp1Z, comp2Z, comp3Z;
   for (register int i=0;i < nxn-1;i++)
     for (register int j=0;j < nyn-1;j++)
	    for (register int k=0;k < nzn-1;k++){
          // X
		  comp1X = .25*(pXX[ns][i+1][j][k] - pXX[ns][i][j][k])*invdx +  .25*(pXX[ns][i+1][j+1][k] - pXX[ns][i][j+1][k])*invdx +  .25*(pXX[ns][i+1][j][k+1] - pXX[ns][i][j][k+1])*invdx +  .25*(pXX[ns][i+1][j+1][k+1] - pXX[ns][i][j+1][k+1])*invdx;
          comp2X = .25*(pXY[ns][i+1][j][k] - pXY[ns][i][j][k])*invdx +  .25*(pXY[ns][i+1][j+1][k] - pXY[ns][i][j+1][k])*invdx +  .25*(pXY[ns][i+1][j][k+1] - pXY[ns][i][j][k+1])*invdx +  .25*(pXY[ns][i+1][j+1][k+1] - pXY[ns][i][j+1][k+1])*invdx;
          comp3X = .25*(pXZ[ns][i+1][j][k] - pXZ[ns][i][j][k])*invdx +  .25*(pXZ[ns][i+1][j+1][k] - pXZ[ns][i][j+1][k])*invdx +  .25*(pXZ[ns][i+1][j][k+1] - pXZ[ns][i][j][k+1])*invdx +  .25*(pXZ[ns][i+1][j+1][k+1] - pXZ[ns][i][j+1][k+1])*invdx;
          // Y
		  comp1Y = .25*(pXY[ns][i][j+1][k] - pXY[ns][i][j][k])*invdy +  .25*(pXY[ns][i+1][j+1][k] - pXY[ns][i+1][j][k])*invdy +  .25*(pXY[ns][i][j+1][k+1] - pXY[ns][i][j][k+1])*invdy +  .25*(pXY[ns][i+1][j+1][k+1] - pXY[ns][i+1][j][k+1])*invdy;
          comp2Y = .25*(pYY[ns][i][j+1][k] - pYY[ns][i][j][k])*invdy +  .25*(pYY[ns][i+1][j+1][k] - pYY[ns][i+1][j][k])*invdy +  .25*(pYY[ns][i][j+1][k+1] - pYY[ns][i][j][k+1])*invdy +  .25*(pYY[ns][i+1][j+1][k+1] - pYY[ns][i+1][j][k+1])*invdy;
          comp3Y = .25*(pYZ[ns][i][j+1][k] - pYZ[ns][i][j][k])*invdy +  .25*(pYZ[ns][i+1][j+1][k] - pYZ[ns][i+1][j][k])*invdy +  .25*(pYZ[ns][i][j+1][k+1] - pYZ[ns][i][j][k+1])*invdy +  .25*(pYZ[ns][i+1][j+1][k+1] - pYZ[ns][i+1][j][k+1])*invdy;
          // Z
		  comp1Z = .25*(pXZ[ns][i][j][k+1] - pXZ[ns][i][j][k])*invdz +  .25*(pXZ[ns][i+1][j][k+1] - pXZ[ns][i+1][j][k])*invdz +  .25*(pXZ[ns][i][j+1][k+1] - pXZ[ns][i][j+1][k])*invdz +  .25*(pXZ[ns][i+1][j+1][k+1] - pXZ[ns][i+1][j+1][k])*invdz;
		  comp2Z = .25*(pYZ[ns][i][j][k+1] - pYZ[ns][i][j][k])*invdz +  .25*(pYZ[ns][i+1][j][k+1] - pYZ[ns][i+1][j][k])*invdz +  .25*(pYZ[ns][i][j+1][k+1] - pYZ[ns][i][j+1][k])*invdz +  .25*(pYZ[ns][i+1][j+1][k+1] - pYZ[ns][i+1][j+1][k])*invdz;
		  comp3Z = .25*(pZZ[ns][i][j][k+1] - pZZ[ns][i][j][k])*invdz +  .25*(pZZ[ns][i+1][j][k+1] - pZZ[ns][i+1][j][k])*invdz +  .25*(pZZ[ns][i][j+1][k+1] - pZZ[ns][i][j+1][k])*invdz +  .25*(pZZ[ns][i+1][j+1][k+1] - pZZ[ns][i+1][j+1][k])*invdz;
		  // calculate the divergences
          divCX[i][j][k] = comp1X + comp1Y + comp1Z;
          divCY[i][j][k] = comp2X + comp2Y + comp2Z;
	      divCZ[i][j][k] = comp3X + comp3Y + comp3Z;
           
     }


}
/** calculate divergence on nodes, given a vector field defined on central points  */
inline void Grid3DCU::divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC, double ***vecFieldZC){
  double compX;
  double compY;
  double compZ;
  for (register int i=1; i < nxc;i++)
	  for (register int j=1; j < nyc;j++)
	   for (register int k=1; k < nzc;k++) {
          compX = .25*(vecFieldXC[i][j][k] - vecFieldXC[i-1][j][k])*invdx + .25*(vecFieldXC[i][j-1][k] - vecFieldXC[i-1][j-1][k])*invdx + .25*(vecFieldXC[i][j][k-1] - vecFieldXC[i-1][j][k-1])*invdx + .25*(vecFieldXC[i][j-1][k-1] - vecFieldXC[i-1][j-1][k-1])*invdx;
          compY = .25*(vecFieldYC[i][j][k] - vecFieldYC[i][j-1][k])*invdy + .25*(vecFieldYC[i-1][j][k] - vecFieldYC[i-1][j-1][k])*invdy + .25*(vecFieldYC[i][j][k-1] - vecFieldYC[i][j-1][k-1])*invdy + .25*(vecFieldYC[i-1][j][k-1] - vecFieldYC[i-1][j-1][k-1])*invdy;
          compZ = .25*(vecFieldZC[i][j][k] - vecFieldZC[i][j][k-1])*invdz + .25*(vecFieldZC[i-1][j][k] - vecFieldZC[i-1][j][k-1])*invdz + .25*(vecFieldZC[i][j-1][k] - vecFieldZC[i][j-1][k-1])*invdz + .25*(vecFieldZC[i-1][j-1][k] - vecFieldZC[i-1][j-1][k-1])*invdz;
		  // calculate the divergence
		  divN[i][j][k] = compX + compY + compZ;
       }
   
}
/** calculate curl on nodes, given a vector field defined on central points  (nuovo)*/
inline void Grid3DCU::curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
    for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
	    for (register int k=1;j < nzn-1;k++){
          // curl - X
          compZDY = .25*(vecFieldZC[i][j][k] - vecFieldZC[i][j-1][k])*invdy +  .25*(vecFieldZC[i-1][j][k] - vecFieldZC[i-1][j-1][k])*invdy +  .25*(vecFieldZC[i][j][k-1] - vecFieldZC[i][j-1][k-1])*invdy +  .25*(vecFieldZC[i-1][j][k-1] - vecFieldZC[i-1][j-1][k-1])*invdy;
          compYDZ = .25*(vecFieldYC[i][j][k] - vecFieldYC[i][j][k-1])*invdz +  .25*(vecFieldYC[i-1][j][k] - vecFieldYC[i-1][j][k-1])*invdz +  .25*(vecFieldYC[i][j-1][k] - vecFieldYC[i][j-1][k-1])*invdz +  .25*(vecFieldYC[i-1][j-1][k] - vecFieldYC[i-1][j-1][k-1])*invdz;
		  // curl - Y
          compXDZ = .25*(vecFieldXC[i][j][k] - vecFieldXC[i][j][k-1])*invdz +  .25*(vecFieldXC[i-1][j][k] - vecFieldXC[i-1][j][k-1])*invdz +  .25*(vecFieldXC[i][j-1][k] - vecFieldXC[i][j-1][k-1])*invdz +  .25*(vecFieldXC[i-1][j-1][k] - vecFieldXC[i-1][j-1][k-1])*invdz;
		  compZDX = .25*(vecFieldZC[i][j][k] - vecFieldZC[i-1][j][k])*invdx +  .25*(vecFieldZC[i][j-1][k] - vecFieldZC[i-1][j-1][k])*invdx +  .25*(vecFieldZC[i][j][k-1] - vecFieldZC[i-1][j][k-1])*invdx +  .25*(vecFieldZC[i][j-1][k-1] - vecFieldZC[i-1][j-1][k-1])*invdx;
          // curl - Z
          compYDX = .25*(vecFieldYC[i][j][k] - vecFieldYC[i-1][j][k])*invdx +  .25*(vecFieldYC[i][j-1][k] - vecFieldYC[i-1][j-1][k])*invdx +  .25*(vecFieldYC[i][j][k-1] - vecFieldYC[i-1][j][k-1])*invdx +  .25*(vecFieldYC[i][j-1][k-1] - vecFieldYC[i-1][j-1][k-1])*invdx;
	      compXDY = .25*(vecFieldXC[i][j][k] - vecFieldXC[i][j-1][k])*invdy +  .25*(vecFieldXC[i-1][j][k] - vecFieldXC[i-1][j-1][k])*invdy +  .25*(vecFieldXC[i][j][k-1] - vecFieldXC[i][j-1][k-1])*invdy +  .25*(vecFieldXC[i-1][j][k-1] - vecFieldXC[i-1][j-1][k-1])*invdy;

          curlXN[i][j][k] = compZDY - compYDZ;
          curlYN[i][j][k] = compXDZ - compZDX;
          curlZN[i][j][k] = compYDX - compXDY;
		}
 

}
/** calculate curl on central points, given a vector field defined on nodes */
inline void Grid3DCU::curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN){
  double compZDY, compYDZ;
  double compXDZ, compZDX;
  double compYDX, compXDY;
    for (register int i=0;i < nxn-1;i++)
     for (register int j=0;j < nyn-1;j++)
	  for (register int k=0;k < nzn-1;k++){
          // curl - X
          compZDY = .25*(vecFieldZN[i][j+1][k] - vecFieldZN[i][j][k])*invdy +  .25*(vecFieldZN[i+1][j+1][k] - vecFieldZN[i+1][j][k])*invdy +  .25*(vecFieldZN[i][j+1][k+1] - vecFieldZN[i][j][k+1])*invdy +  .25*(vecFieldZN[i+1][j+1][k+1] - vecFieldZN[i+1][j][k+1])*invdy;
          compYDZ = .25*(vecFieldYN[i][j][k+1] - vecFieldYN[i][j][k])*invdz +  .25*(vecFieldYN[i+1][j][k+1] - vecFieldYN[i+1][j][k])*invdz + .25*(vecFieldYN[i][j+1][k+1] - vecFieldYN[i][j+1][k])*invdz  +  .25*(vecFieldYN[i+1][j+1][k+1] - vecFieldYN[i+1][j+1][k])*invdz;
		  // curl - Y
          compXDZ = .25*(vecFieldXN[i][j][k+1] - vecFieldXN[i][j][k])*invdz + .25*(vecFieldXN[i+1][j][k+1] - vecFieldXN[i+1][j][k])*invdz  + .25*(vecFieldXN[i][j+1][k+1] - vecFieldXN[i][j+1][k])*invdz  +  .25*(vecFieldXN[i+1][j+1][k+1] - vecFieldXN[i+1][j+1][k])*invdz;
		  compZDX = .25*(vecFieldZN[i+1][j][k] - vecFieldZN[i][j][k])*invdx + .25*(vecFieldZN[i+1][j+1][k] - vecFieldZN[i][j+1][k])*invdx  + .25*(vecFieldZN[i+1][j][k+1] - vecFieldZN[i][j][k+1])*invdx  +  .25*(vecFieldZN[i+1][j+1][k+1] - vecFieldZN[i][j+1][k+1])*invdx;
          // curl - Z
          compYDX = .25*(vecFieldYN[i+1][j][k] - vecFieldYN[i][j][k])*invdx +  .25*(vecFieldYN[i+1][j+1][k] - vecFieldYN[i][j+1][k])*invdx +  .25*(vecFieldYN[i+1][j][k+1] - vecFieldYN[i][j][k+1])*invdx +  .25*(vecFieldYN[i+1][j+1][k+1] - vecFieldYN[i][j+1][k+1])*invdx;
          compXDY = .25*(vecFieldXN[i][j+1][k] - vecFieldXN[i][j][k])*invdy +  .25*(vecFieldXN[i+1][j+1][k] - vecFieldXN[i+1][j][k])*invdy +  .25*(vecFieldXN[i][j+1][k+1] - vecFieldXN[i][j][k+1])*invdy +  .25*(vecFieldXN[i+1][j+1][k+1] - vecFieldXN[i+1][j][k+1])*invdy;


          curlXC[i][j][k] = compZDY - compYDZ;
          curlYC[i][j][k] = compXDZ - compZDX;
          curlZC[i][j][k] = compYDX - compXDY;
      }
 

}
/** calculate laplacian on nodes, given a scalar field defined on nodes */
inline void Grid3DCU::lapN2N(double*** lapN,double ***scFieldN,VirtualTopology *vct){
   // calculate laplacian as divercence of gradient
   // allocate 3 gradients: defined on central points
   double*** gradXC = newArr3(double,nxc,nyc,nzc);
   double*** gradYC = newArr3(double,nxc,nyc,nzc);
   double*** gradZC = newArr3(double,nxc,nyc,nzc);
   gradN2C(gradXC,gradYC,gradZC,scFieldN);
   // communicate centers ?
   
   divC2N(lapN,gradXC,gradYC,gradZC);
   communicateNode(nxn,nyn,nzn,lapN,vct);
    
   // deallocate
   delArr3(gradXC,nxc,nyc);
   delArr3(gradYC,nxc,nyc);
   delArr3(gradZC,nxc,nyc);
}
/** calculate laplacian on central points, given a scalar field defined on central points */
inline void Grid3DCU::lapC2C(double*** lapC,double ***scFieldC,VirtualTopology *vct){
    //  calculate laplacian as divercence of gradient
    // allocate 3 gradients: defined on nodes
    double*** gradXN = newArr3(double,nxn,nyn,nzn);
    double*** gradYN = newArr3(double,nxn,nyn,nzn);
   
    gradC2N(gradXN,gradYN,gradZN,scFieldC);
    
    communicateNode(nxn,nyn,nzn,gradXN,vct); 
    communicateNode(nxn,nyn,nzn,gradYN,vct);
	communicateNode(nxn,nyn,nzn,gradZN,vct);
    // calculate divergence of gradient
    divN2C(lapC,gradXN,gradYN);
    // deallocate
    delArr3(gradXN,nxn,nyn);
    delArr3(gradYN,nxn,nyn);
	delArr3(gradZN,nxn,nyn);

}                              

/** calculate divergence on  boundaries */
inline void Grid3DCU::divBCleft(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dirDER){
	double compX, compY, compZ;
    switch(dirDER){
    // remember this has to calculated on the center

    case 0:    //  DIVERGENCE DIRECTION  X
     for (int j=1; j <  nyn-2;j++)
	  for (int k=1; k <  nzn-2;k++){
         compX = .5*(vectorX[leftActiveNode+1][j][k] - vectorX[leftActiveNode][j][k])*invdx +  .5*(vectorX[leftActiveNode+1][j+1][k] - vectorX[leftActiveNode][j+1][k])*invdx; 
         compY = .5*(vectorY[leftActiveNode][j+1][0] - vectorY[leftActiveNode][j][0])*invdy +  .5*(vectorY[leftActiveNode+1][j+1][0] - vectorY[leftActiveNode+1][j][0])*invdy; // average
         divBC[leftActiveNode][j][0] = compX + compY;
     }
   compX = .5*(vectorX[leftActiveNode+1][nyn-2][0] - vectorX[leftActiveNode][nyn-2][0])*invdx +  .5*(vectorX[leftActiveNode+1][nyn-3][0] - vectorX[leftActiveNode][nyn-3][0])*invdx; // average
   compY = .5*(vectorY[leftActiveNode][nyn-2][0] - vectorY[leftActiveNode][nyn-3][0])*invdy +  .5*(vectorY[leftActiveNode+1][nyn-2][0] - vectorY[leftActiveNode+1][nyn-3][0])*invdy; // average
   divBC[leftActiveNode][nyn-2][0] = compX+ compY;
    
     break;
  

   case 1:    //  DIVERGENCE DIRECTION  Y
     for (int i=1; i < nxn-2;i++){
          compX = .5*(vectorX[i+1][leftActiveNode][0] - vectorX[i][leftActiveNode][0])*invdx +  .5*(vectorX[i+1][leftActiveNode+1][0] - vectorX[i][leftActiveNode+1][0])*invdx;
          compY = .5*(vectorY[i][leftActiveNode+1][0] - vectorY[i][leftActiveNode][0])*invdy +  .5*(vectorY[i+1][leftActiveNode+1][0] - vectorY[i+1][leftActiveNode][0])*invdy;
	  divBC[i][leftActiveNode][0] = compX + compY;
     }
    compX = .5*(vectorX[nxn-2][leftActiveNode][0] - vectorX[nxn-3][leftActiveNode][0])*invdx +  .5*(vectorX[nxn-2][leftActiveNode+1][0] - vectorX[nxn-3][leftActiveNode+1][0])*invdx;
    compY = .5*(vectorY[nxn-2][leftActiveNode+1][0] - vectorY[nxn-2][leftActiveNode][0])*invdy +  .5*(vectorY[nxn-3][leftActiveNode+1][0] - vectorY[nxn-3][leftActiveNode][0])*invdy;
   divBC[nxn-2][leftActiveNode][0] = compX + compY;
    
     break;
     }

}
/** calculate divergence on  boundaries to impose divergence constarint*/
inline void Grid3DCU::divBCright(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ,int rightActiveNode, int dirDER){
   double compX, compY, compZ;
    switch(dirDER){
     case 0:    //  DIVERGENCE DIRECTION  X
     for (int j=1; j <  nyn-2;j++){
          compX = .5*(vectorX[rightActiveNode][j][0] - vectorX[rightActiveNode-1][j][0])*invdx +  .5*(vectorX[rightActiveNode][j+1][0] - vectorX[rightActiveNode-1][j+1][0])*invdx;
          compY = .5*(vectorY[rightActiveNode][j+1][0] - vectorY[rightActiveNode][j][0])*invdy +  .5*(vectorY[rightActiveNode-1][j+1][0] - vectorY[rightActiveNode-1][j][0])*invdy;
          divBC[rightActiveNode][j][0] = compX + compY;
     }
  compX = .5*(vectorX[rightActiveNode][nyn-2][0] - vectorX[rightActiveNode-1][nyn-2][0])*invdx +  .5*(vectorX[rightActiveNode][nyn-3][0] - vectorX[rightActiveNode-1][nyn-3][0])*invdx;
 compY = .5*(vectorY[rightActiveNode][nyn-2][0] - vectorY[rightActiveNode][nyn-3][0])*invdy +  .5*(vectorY[rightActiveNode-1][nyn-2][0] - vectorY[rightActiveNode-1][nyn-3][0])*invdy;
  divBC[rightActiveNode][nyn-2][0] = compX + compY;
    
     break;


     case 1:    //  DIVERGENCE DIRECTION  Y
     
     for (int i=1; i < nxn-2;i++){
          compX = .5*(vectorX[i+1][rightActiveNode][0] - vectorX[i][rightActiveNode-1][0])*invdx +  .5*(vectorX[i+1][rightActiveNode][0] - vectorX[i][rightActiveNode-1][0])*invdx;
          compY = .5*(vectorY[i][rightActiveNode][0]   - vectorY[i][rightActiveNode-1][0])*invdy +  .5*(vectorY[i+1][rightActiveNode-1][0] - vectorY[i+1][rightActiveNode-1][0])*invdy;
	     divBC[i][rightActiveNode][0] =  compX + compY;
     }
   
   compX = .5*(vectorX[nxn-2][rightActiveNode][0] - vectorX[nxn-3][rightActiveNode-1][0])*invdx +  .5*(vectorX[nxn-2][rightActiveNode][0] - vectorX[nxn-3][rightActiveNode-1][0])*invdx;
   compY = .5*(vectorY[nxn-2][rightActiveNode][0]   - vectorY[nxn-2][rightActiveNode-1][0])*invdy +  .5*(vectorY[nxn-3][rightActiveNode-1][0] - vectorY[nxn-3][rightActiveNode-1][0])*invdy;
  divBC[nxn-2][rightActiveNode][0] =  compX + compY;
     
      break;
}
}
/** calculate derivative on boundary */
inline void Grid3DCU::derBC(double ***derBC, double ***vector, int leftActiveNode, int dirDER){
    switch(dirDER){
     case 0:    //  DERIVATIVE DIRECTION  X
     for (register int j=1;j < nyc-1;j++)
            derBC[leftActiveNode][j][0] = .5*(vector[leftActiveNode+1][j][0] - vector[leftActiveNode][j][0])*invdx +  .5*(vector[leftActiveNode+1][j+1][0] - vector[leftActiveNode][j+1][0])*invdx;;
      
     break;
     case 1:    //  DERIVATIVE DIRECTION  Y
     for (register int i=1;i < nxc-1;i++)
           derBC[i][leftActiveNode][0] = .5*(vector[i][leftActiveNode+1][0] - vector[i][leftActiveNode][0])*invdy +  .5*(vector[i+1][leftActiveNode+1][0] - vector[i+1][leftActiveNode][0])*invdy;
       
     break;
}
}
/** interpolate on nodes from central points: do this for the magnetic field (new)*/
inline void Grid3DCU::interpC2N(double ***vecFieldN, double ***vecFieldC){
 for (register int i=1;i < nxn-1;i++)
      for (register int j=1;j < nyn-1;j++)
	    for (register int k=1;k < nzn-1;k++)
         vecFieldN[i][j][k] = .125*(vecFieldC[i][j][k] + vecFieldC[i-1][j][k] + vecFieldC[i][j-1][k] + vecFieldC[i-1][j-1][k] + vecFieldC[i][j][k-1] + vecFieldC[i-1][j][k-1] + vecFieldC[i][j-1][k-1] + vecFieldC[i-1][j-1][k-1]);
}
/** interpolate on central points from nodes*/
inline void Grid3DCU::interpN2C(double ***vecFieldC, double ***vecFieldN){
  for (register int i=0;i < nxc;i++)
	   for (register int j=0;j < nyc;j++)
	      for (register int k=0;j < nzc;k++)
	       vecFieldC[i][j][k] = .125*(vecFieldN[i+1][j][k] + vecFieldN[i][j][k] + vecFieldN[i+1][j+1][k] + vecFieldN[i][j+1][k] + vecFieldN[i+1][j][k+1] + vecFieldN[i][j][k+1] + vecFieldN[i+1][j+1][k+1] + vecFieldN[i][j+1][k+1]);
	   

}
/** interpolate on central points from nodes per species (new)*/
inline void Grid3DCU::interpN2C(double ****vecFieldC, int ns, double ****vecFieldN){
  for (register int i=0;i < nxc;i++)
     for (register int j=0;j < nyc;j++)
	   for (register int k=0;k < nzc;k++)
         vecFieldC[ns][i][j][k] = .125*(vecFieldN[ns][i+1][j][k] + vecFieldN[ns][i][j][k] + vecFieldN[ns][i+1][j+1][k] + vecFieldN[ns][i][j+1][k] + vecFieldN[ns][i+1][j][k+1] + vecFieldN[ns][i][j][k+1] + vecFieldN[ns][i+1][j+1][k+1] + vecFieldN[ns][i][j+1][k+1]);
}


/** get nxc - number of cells in the X-DIRECTION*/
inline int Grid3DCU::getNXC(){
  return(nxc);
}
/** get nxn - number of nodes in the X-DIRECTION*/
inline int Grid3DCU::getNXN(){
  return(nxn);
}
/** get nyc - number of cells in the Y-DIRECTION*/
inline int Grid3DCU::getNYC(){
  return(nyc);
}
/** get nyc - number of cells in the Z-DIRECTION*/
inline int Grid3DCU::getNZC(){
  return(nzc);
}
/** get nyn - number of nodes in the Y-DIRECTION*/
inline int Grid3DCU::getNYN(){
  return(nyn);
}
/** get nzn - number of nodes in the Z-DIRECTION*/
inline int Grid3DCU::getNZN(){
   return(nzn);
}
/** get dx - grid spacing */
inline double Grid3DCU::getDX(){
  return(dx);
}
/** get dy - grid spacing*/
inline double Grid3DCU::getDY(){
  return(dy);
}
/** get dz - grid spacing*/
inline double Grid3DCU::getDZ(){
   return(dz);
}

/** get x-coordinate of the node(X,Y,Z) */
inline double &Grid3DCU::getXN(int indexX, int indexY, int indexZ){
  return(xn[indexX][indexY][indexZ]);
}
/** get y-coordinate of the node(X,Y,Z) */
inline double &Grid3DCU::getYN(int indexX, int indexY, int indexZ){
  return(yn[indexX][indexY][indexZ]);
}
/** get z-coordinate of the node(X,Y,Z) */
inline double &Grid3DCU::getZN(int indexX, int indexY, int indexZ){
   return(zn[indexX][indexY][indexZ]);
}
/** get x-coordinate of the center of cell(X,Y,Z) */
inline double &Grid3DCU::getXC(int indexX, int indexY, int indexZ){
  return(xc[indexX][indexY][indexZ]);
}
/** get y-coordinate of the center of cell(X,Y,Z) */
inline double &Grid3DCU::getYC(int indexX, int indexY, int indexZ){
  return(yc[indexX][indexY][indexZ]);
}
/** get z-coordinate of the center of cell(X,Y,Z) */
inline double &Grid3DCU::getZC(int indexX, int indexY, int indexZ){
   return(zc[indexX][indexY][indexZ]);
}

/** get the whole vector xc with coordinate of center cells*/
inline double*** Grid3DCU::getXC(){
  return(xc);
}
/** get the whole vector yc with coordinate of center cells*/
inline double*** Grid3DCU::getYC(){
  return(yc);
}
/** get the whole vector zc with coordinate of center cells*/
inline double*** Grid3DCU::getZC(){
  return(zc);
}
/** get x coordinate of the first node of Simulation box - X direction */
inline double Grid3DCU::getXstart(){
 return(xStart);
}
/** get x coordinate of the last node of Simulation box - X direction*/
inline double Grid3DCU::getXend(){
 return(xEnd);
}
/** get y coordinate of the first node of Simulation box - Y direction */
inline double Grid3DCU::getYstart(){
 return(yStart);
}
/** get y coordinate of the last node of Simulation box - Y direction */
inline double Grid3DCU::getYend(){
 return(yEnd);
}
/** get y coordinate of the first node of Simulation box - Z direction */
inline double Grid3DCU::getZstart(){
  return(zStart);
}
/** get y coordinate of the last node of Simulation box - Z direction */
inline double Grid3DCU::getZend(){
  return(zEnd);
}
/** get the inverse of volume */
inline double Grid3DCU::getInvVOL(){
 return(invVOL);
}

/** calculate divergence on nodes, given a vector field defined on central points  */
inline void Grid3DCU::divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC){
  cout << "2D divC2N called from a 3D grid!" << endl;
}
/** calculate gradient on nodes, given a scalar field defined on central points  */
inline void Grid3DCU::gradN2C(double ***gradXC, double ***gradYC, double ***scFieldN){
   cout << "GradN2C in 2D called. Call the grad in 3D! " << endl;
}

/** calculate gradient on nodes, given a scalar field defined on central points  (nuovo)*/
inline void Grid3DCU::gradC2N(double***  gradXN, double***  gradYN, double*** scFieldC){
    cout << "GradC2N in 2D called. Call the grad in 3D! " << endl;
}
/** calculate divergence on central points, given a Tensor field defined on nodes */
inline void Grid3DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ****pXX, double**** pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
  cout << "Call to divSymmTensorN2C in 2D from 3D geometry!" << endl;
}








#endif

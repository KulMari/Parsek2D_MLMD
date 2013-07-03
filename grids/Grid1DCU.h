 /*******************************************************************************************************
 Grid1DCU.h  -  uniform cartesian 1D local grid for each process, including che guard cells
                             -------------------
 developers: Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
 *******************************************************************************************************/

#ifndef Grid1DCU_H
#define Grid1DCU_H

#include <iostream>

#include "Grid.h"
#include "../communication/ComInterpNodes.h"
#include "../communication/ComNodes.h"
#include "../utility/Alloc.h"

using std::cout;
using std::endl;

/**
* 1D Uniform cartesian local(each processor has its own grid) grid 
*
* Note that GRID1DCU is implementing the abstract class Grid.h, therefore
* it must implements all the virtual methods in Grid.h
*
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Grid1DCU : public Grid {
  public:
      /** constructor */
      Grid1DCU(CollectiveIO *col, VirtualTopology *vct);
      /** destructor */
      ~Grid1DCU();
      /** allocate grid arrays for this domain */
      void allocate(CollectiveIO *ptC, VirtualTopology *ptVCT);
      /** deallocate grid arrays for this domain */
      void deallocate();
      /** print grid info */
      void print(VirtualTopology* ptVCT);
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradC2N(double*** gradXN,double*** gradYN,double***scFieldC );
      /** calculate gradient on nodes, given a scalar field defined on central points  */
      void gradN2C(double ***gradXC, double ***gradYC, double ***scFieldN);
      /** calculate divergence on central points, given a vector field defined on nodes  */
      void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN);
      /** calculate divergence on nodes, given a vector field defined on central points  */
      void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC);
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
     /** number of nodes - Y direction, including + 2 extra nodes for guard cells */
     int nyn;
     /** dx = grid spacing - X direction */
     double dx;
     /** invdx = 1/dx */
     double invdx;
     /** invol = inverse of volume*/
     double invVOL;
     /** node - X coordinate  */
     double ***xn;
     /** centre of cell - X coordinate */
     double ***xc;
     /** local grid boundaries coordinates  */
     double xStart, xEnd;
     
};
/** constructor */
inline Grid1DCU::Grid1DCU(CollectiveIO* col, VirtualTopology* vct){
   // add 2 for the guard cells
   nxc = (col->getNxc())/(vct->getXLEN()) + 2;
   nxn = nxc + 1;
   dx  = col->getLx()/col->getNxc();
   invVOL = 1.0/dx;
   invdx = 1.0/dx;
   // local grid dimensions and boundaries of active nodes
   xStart = vct->getCoordinates(0)*(col->getLx()/(double) vct->getXLEN());
   xEnd   = xStart + (col->getLx()/(double) vct->getXLEN());
   // arrays allocation: nodes ---> the first node has index 1, the last has index nxn-2!
   xn = newArr3(double,nxn,1,1);
   for (int i=0; i < nxn; i++)
     xn[i][0][0] = xStart + (i-1)*dx;
   // arrays allocation: cells ---> the first cell has index 1, the last has index ncn-2!
   xc = newArr3(double,nxc,1,1);
   for (int i=0; i < nxc; i++)
    xc[i][0][0] = .5*(xn[i][0][0] + xn[i+1][0][0]);
       
}

/** deallocate the local grid */
inline Grid1DCU::~Grid1DCU(){
   // deallocate nodes
   delArr3(xn,nxn,1);
   // centers cells
   delArr3(xc,nxc,1);
}
/** print the local grid info */
inline void Grid1DCU::print(VirtualTopology* ptVCT){
    cout << endl;
    cout <<  "Subgrid (" << ptVCT->getCoordinates(0) << "," << ptVCT->getCoordinates(1) << ")"<< endl;
    cout <<  "Number of cell: -X=" << nxc-2 << endl;
    cout <<  "Xin = " << xn[1][1][0] << "; Xfin = " << xn[nxn-2][1][0] << endl;
    cout << endl;
    
}
/** calculate gradient on nodes, given a scalar field defined on central points  (nuovo)*/
inline void Grid1DCU::gradC2N(double***  gradXN, double***  gradYN, double*** scFieldC){
     for (register int i=1;i < nxc;i++)
	    gradXN[i][0][0] = (scFieldC[i][0][0] - scFieldC[i-1][0][0])*invdx;
}
/** calculate gradient on nodes, given a scalar field defined on central points  */
inline void Grid1DCU::gradN2C(double ***gradXC, double ***gradYC, double ***scFieldN){
  for (register int i=0;i < nxn-1; i++)
          gradXC[i][0][0] = (scFieldN[i+1][0][0] - scFieldN[i][0][0])*invdx;
         
}
/** calculate laplacian on nodes, given a scalar field defined on nodes */
inline void Grid1DCU::lapN2N(double*** lapN,double ***scFieldN,VirtualTopology *vct){
         // be sure it has ghost cells
	 for (int i=1; i < nxn-1;i++) 
	   lapN[i][0][0] = (scFieldN[i-1][0][0] - 2*scFieldN[i][0][0] + scFieldN[i+1][0][0])*invdx*invdx;
	 // communicate
}
/** calculate laplacian on central points, given a scalar field defined on central points */
inline void Grid1DCU::lapC2C(double*** lapC,double ***scFieldC,VirtualTopology *vct){
         // be sure it has ghost cells
	 for (int i=1; i < nxc-1;i++)
	   lapC[i][0][0] = (scFieldC[i-1][0][0] - 2*scFieldC[i][0][0] + scFieldC[i+1][0][0])*invdx*invdx;
	 // communicate
	 
}                              
/** interpolate on nodes from central points: do this for the magnetic field (new)*/
inline void Grid1DCU::interpC2N(double ***vecFieldN, double ***vecFieldC){
 // be sure that you comunicate here
 for (register int i=1;i < nxc;i++)
	  vecFieldN[i][0][0] = .5*(vecFieldC[i-1][0][0] + vecFieldC[i][0][0]); 
 // be sure you communicate after

}
/** interpolate on central points from nodes (new)*/
inline void Grid1DCU::interpN2C(double ***vecFieldC, double ***vecFieldN){
  // be sure that you comunicate here
   for (register int i=0;i < nxn-1; i++)
	  vecFieldC[i][0][0] = .5*(vecFieldN[i+1][0][0] + vecFieldN[i][0][0]);
  // be sure you communicate after

}
/** interpolate on central points from nodes per species (new)*/
inline void Grid1DCU::interpN2C(double ****vecFieldC, int ns, double ****vecFieldN){
    for (register int i=0;i < nxc-1;i++)
       vecFieldC[ns][i][0][0] = .5*(vecFieldN[ns][i+1][0][0] + vecFieldN[ns][i][0][0]);
}
/** get nxc - number of cells in the X-DIRECTION*/
inline int Grid1DCU::getNXC(){
  return(nxc);
}
/** get nxn - number of nodes in the X-DIRECTION*/
inline int Grid1DCU::getNXN(){
  return(nxn);
}
/** get dx - grid spacing */
inline double Grid1DCU::getDX(){
  return(dx);
}
/** get x-coordinate of the node(X,Y,Z) */
inline double &Grid1DCU::getXN(int indexX, int indexY, int indexZ){
  return(xn[indexX][0][0]);
}

/** get x-coordinate of the center of cell(X,Y,Z) */
inline double &Grid1DCU::getXC(int indexX, int indexY, int indexZ){
  return(xc[indexX][0][0]);
}
/** get the whole vector xc with coordinate of center cells*/
inline double*** Grid1DCU::getXC(){
  return(xc);
}

/** get x coordinate of the first node of Simulation box - X direction */
inline double Grid1DCU::getXstart(){
 return(xStart);
}
/** get x coordinate of the last node of Simulation box - X direction*/
inline double Grid1DCU::getXend(){
 return(xEnd);
}

/** get the inverse of volume */
inline double Grid1DCU::getInvVOL(){
 return(invVOL);
}



///////  ***** METHODS NOT IMPLEMENTED ********* ////////////
/** calculate divergence on central points, given a vector field defined on nodes  */
inline void Grid1DCU::divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN){
      cout << "method not implemented in 1D" << endl;
}
/** calculate divergence on central points, given a Tensor field defined on nodes */
// THIS NEED TO CHECKED!!!!!!!!!!!!!!!
inline void Grid1DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ****pXX, double**** pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
  cout << "method not implemented in 1D" << endl;

}
/** calculate divergence on central points, given a Tensor field defined on nodes */
inline void Grid1DCU::divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double**** pXY, double**** pXZ, double**** pYY, double**** pYZ, double**** pZZ,int ns){
   cout << "method not implemented in 1D" << endl;
}

/** calculate divergence on nodes, given a vector field defined on central points  */
inline void Grid1DCU::divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC){
   cout << "method not implemented in 1D" << endl;
   
}
/** calculate curl on nodes, given a vector field defined on central points  (nuovo)*/
inline void Grid1DCU::curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC){
   cout << "method not implemented in 1D" << endl;

}
/** calculate curl on central points, given a vector field defined on nodes  (nuovo)*/
inline void Grid1DCU::curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN){
    cout << "method not implemented in 1D" << endl;
}
/** calculate divergence on  boundaries */
inline void Grid1DCU::divBCleft(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ, int leftActiveNode, int dirDER){
     cout << "method not implemented in 1D" << endl;

}
/** calculate divergence on  boundaries to impose divergence constarint*/
inline void Grid1DCU::divBCright(double ***divBC,double ***vectorX, double ***vectorY, double ***vectorZ,int rightActiveNode, int dirDER){
     cout << "method not implemented in 1D" << endl;
}
/** calculate derivative on boundary */
inline void Grid1DCU::derBC(double ***derBC, double ***vector, int leftActiveNode, int dirDER){
     cout << "method not implemented in 1D" << endl;
}

/** get nyn - number of nodes in the Y-DIRECTION*/
inline int Grid1DCU::getNYN(){
  //cout << "No need to call getNYN(). You are using 1D grid X" << endl;
  return(3);
}
/** get nyc - number of cells in the Y-DIRECTION*/
inline int Grid1DCU::getNYC(){
  // cout << "No need to call getNYC(). You are using 1D grid X" << endl;
  return(3);
}
/** get dy - grid spacing*/
inline double Grid1DCU::getDY(){
  //cout << "No need to call getDY(). You are using 1D grid X" << endl;
  return(0.0);
}
/** get y-coordinate of the node(X,Y,Z) */
inline double &Grid1DCU::getYN(int indexX, int indexY, int indexZ){
  //cout << "No need to call getYN(). You are using 1D grid X" << endl;
  return(xc[indexX][0][0]);
}
/** get y-coordinate of the center of cell(X,Y,Z) */
inline double &Grid1DCU::getYC(int indexX, int indexY, int indexZ){
  //cout << "No need to call getYC(). You are using 1D grid X" << endl;
  return(xc[indexX][0][0]);
}
/** get the whole vector yc with coordinate of center cells*/
inline double*** Grid1DCU::getYC(){
  //cout << "No need to call getYC(). You are using 1D grid X" << endl;
    return(NULL);
}
/** get y coordinate of the first node of Simulation box - Y direction */
inline double Grid1DCU::getYstart(){
 cout << "No need to call getYstart(). You are using 1D grid X" << endl;
 return(0.0);
}
/** get y coordinate of the last node of Simulation box - Y direction */
inline double Grid1DCU::getYend(){
  cout << "No need to call getYend(). You are using 1D grid X" << endl;
 return(0.0);
}
/** get nyc - number of cells in the Z-DIRECTION*/
inline int Grid1DCU::getNZC(){
  cout << "No need to call getNZC(). You are using 1D grid X" << endl;
  return(1);
}
/** get nzn - number of nodes in the Z-DIRECTION*/
inline int Grid1DCU::getNZN(){
   cout << "No need to call getNZN(). You are using 1D grid X" << endl;
   return(1);
}
/** get dz - grid spacing*/
inline double Grid1DCU::getDZ(){
   cout << "No need to call getDZ(). You are using 1D grid X" << endl;
   return(1.0);
}
/** get z-coordinate of the center of cell(X,Y,Z) */
inline double &Grid1DCU::getZC(int indexX, int indexY, int indexZ){
   cout << "No need to call getZC(int indexX, int indexY, int indexZ). You are using 1D grid X" << endl;
   return(xc[0][0][0]);
}
/** get z-coordinate of the node(X,Y,Z) */
inline double &Grid1DCU::getZN(int indexX, int indexY, int indexZ){
    cout << "No need to call getZN(int indexX, int indexY, int indexZ). You are using 1D grid X" << endl;
   return(xn[0][0][0]);
}
/** get the whole vector zc with coordinate of center cells*/
inline double*** Grid1DCU::getZC(){
  cout << "No need to call getZC(). You are using 1D grid X" << endl;
  return(NULL);
}

/** get y coordinate of the first node of Simulation box - Z direction */
inline double Grid1DCU::getZstart(){
  cout << "No need to call getZstart(). You are using 1D grid X" << endl;
  return(0.0);
}
/** get y coordinate of the last node of Simulation box - Z direction */
inline double Grid1DCU::getZend(){
  cout << "No need to call getZend(). You are using 1D grid X" << endl;
  return(0.0);
}
#endif

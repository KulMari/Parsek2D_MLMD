/***************************************************************************
Grid.h  -  Abstract Base class for different GRIDS
                         
developers: Stefano Markidis, Giovanni Lapenta
   
 ***************************************************************************/

#ifndef Grid_H
#define Grid_H
/**
*  Abstract base class for different kinds of grid. Each subclass
*  must implement these methods
*
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta
* @version 2.0
*/
class Grid {
  public:

  // put them in Grid.h, to make them publicly accessible                                                                   
  /** node - X coordinate (indexX, indexY)   */
  double ***xn;                                                                                                           
  /** node - Y coordinate (indexX, indexY)   */
  double ***yn;              

    /** print grid info */
    virtual void print(VirtualTopology* ptVCT)=0;
    /** calculate gradient on nodes, given a scalar field defined on central points  */
    virtual   void gradC2N(double ***gradXN, double ***gradYN, double ***scFieldC,VirtualTopology *vct)=0;
    /** calculate gradient on nodes, given a scalar field defined on central points with one sided derivative on the BC  */
    virtual   void gradC2N_onesided_derBC(double ***gradXN, double ***gradYN, double ***scFieldC,VirtualTopology *vct)=0;
    /** calculate gradient on nodes, given a scalar field defined on central points  */
    virtual  void gradN2C(double ***gradXC, double ***gradYC, double ***scFieldN)=0;
    /** calculate gradient on nodes, given a scalar field defined on central points, includes ghost centers calculated using ghost nodes  */
    virtual  void gradN2C_plusghost(double ***gradXC, double ***gradYC, double ***scFieldN)=0;
    /** calculate divergence on central points, given a vector field defined on nodes  */
    virtual void divN2C(double ***divC, double ***vecFieldXN, double ***vecFieldYN)=0;
    /** calculate divergence on central points, given a vector field defined on nodes, includes ghost centers calculated using ghost nodes  */
    virtual void divN2C_plusghost(double ***divC, double ***vecFieldXN, double ***vecFieldYN)=0;
    /** calculate divergence on nodes, given a vector field defined on central points  */
    virtual  void divC2N(double ***divN, double ***vecFieldXC, double ***vecFieldYC,VirtualTopology *vct)=0;
	/** calculate divergence on nodes, given a vector field defined on central points with one sided derivative on the BC */
    virtual  void divC2N_onesided_derBC(double ***divN, double ***vecFieldXC, double ***vecFieldYC,VirtualTopology *vct)=0;
    /** calculate curl on nodes, given a vector field defined on central points  */
    virtual  void curlC2N(double ***curlXN, double ***curlYN, double ***curlZN,double ***vecFieldXC, double ***vecFieldYC, double*** vecFieldZC,VirtualTopology *vct)=0;
    /** calculate curl on central points, given a vector field defined on nodes  */
    virtual  void curlN2C(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN)=0;
    /** calculate curl on central points, given a vector field defined on nodes including ghost cells  */
    virtual  void curlN2C_withghost(double ***curlXC, double ***curlYC, double ***curlZC,double ***vecFieldXN, double ***vecFieldYN, double*** vecFieldZN)=0;
    /** calculate laplacian on nodes, given a scalar field defined on nodes */
    virtual  void lapN2N(double ***lapN,double ***scFieldN,VirtualTopology *vct)=0;
    /** calculate laplacian on active nodes, given a scalar field defined on nodes including ghost nodes */
    virtual  void lapN2N_plusghost(double ***lapN,double ***scFieldN,VirtualTopology *vct)=0;
    /** calculate laplacian on central points, given a scalar field defined on central points */
    virtual  void lapC2C(double*** lapC,double ***scFieldC,VirtualTopology *vct)=0;
    /** calculate laplacian on central points, given a scalar field defined on central points; patch for ghost centers */
    virtual  void lapC2C_plusghost(double*** lapC,double ***scFieldC,VirtualTopology *vct)=0;
    /** calculate divergence on central points, given a Tensor field defined on nodes  */
    virtual  void divSymmTensorN2C(double ***divCX, double ***divCY, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ, int ns)=0;
    virtual  void divSymmTensorN2C(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ, int ns)=0;

    virtual void divSymmTensorN2C_alsoGC(double ***divCX, double ***divCY, double ***divCZ, double ****pXX, double ****pXY, double ****pXZ, double ****pYY, double ****pYZ, double ****pZZ,int ns)=0;

    /** interpolate on nodes from central points */
    virtual   void interpC2N(double ***vecFieldN, double ***vecFieldC,VirtualTopology *vct)=0;
	/** interpolate on nodes from central points */
    virtual   void interpC2N_BC(double ***vecFieldN, double ***vecFieldC,VirtualTopology *vct)=0;
    // interpolate on nodes from central points, also ghost nodes
    virtual void interpC2N_BC_alsoGN(double ***vecFieldN, double ***vecFieldC, VirtualTopology *vct)=0;
    /** interpolate on central points from nodes */
    virtual  void interpN2C(double ***vecFieldC, double ***vecFieldN)=0;
    /** interpolate on central points from nodes */
    virtual  void interpN2C(double ****vecFieldC, int ns, double ****vecFieldN)=0;
    /** interpolate on central points from nodes including ghost cells*/
    virtual  void interpN2C_alsoGC(double ****vecFieldC, int ns, double ****vecFieldN)=0;
    virtual  void interpN2C_alsoGC(double ***vecFieldC, double ***vecFieldN)=0;

    /** return grid level */
    virtual int getLevel()=0;
    /** return ratio */
    virtual double getRatio()=0;
    /** return origin in x **/
    virtual double &getOx(int level)=0;
    /** return origin in y **/
    virtual double &getOy(int level)=0;
    /** return number of cells in the X-DIRECTION */
    virtual int getNXC()=0;
    /** return number of nodes in the X-DIRECTION */
    virtual int getNXN()=0;
    /** return number of cells in the Y-DIRECTION*/
    virtual int getNYC()=0;
    /** return number of nodes in the Y-DIRECTION */
    virtual int getNYN()=0;
    /** return number of cells in the Z-DIRECTION */
    virtual int getNZC()=0;
    /** return number of nodes in the Z-DIRECTION */
    virtual int getNZN()=0;
    /** return the grid spacing dx*/
    virtual double getDX()=0;
    /** return the grid spacing dy*/
    virtual double getDY()=0;
    /** return the grid spacing dz*/
    virtual double getDZ()=0;
    /** get x-coordinate of the node(X,Y,Z) */
    virtual double &getXN(int indexX, int indexY,int indexZ)=0;
    /** get y-coordinate of the node(X,Y,Z) */
    virtual double &getYN(int indexX, int indexY,int indexZ)=0;
    /** get z-coordinate of the node(X,Y,Z) */
    virtual double getModifiedXN(int indexX, int indexY,int indexZ)=0;
    /** get y-coordinate of the node(X,Y,Z) */
    virtual double getModifiedYN(int indexX, int indexY,int indexZ)=0;
    /** get z-coordinate of the node(X,Y,Z) */
    virtual double &getZN(int indexX, int indexY,int indexZ)=0;
    /** get x-coordinate of the center of cell(X,Y,Z) */
    virtual double &getXC(int indexX, int indexY,int indexZ)=0;
    /** get y-coordinate of the center of cell(X,Y,Z) */
    virtual double &getYC(int indexX, int indexY,int indexZ)=0;
    /** get z-coordinate of the center of cell(X,Y,Z) */
    virtual double &getZC(int indexX, int indexY,int indexZ)=0;
    /** get the whole vector xc with coordinate of center cells*/
    virtual double*** getXC()=0;
    /** get the whole vector yc with coordinate of center cells*/
    virtual double*** getYC()=0;
    /** get the whole vector zc with coordinate of center cells*/
    virtual double*** getZC()=0;
    /** get x coordinate of the first node of Simulation box - X direction */
    virtual double getXstart()=0;
    /** get x coordinate of the first node of Simulation box - X direction or of node 0 if this proc is on the x=0 border */
    virtual double getmodifiedXstart(VirtualTopology *vct)=0;
    /** get x coordinate of the last node of Simulation box - X direction*/
    virtual double getXend()=0;
    /** get x coordinate of the last node of Simulation box - X direction or of the VERY last node if this proc is on the x=Lx border*/
    virtual double getmodifiedXend(VirtualTopology *vct)=0;
    /** get y coordinate of the first node of Simulation box - Y direction */
    virtual double getYstart()=0;
    /** get y coordinate of the first node of Simulation box - Y direction or of node 0 if this proc is on the y=0 border */
    virtual double getmodifiedYstart(VirtualTopology *vct)=0;
    /** get y coordinate of the last node of Simulation box - Y direction */
    virtual double getYend()=0;
    /** get y coordinate of the first node of Simulation box - Y direction or of the VERY last node if this proc is on the y=Ly border */
    virtual double getmodifiedYend(VirtualTopology *vct)=0;
    /** get z coordinate of the first node of Simulation box - Z direction */
    virtual double getZstart()=0;
     /** get z coordinate of the last node of Simulation box - Z direction */
    virtual double getZend()=0;
    /** get the inverse of Volume */
    virtual double getInvVOL()=0;
};
#endif

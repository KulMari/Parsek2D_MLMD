/*******************************************************************************************
Field.h  -  Abstract class for fields
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef Field_H
#define Field_H
/**
* 
* Abstract class for fields
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Field {
public:
     /** print field information */
     virtual void print(void)const=0;
     
     
     /** get Potential array */
     virtual double*** getPHI() =0;
     /** get Electric field X component array */
     virtual double*** getEx() =0;
     /** get Electric field Y component array */
     virtual double*** getEy() =0;
     /** get Electric field Z component array */
     virtual double*** getEz() =0;
     /** get Magnetic field X component array */
     virtual double*** getBx() =0;
     /** get Magnetic field Y component array */
     virtual double*** getBy() =0;
     /** get Magnetic field Z component array */
     virtual double*** getBz() =0;
     /** get Electric Field component X defined on node(indexX,indexY,indexZ) */
     virtual double &getEx(int indexX, int indexY, int indexZ)const=0;
     /** get Electric Field component Y defined on node(indexX,indexY,indexZ) */
     virtual double &getEy(int indexX, int indexY, int indexZ)const=0;
     /** get Electric Field component Z defined on node(indexX,indexY,indexZ) */
     virtual double &getEz(int indexX, int indexY, int indexZ)const=0;
     /** get Magnetic Field component X defined on node(indexX,indexY,indexZ) */
     virtual double &getBx(int indexX, int indexY, int indexZ)const=0;
     /** get Magnetic Field component Y defined on node(indexX,indexY,indexZ) */
     virtual double &getBy(int indexX, int indexY, int indexZ)const=0;
     /** get Magnetic Field component Z defined on node(indexX,indexY,indexZ) */
     virtual double &getBz(int indexX, int indexY, int indexZ)const=0;
     /** get density on cell(indexX,indexY,indexZ)  */
     virtual double &getRHOc(int indexX, int indexY, int indexZ)const=0;
     /** get density on node(indexX,indexY,indexZ)  */
     virtual double &getRHOn(int indexX, int indexY, int indexZ)const=0;
     /** SPECIES: get density defined on center of cells*/
     virtual double &getRHOcs(int indexX, int indexY, int indexZ, int ns)const=0;
     /** SPECIES: get density defined on nodes */
     virtual double &getRHOns(int indexX, int indexY, int indexZ, int ns)const=0;
     /** get current -Direction X */
     virtual double &getJx(int indexX, int indexY, int indexZ) const=0;
     /** get current -Direction Y  */
     virtual double &getJy(int indexX, int indexY, int indexZ) const=0;
     /** get current -Direction Z */
     virtual double &getJz(int indexX, int indexY, int indexZ) const=0;
     /** get current -Direction X */
     virtual double &getJxs(int indexX, int indexY, int indexZ, int is) const=0;
     /** get current -Direction Y  */
     virtual double &getJys(int indexX, int indexY, int indexZ, int is) const=0;
     /** get current -Direction Z */
     virtual double &getJzs(int indexX, int indexY, int indexZ, int is) const=0;
     
     /** get pressure -Direction XX */
     virtual double &getpXXns(int indexX, int indexY, int indexZ, int is) const=0;
     /** get pressure -Direction XY  */
     virtual double &getpXYns(int indexX, int indexY, int indexZ, int is) const=0;
     /** get pressure -Direction XZ */
     virtual double &getpXZns(int indexX, int indexY, int indexZ, int is) const=0;
     /** get pressure -Direction YY */
     virtual double &getpYYns(int indexX, int indexY, int indexZ, int is) const=0;
     /** get pressure -Direction YY  */
     virtual double &getpYZns(int indexX, int indexY, int indexZ, int is) const=0;
     /** get pressure -Direction ZZ */
     virtual double &getpZZns(int indexX, int indexY, int indexZ, int is) const=0;

     /** get density array defined on centers cells */
     virtual double*** getRHOc() =0;
     /** get density array defined on nodes*/
     virtual double*** getRHOn() =0;
     /** SPECIES: get density array defined on nodes  */
     virtual double**** getRHOns() =0;
     /** get current array X component */
     virtual double*** getJx() =0;
     /** get current array Y component */
     virtual double*** getJy() =0;
     /** get current array Z component */
     virtual double*** getJz() =0;
     /** SPECIES: get current array X component */
     virtual double**** getJxs() =0;
     /** SPECIES: get current array Y component */
     virtual double**** getJys() =0;     
     /** SPECIES: get current array X component */
     virtual double**** getJzs() =0;
     /** SPECIES: get pressure tensor component XX defined on nodes */
     virtual double**** getpXXsn() =0;
     /** SPECIES: get pressure tensor component XY defined on nodes */
     virtual double**** getpXYsn() =0;
     /** SPECIES: get pressure tensor component XZ defined on nodes */
     virtual double**** getpXZsn() =0;
     /** SPECIES: get pressure tensor component YY defined on nodes */
     virtual double**** getpYYsn() =0;
     /** SPECIES: get pressure tensor component YZ defined on nodes */
     virtual double**** getpYZsn() =0;
     /** SPECIES: get pressure tensor component ZZ defined on nodes */
     virtual double**** getpZZsn() =0;
     
     // virtual sets
     virtual void setRHOns(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setJxs(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setJys(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setJzs(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setpXXns(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setpXYns(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setpXZns(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setpYYns(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setpYZns(double value, int indexX, int indexY,int indexZ,int is)=0;
     virtual void setpZZns(double value, int indexX, int indexY,int indexZ,int is)=0;
    
    

     ////////////////////////// INTERPOLATION ///////////////////////////////
     /** set to 0 all the densities fields */
     virtual void setZeroDensities()=0;
     
     /** add an amount of charge density to charge density field at node X,Y,Z */
     virtual void addRho(double ***weight, int X, int Y, int Z, int ns)=0;
     /** add an amount of current density - direction X to current density field at node X,Y,Z */
     virtual void addJx(double ***weight, int X, int Y, int Z, int ns)=0;
     /** add an amount of current density - direction Y to current density field at node X,Y,Z */
     virtual void addJy(double ***weight, int X, int Y,int Z, int ns)=0;
     /** add an amount of current density - direction Z to current density field at node X,Y,Z */
     virtual void addJz(double ***weight, int X, int Y,int Z, int ns)=0;
     /** add an amount of pressure density - direction XX to current density field at node X,Y,Z */
     virtual void addPxx(double ***weight, int X, int Y,int Z, int ns)=0;
     /** add an amount of pressure density - direction XY to current density field at node X,Y,Z */
     virtual void addPxy(double ***weight, int X, int Y, int Z,int ns)=0;
     /** add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
     virtual void addPxz(double ***weight, int X, int Y, int Z,int ns)=0;
     /** add an amount of pressure density - direction YY to current density field at node X,Y,Z */
     virtual void addPyy(double ***weight, int X, int Y, int Z,int ns)=0;
     /** add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
     virtual void addPyz(double ***weight, int X, int Y, int Z,int ns)=0;
     /** add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
     virtual void addPzz(double ***weight, int X, int Y, int Z,int ns)=0;

     /** add an amount of charge density to charge density field at node X,Y,Z */
     virtual void addRho_OS(double ***weight, int X, int Y, int Z, int ns, Grid *g)=0;
     /** add an amount of current density - direction X to current density field at node X,Y,Z */
     virtual void addJx_OS(double ***weight, int X, int Y, int Z, int ns, Grid *g)=0;
     /** add an amount of current density - direction Y to current density field at node X,Y,Z */
     virtual void addJy_OS(double ***weight, int X, int Y,int Z, int ns, Grid *g)=0;
     /** add an amount of current density - direction Z to current density field at node X,Y,Z */
     virtual void addJz_OS(double ***weight, int X, int Y,int Z, int ns, Grid *g)=0;
     /** add an amount of pressure density - direction XX to current density field at node X,Y,Z */
     virtual void addPxx_OS(double ***weight, int X, int Y,int Z, int ns, Grid *g)=0;
     /** add an amount of pressure density - direction XY to current density field at node X,Y,Z */
     virtual void addPxy_OS(double ***weight, int X, int Y, int Z,int ns, Grid *g)=0;
     /** add an amount of pressure density - direction XZ to current density field at node X,Y,Z */
     virtual void addPxz_OS(double ***weight, int X, int Y, int Z,int ns, Grid *g)=0;
     /** add an amount of pressure density - direction YY to current density field at node X,Y,Z */
     virtual void addPyy_OS(double ***weight, int X, int Y, int Z,int ns, Grid *g)=0;
     /** add an amount of pressure density - direction YZ to current density field at node X,Y,Z */
     virtual void addPyz_OS(double ***weight, int X, int Y, int Z,int ns, Grid *g)=0;
     /** add an amount of pressure density - direction ZZ to current density field at node X,Y,Z */
     virtual void addPzz_OS(double ***weight, int X, int Y, int Z,int ns, Grid *g)=0;
     
     /** communicate ghost for grid -> Particles interpolation */
     virtual void communicateGhostP2G(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct)=0;
     /** communicate ghost for grid -> Particles interpolation for OS */
     virtual void communicateGhostP2GOS(int ns, int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft, VirtualTopology *vct)=0;
     /** Sum density over different species */
     virtual void sumOverSpecies(VirtualTopology *vct)=0;
     /** Sum current over different species */
     virtual void sumOverSpeciesJ()=0;
     /** communicate ghost for densities and interp rho from node to center */
     virtual void interpDensitiesN2C(VirtualTopology *vct, Grid *grid)=0;
     

     //////////////////   FIELD SOLUTION  /////////////////////////////////////////
     /** Image of Maxwell solver */
     virtual void MaxwellImage(double *im, double *vector, Grid *grid,VirtualTopology *vct)=0;     
     /** maxwell Source */
     virtual void MaxwellSource(double *bkrylov, Grid *grid, VirtualTopology *vct)=0;
     /** Image of Poisson Solver */
     virtual void PoissonImage(double *image, double *vector, Grid *grid, VirtualTopology *vct)=0;

     // AMR, ME
     virtual int getNmessageBC()=0;
     virtual int *getTargetBC()=0;
     virtual int *getBCSide()=0;
     virtual int getTargetBOTTOM()=0;
     virtual int getTargetTOP()=0;
     virtual int getTargetLEFT()=0;
     virtual int getTargetRIGHT()=0;
     
     virtual int getNmessagerecuBC()=0;

     virtual void printRhons(int ns, VirtualTopology *vct)=0;
     virtual void printRhoc(VirtualTopology *vct)=0;
     virtual void setRhoc(double val)=0;  // for debugging purposes                                    
     virtual void printVecN(double ***vec, VirtualTopology *vct)=0;
     virtual void printVecC(double ***vec, VirtualTopology *vct)=0;
     virtual void printVecN(double ****vec, VirtualTopology *vct, int is)=0;
     virtual void printVecC(double ****vec, VirtualTopology *vct, int is)=0;
     virtual void printJxs(int ns, VirtualTopology *vct)=0;
     virtual void printPxxsn(int ns, VirtualTopology *vct)=0;
     // end, AMR ME
};
#endif

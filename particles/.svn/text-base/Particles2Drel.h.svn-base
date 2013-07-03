/*******************************************************************************************
Particles2Drel.h  -  Class for particles of the same species, for relativstic case and ES field, in a 2D space and 3component velocity
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef Part2Drel_H
#define Part2Drel_H


#include "Particles2Dcomm.h"
/**
* 
* Class for particles of the same species, for relativstic case and ES field, in a 2D space and 3component velocity
* @date Fri Jun 4 2007
* @author Stefano Markidis, Enrico Camporeale, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Particles2Drel : public Particles2Dcomm {

  public:
   /** default constructor */
   Particles2Drel();
   /** overload the interpolation because we use just ES field */
   void interpP2G(Field* ESf, Grid *grid, VirtualTopology* vct);
   /** overload the mover with an explicit scheme we use relaticistic EOM*/
   void mover_explicit(Grid* grid,VirtualTopology* vct, Field* EMf);
   /** initialize particle distribution for the relativistic two stream instability */
   void two_beams(Grid* grid,Field* EMf,VirtualTopology* vct);

};



#endif






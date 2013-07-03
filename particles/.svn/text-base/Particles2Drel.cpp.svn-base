/*******************************************************************************************
Particles2Drel.cpp  -  Class for particles of the same species, for relativstic case and ES field, in a 2D space and 3component velocity
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/


#include "../processtopology/VirtualTopology.h"
#include "../processtopology/VCtopology.h"
#include "../inputoutput/CollectiveIO.h"
#include "../inputoutput/Collective.h"
#include "../grids/Grid.h"
#include "../grids/Grid2DCU.h"
#include "../fields/Field.h"
#include "../mathlib/Basic.h"





#include "Particles2Drel.h"

/** the relativistic particles inherits the constructor from the base class */
Particles2Drel::Particles2Drel():Particles2Dcomm(){}

/** overload the mover with an explicit relativistic mover scheme in the ES limit*/
void Particles2Drel::mover_explicit(Grid* grid,VirtualTopology* vct, Field* EMf){
  
 	int innter,temp1,temp2,ix,iy;
  	double  qomdt = qom*dt,weight11, weight00, weight10, weight01,Exl, Eyl, Ezl;
  	// relativistic factor
        double gamma;
        // momenta
        double px, py;
        for (register int i=0; i < nop; i++){
	 	ix = 2 +  int(floor((x[i]-grid->getXstart())/grid->getDX()));
        	iy = 2 +  int(floor((y[i]-grid->getYstart())/grid->getDY()));
        	temp1 = (int) min(ix, nxn-2);
        	temp2 = (int) min(iy, nyn-2);
        	ix = (int) max(temp1,2);
        	iy = (int) max(temp2,2);
        	weight11 = ((x[i] - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((y[i] - grid->getYN(ix-1,iy-1,0))/grid->getDY());
        	weight10 = ((x[i] - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - y[i])/grid->getDY());
        	weight01 = ((grid->getXN(ix,iy-1,0) - x[i])/grid->getDX())*((y[i] - grid->getYN(ix,iy-1,0))/grid->getDY());
        	weight00 = ((grid->getXN(ix,iy,0) - x[i])/grid->getDX())*((grid->getYN(ix,iy,0) - y[i])/grid->getDY());
        	Exl = weight00*EMf->getEx(ix-1,iy-1,0) + weight01*EMf->getEx(ix-1,iy,0) + weight10*EMf->getEx(ix,iy-1,0) +
		weight11*EMf->getEx(ix,iy,0);
        	Eyl = weight00*EMf->getEy(ix-1,iy-1,0) + weight01*EMf->getEy(ix-1,iy,0) + weight10*EMf->getEy(ix,iy-1,0) +
		weight11*EMf->getEy(ix,iy,0);
        	// calculate gamma for this particle
        	gamma = 1/(sqrt(1 - (u[i]*u[i] + v[i]*v[i])/(c*c)));
		px = gamma*u[i];
                py = gamma*v[i];
                // calculate the new momenta using equation of motion
                px += qomdt*Exl;
                py += qomdt*Eyl;
                // calculate the new gamma
                gamma = sqrt(1 + px*px + py*py);
                // update velocity and position
                u[i] = px/gamma;  
                v[i] = py/gamma;
                x[i] += u[i]*dt;
		y[i] += v[i]*dt;
		
		
	}
	communicate(vct);
        MPI_Barrier(MPI_COMM_WORLD);
	
 
  	
}
/** overload the interpolation, using just the electric field (ES limit) */
void Particles2Drel::interpP2G(Field* ESf, Grid *grid, VirtualTopology* vct){
     double*** weight  = newArr3(double,2,2,1);
     double*** temp    = newArr3(double,2,2,1);
     int ix,iy, temp2,temp1;
     for (register int i=0; i < nop; i++){
         ix = 2 +  int(floor((x[i]-grid->getXstart())/grid->getDX()))  ;
         iy = 2 +  int(floor((y[i]-grid->getYstart())/grid->getDY()))  ;
         temp1 = (int) min(ix, nxn-2);
         temp2 = (int) min(iy, nyn-2);
         ix = (int) max(temp1,2);
         iy = (int) max(temp2,2);
         weight[1][1][0] = ((x[i] - grid->getXN(ix-1,iy-1,0))/grid->getDX())*((y[i] - grid->getYN(ix-1,iy-1,0))/grid->getDY());
         weight[1][0][0] = ((x[i] - grid->getXN(ix-1,iy,0))/grid->getDX())*((grid->getYN(ix-1,iy,0) - y[i])/grid->getDY());
         weight[0][1][0] = ((grid->getXN(ix,iy-1,0) - x[i])/grid->getDX())*((y[i] - grid->getYN(ix,iy-1,0))/grid->getDY());
         weight[0][0][0] = ((grid->getXN(ix,iy,0) - x[i])/grid->getDX())*((grid->getYN(ix,iy,0) - y[i])/grid->getDY());

         scale(weight,q[i],2,2);
         
	 // rho
         ESf->addRho(weight,ix,iy,1,ns);
         
    }
    // communicate contribution from ghost cells     
    ESf->communicateGhostP2G(ns,0,0,0,0,vct);
    delArr3(weight,2,2);
    delArr3(temp,2,2);
	 
}

/** initialize particle distribution for the relativistic two stream instability
 */
void Particles2Drel::two_beams(Grid* grid,Field* EMf,VirtualTopology* vct){
 
 double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY();
 double beam_rnd;
 int counter=0;
 /* initialize random generator */
 srand (vct->getCartesian_rank()+1+ns);
 for (int i=1; i< grid->getNXC()-1;i++)
    for (int j=1; j< grid->getNYC()-1;j++)
          for (int ii=0; ii < npcelx; ii++)
               for (int jj=0; jj < npcely; jj++){
                     x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,0);   
                     y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
                     // q = charge
                     q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,0,ns)/npcel)*(1.0/invVOL);

                     // u
	             u[counter] = c; v[counter] = c; w[counter] = c;
		     while ((fabs(u[counter])>=c) | (fabs(v[counter])>=c) | (fabs(w[counter])>=c)){
                     harvest =   rand()/(double)RAND_MAX;
                     prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                     harvest =   rand()/(double)RAND_MAX;
                     theta = 2.0*M_PI*harvest;
                     // the motion of the two beams are in the X-direction
                     beam_rnd = rand()/(double)RAND_MAX;
                     if (beam_rnd <= .5)
                       u[counter] = u0 + uth*prob*cos(theta);
                     else
                       u[counter] = -u0 + uth*prob*cos(theta);

                     // v
                     v[counter] = v0 + vth*prob*sin(theta);
                     // w
                     harvest =   rand()/(double)RAND_MAX;
                     prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                     harvest =   rand()/(double)RAND_MAX;
                     theta = 2.0*M_PI*harvest;
                     w[counter] = w0 + wth*prob*cos(theta);}
		     if (TrackParticleID)	      
			     ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                     counter++ ;
            }

 }









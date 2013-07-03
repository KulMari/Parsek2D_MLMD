/*******************************************************************************************
Particles1D.cpp  -  Class for particles of the same species in 1D
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#include <iostream>
#include <math.h>
#include "../processtopology/VirtualTopology.h"
#include "../processtopology/VCtopology.h"
#include "../inputoutput/CollectiveIO.h"
#include "../inputoutput/Collective.h"
#include "../mathlib/Basic.h"
#include "../bc/BcParticles.h"
#include "../grids/Grid.h"
#include "../grids/Grid2DCU.h"
#include "../fields/Field.h"

#include "Particles1D.h"

#include "hdf5.h"
#include <complex>

using std::cout;
using std::cerr;
using std::endl;





/**
* 
* Class for particles of the same species 1D geometry
* @date Fri Jun 4 2007
* @author Stefano Markidis, Enrico Camporeale, Enrico Camporeale, David Burgess
* @version 2.0
*
*/

/** constructor */
Particles1D::Particles1D(){
   // see allocate(int species, CollectiveIO* col, VirtualTopology* vct, Grid* grid)

}
/** deallocate particles */
Particles1D::~Particles1D(){
    delete[] x;
    delete[] u;
    delete[] q;
}

/** particles are uniformly distributed with zero velocity   */
void Particles1D::uniform_background(Grid* grid,Field* EMf){
 double dx = grid->getDX(), dy = grid->getDY();
 int counter=0;
 for (int i=1; i< grid->getNXC()-1;i++)
          for (int ii=0; ii < npcelx; ii++){
                     x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,0,0);
                     u[counter] = 0.0;
		     q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,0,0,ns)/npcel)*(1.0/invVOL);
		     if (TrackParticleID)	      
			     ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
                     counter++;
                  }
      
      cout << "Uniform Distribution with 0 velocity " << endl;
}
/** neutralizing background ions: used to neutralize completly the electrons beam  */
void Particles1D::uniform_background_ions(Grid* grid,Field* EMf){
 double dx = grid->getDX();
 int counter=0;
 for (int i=1; i< grid->getNXC()-1;i++)
          for (int ii=0; ii < npcelx; ii++){
                     x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,0,0);
                     u[counter] = 0.0;
		             q[counter] =  -(EMf->getRHOc(i,0,0)/((double)npcel))*(1.0/invVOL); // get electron charge density!
		             if (TrackParticleID)	      
			           ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
                     counter++;
                  }
}
/** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
<ul>
<li> dim = 0 --> constant velocity on X direction </li>
<li> dim = 1 --> constant velocity on Y direction </li>
<li> dim = 2 --> constant velocity on Z direction </li>
</ul>

*/
void Particles1D::constantVelocity(double vel, int dim,Grid* grid,Field* EMf){
 double harvest, prob, theta, dx = grid->getDX();
 double theta_pert;
 double mode_pert = 1.0;
 int counter=0;
    /* initialize random generator */
 srand (1+ns);
 for (int i=1; i< grid->getNXC()-1;i++)
    for (int ii=0; ii < npcelx; ii++)  {
                     // try to initialize insde the simulation box
                     x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,0,0);   
                     // q = charge
                     q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,0,0,ns)/npcel)*(1.0/invVOL);
                     // u
	                 u[counter] = c;
		             while ((fabs(u[counter])>=c)){
                       harvest =   rand()/(double)RAND_MAX;
                       prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                       harvest =   rand()/(double)RAND_MAX;
                       theta = 2.0*M_PI*harvest;
                       // the motion of the two beams are in the X-direction
                       u[counter] = u0 + uth*prob*cos(theta);
					  }
					  // perturbation
					  theta_pert = (6.2831853*mode_pert*x[counter])/Lx;
                       //x[counter] += .001*cos(theta_pert); OKKIO ALLA PERTURBAZIONE CHE PUO FARE USCIRE PARTICELLE
					 u[counter] += .00*sin(theta_pert);
                     if (TrackParticleID)	      
			         ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

					counter++ ;
            }
		

 }
/** Maxellian random velocity and uniform spatial distribution */
void Particles1D::maxwellian(Grid* grid,Field* EMf,VirtualTopology* vct){
 double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY();
 int counter=0;
/* initialize random generator */
 srand (vct->getCartesian_rank()+1+ns);
 for (int i=1; i< grid->getNXC()-1;i++)
          for (int ii=0; ii < npcelx; ii++){
                     x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,0,0);   
		     // q = charge
                     q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,0,0,ns)/npcel)*(1.0/invVOL);
                     // u
	             u[counter] = c; 
		     while (fabs(u[counter])>=c){
                               harvest =   rand()/(double)RAND_MAX;
                               prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                               harvest =   rand()/(double)RAND_MAX;
                               theta = 2.0*M_PI*harvest;
                               u[counter] = u0 + uth*prob*cos(theta);
                      }
		      if (TrackParticleID)	      
			     ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
			     counter++ ;
                       }

 }
 /** initialize particle distribution for the relativistic two stream instability
 Note that I'm using SI*/
void Particles1D::two_beamsWARP(Grid* grid,Field* EMf,VirtualTopology* vct){
 //      --- non-relativistic
 //       --- Note that in the expression for vbeam, amu is outside of the dvnz
 //        --- macro since it causes a loss of accuracy since amu is so small.
 //       --- This makes the assumption that amu would never be set to zero
 //       --- (it should never even be changed and should in fact be truly a
 //       --- constant).
 // take the u0 from the kinetic energy, override the value from the inputfile
 
 double harvest, prob, theta, dx = grid->getDX();
 double beam_rnd;
 int counter=0;
 double theta_pert;
 double mode_pert = 1.0;
 /* initialize random generator */
 srand (vct->getCartesian_rank()+1+ns);
 for (int i=1; i< grid->getNXC()-1;i++)
    for (int ii=0; ii < npcelx; ii++)  {
			// try to initialize insde the simulation box
                     x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,0,0);   
                     // q = charge
                     q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,0,0,ns)/npcel)*(1.0/invVOL);
			
	                 u[counter] = c; 
		             while ((fabs(u[counter])>=c)){
                       harvest =   rand()/(double)RAND_MAX;
                       prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                       harvest =   rand()/(double)RAND_MAX;
                       theta = 2.0*M_PI*harvest;
                       // the motion of the two beams are in the X-direction
                       beam_rnd = rand()/(double)RAND_MAX;
                       if (beam_rnd <= .5)
                         u[counter] = u0 + uth*prob*cos(theta);
                       else{
                         u[counter] = -u0 + uth*prob*cos(theta);
					
					   }
                       // perturbation
					   theta_pert = (6.2831853*mode_pert*x[counter])/Lx;
                       //x[counter] += .001*cos(theta_pert); 
                       u[counter] += .00*sin(theta_pert);
                     }
					 
					        		     if (TrackParticleID)	      
			     ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

                     counter++ ;
            }

 }
 /** two beams of particles of the same species moving in opposite direction */
 void Particles1D::two_beams(Grid* grid,Field* EMf,VirtualTopology* vct){
 double harvest, prob, theta, dx = grid->getDX();
 double beam_rnd;
 int counter=0;
 // PERTURBATION
 double theta_pert;
 double mode_pert = 1.0;
 double AmpPert = 0.01;
 
 srand (vct->getCartesian_rank()+1+ns); // particles are sampled with the same seed
 for (int i=1; i< grid->getNXC()-1;i++)
    for (int ii=0; ii < npcelx; ii++)  {
                     x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,0,0);  
			         q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,0,0,ns)/npcel)*(1.0/invVOL);  // q = charge
					 
	                 u[counter] = c;
		             while ((fabs(u[counter])>=c)){
                       harvest =   rand()/(double)RAND_MAX;
                       prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                       harvest =   rand()/(double)RAND_MAX;
                       theta = 2.0*M_PI*harvest;
                       // the motion of the two beams are in the X-direction
                       beam_rnd = rand()/(double)RAND_MAX;
                       if (beam_rnd <= .5)
                         u[counter] = u0 + uth*prob*cos(theta);
                       else{
                         u[counter] = -u0 + uth*prob*cos(theta);
					
					   }
                       // perturbation
					   theta_pert = (6.2831853*mode_pert*x[counter])/Lx;
                       x[counter] += AmpPert*sin(theta_pert); // perturbation can kick out the particle from the proc domain
                       u[counter] += .00*cos(theta_pert);
                     }
		             if (TrackParticleID)	      
			             ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
					counter++ ;
            }
			// communicate if particles are out of the sim box because of the perturbation
			communicate(vct); // communicate particles after the perturbation
	        MPI_Barrier(MPI_COMM_WORLD); // sync

 }
 /** explicit mover */
 void Particles1D::mover_explicit(Grid* grid,VirtualTopology* vct, Field* EMf){
 	int ix;
  	double  qomdt = qom*dt,weight1, weight0,Exl;
  	for (register int i=0; i < nop; i++){
	 	ix = 1 +  int(floor((x[i]-grid->getXstart())/grid->getDX()));
        	weight1 = ((x[i] - grid->getXN(ix,0,0))/grid->getDX());
		    weight0 = ((grid->getXN(ix+1,0,0) - x[i])/grid->getDX());
        	Exl = weight1*EMf->getEx(ix+1,0,0) + weight0*EMf->getEx(ix,0,0); 
        	u[i] += qomdt*Exl; 
			x[i] += u[i]*dt;
	}
        
	
	communicate(vct); // communicate particles
	MPI_Barrier(MPI_COMM_WORLD); // sync
	
 }
  /** explicit relativistic mover */
 void Particles1D::mover_relativistic(Grid* grid,VirtualTopology* vct, Field* EMf){
 	int ix;
  	double  qomdt = qom*dt,weight1, weight0,Exl;
	double gamma, px;
  	for (register int i=0; i < nop; i++){
	 	ix = 1 +  int(floor((x[i]-grid->getXstart())/grid->getDX()));
		weight1 = ((x[i] - grid->getXN(ix,0,0))/grid->getDX());
		weight0 = ((grid->getXN(ix+1,0,0) - x[i])/grid->getDX());
		Exl = weight1*EMf->getEx(ix+1,0,0) + weight0*EMf->getEx(ix,0,0); 
		// calculate gamma and normalized momentum over m0*c^2
		gamma = 1.0/sqrt(1 - (u[i]*u[i]));
		px = gamma*u[i];
		// update the momentum
		px += qomdt*Exl;
		// calculate the new gamma
		gamma = sqrt(1 + px*px);
		// calculate the new position
		u[i]  = px/gamma;  // new velocity
		x[i] += u[i]*dt;   // new position
	}
        
	// communicate particles
	communicate(vct);
	// sync
	MPI_Barrier(MPI_COMM_WORLD);
	
 }

/** print beam information */
void Particles1D::Print(VirtualTopology* ptVCT)const{
   cout << "*********************************" << endl;
   cout << "* Species = " << ns << endl;
   cout << "* Beam Energy = " << ekin << " eV" << endl;
   cout << "* Beam velocity = " << u0 << " m/s" << endl;
   cout << "* Beam current = " << ibeam << " A" << endl;
   cout << "* A ion        = " << aion << endl;
   cout << "* Z ion        = " << zion << endl;
   cout << "* weight       = " << w << endl;
   cout << "* qom          = " << qom << " C/Kg" << endl;
   cout << "* charge of particle = " << q[0] << endl;
   cout << "* time step    = " << dt << " s" << endl;
   cout << "* eps_0        = " << eps0 << endl;
   // calculate the plasma frequency for this species : length 1m
   double beam_duration = 5e-9;
   double beam_length  = beam_duration*u0;
   double part_density = ibeam/(u0*echarge*beam_length);
   double omega_p = sqrt(part_density*echarge*echarge/(emass*eps0));
   cout << "* particle density = " << part_density << " particles / m^3" << endl;
   cout << "* plasma frequency = " << omega_p << " Hz " << endl;
   cout << "* plasma period    = " << (2*pi)/omega_p << " s" << endl;
   cout << "* BEAM v/c         = " << u0/c << endl;
   cout << "* Space covered in dt by beam = " << u0*dt << " m " << endl;
   cout << "* Length for 2 stream 1 mode  = " << (2*pi)*u0/(1.4142*omega_p) << " m" <<  endl;
   cout << "* Beam permance in system     = " << 1/u0 <<  "s" << endl;
   cout << "*********************************" << endl;
}

/** alternative routine maxellian random velocity and uniform spatial distribution */
void Particles1D::alt_maxwellian(Grid* grid,Field* EMf,VirtualTopology* vct){
   
   cout << " Need to be implemented = " << endl;
}
/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void Particles1D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid* grid){
  cout << " Need to be implemented" << endl;
}
/** mover with a Predictor-Corrector scheme */
void Particles1D::mover_PC(Grid* grid,VirtualTopology* vct, Field* EMf){
  cout << "Predictor corrector scheme not implemented in the EOM" << endl;
}
/** interpolation Particle->Grid only for pressure tensor */    
void Particles1D::interpP2G_onlyP(Field* EMf, Grid *grid, VirtualTopology* vct){
  cout << "interpP2G_only not implemented in Particles 1D" << endl;
}
/** interpolation Particle->Grid only charge density, current */    
void Particles1D::interpP2G_notP(Field* EMf, Grid *grid, VirtualTopology* vct){
   cout << "interpP2G_notP not implemented in Particles 1D" << endl; 
}
/** apply a linear perturbation to particle distribution */
void Particles1D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r,     double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid* grid,Field* EMf,VirtualTopology* vct){
   cout << "linear_perturbation not implemented in Particles 1D" << endl;
}

/** Linear delta f for bi-maxwellian plasma */
double Particles1D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase,double Ez_mod, double Ez_phase, double theta, Field* EMf){
   cout << "delta_f not implemented in Particles 1D" << endl;
}





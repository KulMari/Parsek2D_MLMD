/*******************************************************************************************
Particles1D.h  -  Class for particles of the same species, in 1D
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef Part1D_H
#define Part1D_H

#include "Particles1Dcomm.h"
/**
* 
* Class for particles of the same species, in a 2D space and 3 component velocity
* 
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Particles1D : public Particles1Dcomm {
  public:
     /** constructor */
     Particles1D();
     /** destructor */
     ~Particles1D();
     /** Initial condition: uniform in space and motionless */
     void uniform_background(Grid* grid,Field* EMf);
     /** neutralizing background ion */ 
	 void uniform_background_ions(Grid* grid,Field* EMf);
	 /** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
	<ul>
		<li> dim = 0 --> constant velocity on X direction </li>
		<li> dim = 1 --> constant velocity on Y direction </li>
	</ul>
      */
	 void constantVelocity(double vel, int dim,Grid* grid,Field* EMf);
     /** Initial condition: uniform in space and maxwellian in velocity */
     void maxwellian(Grid* grid,Field* EMf, VirtualTopology* vct);
     /** two beam moving in opposite directions */
     void two_beams(Grid* grid,Field* EMf,VirtualTopology* vct);
	 /** two beam moving in opposite directions */
     void two_beamsWARP(Grid* grid,Field* EMf,VirtualTopology* vct);
     /** Initial condition: uniform in space and maxwellian in velocity */
     void alt_maxwellian(Grid* grid,Field* EMf, VirtualTopology* vct);
     /** Linear_perturbation */
     void linear_perturbation(double deltaBX, double kx, double ky, double theta, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid* grid,Field* EMf,VirtualTopology* vct);
     /**Add a periodic perturbation in velocity exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
     void AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid* grid);
/** Linear delta f for bi-maxwellian plasma */
     double delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_ampl, double Ex_phase,double Ey_ampl, double Ey_phase,double Ez_ampl, double Ez_phase,double theta, Field* EMf);
     /** Derivative of f0 wrt vpar */
     double df0_dvpar(double vpar,double vperp);
     /** Derivative of f0 wrt vperp */
     double df0_dvperp(double vpar,double vperp);
     /** Equilibrium bi-maxwellian f0 */
     double f0(double vpar,double vperp);
     /** Rotate velocities in plane XY of angle theta */
     void RotatePlaneXY(double theta);
     /** mover with the esplicit non relativistic scheme */
     void mover_explicit(Grid* grid,VirtualTopology* vct, Field* EMf);
	  /** mover with the esplicit non relativistic scheme */
     void mover_relativistic(Grid* grid,VirtualTopology* vct, Field* EMf);
     /** mover with a Predictor-Corrector Scheme */
     void mover_PC(Grid* grid,VirtualTopology* vct, Field* EMf);
     /** interpolation Particle->Grid only charge density, current */  
     void interpP2G_notP(Field* EMf, Grid *grid, VirtualTopology* vct);
     /** interpolation Particle->Grid only for pressure tensor */  
     void interpP2G_onlyP(Field* EMf, Grid *grid, VirtualTopology* vct);
     /** print beam information */
	 void Print(VirtualTopology* ptVCT)const;
};


#endif






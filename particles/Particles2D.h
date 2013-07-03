/*******************************************************************************************
Particles2D.h  -  Class for particles of the same species, in a 2D space and 3 component velocity
                            -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/

#ifndef Part2D_H
#define Part2D_H

#include "Particles2Dcomm.h"
#include "../communication/ComParticles.h"

/**
* 
* Class for particles of the same species, in a 2D space and 3 component velocity
* 
* @date Fri Jun 4 2007
* @author Stefano Markidis, Giovanni Lapenta, Enrico Camporeale, David Burgess
* @version 2.0
*
*/
class Particles2D : public Particles2Dcomm {
  public:
     /** constructor */
     Particles2D();
     /** destructor */
     ~Particles2D();
     /** Initial condition: uniform in space and motionless */
     void uniform_background(Grid* grid,Field* EMf);
     /** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
	<ul>
		<li> dim = 0 --> constant velocity on X direction </li>
		<li> dim = 1 --> constant velocity on Y direction </li>
	</ul>
      */
     void constantVelocity(double vel, int dim,Grid* grid,Field* EMf);
	 /** Delete the particles inside the sphere with radius R and center x_center y_center */
     int deleteParticlesInsideSphere(double radius, double x_center, double y_center);
     /** initialize particles for double harris **/
     void DoubleHarris(Grid* grid,Field* EMf,VirtualTopology* vct);
     /** Initial condition: RANDOM in space and maxwellian in velocity */
     void maxwellian(Grid* grid,Field* EMf, VirtualTopology* vct);
     /** Initial condition: RANDOM in space and maxwellian in velocity */
     void maxwellian_ball(Grid* grid,Field* EMf, VirtualTopology* vct);
     /** Initial condition: RANDOM in space and maxwellian in velocity: test routine to have the same particle init condition with different topologies */
     int maxwellian_sameParticleInit(Grid* grid,Field* EMf, VirtualTopology* vct);
     /** force free */
     void force_free(Grid *grid, Field* EMf, VirtualTopology* vct);
     /** Initial condition: uniform in space and maxwellian in velocity */
     void alt_maxwellian(Grid* grid,Field* EMf, VirtualTopology* vct);
	 /** Injection of a beam at a certain position and velocity*/
     int inject_beam(double x_source, double y_source, double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, double radius_beam, int np_inj, double qom_beam, double rho_init, VirtualTopology* vct);
     /** inject particles with uniform distribution from the wall */
	 int inject_Xleft_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind);
	 /** inject particles with uniform distribution from the wall */
	 int inject_Xright_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind);
	  /** inject particles with uniform distribution from the wall */
	 int inject_Yleft_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind);
	 /** inject particles with uniform distribution from the wall */
	 int inject_Yright_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind);
	 /** generate the vleocity of the particle with the rejection method */
     double generate_v(double v00, double vvth);
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
     /** mover with a Predictor-Corrector Scheme */
     int mover_PC(Grid* grid,VirtualTopology* vct, Field* EMf);
	 /** relativistic mover with a Predictor-Corrector scheme */
     int mover_relativistic(Grid* grid,VirtualTopology* vct, Field* EMf);
	 /** mover with a Predictor-Corrector scheme wit dipole and earth */
     int moverDipole(double Bdipole, double x_dipole, double y_dipole, double radius_Earth, Grid* grid,VirtualTopology* vct, Field* EMf);
     /** interpolation Particle->Grid only charge density, current */  
     void interpP2G_notP(Field* EMf, Grid *grid, VirtualTopology* vct);
     /** interpolation Particle->Grid only for pressure tensor */  
     void interpP2G_onlyP(Field* EMf, Grid *grid, VirtualTopology* vct);
     /** get the total number of particles  in the system */
     int getNPsystem(VirtualTopology* vct);
	 /** get the total number of particles eliminated at the boundary */
     int getNPdeletedBoundary(VirtualTopology* vct);
	 /** get the total number of particles eliminated from the dipole */
     int getNPdeletedDipole(VirtualTopology* vct);
     // ME
     // repopulate the ghost cells of the coarse level
     int RepopulateCoarseLevel(Grid *grid, VirtualTopology* vct);

};


#endif






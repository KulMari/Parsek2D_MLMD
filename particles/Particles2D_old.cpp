/*******************************************************************************************
Particles2D.cpp  -  Class for particles of the same species, in a 2D space and 3component velocity
-------------------
developers: Stefano Markidis, Giovanni Lapenta
	********************************************************************************************/
#include "../processtopology/VirtualTopology.h"
#include "../processtopology/VCtopology.h"
#include <iostream>
#include <math.h>
#include "../inputoutput/CollectiveIO.h"
#include "../inputoutput/Collective.h"
#include "../mathlib/Basic.h"
#include "../bc/BcParticles.h"
#include "../grids/Grid.h"
#include "../grids/Grid2DCU.h"
#include "../fields/Field.h"

#include "Particles2D.h"

#include "hdf5.h"
#include <complex>

using std::cout;
using std::cerr;
using std::endl;

#define min(a,b) (((a)<(b))?(a):(b));
#define max(a,b) (((a)>(b))?(a):(b));
#define MIN_VAL   1E-16
/**
* 
 * Class for particles of the same species, in a 2D space and 3component velocity
 * @date Fri Jun 4 2007
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 2.0
 *
 */

/** constructor */
Particles2D::Particles2D(){
	// see allocate(int species, CollectiveIO* col, VirtualTopology* vct, Grid* grid)
	
}
/** deallocate particles */
Particles2D::~Particles2D(){
    delete[] x;
    delete[] y;
    delete[] u;
    delete[] v;
    delete[] w;
    delete[] q;
    delete[] xptilde;
    delete[] yptilde;
    delete[] uptilde;
    delete[] vptilde;
    delete[] wptilde;
}
/** Delete the particles inside the sphere with radius R and center x_center y_center */
int Particles2D::deleteParticlesInsideSphere(double radius, double x_center, double y_center){
	// take track how many particles are eliminated because inside the earth
	npDeletedDipole = 0;
	
	int np_current = 0, nplast = nop-1;
	double distance_sq=0.0;
	while (np_current < nplast+1){
		distance_sq = (x[np_current] - x_center)*(x[np_current] - x_center) + (y[np_current] - y_center)*(y[np_current] - y_center);
		if (distance_sq < radius*radius){
			// delete the particle and pack the particle array, the value of nplast changes
			del_pack(np_current,&nplast);
			npDeletedDipole++;
		} else {
			// particle is still in the domain, procede with the next particle
			np_current++;
		}
	}
	nop = nplast +1;
	return(npDeletedDipole);
}
/** Double Harris **/
void Particles2D::DoubleHarris(Grid* grid,Field* EMf,VirtualTopology* vct){
	double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY();
        double Ox,Oy,coarsedx,coarsedy;
	int counter=0;
	double shaperx, shapery, shaperz;
	double flvx=1.0,flvy=1.0,flvz=1.0;
	double segno = (qom/fabs(qom));
        int xcS, xcE, ycS, ycE;

        //particles also in the GC of the coarse grid
        //if (vct->getXleft_neighbor()== MPI_PROC_NULL) {xcS=0;}
        if (vct->getXleft_neighbor()== MPI_PROC_NULL) {
            xcS=0;
            cout << "Proc " << vct->getCartesian_rank_COMMTOTAL() << "has MPI_PROC_NULL left X neighbor" << endl;
        } //if periodic boundary conditions, no particles in GC
        else {xcS=1;}
        //if (vct->getXright_neighbor()== MPI_PROC_NULL) {xcE=grid->getNXC();}
        if (vct->getXright_neighbor()== MPI_PROC_NULL ) {
            xcE=grid->getNXC();
            cout << "Proc " << vct->getCartesian_rank_COMMTOTAL() << "has MPI_PROC_NULL right X neighbor" << endl;
        } 
        else {xcE=grid->getNXC()-1;}
        //if (vct->getYleft_neighbor()== MPI_PROC_NULL) {ycS=0;}
        if (vct->getYleft_neighbor()== MPI_PROC_NULL ) {
            ycS=0;
            cout << "Proc " << vct->getCartesian_rank_COMMTOTAL() << "has MPI_PROC_NULL left Y neighbor" << endl;
            } 
        else {ycS=1;}
        //if (vct->getYright_neighbor()== MPI_PROC_NULL) {ycE=grid->getNYC();}
        if (vct->getYright_neighbor()== MPI_PROC_NULL ) {
            ycE=grid->getNYC();
            cout << "Proc " << vct->getCartesian_rank_COMMTOTAL() << "has MPI_PROC_NULL right Y neighbor" << endl;
        } 
        else {ycE=grid->getNYC()-1;}
        /*xcS = 1;
        ycS = 1;
        xcE = grid->getNXC()-1;
        ycE = grid->getNYC()-1;*/

        Ox = 0;
        Oy = 0;
        for (int i=1; i < grid->getLevel()+1; i++){
            Ox += grid->getOx(i);   
            Oy += grid->getOy(i);   
        }
        coarsedx = grid->getDX()*pow(grid->getRatio(),grid->getLevel()); 
        coarsedy = grid->getDY()*pow(grid->getRatio(),grid->getLevel()); 

	/* initialize random generator */
	srand (vct->getCartesian_rank()+1+ns);
	for (int i=xcS; i< xcE;i++)
		for (int j=ycS; j< ycE;j++)
			for (int ii=0; ii < npcelx; ii++)
				for (int jj=0; jj < npcely; jj++){
					x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,0);
					y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
					// q = charge
                                        /*/if ( vct->getCartesian_rank_COMMTOTAL()==72){
                                            cout << i << " " << j <<" "<<EMf->getRHOcs(i,j,0,ns)<< endl;
                                        }*/
					q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,0,ns)/npcel)*(1.0/invVOL);
					//shaperx = tanh((y[counter] - Ly/2)/delta)/cosh((y[counter] - Ly/2)/delta) + 0.1*(M_PI/Ly)*cos(2*M_PI*x[counter]/Lx)*sin(M_PI*(y[counter]- Ly/2)/Ly);
                    //shaperz = 1.0/(cosh((y[counter] - Ly/2)/delta)*cosh((y[counter] - Ly/2)/delta));
					//shapery = shapery - 0.1*(2*M_PI/Lx)*sin(2*M_PI*x[counter]/Lx)*cos((M_PI*y[counter]- Ly/2)/Ly);
					//shaperz = -tanh((y[counter] - Ly/2)/delta) ;
					if (grid->getLevel()==0)
					  shaperz = -tanh((y[counter] - Ly/2)/delta) ; 
					else if (grid->getLevel()==1)
					  shaperz = tanh((y[counter] - Ly)/delta) ;
					else
					  cout <<"This init will work only for 2 grids, with refined grid in the uèpper half\n";
					  
					shaperx = 1.0;
					shapery = 1.0;
					// new drift velocity to satisfy JxB=0
					//flvx =u0*flvx*shaperx;
					//flvz =w0*flvz*shaperz;
					//flvy =v0*flvy*shapery;
					flvx =u0*shaperx;
					flvz =w0*shaperz;
					flvy =v0*shapery;
					u[counter] = c; v[counter] = c; w[counter] = c;
					while ((fabs(u[counter])>=c) | (fabs(v[counter])>=c) | (fabs(w[counter])>=c)){
					    harvest =   rand()/(double)RAND_MAX;
						prob  = sqrt(-2.0*log(1.0-.999999*harvest));
						harvest =   rand()/(double)RAND_MAX;
						theta = 2.0*M_PI*harvest;
						u[counter] = flvx + uth*prob*cos(theta);
						// v
						v[counter] = flvy + vth*prob*sin(theta);
						// w
						harvest =   rand()/(double)RAND_MAX;
						prob  = sqrt(-2.0*log(1.0-.999999*harvest));
						harvest =   rand()/(double)RAND_MAX;
						theta = 2.0*M_PI*harvest;
						w[counter] = flvz + wth*prob*cos(theta);}
					if (TrackParticleID)
						ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];

					counter++ ;
				}
  nop= counter;

}

/** particles are uniformly distributed with zero velocity   */
void Particles2D::uniform_background(Grid* grid,Field* EMf){
	int counter=0;
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int ii=0; ii < npcelx; ii++)
				for (int jj=0; jj < npcely; jj++){
					x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,0);
					y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
					u[counter] = 0.0;
					v[counter] = 0.0;
					w[counter] = 0.0;
					q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,0,ns)/npcel)*(1.0/invVOL);
					if (TrackParticleID)	      
						ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
					counter++;
				}
					
					cout << "Velocity Maxwellian Distribution " << endl;
}
/** Initialize particles with a constant velocity in dim direction. Depending on the value of dim:
<ul>
<li> dim = 0 --> constant velocity on X direction </li>
<li> dim = 1 --> constant velocity on Y direction </li>
<li> dim = 2 --> constant velocity on Z direction </li>
</ul>

*/
void Particles2D::constantVelocity(double vel, int dim,Grid* grid,Field* EMf){
	switch(dim)
	{
		case 0:
			for (int i=0; i < nop; i++)
				u[i] = vel, v[i] = 0.0, w[i] = 0.0;
			break;
		case 1:
			for (register int i=0; i < nop; i++)
				u[i] = 0.0,v[i] = vel,w[i] = 0.0;
			break;
		case 2:
			for (register int i=0; i < nop; i++)
				u[i] = 0.0, v[i] = 0.0, w[i] = vel;
			break;
			
	}
	
}

/** alternative routine maxellian random velocity and uniform spatial distribution */
void Particles2D::alt_maxwellian(Grid* grid,Field* EMf,VirtualTopology* vct){
	double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY(), U,V,W,X,Y, weight00,weight01,weight10,weight11;
	int counter=0;
	/* initialize random generator */
	srand (vct->getCartesian_rank()+1+ns);
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int ii=0; ii < npcelx/2; ii++)
				for (int jj=0; jj < npcely/4; jj++){
					// u
					harvest =   rand()/(double)RAND_MAX;
					prob  = sqrt(-2.0*log(1.0-.999999*harvest));
					harvest =   rand()/(double)RAND_MAX;
					theta = 2.0*M_PI*harvest;
					U = u0 + uth*prob*cos(theta);
					// v
					V = v0*(qom/fabs(qom)) + vth*prob*sin(theta);
					// w
					harvest =   rand()/(double)RAND_MAX;
					prob  = sqrt(-2.0*log(1.0-.999999*harvest));
					harvest =   rand()/(double)RAND_MAX;
					theta = 2.0*M_PI*harvest;
					W = w0 + wth*prob*cos(theta);
					X = (ii + .5)*(2*dx/npcelx) + grid->getXN(i,j,0);   
					Y = (jj + .5)*(4*dy/npcely) + grid->getYN(i,j,0);      		     
					
					for (register int p=0; p<8; p++){
						x[counter] = X;
						y[counter] = Y;
						weight11 = ((x[counter] - grid->getXN(i,j,0))/dx)*((y[counter] - grid->getYN(i,j,0))/dy);
						weight10 = ((x[counter] - grid->getXN(i,j+1,0))/dx)*((grid->getYN(i,j+1,0) - y[counter])/dy);
						weight01 = ((grid->getXN(i+1,j,0) - x[counter])/dx)*((y[counter] - grid->getYN(i+1,j,0))/dy);
						weight00 = ((grid->getXN(i+1,j+1,0) - x[counter])/dx)*((grid->getYN(i+1,j+1,0) - y[counter])/dy);
						u[counter] = U*pow(-1.0,p) ;
						v[counter] = V*pow(-1.0,p/2);
						w[counter] = W*pow(-1.0,p/4);
						q[counter] =  (qom/fabs(qom))*(weight00*EMf->getRHOns(i,j,0,ns)+weight01*EMf->getRHOns(i,j+1,0,ns)+weight10*EMf->getRHOns(i+1,j,0,ns)+weight11*EMf->getRHOns(i+1,j+1,0,ns))*(1.0/invVOL/npcel);
						
						if (TrackParticleID)	      
							ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
						counter++ ;}
				}
					
}

/** Maxellian random velocity and RANDOM spatial distribution */
void Particles2D::maxwellian(Grid* grid,Field* EMf,VirtualTopology* vct){
  //cout << "Maxwellian also in ghost cells\n"; 
  //AMR, ME: also first and last ghost cell with particles in the refined grids
  int xcS, xcE, ycS, ycE;

  //particles also in the GC of the coarse grid

  if (vct->getXleft_neighbor()== MPI_PROC_NULL) {xcS=0;} //if periodic boundary conditions, no particles in GC
  else {xcS=1;}

  if (vct->getXright_neighbor()== MPI_PROC_NULL) {xcE=grid->getNXC();} 
  else {xcE=grid->getNXC()-1;}

  if (vct->getYleft_neighbor()== MPI_PROC_NULL) {ycS=0;} 
  else {ycS=1;}

  if (vct->getYright_neighbor()== MPI_PROC_NULL) {ycE=grid->getNYC();} 
  else {ycE=grid->getNYC()-1;}
  
  double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY();
  int counter=0;
  /* initialize random generator */
  srand (vct->getCartesian_rank()+1+ns);  //different seeding if different XLEN, YLEN
  //srand(5);
  for (int i=xcS; i< xcE;i++)
    //for (int i=1; i< grid->getNXC()-1;i++)
    for (int j=ycS; j< ycE;j++)
      //for (int j=1; j< grid->getNYC()-1;j++)
      for (int ii=0; ii < npcelx; ii++)
	for (int jj=0; jj < npcely; jj++){
	  // in fixed position
	  //x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,0);   
	  //y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
	  // random position in the cell
	  x[counter]=grid->getXN(i,j,0) + (double)(rand()/(double)RAND_MAX * dx);
	  y[counter]=grid->getYN(i,j,0) + (double)(rand()/(double)RAND_MAX * dy);
	  // q = charge
	  q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,0,ns)/npcel)*(1.0/invVOL);
	  // u
	  u[counter] = c; v[counter] = c; w[counter] = c;
	  while ((fabs(u[counter])>=c) | (fabs(v[counter])>=c) | (fabs(w[counter])>=c)){
	    harvest =   rand()/(double)RAND_MAX;
	    prob  = sqrt(-2.0*log(1.0-.999999*harvest));
	    harvest =   rand()/(double)RAND_MAX;
	    theta = 2.0*M_PI*harvest;
	    u[counter] = u0 + uth*prob*cos(theta);
	    // v
	    v[counter] = v0 + vth*prob*sin(theta);
	    // w
	    harvest =   rand()/(double)RAND_MAX;
	    prob  = sqrt(-2.0*log(1.0-.999999*harvest));
	    harvest =   rand()/(double)RAND_MAX;
	    theta = 2.0*M_PI*harvest;
	    w[counter] = w0 + wth*prob*cos(theta);}
	  if (TrackParticleID){	      
	    ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
	  }
	  counter++ ;
	}
 
  nop= counter;
 
}

/** Maxellian random velocity and RANDOM spatial distribution within a ball */
/*void Particles2D::maxwellian_ball(Grid* grid,Field* EMf,VirtualTopology* vct){
  //cout << "Maxwellian also in ghost cells\n"; 
  //AMR, ME: also first and last ghost cell with particles in the refined grids
  int xcS, xcE, ycS, ycE,i;
  double Rsquared, xcenter, ycenter;
  Rsquared=0.01;
  xcenter = 0.5;//X Coordinate of the ball center in the coarsest grid
  ycenter = 0.5;//Y Coordinate of the ball center in the coarsest grid

  for (i=grid->getLevel();i>0;i-=1){
      xcenter-= grid->getOx(i);
      ycenter-= grid->getOy(i);
  }

  //particles also in the GC of the coarse grid
  if (vct->getXleft_neighbor()== MPI_PROC_NULL) {xcS=0;}
  else {xcS=1;}
  if (vct->getXright_neighbor()== MPI_PROC_NULL) {xcE=grid->getNXC();}
  else {xcE=grid->getNXC()-1;}
  if (vct->getYleft_neighbor()== MPI_PROC_NULL) {ycS=0;}
  else {ycS=1;}
  if (vct->getYright_neighbor()== MPI_PROC_NULL) {ycE=grid->getNYC();}
  else {ycE=grid->getNYC()-1;}
  
  double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY();
  int counter=0;
  // initialize random generator 
  srand (vct->getCartesian_rank()+1+ns);  //different seeding if different XLEN, YLEN
  //srand(5);
  if (grid->getXstart()< xcenter + sqrt(Rsquared) && grid->getXend()> xcenter - sqrt(Rsquared) && grid->getYstart()< ycenter + sqrt(Rsquared) && grid->getYend()> ycenter - sqrt(Rsquared)){   
      for (int i=xcS; i< xcE;i++)
        //for (int i=1; i< grid->getNXC()-1;i++)
        for (int j=ycS; j< ycE;j++)
          //for (int j=1; j< grid->getNYC()-1;j++)
          for (int ii=0; ii < npcelx; ii++)
            for (int jj=0; jj < npcely; jj++){
              // in fixed position
              //x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,0);   
              //y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
              // random position in the cell
              x[counter]=grid->getXN(i,j,0) + (double)(rand()/(double)RAND_MAX * dx);
              y[counter]=grid->getYN(i,j,0) + (double)(rand()/(double)RAND_MAX * dy);
              // q = charge
              if ((x[counter]-xcenter)*(x[counter]-xcenter)+(y[counter]-ycenter)*(y[counter]-ycenter) >Rsquared){
                  q[counter] =  0;
              } else {
                  q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,0,ns)/npcel)*(1.0/invVOL);
              }
              // u
              u[counter] = c; v[counter] = c; w[counter] = c;
              while ((fabs(u[counter])>=c) | (fabs(v[counter])>=c) | (fabs(w[counter])>=c)){
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                u[counter] = u0 + uth*prob*cos(theta);
                // v
                v[counter] = v0 + vth*prob*sin(theta);
                // w
                harvest =   rand()/(double)RAND_MAX;
                prob  = sqrt(-2.0*log(1.0-.999999*harvest));
                harvest =   rand()/(double)RAND_MAX;
                theta = 2.0*M_PI*harvest;
                w[counter] = w0 + wth*prob*cos(theta);}
              if (TrackParticleID){	      
                ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
              }
              counter++ ;
            }
   }
      nop= counter;
}*/

/** Maxellian random velocity and uniform spatial distribution */
int Particles2D::maxwellian_sameParticleInit(Grid* grid,Field* EMf,VirtualTopology* vct){
  //AMR, ME: also first and last ghost cell with particles in the refined grids
  // only works with uniform density
  int xcS, xcE, ycS, ycE;

  //particles also in the GC of the coarse grid
  xcS=0;
  xcE=(grid->getNXC()-2)*vct->getXLEN()+2;
  ycS=0;
  ycE=(grid->getNYC()-2)*vct->getYLEN()+2;
  
  double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY();
  int counter=0;
  /* initialize random generator */
  //srand (vct->getCartesian_rank()+1+ns);  //different seeding if different XLEN, YLEN
  srand(5);
  if ( vct->getCartesian_rank_COMMTOTAL() % (vct->getXLEN()*vct->getYLEN())==0  )
    {
      cout << "R" <<vct->getCartesian_rank_COMMTOTAL() <<" generating particles"<<endl;
      cout << "xcS " << xcS <<" xcE " <<xcE <<endl;
      
      for (int i=xcS; i< xcE;i++)
	for (int j=ycS; j< ycE;j++)
	  for (int ii=0; ii < npcelx; ii++)
	    for (int jj=0; jj < npcely; jj++){
	      // in fixed position
	      //x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,0);   
	      //y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
	      // random position in the cell
	      //x[counter]=grid->getXN(i,j,0) + (double)(rand()/(double)RAND_MAX * dx);
	      //y[counter]=grid->getYN(i,j,0) + (double)(rand()/(double)RAND_MAX * dy);
	      // with respect to the normal generation, modidy the starting point because here not all the grid is visible
	      x[counter]=grid->getDX()*(i-1) + (double)(rand()/(double)RAND_MAX * dx);  
              y[counter]=grid->getDY()*(j-1) + (double)(rand()/(double)RAND_MAX * dy);
	      // q = charge
	      //q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,0,ns)/npcel)*(1.0/invVOL);
	      q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(1,1,0,ns)/npcel)*(1.0/invVOL); 
	      // u
	      u[counter] = c; v[counter] = c; w[counter] = c;
	      while ((fabs(u[counter])>=c) | (fabs(v[counter])>=c) | (fabs(w[counter])>=c)){
		harvest =   rand()/(double)RAND_MAX;
		prob  = sqrt(-2.0*log(1.0-.999999*harvest));
		harvest =   rand()/(double)RAND_MAX;
		theta = 2.0*M_PI*harvest;
		u[counter] = u0 + uth*prob*cos(theta);
		// v
		v[counter] = v0 + vth*prob*sin(theta);
		// w
		harvest =   rand()/(double)RAND_MAX;
		prob  = sqrt(-2.0*log(1.0-.999999*harvest));
		harvest =   rand()/(double)RAND_MAX;
		theta = 2.0*M_PI*harvest;
		w[counter] = w0 + wth*prob*cos(theta);}
	      if (TrackParticleID){      
		ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
	      }
	      counter++ ;
	    }
      
      nop= counter;
      
    }// end particle generation in the proc with lower rank in the grid
  
  int BC_partCommunicate=0;
  
  int avail = communicate(vct, grid, BC_partCommunicate);

  if (avail < 0)
    return(-1);
  MPI_Barrier(vct->getCART_COMM());
  while(isMessagingDone(vct) >0)
    {
      avail = communicate(vct, grid, BC_partCommunicate);
    
      if (avail < 0)
	return(-1);
      MPI_Barrier(vct->getCART_COMM());
    }
}

/** Force Free initialization (JxB=0) for particles */
void Particles2D::force_free(Grid* grid,Field* EMf,VirtualTopology* vct){
	double harvest, prob, theta, dx = grid->getDX(),dy = grid->getDY();
	int counter=0;
	double shaperx, shapery, shaperz;
	double flvx=1.0,flvy=1.0,flvz=1.0;
	double segno = (qom/fabs(qom));
	
	/* initialize random generator */
	srand (vct->getCartesian_rank()+1+ns);
	for (int i=1; i< grid->getNXC()-1;i++)
		for (int j=1; j< grid->getNYC()-1;j++)
			for (int ii=0; ii < npcelx; ii++)
				for (int jj=0; jj < npcely; jj++){
					flvx=1.0;
                                        flvy=1.0;
                                        flvz=1.0;
                                        x[counter] = (ii + .5)*(dx/npcelx) + grid->getXN(i,j,0);   
					y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
					// q = charge
					q[counter] =  (qom/fabs(qom))*(EMf->getRHOcs(i,j,0,ns)/npcel)*(1.0/invVOL);
					// delta_v = 10 * delta_B
					shaperx = 1/(cosh((y[counter] - Ly/2)/(10*delta))*cosh((y[counter] - Ly/2)/(10*delta)))/(10*delta);
					//shaperx = tanh((y[counter] - Ly/2)/delta)/cosh((y[counter] - Ly/2)/delta)/delta;
                    //shaperz = 1.0/(cosh((y[counter] - Ly/2)/delta)*cosh((y[counter] - Ly/2)/delta))/delta;
					shaperz = 0.0;
					shapery = 0.0;
					// new drift velocity to satisfy JxB=0
					flvx =u0*flvx*shaperx;
					flvz =w0*flvz*shaperz;
					flvy =v0*flvy*shapery;
					u[counter] = c; v[counter] = c; w[counter] = c;
					while ((fabs(u[counter])>=c) | (fabs(v[counter])>=c) | (fabs(w[counter])>=c)){
					    harvest =   rand()/(double)RAND_MAX;
						prob  = sqrt(-2.0*log(1.0-.999999*harvest));
						harvest =   rand()/(double)RAND_MAX;
						theta = 2.0*M_PI*harvest;
						u[counter] = flvx + uth*prob*cos(theta);
						// v
						v[counter] = flvy + vth*prob*sin(theta);
						// w
						harvest =   rand()/(double)RAND_MAX;
						prob  = sqrt(-2.0*log(1.0-.999999*harvest));
						harvest =   rand()/(double)RAND_MAX;
						theta = 2.0*M_PI*harvest;
						w[counter] = flvz + wth*prob*cos(theta);}
					if (TrackParticleID)	      
						ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
					
					counter++ ;
				}
					
}

/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
inline void Particles2D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0, Grid* grid){
	
	// rescaling of amplitudes according to deltaBoB //
	double alpha;
	alpha=deltaBoB*B0/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);
	jx_mod *=alpha;
	jy_mod *=alpha;
	jz_mod *=alpha;
	for (register int i=0; i<nop; i++){
		u[i] += jx_mod/q[i]/npcel/invVOL* cos(kx*x[i] + ky*y[i] + jx_phase);
		v[i] += jy_mod/q[i]/npcel/invVOL* cos(kx*x[i] + ky*y[i] + jy_phase);
		w[i] += jz_mod/q[i]/npcel/invVOL* cos(kx*x[i] + ky*y[i] + jz_phase);}
}


/** explicit mover */
void Particles2D::mover_explicit(Grid* grid,VirtualTopology* vct, Field* EMf){
 	int innter,temp1,temp2,ix,iy;
  	double  qomdt = qom*dt,weight11, weight00, weight10, weight01,Exl, Eyl, Ezl, Bxl, Byl, Bzl;
  	for (register int i=0; i < nop; i++){
	 	ix = 2 +  int(floor((x[i]-xstart)/dx));
		iy = 2 +  int(floor((y[i]-ystart)/dy));
		weight11 = ((x[i] - grid->getXN(ix-1,iy-1,0))/dx)*((y[i] - grid->getYN(ix-1,iy-1,0))/dy);
		weight10 = ((x[i] - grid->getXN(ix-1,iy,0))/dx)*((grid->getYN(ix-1,iy,0) - y[i])/dy);
		weight01 = ((grid->getXN(ix,iy-1,0) - x[i])/dx)*((y[i] - grid->getYN(ix,iy-1,0))/dy);
		weight00 = ((grid->getXN(ix,iy,0) - x[i])/dx)*((grid->getYN(ix,iy,0) - y[i])/dy);
		Exl = weight00*EMf->getEx(ix-1,iy-1,0) + weight01*EMf->getEx(ix-1,iy,0) + weight10*EMf->getEx(ix,iy-1,0) + weight11*EMf->getEx(ix,iy,0);
		Eyl = weight00*EMf->getEy(ix-1,iy-1,0) + weight01*EMf->getEy(ix-1,iy,0) + weight10*EMf->getEy(ix,iy-1,0) + weight11*EMf->getEy(ix,iy,0);
		Ezl = weight00*EMf->getEz(ix-1,iy-1,0) + weight01*EMf->getEz(ix-1,iy,0) + weight10*EMf->getEz(ix,iy-1,0) + weight11*EMf->getEz(ix,iy,0);
		Bxl = weight00*EMf->getBx(ix-1,iy-1,0) + weight01*EMf->getBx(ix-1,iy,0) + weight10*EMf->getBx(ix,iy-1,0) + weight11*EMf->getBx(ix,iy,0);
		Byl = weight00*EMf->getBy(ix-1,iy-1,0) + weight01*EMf->getBy(ix-1,iy,0) + weight10*EMf->getBy(ix,iy-1,0) + weight11*EMf->getBy(ix,iy,0);
		Bzl = weight00*EMf->getBz(ix-1,iy-1,0) + weight01*EMf->getBz(ix-1,iy,0) + weight10*EMf->getBz(ix,iy-1,0) + weight11*EMf->getBz(ix,iy,0);
		x[i] += u[i]*dt;
		y[i] += v[i]*dt;
		u[i] += qomdt*(Exl + v[i]*Bzl -w[i]*Byl);  
		v[i] += qomdt*(Eyl + w[i]*Bxl -u[i]*Bzl);
		
	}
	communicate(vct);
	MPI_Barrier(vct->getCART_COMM());
	
}
/** mover with a Predictor-Corrector scheme */
int Particles2D::mover_PC(Grid* grid,VirtualTopology* vct, Field* EMf){
  int BC_partCommunicate;
  int cartesian_rank_total;
  MPI_Comm_rank(vct->getCART_COMM_TOTAL(), &cartesian_rank_total);     
  /*if (0)
    {
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<<"qom" <<qom << "Nop Before Mover: " << nop << endl;
      cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<<"qom" <<qom << " np_REPOP_b_BOTTOM " << np_REPOP_b_BOTTOM << " np_REPOP_b_TOP " << np_REPOP_b_TOP << " np_REPOP_b_LEFT " << np_REPOP_b_LEFT << " np_REPOP_b_RIGHT " << np_REPOP_b_RIGHT << endl;
      }*/

  int innter,temp1,temp2;
  int avail;
  double dto2 = .5*dt, qomdt = qom*dto2, omdtsq, denom, ut, vt, wt, udotb, weight00, weight01,weight10,weight11;
  double Exl=0.0, Eyl=0.0, Ezl=0.0, Bxl=0.0, Byl=0.0, Bzl=0.0;
  double inv_dx = 1.0/dx;
  double inv_dy = 1.0/dy;
  int ix,iy;
  for (int i=0; i < nop; i++){
    xptilde[i] = x[i];
    yptilde[i] = y[i];

    //for MLDM ops
    if (FinerLevels_PRAOps==1 && PRACollectionMethod==0)
      {
        AlreadyAccumulated[i]= false;
      }
  }
  innter=0;
  while (innter< NiterMover)
  {	  
    // move each particle with new fields
    for (int i=0; i <  nop; i++)
    {
      // interpolation G-->P
      ix = 2 +  int(floor((x[i]-xstart)/dx));
      iy = 2 +  int(floor((y[i]-ystart)/dy));
                  
      // ME, AMR
      /*if (ix-1 <0 || iy-1 <0 || ix> grid->getNXN()-1 || iy> grid->getNYN()-1)
      {
	cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()  << "Particle mess in the mover, exiting " << endl;
	cout << "x: " << x[i] << ", y: " << y[i] << ", qom: " << qom << ", ID: " << ParticleID[i]<< ", PRA_oxStartLeft: " <<  PRA_oxStartLeft <<", PRA_oyStartLeft: " <<  PRA_oyStartLeft<< ", xstart: " << xstart << ", xend: " << xend << ", ystart: " << ystart << ", yend: " << yend<<", xstart-dx: " << xstart-dx <<",ystart-dy: " << ystart-dy<<", Lx: " <<Lx <<", Ly: " <<Ly  << ", i: " <<i << ", nop: "<<nop <<", innter: " <<innter <<  endl;
	return -1;
	}*/
      // end ME, AMR
      weight11 = ((x[i] - grid->getXN(ix-1,iy-1,0))*inv_dx)*((y[i] - grid->getYN(ix-1,iy-1,0))*inv_dy);
      weight10 = ((x[i] - grid->getXN(ix-1,iy,0))*inv_dx)*((grid->getYN(ix-1,iy,0) - y[i])*inv_dy);
      weight01 = ((grid->getXN(ix,iy-1,0) - x[i])*inv_dx)*((y[i] - grid->getYN(ix,iy-1,0))*inv_dy);
      weight00 = ((grid->getXN(ix,iy,0) - x[i])*inv_dx)*((grid->getYN(ix,iy,0) - y[i])*inv_dy);
      Exl = weight00*EMf->getEx(ix-1,iy-1,0) + weight01*EMf->getEx(ix-1,iy,0) + weight10*EMf->getEx(ix,iy-1,0) + weight11*EMf->getEx(ix,iy,0);
      Eyl = weight00*EMf->getEy(ix-1,iy-1,0) + weight01*EMf->getEy(ix-1,iy,0) + weight10*EMf->getEy(ix,iy-1,0) + weight11*EMf->getEy(ix,iy,0);
      Ezl = weight00*EMf->getEz(ix-1,iy-1,0) + weight01*EMf->getEz(ix-1,iy,0) + weight10*EMf->getEz(ix,iy-1,0) + weight11*EMf->getEz(ix,iy,0);
      Bxl = weight00*EMf->getBx(ix-1,iy-1,0) + weight01*EMf->getBx(ix-1,iy,0) + weight10*EMf->getBx(ix,iy-1,0) + weight11*EMf->getBx(ix,iy,0);
      Byl = weight00*EMf->getBy(ix-1,iy-1,0) + weight01*EMf->getBy(ix-1,iy,0) + weight10*EMf->getBy(ix,iy-1,0) + weight11*EMf->getBy(ix,iy,0);
      Bzl = weight00*EMf->getBz(ix-1,iy-1,0) + weight01*EMf->getBz(ix-1,iy,0) + weight10*EMf->getBz(ix,iy-1,0) + weight11*EMf->getBz(ix,iy,0);
      
      // end interpolation
      omdtsq = qomdt*qomdt/c/c*(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
      denom = 1.0/(1.0 + omdtsq);
      // solve the position equation
      ut= u[i] + qomdt*Exl;
      vt= v[i] + qomdt*Eyl;
      wt= w[i] + qomdt*Ezl;
      udotb = ut*Bxl + vt*Byl + wt*Bzl;
      // solve the velocity equation 
      uptilde[i] = (ut+qomdt/c*(vt*Bzl -wt*Byl + qomdt/c*udotb*Bxl))*denom; 
      vptilde[i] = (vt+qomdt/c*(wt*Bxl -ut*Bzl + qomdt/c*udotb*Byl))*denom;
      wptilde[i] = (wt+qomdt/c*(ut*Byl -vt*Bxl + qomdt/c*udotb*Bzl))*denom;
      // update position
      x[i] = xptilde[i] + uptilde[i]*dto2;
      y[i] = yptilde[i] + vptilde[i]*dto2;
      
    }// end index on particle
   
    LastCommunicate=0;  // this are just the communication inside the inner mover:
    // do not save parts for finer grids at this stage

    // 0: coarse grid
    // 1: refined grid, inner mover
    // 2: refined grid, final position
    if (cartesian_rank_total < vct->getXLEN()*vct->getYLEN())// coarse grid
    {
    	BC_partCommunicate=0;
      }
    else
      {
    BC_partCommunicate=1; // after May4: 1: inner mover 2: final position for both grids
    }
    avail = communicate(vct, grid, BC_partCommunicate);
    if (avail < 0)
      return(-1);
    MPI_Barrier(vct->getCART_COMM());
    while(isMessagingDone(vct) >0)
    {
      avail = communicate(vct, grid, BC_partCommunicate);
      if (avail < 0)
	return(-1);
      MPI_Barrier(vct->getCART_COMM());
    }
    innter++;
  }  // end inner iterations
  for (int i=0; i < nop; i++)
  {
    u[i]= 2.0*uptilde[i] - u[i];
    v[i]= 2.0*vptilde[i] - v[i];
    w[i]= 2.0*wptilde[i] - w[i];
    x[i] = xptilde[i] + uptilde[i]*dt;
    y[i] = yptilde[i] + vptilde[i]*dt;
  }
  // COMMUNICATION
  // AMR, ME: different communicate for finer grids (different particle BC)
  /*if (0 && grid->getLevel()==0)
  {
    cout<< endl << endl <<"R" <<vct->getCartesian_rank_COMMTOTAL() << " np_REPOP_b_BOTTOM " << np_REPOP_b_BOTTOM << " np_REPOP_b_TOP " <<np_REPOP_b_TOP << " np_REPOP_b_LEFT " <<np_REPOP_b_LEFT << " np_REPOP_b_RIGHT " <<np_REPOP_b_RIGHT << "first entry BOTTOM " <<REPOP_b_BOTTOM[0] << "first entry TOP "<<REPOP_b_TOP[0] << "first entry LEFT "<<REPOP_b_LEFT[0] << "first entry RIGHT "<<REPOP_b_RIGHT[0] << endl << endl << endl;
    } */
  
  if (PRACollectionMethod==1)
  {
    LastCommunicate=0; // with this switch, the particle collection for repopulation is done for the entire vector after the mover
  }
  else if (PRACollectionMethod==0) // with this switch, the particle collection for repopulation is dome one by one
  {
    LastCommunicate=1;    //these are the communicate ops after the definitive positions:
    // decide which PRA particles to send to the finer grids now
  }
  else 
  {
    cout << "Problem with the choice of the repopulation method, exiting...\n";
    return -1;
  }
  
  // 0: coarse grid                                                                                                                                   // 1: refined grid, inner mover                                                                                                                  
  // 2: refined grid, final position  
  //MPI_Comm_rank(vct->getCART_COMM_TOTAL(), &cartesian_rank_total);
  if (cartesian_rank_total < vct->getXLEN()*vct->getYLEN())// coarse grid                                                                       
    {
      BC_partCommunicate=0;
    }
  else
    {
      BC_partCommunicate=2;
  }

  avail = communicate(vct, grid, BC_partCommunicate);
  //  cout<<"R" <<vct->getCartesian_rank_COMMTOTAL() << ": (last) communicate, qom " <<qom <<"\n";
  if (avail < 0)
    return(-1);
  MPI_Barrier(vct->getCART_COMM());
  while(isMessagingDone(vct) >0)
  {
    avail = communicate(vct, grid, BC_partCommunicate);
    //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL() << ": (last) communicate, qom " <<qom <<"\n";
    if (avail < 0)
      return(-1);
    MPI_Barrier(vct->getCART_COMM());
  }
  //cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()  << "L" << grid->getLevel() <<"QOM" << qom<<": end communicate after definitive params\n";

  /*if (0)
  {
    cout <<"R" <<vct->getCartesian_rank_COMMTOTAL()<<"qom" <<qom << "Nop After Mover: " << nop <<endl;
  }
  if (0) // debug
  {
    if (grid->getLevel()<  vct->getNgrids()-1)
    {
      // print the number of particles in the repopulation buffers to send to the finer grid
      cout  <<"R" <<vct->getCartesian_rank_COMMTOTAL() <<"qom" <<qom << "nop in repop buffers: BOTTOM: " << np_REPOP_b_BOTTOM << " TOP " << np_REPOP_b_TOP << " LEFT " << np_REPOP_b_LEFT << " RIGHT " << np_REPOP_b_RIGHT << " out of "<< MAX_NP_REPOP_SIZE  << endl;
    }
    }*/
  LastCommunicate=0;  // because Last Communicate is checked also in unbuffer, which is used again when distributing splitted particles
  
  return(0); // exit succcesfully	
}
/** relativistic mover with a Predictor-Corrector scheme */
int Particles2D::mover_relativistic(Grid* grid,VirtualTopology* vct, Field* EMf){
	int innter,temp1,temp2;
	int avail;
	double dto2 = .5*dt, qomdt = qom*dto2, omdtsq, denom, ut, vt, wt, udotb, weight00, weight01,weight10,weight11;
	double Exl=0.0, Eyl=0.0, Ezl=0.0, Bxl=0.0, Byl=0.0, Bzl=0.0;
	double inv_dx = 1.0/dx;
	double inv_dy = 1.0/dy;
	int ix,iy;
	double gamma, gamma0, u02, v2, vdu, cfa, cfb, cfc, delta_rel;
	for (int i=0; i < nop; i++){
		xptilde[i] = x[i];
		yptilde[i] = y[i];
	}
	innter=0;
	while (innter< NiterMover){	  
		// move each particle with new fields
		for (int i=0; i <  nop; i++){
			// interpolation G-->P
			ix = 2 +  int(floor((x[i]-xstart)/dx));
			iy = 2 +  int(floor((y[i]-ystart)/dy));
			weight11 = ((x[i] - grid->getXN(ix-1,iy-1,0))*inv_dx)*((y[i] - grid->getYN(ix-1,iy-1,0))*inv_dy);
			weight10 = ((x[i] - grid->getXN(ix-1,iy,0))*inv_dx)*((grid->getYN(ix-1,iy,0) - y[i])*inv_dy);
			weight01 = ((grid->getXN(ix,iy-1,0) - x[i])*inv_dx)*((y[i] - grid->getYN(ix,iy-1,0))*inv_dy);
			weight00 = ((grid->getXN(ix,iy,0) - x[i])*inv_dx)*((grid->getYN(ix,iy,0) - y[i])*inv_dy);
			Exl = weight00*EMf->getEx(ix-1,iy-1,0) + weight01*EMf->getEx(ix-1,iy,0) + weight10*EMf->getEx(ix,iy-1,0) + weight11*EMf->getEx(ix,iy,0);
			Eyl = weight00*EMf->getEy(ix-1,iy-1,0) + weight01*EMf->getEy(ix-1,iy,0) + weight10*EMf->getEy(ix,iy-1,0) + weight11*EMf->getEy(ix,iy,0);
			Ezl = weight00*EMf->getEz(ix-1,iy-1,0) + weight01*EMf->getEz(ix-1,iy,0) + weight10*EMf->getEz(ix,iy-1,0) + weight11*EMf->getEz(ix,iy,0);
			Bxl = weight00*EMf->getBx(ix-1,iy-1,0) + weight01*EMf->getBx(ix-1,iy,0) + weight10*EMf->getBx(ix,iy-1,0) + weight11*EMf->getBx(ix,iy,0);
			Byl = weight00*EMf->getBy(ix-1,iy-1,0) + weight01*EMf->getBy(ix-1,iy,0) + weight10*EMf->getBy(ix,iy-1,0) + weight11*EMf->getBy(ix,iy,0);
			Bzl = weight00*EMf->getBz(ix-1,iy-1,0) + weight01*EMf->getBz(ix-1,iy,0) + weight10*EMf->getBz(ix,iy-1,0) + weight11*EMf->getBz(ix,iy,0);
			
			// end interpolation
			omdtsq = qomdt*qomdt/c/c*(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
			denom = 1.0/(1.0 + omdtsq);
			// solve the momentum equation
			
			/////////////////////////
			//relativistic part
			gamma0 = 1.0/(sqrt(1.0 - u[i]*u[i] - v[i]*v[i] - w[i]*w[i]));
			
			ut= gamma0*u[i] + qomdt*Exl;
			vt= gamma0*v[i] + qomdt*Eyl;
			wt= gamma0*w[i] + qomdt*Ezl;
			
			//gamma = sqrt(1.0 + ut*ut + vt*vt + wt*wt);
			gamma = gamma0;
                        Bxl /=gamma;
			Byl /=gamma;
			Bzl /=gamma;
			denom /=gamma;
			//////////////
			/////////////
			
			udotb = ut*Bxl + vt*Byl + wt*Bzl;
			// solve the velocity equation 
			uptilde[i] = (ut+qomdt/c*(vt*Bzl -wt*Byl + qomdt/c*udotb*Bxl))*denom; 
			vptilde[i] = (vt+qomdt/c*(wt*Bxl -ut*Bzl + qomdt/c*udotb*Byl))*denom;
			wptilde[i] = (wt+qomdt/c*(ut*Byl -vt*Bxl + qomdt/c*udotb*Bzl))*denom;
			// update position
			x[i] = xptilde[i] + uptilde[i]*dto2;
		    y[i] = yptilde[i] + vptilde[i]*dto2;
			
		}
		// Communication
		avail = communicate(vct);
		if (avail < 0)
			return(-1);
		MPI_Barrier(vct->getCART_COMM());
		while(isMessagingDone(vct) >0){
			// communication particles
			avail = communicate(vct);
			if (avail < 0)
				return(-1);
			MPI_Barrier(vct->getCART_COMM());
		}
		innter++;
	}
	// Advance solution to final values
	for (int i=0; i < nop; i++){
		 // relativistic velocity update
		 gamma0 = 1.0/(sqrt(1.0 - u[i]*u[i] - v[i]*v[i] - w[i]*w[i]));
		 ut = u[i]*gamma0;
		 vt = v[i]*gamma0;
		 wt = w[i]*gamma0;
		 u02 = ut*ut + vt*vt + wt*wt;
		 v2 = uptilde[i]*uptilde[i] + vptilde[i]*vptilde[i] + wptilde[i]*wptilde[i];
		 vdu = ut*uptilde[i] + vt*vptilde[i] + wt*wptilde[i];
		 
		 cfa = 1.0 -v2;
		 cfb =-2.0*(-vdu+gamma0*v2);
		 cfc =-1.0-gamma0*gamma0*v2+2.0*gamma0*vdu - u02; 
		 
		 delta_rel=cfb*cfb -4.0*cfa*cfc;
		 // update velocity
		 if (delta_rel < 0.0){
		     cout << "Relativity violated: gamma0=" << gamma0 << ",  v2=" << v2;
		     u[i] = (2.0*gamma)*uptilde[i] - u[i]*gamma0;
			 v[i] = (2.0*gamma)*vptilde[i] - v[i]*gamma0;
			 w[i] = (2.0*gamma)*wptilde[i] - w[i]*gamma0;
		  } else {
		     gamma = (-cfb+sqrt(delta_rel))/2.0/cfa;
			 u[i] = (gamma + gamma0)*uptilde[i] - ut;
			 v[i] = (gamma + gamma0)*vptilde[i] - vt;
			 w[i] = (gamma + gamma0)*wptilde[i] - wt;
			 u[i] /=gamma;
			 v[i] /=gamma;
			 w[i] /=gamma;
		 }
		 // update position
		 x[i] = xptilde[i] + uptilde[i]*dt;
		 y[i] = yptilde[i] + vptilde[i]*dt;
	}
	// COMMUNICATION
	avail = communicate(vct);
	if (avail < 0)
		return(-1);
	MPI_Barrier(vct->getCART_COMM());
	while(isMessagingDone(vct) >0){
		// COMMUNICATION
		avail = communicate(vct);
		if (avail < 0)
	        return(-1);
		MPI_Barrier(vct->getCART_COMM());
	}
	
	return(0); // exit succcesfully	
}

/** mover with a Predictor-Corrector scheme wit dipole and earth */
int Particles2D::moverDipole(double B_dipole, double x_dipole, double y_dipole, double radius_Earth, Grid* grid,VirtualTopology* vct, Field* EMf){
	int innter,temp1,temp2;
	int avail;
	double distance_sq, x_displ, y_displ;
	double dto2 = .5*dt, qomdt = qom*dto2, omdtsq, denom, ut, vt, wt, udotb, weight00, weight01,weight10,weight11;
	double Exl=0.0, Eyl=0.0, Ezl=0.0, Bxl=0.0, Byl=0.0, Bzl=0.0;
	double inv_dx = 1.0/dx;
	double inv_dy = 1.0/dy;
	int ix,iy;
	for (int i=0; i < nop; i++){
		xptilde[i] = x[i];
		yptilde[i] = y[i];
	}
	innter=0;
	while (innter< NiterMover){	  
		// move each particle with new fields
		for (int i=0; i <  nop; i++){
			
			// interpolation G-->P
			ix = 2 +  int(floor((x[i]-xstart)/dx));
			iy = 2 +  int(floor((y[i]-ystart)/dy));
			weight11 = ((x[i] - grid->getXN(ix-1,iy-1,0))*inv_dx)*((y[i] - grid->getYN(ix-1,iy-1,0))*inv_dy);
			weight10 = ((x[i] - grid->getXN(ix-1,iy,0))*inv_dx)*((grid->getYN(ix-1,iy,0) - y[i])*inv_dy);
			weight01 = ((grid->getXN(ix,iy-1,0) - x[i])*inv_dx)*((y[i] - grid->getYN(ix,iy-1,0))*inv_dy);
			weight00 = ((grid->getXN(ix,iy,0) - x[i])*inv_dx)*((grid->getYN(ix,iy,0) - y[i])*inv_dy);
			// el field
			Exl = weight00*EMf->getEx(ix-1,iy-1,0) + weight01*EMf->getEx(ix-1,iy,0) + weight10*EMf->getEx(ix,iy-1,0) + weight11*EMf->getEx(ix,iy,0);
			Eyl = weight00*EMf->getEy(ix-1,iy-1,0) + weight01*EMf->getEy(ix-1,iy,0) + weight10*EMf->getEy(ix,iy-1,0) + weight11*EMf->getEy(ix,iy,0);
			Ezl = weight00*EMf->getEz(ix-1,iy-1,0) + weight01*EMf->getEz(ix-1,iy,0) + weight10*EMf->getEz(ix,iy-1,0) + weight11*EMf->getEz(ix,iy,0);
			// magn field
			Bxl = weight00*EMf->getBx(ix-1,iy-1,0) + weight01*EMf->getBx(ix-1,iy,0) + weight10*EMf->getBx(ix,iy-1,0) + weight11*EMf->getBx(ix,iy,0);
			Byl = weight00*EMf->getBy(ix-1,iy-1,0) + weight01*EMf->getBy(ix-1,iy,0) + weight10*EMf->getBy(ix,iy-1,0) + weight11*EMf->getBy(ix,iy,0);
			Bzl = weight00*EMf->getBz(ix-1,iy-1,0) + weight01*EMf->getBz(ix-1,iy,0) + weight10*EMf->getBz(ix,iy-1,0) + weight11*EMf->getBz(ix,iy,0);
			// add dipole contribution only outside the sphere
			distance_sq = (x[i] - x_dipole)*(x[i] - x_dipole) + (y[i] - y_dipole)*(y[i] - y_dipole);
			if (distance_sq > radius_Earth*radius_Earth){ 
				x_displ = x[i]-  x_dipole;
				y_displ = y[i] - y_dipole;
				Bxl += B_dipole*(x_displ*y_displ)/pow(x_displ*x_displ + y_displ*y_displ,2.5);
				Byl += B_dipole*(-(1/3)*(x_displ*x_displ+ y_displ*y_displ))/pow(x_displ*x_displ + y_displ*y_displ,2.5);
			}
			// end interpolation
			omdtsq = qomdt*qomdt/c/c*(Bxl*Bxl+Byl*Byl+Bzl*Bzl);
			denom = 1.0/(1.0 + omdtsq);
			// solve the position equation
			ut= u[i] + qomdt*Exl;
			vt= v[i] + qomdt*Eyl;
			wt= w[i] + qomdt*Ezl;
			udotb = ut*Bxl + vt*Byl + wt*Bzl;
			// solve the velocity equation 
			uptilde[i] = (ut+qomdt/c*(vt*Bzl -wt*Byl + qomdt/c*udotb*Bxl))*denom; 
			vptilde[i] = (vt+qomdt/c*(wt*Bxl -ut*Bzl + qomdt/c*udotb*Byl))*denom;
			wptilde[i] = (wt+qomdt/c*(ut*Byl -vt*Bxl + qomdt/c*udotb*Bzl))*denom;
			// update position
			x[i] = xptilde[i] + uptilde[i]*dt;
			y[i] = yptilde[i] + vptilde[i]*dt;
			
		}
		// Communication
		avail = communicate(vct);
		if (avail < 0)
			return(-1);
		MPI_Barrier(vct->getCART_COMM());
		while(isMessagingDone(vct) >0){
			// communication particles
			avail = communicate(vct);
			if (avail < 0)
				return(-1);
			MPI_Barrier(vct->getCART_COMM());
		}
		innter++;
	}
	for (int i=0; i < nop; i++){
        u[i]= 2.0*uptilde[i] - u[i];
        v[i]= 2.0*vptilde[i] - v[i];
        w[i]= 2.0*wptilde[i] - w[i];
        x[i] = xptilde[i] + uptilde[i]*dt;
        y[i] = yptilde[i] + vptilde[i]*dt;
	}
	
	// COMMUNICATION
	avail = communicate(vct);
	if (avail < 0)
		return(-1);
	MPI_Barrier(vct->getCART_COMM());
	while(isMessagingDone(vct) >0){
		// COMMUNICATION
		avail = communicate(vct);
		if (avail < 0)
	        return(-1);
		MPI_Barrier(vct->getCART_COMM());
	}
	
	return(0); // exit succcesfully	
}
/** Injection of a beam at a certain position and velocity*/
int Particles2D::inject_beam(double x_source, double y_source, double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, double radius_beam, int np_inj, double qom_beam, double rho_init, VirtualTopology* vct){
	double rho_1, rho_2, rho_3, r, theta, harvest, prob;
	double const_density = (rho_init/4*M_PI); 
	// this ensures that all the particles when generated
	// are all on the same processor
	// this is done to ensure that the number of particles injected is np_inj
	if (  x_source > xstart && x_source < xend && y_source > ystart && y_source < yend){	 
		for (int i=0; i < np_inj;i++){
			// sample three random numbers for the position
			rho_1 =  rand()/(double)RAND_MAX;
			rho_2 =  rand()/(double)RAND_MAX;
			rho_3 =  rand()/(double)RAND_MAX;
			r = radius_beam*sqrt(rho_1);
			theta = 2*M_PI*rho_2;
			if (u_beam != 0.0){  // beam moving in X-direction
				y[nop] = y_source + r*sin(theta);
				q[nop] = (qom/fabs(qom))*2*radius_beam*const_density*dt/np_inj;
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				u[nop] = generate_v(u_beam,uth_beam);
				x[nop] = x_source + u[nop]*dt*rho_3;
				// v
				v[nop] = v_beam + vth_beam*prob*sin(theta);
				// w
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				w[nop] = w_beam + wth_beam*prob*cos(theta);
				if (TrackParticleID)	      
					ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
			} else if (v_beam != 0.0){  // beam moving in Y-direction	   
				
				x[nop] = x_source + r*sin(theta);
				q[nop] = (qom_beam/fabs(qom_beam))*2*radius_beam*const_density*dt/np_inj;
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				u[nop] = u_beam + uth_beam*prob*cos(theta);
				// v
				v[nop] =  generate_v(v_beam,vth_beam);
				y[nop] = y_source + v[nop]*dt*rho_3;
				// w
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				w[nop] = w_beam + wth_beam*prob*cos(theta);
				if (TrackParticleID)	      
					ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
				
			}
			// update the number of particles
			nop++; 
			if (nop > npmax*.95)
				return(-1);
		}
	}
	// since the the particles are all generated on the same processor
	// they need to be placed on the right processor
	// COMMUNICATION
	int  avail = communicate(vct);
	if (avail < 0)
	    return(-1);
	MPI_Barrier(vct->getCART_COMM());
	while(isMessagingDone(vct) >0){
		// communication particles
		avail = communicate(vct);
		if (avail < 0)
	        return(-1);
		MPI_Barrier(vct->getCART_COMM());
	} 
	
	
	return(0); // exit succcesfully
	
}
/** generate the vleocity of the particle with the rejection method */
double Particles2D::generate_v(double v00, double vvth){
	double sotto = v00*vvth*sqrt(M_PI)/2.0*(1-erf(-v00/vvth)) + vvth*vvth/2.0*exp(-v00*v00/(vvth*vvth));
	double vmax = .5*(v00+sqrt(2.0*vvth*vvth+v00*v00));
	double fmax = vmax*exp(-(vmax-v00)*(vmax-v00)/(vvth*vvth))/sotto;
	double max_v;
	if (v00 > 0.0)
		max_v = v00;
	else
		max_v = 0.0;
	double vv = (2.0*vvth + max_v)*rand()/(double)RAND_MAX;
	double r = rand()/(double)RAND_MAX;
	double ff = vv*exp(-(vv - v00)*(vv - v00)/(vvth*vvth))/sotto/fmax;
	while(r > ff){
		vv = (2.0*vvth + max_v)*rand()/(double)RAND_MAX;
		r = rand()/(double)RAND_MAX;
		ff = vv*exp(-(vv - v00)*(vv - v00)/(vvth*vvth))/sotto/fmax;
	}
	return(vv);
}
/** inject particles from the X leftwall */
int Particles2D::inject_Xleft_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind){
	double harvest, prob, theta;
	int nyc = grid->getNYC();
	if(vct->getXleft_neighbor()==MPI_PROC_NULL){ // only these processors injects
		double flxrnd=exp(-u_beam*u_beam/(uth_beam*uth_beam))*uth_beam/2.0/sqrt(M_PI);
		double flxdir=.5*u_beam*(1.+erf(u_beam/uth_beam));
		double fluxin= (flxrnd+flxdir);
		double* density_profile = new double[nyc];
		switch(dens_profile_kind){
			case 0:   // constant density
				for (int i=0; i < nyc; i++)
					density_profile[i] = rho_init/(4*M_PI); 
				break;
			case 1:    // newton challange density profile
				
				break;
		}
		for (int i=1; i < nyc-1; i++){
			for (int ii=0; ii < npcely; ii++){
				u[nop] = generate_v(u_beam,uth_beam);
				x[nop] = u[nop]*dt;
				q[nop] = (qom/fabs(qom))*fluxin*(density_profile[i]*dt/npcely)*dy;
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				v[nop] = v_beam + vth_beam*prob*sin(theta);
				y[nop] = (ii + .5)*(dy/npcely) + grid->getYN(1,i,0);
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				w[nop] = w_beam + wth_beam*prob*cos(theta);
				if (TrackParticleID)	      
					ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
				// update the number of particles
				nop++;
				if (nop > npmax*.95)
					return(-1);
				
			}
		}
		
		
		
		delete[] density_profile;
	}
	int  avail = communicate(vct);
	if (avail < 0)
	    return(-1);
	MPI_Barrier(vct->getCART_COMM());
	while(isMessagingDone(vct) >0){
		// communication particles
		avail = communicate(vct);
		if (avail < 0)
	        return(-1);
		MPI_Barrier(vct->getCART_COMM());
	} 
	
	return(0);
}
/** inject particles from the X right wall */
int Particles2D::inject_Xright_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind){
	double harvest, prob, theta;
	int nyc = grid->getNYC();
	if(vct->getXright_neighbor()==MPI_PROC_NULL){ // only these processors injects
		double flxrnd=exp(-u_beam*u_beam/(uth_beam*uth_beam))*uth_beam/2.0/sqrt(M_PI);
		double flxdir=.5*u_beam*(1.+erf(u_beam/uth_beam));
		double fluxin= (flxrnd+flxdir);
		double* density_profile = new double[nyc];
		switch(dens_profile_kind){
			case 0:   // constant density
				for (int i=0; i < nyc; i++)
					density_profile[i] = rho_init/(4*M_PI); 
				break;
			case 1:    // newton challenge density profile ?
				
				break;
		}
		for (int i=1; i < nyc-1; i++){
			for (int ii=0; ii < npcely; ii++){
				u[nop] = - generate_v(u_beam,uth_beam);
				x[nop] = Lx + u[nop]*dt;
				q[nop] = (qom/fabs(qom))*fluxin*(density_profile[i]*dt/npcely)*dy;
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				v[nop] = v_beam + vth_beam*prob*sin(theta);
				y[nop] = (ii + .5)*(dy/npcely) + grid->getYN(nxn-2,i,0);
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				w[nop] = w_beam + wth_beam*prob*cos(theta);
				if (TrackParticleID)	      
					ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
				// update the number of particles
				nop++;
				if (nop > npmax*.95)
					return(-1);
				
			}
		}
		
		
		
		delete[] density_profile;
	}
	int  avail = communicate(vct);
	if (avail < 0)
	    return(-1);
	MPI_Barrier(vct->getCART_COMM());
	while(isMessagingDone(vct) >0){
		// communication particles
		avail = communicate(vct);
		if (avail < 0)
	        return(-1);
		MPI_Barrier(vct->getCART_COMM());
	} 
	
	return(0);
}
/** inject particles from the Y left wall */
int Particles2D::inject_Yleft_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind){
	double harvest, prob, theta;
	int nxc = grid->getNXC();
	if(vct->getYleft_neighbor()==MPI_PROC_NULL){ // only these processors injects
		double* density_profile = new double[nxc];
		double flxrnd=exp(-v_beam*v_beam/(vth_beam*vth_beam))*vth_beam/2.0/sqrt(M_PI);
		double flxdir=.5*v_beam*(1.+erf(v_beam/vth_beam));
		double fluxin= (flxrnd+flxdir);
		switch(dens_profile_kind){
			case 0:   // constant density
				for (int i=0; i < nxc; i++)
					density_profile[i] = rho_init/(4*M_PI); 
				break;
			case 1:    // newton challange density profile
				
				break;
		}
		for (int i=1; i < nxc-1; i++){
			for (int ii=0; ii < npcelx; ii++){
				v[nop] = generate_v(v_beam,vth_beam);
				y[nop] = v[nop]*dt;
				q[nop] = (qom/fabs(qom))*fluxin*(density_profile[i]*dt/npcelx)*dx;
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				u[nop] = u_beam + uth_beam*prob*sin(theta);
				x[nop] = (ii + .5)*(dx/npcelx) + grid->getXN(i,1,0);
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				w[nop] = w_beam + wth_beam*prob*cos(theta);
				if (TrackParticleID)	      
					ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
				// update the number of particles
				nop++;
				if (nop > npmax*.95)
					return(-1);
				
			}
		}
		delete[] density_profile;
	}
	int  avail = communicate(vct);
	if (avail < 0)
	    return(-1);
	MPI_Barrier(vct->getCART_COMM());
	while(isMessagingDone(vct) >0){
		// communication particles
		avail = communicate(vct);
		if (avail < 0)
	        return(-1);
		MPI_Barrier(vct->getCART_COMM());
	} 
	
	return(0);
}
/** inject particles from the Y left wall */
int Particles2D::inject_Yright_Wall(double u_beam, double v_beam, double w_beam, double uth_beam, double vth_beam, double wth_beam, Grid *grid,VirtualTopology* vct, double rho_init, int dens_profile_kind){
	double harvest, prob, theta;
	int nxc = grid->getNXC();
	if(vct->getYright_neighbor()==MPI_PROC_NULL){ // only these processors injects
		double flxrnd=exp(-v_beam*v_beam/(vth_beam*vth_beam))*vth_beam/2.0/sqrt(M_PI);
		double flxdir=.5*v_beam*(1.+erf(v_beam/vth_beam));
		double fluxin= (flxrnd+flxdir);
		double* density_profile = new double[nxc];
		switch(dens_profile_kind){
			case 0:   // constant density
				for (int i=0; i < nxc; i++)
					density_profile[i] = rho_init/(4*M_PI); 
				break;
			case 1:    // newton challenge density profile
				
				break;
		}
		for (int i=1; i < nxc-1; i++){
			for (int ii=0; ii < npcelx; ii++){
				v[nop] = generate_v(v_beam,vth_beam);
				y[nop] = Ly + v[nop]*dt;
				q[nop] = (qom/fabs(qom))*fluxin*(density_profile[i]*dt/npcelx)*dx;
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				u[nop] = u_beam + uth_beam*prob*sin(theta);
				x[nop] = (ii + .5)*(dx/npcelx) + grid->getXN(i,nyn-2,0);
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				w[nop] = w_beam + wth_beam*prob*cos(theta);
				if (TrackParticleID)	      
					ParticleID[nop]= nop*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
				// update the number of particles
				nop++;
				if (nop > npmax*.95)
					return(-1);
				
			}
		}
		delete[] density_profile;
	}
	int  avail = communicate(vct);
	if (avail < 0)
	    return(-1);
	MPI_Barrier(vct->getCART_COMM());
	while(isMessagingDone(vct) >0){
		// communication particles
		avail = communicate(vct);
		if (avail < 0)
	        return(-1);
		MPI_Barrier(vct->getCART_COMM());
	} 
	return(0);
}
/** interpolation Particle->Grid only for pressure tensor */    
void Particles2D::interpP2G_onlyP(Field* EMf, Grid *grid, VirtualTopology* vct){
	double*** weight = newArr3(double,2,2,1);
	double*** temp = newArr3(double,2,2,1);
	int ix,iy, temp2,temp1;
	for (register int i=0; i < nop; i++){
		ix = 2 +  int(floor((x[i]-xstart)/dx))  ;
		iy = 2 +  int(floor((y[i]-ystart)/dy))  ;
		weight[1][1][0] = ((x[i] - grid->getXN(ix-1,iy-1,0))/dx)*((y[i] - grid->getYN(ix-1,iy-1,0))/dy);
		weight[1][0][0] = ((x[i] - grid->getXN(ix-1,iy,0))/dx)*((grid->getYN(ix-1,iy,0) - y[i])/dy);
		weight[0][1][0] = ((grid->getXN(ix,iy-1,0) - x[i])/dx)*((y[i] - grid->getYN(ix,iy-1,0))/dy);
		weight[0][0][0] = ((grid->getXN(ix,iy,0) - x[i])/dx)*((grid->getYN(ix,iy,0) - y[i])/dy);
		
		scale(weight,q[i],2,2);
		//Pxx
		eqValue(0.0,temp,2,2);
		addscale(u[i]*u[i],temp,weight,2,2);
		EMf->addPxx(temp,ix,iy,0,ns);
		// Pxy
		eqValue(0.0,temp,2,2);
		addscale(u[i]*v[i],temp,weight,2,2);
		EMf->addPxy(temp,ix,iy,0,ns);
		// Pxz
		eqValue(0.0,temp,2,2);
		addscale(u[i]*w[i],temp,weight,2,2);
		EMf->addPxz(temp,ix,iy,0,ns);
		// Pyy
		eqValue(0.0,temp,2,2);
		addscale(v[i]*v[i],temp,weight,2,2);
		EMf->addPyy(temp,ix,iy,0,ns);
		// Pyz
		eqValue(0.0,temp,2,2);
		addscale(v[i]*w[i],temp,weight,2,2);
		EMf->addPyz(temp,ix,iy,0,ns);
		// Pzz
		eqValue(0.0,temp,2,2);
		addscale(w[i]*w[i],temp,weight,2,2);
		EMf->addPzz(temp,ix,iy,0,ns);
	}
	delArr3(weight,2,2);
	delArr3(temp,2,2);
}
/** interpolation Particle->Grid only charge density, current */    
void Particles2D::interpP2G_notP(Field* EMf, Grid *grid, VirtualTopology* vct){
	double*** weight = newArr3(double,2,2,1);
	double*** temp = newArr3(double,2,2,1);
	int ix,iy, temp2,temp1;
	for (register int i=0; i < nop; i++){
		ix = 2 +  int(floor((x[i]-xstart)/dx));
		iy = 2 +  int(floor((y[i]-ystart)/dy));
		weight[1][1][0] = ((x[i] - grid->getXN(ix-1,iy-1,0))/dx)*((y[i] - grid->getYN(ix-1,iy-1,0))/dy);
		weight[1][0][0] = ((x[i] - grid->getXN(ix-1,iy,0))/dx)*((grid->getYN(ix-1,iy,0) - y[i])/dy);
		weight[0][1][0] = ((grid->getXN(ix,iy-1,0) - x[i])/dx)*((y[i] - grid->getYN(ix,iy-1,0))/dy);
		weight[0][0][0] = ((grid->getXN(ix,iy,0) - x[i])/dx)*((grid->getYN(ix,iy,0) - y[i])/dy);
		
		scale(weight,q[i],2,2);
		
		// rho
		EMf->addRho(weight,ix,iy,0,ns);
		// Jx
		eqValue(0.0,temp,2,2);
		addscale(u[i],temp,weight,2,2);
		EMf->addJx(temp,ix,iy,0,ns);
		// Jy
		eqValue(0.0,temp,2,2);
		addscale(v[i],temp,weight,2,2);
		EMf->addJy(temp,ix,iy,0,ns);
		// Jz
		eqValue(0.0,temp,2,2);
		addscale(w[i],temp,weight,2,2);
		EMf->addJz(temp,ix,iy,0,ns);
		
	}
	// communicate contribution from ghost cells     
	EMf->communicateGhostP2G(ns,0,0,0,0,vct);
	delArr3(weight,2,2);
	delArr3(temp,2,2);
	
}
/** apply a linear perturbation to particle distribution */
void Particles2D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r,     double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, Grid* grid,Field* EMf,VirtualTopology* vct){
	
	double value1=0.0,value2=0.0,max_value=0.0,min_value=0.0,phi,n;
	int counter=0, total_generated=0;
	bool rejected;
	double harvest, prob, theta;
	// rescaling of amplitudes according to deltaBoB //
	double alpha;
	double integral=0.0;
	
	alpha=deltaBoB*sqrt(EMf->getBx(1,1,0)*EMf->getBx(1,1,0)+EMf->getBy(1,1,0)*EMf->getBy(1,1,0)+EMf->getBz(1,1,0)*EMf->getBz(1,1,0))/sqrt(Bx_mod*Bx_mod+By_mod*By_mod+Bz_mod*Bz_mod);
	
	Ex_mod *= alpha;
	Ey_mod *= alpha;
	Ez_mod *= alpha;
	Bx_mod *= alpha;
	By_mod *= alpha;
	Bz_mod *= alpha;
	
	
	
	// find the maximum value of f=1+delta_f/f0
	for (register double vpar=-2*uth; vpar<=2*uth; vpar += 0.0005)
		for (register double vperp=1e-10; vperp<=2*vth; vperp += 0.0005)
			for (register double X=xstart; X<=xend; X += 2*grid->getDX())
				for (register double Y=ystart; Y<=yend; Y += 2*grid->getDY()){
					value1=1+delta_f(vpar,vperp,0.0, X, Y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle,  EMf)/f0(vpar,vperp);
					
					if (value1>max_value)
						max_value=value1;
					
					
				}
					
					
					
					max_value *=3.2;phi=1.48409;n=2.948687; // security factor...
					if (ns==1){
						max_value *=3.0;phi=-1.65858;n=2.917946;} // security factor...
					cout<<"max-value="<<max_value<<" min-value="<<min_value<<endl;
					
					/* initialize random generator */
					srand (vct->getCartesian_rank()+2);
					
					for (int i=1; i< grid->getNXC()-1;i++)
						for (int j=1; j< grid->getNYC()-1;j++)
							for (int ii=0; ii < npcelx+round(2*n*(cos(2*M_PI*0.4125*grid->getXN(i,j,0)+phi))); ii++)
								for (int jj=0; jj < npcely; jj++){
									x[counter] = (ii + .5)*(dx/(npcelx+round(2*n*(cos(2*M_PI*0.4125*grid->getXN(i,j,0)+phi))))) + grid->getXN(i,j,0);   
									y[counter] = (jj + .5)*(dy/npcely) + grid->getYN(i,j,0);
									q[counter] =  (qom/fabs(qom))*((0.19635)/npcel)*(1.0/invVOL);
									
									// apply rejection method in velocity space
									rejected=true;
									while (rejected){
										total_generated ++;
										harvest =   rand()/(double)RAND_MAX;
										prob  = sqrt(-2.0*log(1.0-.999999*harvest));
										harvest =   rand()/(double)RAND_MAX;
										theta = 2.0*M_PI*harvest;
										// u
										u[counter] = u0 + uth*prob*cos(theta);
										// v
										v[counter] = v0 + vth*prob*sin(theta);
										// w
										harvest =   rand()/(double)RAND_MAX;
										prob  = sqrt(-2.0*log(1.0-.999999*harvest));
										harvest =   rand()/(double)RAND_MAX;
										theta = 2.0*M_PI*harvest;
										w[counter] = w0 + wth*prob*cos(theta);
										
										// test: if rand < (1+delta_f/f0)/max_value --> accepted
										if ( rand()/(double)RAND_MAX <= (1+ delta_f(u[counter],v[counter],w[counter],x[counter], y[counter], kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase,Ez_mod, Ez_phase, angle,  EMf)/f0(u[counter],sqrt(v[counter]*v[counter]+w[counter]*w[counter])))/max_value)
											rejected=false;
										
									}
									if (TrackParticleID)	      
										ParticleID[counter]= counter*(unsigned long)pow(10.0,BirthRank[1])+BirthRank[0];
									counter++ ;
								}
									nop=counter+1;
					//		     if (vct->getCartesian_rank()==0)
					cout<<"Rejection method: "<<(counter+1)/ double(total_generated) * 100<< " % of particles are accepted for species "<<ns<<" counter="<<counter<<endl; 
}

/** Linear delta f for bi-maxwellian plasma */
double Particles2D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase,double Ez_mod, double Ez_phase, double theta, Field* EMf){
	const complex<double> I(0.0,1.0);
	const double vperp = sqrt(v*v+w*w);
	const double vpar = u;
	const double kpar= kx;
	double kperp;
	if (ky==0.0) // because this formula is not valid for exactly parallel
		kperp = 1e-9;
	else kperp = ky;
	const double om_c = qom/c*sqrt(EMf->getBx(1,1,0)*EMf->getBx(1,1,0)+EMf->getBy(1,1,0)*EMf->getBy(1,1,0))/2/M_PI;
	const double phi = atan2(w,v);
	const double lambda = kperp*vperp/om_c;
	const complex<double> omega (omega_re,omega_i);
	
	const int lmax=5; // sum from -lmax to lmax
	
	double bessel_Jn_array[ lmax + 2 ];
	double bessel_Jn_prime_array[ lmax +1 ];
	complex <double> a1[2*lmax+1], a2[2*lmax+1],a3[2*lmax+1] ;
	complex <double> factor, deltaf; 
	
	// rotation of x,y
	double temp;
	temp=x;
	x=x*cos(theta)-y*sin(theta);
	y=temp*sin(theta)+y*cos(theta);
	
	
	/** for compilation issues comment this part: PUT in the math stuff */
	//calc_bessel_Jn_seq(lambda, lmax, bessel_Jn_array, bessel_Jn_prime_array);
	factor = (kpar*vperp/omega*df0_dvpar(vpar,vperp) + (1.0-(kpar*vpar/omega))*df0_dvperp(vpar,vperp) );
	for (register int l=-lmax; l<0; l++){ // negative index
		a1[l+lmax] = factor/lambda*pow(-1.0,-l)*bessel_Jn_array[-l];
		a1[l+lmax] *= (double) l;
		a2[l+lmax] = factor* I* 0.5*pow(-1.0,-l)*(bessel_Jn_array[-l-1]-bessel_Jn_array[-l+1]);
		a3[l+lmax] = kperp/omega*(vpar*df0_dvperp(vpar,vperp) - vperp*df0_dvpar(vpar,vperp))/lambda*pow(-1.0,-l)*bessel_Jn_array[-l];
		a3[l+lmax] *= (double) l;
		a3[l+lmax] += df0_dvpar(vpar,vperp)*pow(-1.0,-l)*bessel_Jn_array[-l];
	}
	
	for (register int l=0; l<lmax+1; l++){ //positive index
		a1[l+lmax] = factor/lambda*bessel_Jn_array[l];
		a1[l+lmax] *= (double) l;
		a2[l+lmax] = factor* I* bessel_Jn_prime_array[l];
		a3[l+lmax] = kperp/omega*(vpar*df0_dvperp(vpar,vperp) - vperp*df0_dvpar(vpar,vperp))/lambda*bessel_Jn_array[l];
		a3[l+lmax] *= (double) l;
		a3[l+lmax] += df0_dvpar(vpar,vperp)*bessel_Jn_array[l];
	}
	
	deltaf=(0.0,0.0);
	for (register int l=-lmax; l<lmax+1; l++){
		deltaf += (a3[l+lmax]*Ex_mod*exp(I*Ex_phase) + a1[l+lmax]*Ey_mod*exp(I*Ey_phase) + a2[l+lmax]*Ez_mod*exp(I*Ez_phase))/(kpar*vpar+l*om_c-omega)*exp(-I*phi*(double)l);}
	deltaf *= I*qom*exp(I*lambda*sin(phi))*exp(I*(2*M_PI*kx*x+2*M_PI*ky*y));
	
	return(real(deltaf));
}

double Particles2D::df0_dvpar(double vpar,double vperp){
	double result;
	result = -2*(vpar-u0)/uth/uth*exp(-(vperp*vperp/vth/vth + (vpar-u0) * (vpar-u0) /uth/uth));
	result *= 3.92e6/pow(M_PI,3/2)/vth/vth/uth;
	return(result);
}

double Particles2D::df0_dvperp(double vpar,double vperp){
	double result;
	result = -2*(vperp)/vth/vth*exp(-(vperp*vperp/vth/vth + (vpar-u0) * (vpar-u0) /uth/uth));
	result *= 3.92e6/pow(M_PI,3/2)/vth/vth/uth;
	return(result);
}

double Particles2D::f0(double vpar,double vperp){
	double result;
	result = exp(-(vperp*vperp/vth/vth + (vpar-u0) * (vpar-u0) /uth/uth));
	result *= 3.92e6/pow(M_PI,3/2)/vth/vth/uth;
	return(result);
}

void Particles2D::RotatePlaneXY(double theta){
	double temp,temp2;
	for (register int s=0; s < nop; s++){
		temp=u[s];temp2=v[s];
		u[s]=temp*cos(theta)+v[s]*sin(theta);
		v[s]=-temp*sin(theta)+temp2*cos(theta);
	}
}
/** get the total number of particles  in the system */
int Particles2D::getNPsystem(VirtualTopology* vct){
	int result = reduceNumberParticles(vct->getCART_COMM(),nop);
	return(result);
}
/** get the total number of particles eliminated at the boundary */
int Particles2D::getNPdeletedBoundary(VirtualTopology* vct){
	int result = reduceNumberParticles(vct->getCART_COMM(),npDeletedBoundary);
	return(result);
}
/** get the total number of particles eliminated from the dipole */
int Particles2D::getNPdeletedDipole(VirtualTopology* vct){
	int result = reduceNumberParticles(vct->getCART_COMM(),npDeletedDipole);
	return(result);
}

// ME: repopulate the ghost cells of the coarse level
int Particles2D::RepopulateCoarseLevel(Grid *grid, VirtualTopology* vct)
{
  
}




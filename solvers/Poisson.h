/***************************************************************************
                          Poisson.h  -  Poisson solver
                             -------------------
    begin             : Thu Nov 09 2006
    developers        : Enrico Camporeale
 ***************************************************************************/


#ifndef Poisson_H
#define Poisson_H

#include <iostream>
#include <fftw3.h>
#include <math.h>
#include "../mathlib/Basic.h"
#include "../utility/TransArraySpace.h"

typedef  void (Field::*FIELD_IMAGE)(double*, double*, Grid*, VirtualTopology*);
typedef  void (*GENERIC_IMAGE)(double*, double*, Grid*, VirtualTopology*);
using std::cout;
using std::cerr;
using std::endl;

#define PI 4*atan(1) 

/** Poisson solver
*
* @date Thu Nov 9 2006
* @author Enrico Camporeale
* @version 1.0
*/

void applyBCpoissonRight(double** vector,double** x, int bcPHI ,int DIM,int nx, int ny, double dx, double dy);
void  applyBCpoissonLeft(double** vector,double** x, int bcPHI ,int DIM,int nx, int ny, double dx, double dy);

// attention: bPoisson e x are with ghosts !!!			
inline void PoissonSolver(int nx, int ny, double dx, double dy, double** x, double** bPoisson,  double tol, int maxiter, int bcPHIfaceXright, int bcPHIfaceYright,int bcPHIfaceXleft, int bcPHIfaceYleft, FIELD_IMAGE FunctionImage, Grid *grid, VirtualTopology *vct, Field *field){

// b_fftw,phi defined only on active grid
double* b_fftw = (double*) fftw_malloc( sizeof (double) * (nx-2)*(ny-2));	
double* phi_fftw = (double*) fftw_malloc(sizeof (double) * (nx-2)*(ny-2));	
double *r   = new double[(nx-2)*(ny-2)];
double *b_solver = new double[(nx-2)*(ny-2)];
double *image = new double[(nx-2)*(ny-2)];
double **im  = newArr(double,nx,ny);

double **b = newArr(double,nx-2,ny-2);
double **b_initial = newArr(double,nx-2,ny-2);
double **laplacian = newArr(double,nx,ny);
double **gradx = newArr(double,nx,ny);

double **prova=newArr(double,nx,ny);
double error, normb;
int counter;

bool PoissonVerbose = true;

//eqValue(0.0,x,nx,ny);
//eqValue(0.0,b_solver,(nx-2)*(ny-2));

phys2solver(b_solver, bPoisson, nx, ny); // transform the known term in solver format
getRidGhost(b,bPoisson, nx, ny);
eq(b_initial,b,nx-2,ny-2);
counter=0;error=1;

 fftw_plan plan_b, plan_phi,plan_b_inv;
 plan_b   = fftw_plan_r2r_2d(nx-2, ny-2, b_fftw, b_fftw, FFTW_RODFT00, FFTW_RODFT00,0);
 plan_phi = fftw_plan_r2r_2d(nx-2, ny-2, phi_fftw, phi_fftw, FFTW_RODFT00, FFTW_RODFT00,0);

 normb = dotP(vct->getCART_COMM(), b_solver,b_solver,(nx-2)*(ny-2));
cout<<"normb = "<<normb<<endl;
 
// start iteration...
//while (counter<maxiter && error/normb> tol){
	while (counter<maxiter ){
	
// apply BC 
if (vct->getXright_neighbor() != MPI_PROC_NULL){
	for (int j=0; j<ny-2; j++) 	
		b[nx-3][j] -= x[nx-1][j+1]/dx/dx;}
		
else
        applyBCpoissonRight(b,x,bcPHIfaceXright,0,nx,ny,dx,dy);
		

if (vct->getXleft_neighbor() != MPI_PROC_NULL){
	for (int j=0; j<ny-2; j++)
		b[0][j] -= x[0][j+1]/dx/dx;
//		b[0][j] -= 10/dx/dx;
}
else
	applyBCpoissonLeft(b,x,bcPHIfaceXleft,0,nx,ny,dx,dy);

if (vct->getYright_neighbor() != MPI_PROC_NULL){
	for (int i=0; i<nx-2; i++)
		b[i][ny-3] -= x[i+1][ny-1]/dy/dy;}
else
	applyBCpoissonRight(b,x,bcPHIfaceYright,1,nx,ny,dx,dy);

if (vct->getYleft_neighbor() != MPI_PROC_NULL){
	for (int i=0; i<nx-2; i++)
		b[i][0] -= x[i+1][0]/dy/dy;}
else
	applyBCpoissonLeft(b,x,bcPHIfaceYleft,1,nx,ny,dx,dy);
	
// end BC
	
// transform the vector b in FFTW-format
 for (register int i=0; i < nx-2; i++)
	 for (register int j=0; j < ny-2; j++){
    b_fftw[ j + (ny-2)*i]= b[i][j];
//	cout<<"  bfftw["<<i<<"]["<<j<<"] ="<<b_fftw[j + (ny-2)*i]<<endl;
} 
fftw_execute(plan_b); // apply forward transform
	 
for (int i=0; i<nx-2; i++)
	for (int j=0; j<ny-2; j++){
	b_fftw[ j + (ny-2)*i] = b_fftw[ j + (ny-2)*i]/2/sqrt((nx-1)*(ny-1)); // normalization
        phi_fftw[ j + (ny-2)*i] = dx*dx*dy*dy*b_fftw[j + (ny-2)*i]/2/(dy*dy*cos(PI*(i+1)/(nx-1)) + dx*dx*cos(PI*(j+1)/(ny-1)) -(dx*dx+dy*dy));
}					
fftw_execute(plan_phi); // backward transform

for (int i=0; i<nx-2; i++)
	for (int j=0; j<ny-2; j++){
	phi_fftw[ j + (ny-2)*i] = phi_fftw[ j + (ny-2)*i]/2/sqrt((nx-1)*(ny-1)); // normalization
}

// Compute r = b -A*PHI

for (register int i=0; i < nx-2; i++)
for (register int j=0; j < ny-2; j++){
	x[i+1][j+1] = phi_fftw[j + (ny-2)*i];
	
//	for (int j=0; j<ny-1; j++)
//		x[0][j+1]=10;
	//-phi_fftw[0]; // phi_fftw[0] is just an additive constant

//	cout<<"  phi["<<i<<"]["<<j<<"] ="<<phi_fftw[j + (ny-2)*i]<<endl;
	}
/*// prova media 1/2
eq(prova,x,nx,ny);	
	communicateNode(nx,ny,x,vct);
	for (register int i=0; i < ny; i++){
		x[0][i]=(0.01*x[0][i]+0.99*prova[0][i]);
	x[nx-1][i]=(0.01*x[nx-1][i]+0.99*prova[nx-1][i]);}
	for (register int j=0; j < nx; j++){
		x[j][0]=(0.01*x[j][0]+0.99*prova[j][0]);
	x[j][ny-1]=(0.01*x[j][ny-1]+0.99*prova[j][ny-1]);}
*/		
	communicateNode(nx,ny,x,vct);
//	communicateGhostBC(nx,ny,x,bcPHIfaceXright,bcPHIfaceXleft,bcPHIfaceYright,bcPHIfaceYleft,vct);
	
	for (register int i=0; i < nx; i++)
		for (register int j=0; j < ny; j++)
//		cout<<"x["<<i<<"]["<<j<<"] ="<<x[i][j]<<endl;
	

	for (register int i=1; i < nx-1; i++)
		for (register int j=1; j < ny-1; j++){
			laplacian[i][j]=(x[i-1][j]+x[i+1][j]-2*x[i][j])/dx/dx+(x[i][j-1]+x[i][j+1]-2*x[i][j])/dy/dy;
 			gradx[i][j]=0.5*(x[i+1][j]-x[i-1][j])/dx;
//			cout<<"laplacian phi["<<i<<"]["<<j<<"]  = "<<laplacian[i][j]<<" bPoisson= "<<bPoisson[i][j]<<" gradx ="<<gradx[i][j]<<endl;
			
		}
// calculate the laplacian
	grid->lapC2C(im,x,vct);
    // move from physical space to krylov space, in this case ghost cells don't count
//	phys2solver(image,im,nx,ny);
	phys2solver(image,laplacian,nx,ny);
sub(r,b_solver,image,(nx-2)*(ny-2));

for (register int i=0; i < (nx-2)*(ny-2); i++)
//	cout<<"r["<<i<<"] ="<<r[i]<<endl; 
error = dotP(vct->getCART_COMM(),r,r,(nx-2)*(ny-2));

eq(b,b_initial,nx-2,ny-2);


if (PoissonVerbose && vct->getCartesian_rank()==0)
	cout<<"Poisson solver (FFTW): iteration #" <<counter<< " ; relative error:"<<error/normb<<endl;

counter++;

};

communicateNode(nx,ny,x,vct);
	for (register int i=1; i < nx-1; i++)
		for (register int j=1; j < ny-1; j++){
		laplacian[i][j]=(x[i-1][j]+x[i+1][j]-2*x[i][j])/dx/dx+(x[i][j-1]+x[i][j+1]-2*x[i][j])/dy/dy;
		gradx[i][j]=0.5*(x[i+1][j]-x[i-1][j])/dx;
//			cout<<"laplacian phi["<<i<<"]["<<j<<"]  = "<<laplacian[i][j]<<" bPoisson= "<<bPoisson[i][j]<<" gradx ="<<gradx[i][j]<<endl;
		}
		phys2solver(image,laplacian,nx,ny);
		sub(r,b_solver,image,(nx-2)*(ny-2));
		error = dotP(vct->getCART_COMM(),r,r,(nx-2)*(ny-2));


if (PoissonVerbose && vct->getCartesian_rank()==0)
	cout<<"Poisson solver (FFTW) has reached iteration #" <<counter<< " with relative error:"<<error/normb<<endl;

fftw_destroy_plan(plan_b);
fftw_destroy_plan(plan_phi);

}

void applyBCpoissonRight(double** b,double** x, int bcPHI ,int DIM, int nx, int ny, double dx, double dy){
switch(DIM){
	case 0: // along x
		switch(bcPHI){
		 case 0: // don't get the difference with case 1...just put the same !!
		 for (int j=1; j<ny-3; j++) 			//without corners
		 b[nx-3][j] -= -x[nx-3][j+1]/dx/dx;
		 b[nx-3][0] -= x[nx-3][1]/dx/dx ;	// adjustment for corner (it will take another contribution)	
		 b[nx-3][ny-3] -= x[nx-3][ny-2]/dx/dx ;
		 break;
	 
		 case 1: //Dirichlet --> x[nx-1]= -x[nx-3]
 		 for (int j=1; j<ny-3; j++) 			//without corners
		   b[nx-3][j] -= -x[nx-3][j+1]/dx/dx;
		   b[nx-3][0] -= x[nx-3][1]/dx/dx ;	// adjustment for corner (it will take another contribution)	
	 	   b[nx-3][ny-3] -= x[nx-3][ny-2]/dx/dx ;
		  break;
		 
		 case 2: // Neumann --> x[nx-1]=x[nx-3]
		 for (int j=1; j<ny-3; j++) 			//without corners
 		  b[nx-3][j] -= x[nx-3][j+1]/dx/dx;
		  b[nx-3][0] -= x[nx-3][1]/dx/dx; 		// adjustment for corner	
		  b[nx-3][ny-3] -= x[nx-3][ny-2]/dx/dx; 
		 break;
		}
	break;
	case 1: //along y
		switch(bcPHI){
 		 case 0: // don't get the difference with case 1...just put the same !! 
		  for (int i=1; i<nx-3; i++)
		  b[i][ny-3] -= -x[i+1][ny-3]/dy/dy;
		  b[0][ny-3] -= x[1][ny-3]/dy/dy;
		  b[nx-3][ny-3] -= x[nx-2][ny-3]/dy/dy;
		  break;
	 	 case 1: //Dirichlet --> x[ny-1]= -x[ny-3]
  		  for (int i=1; i<nx-3; i++)
		  b[i][ny-3] -= -x[i+1][ny-3]/dy/dy;
		  b[0][ny-3] -= x[1][ny-3]/dy/dy;
		  b[nx-3][ny-3] -= x[nx-2][ny-3]/dy/dy;
		  break;
		  case 2: // Neumann --> x[ny-1]=x[ny-3]
		  for (int i=1; i<nx-3; i++)
		  b[i][ny-3] -= x[i+1][ny-3]/dy/dy;
		  b[0][ny-3] -= x[1][ny-3]/dy/dy;
		  b[nx-3][ny-3] -= x[nx-2][ny-3]/dy/dy;
		  break;
		}
	break;
	}
}

void applyBCpoissonLeft(double** b,double** x, int bcPHI ,int DIM, int nx, int ny, double dx, double dy){
switch(DIM){
	case 0: // along x
	 switch(bcPHI){
	  case 0: // don't get the difference with case 1...just put the same !!
	   for (int j=1; j<ny-3; j++)
	   b[0][j] -= -x[2][j+1]/dx/dx;
	   b[0][0] -= x[2][1]/dx/dx;
	   b[0][ny-3] -= x[2][ny-2]/dx/dx;
	   break;	
	  case 1: //Dirichlet --> x[0]= -x[2]
	   for (int j=1; j<ny-3; j++)
	   b[0][j] -= -x[2][j+1]/dx/dx;
	   b[0][0] -= x[2][1]/dx/dx;
	   b[0][ny-3] -= x[2][ny-2]/dx/dx;
	   break;	
	  case 2: // Neumann --> x[0]=x[2]
	   for (int j=1; j<ny-3; j++)
	   b[0][j] -= x[2][j+1]/dx/dx;
	   b[0][0] -= x[2][1]/dx/dx;
	   b[0][ny-3] -= x[2][ny-2]/dx/dx;
	   break;
		}
	break;
	  case 1: // along y
	   switch(bcPHI){
	    case 0: // don't get the difference with case 1...just put the same !!
	     for (int i=1; i<nx-3; i++)
	     b[i][0] -= -x[i+1][2]/dy/dy;
	     b[0][0] -= x[1][2]/dy/dy;	
	     b[nx-3][0] -= x[nx-2][2]/dy/dy;
	     break;
	    case 1: //Dirichlet --> x[0]= -x[2]
	     for (int i=1; i<nx-3; i++)
	     b[i][0] -= -x[i+1][2]/dy/dy;
	     b[0][0] -= x[1][2]/dy/dy;	
	     b[nx-3][0] -= x[nx-2][2]/dy/dy;
	     break;
	    case 2: // Neumann --> x[0]=x[2]
	     for (int i=1; i<nx-3; i++)
	     b[i][0] -= x[i+1][2]/dy/dy;
	     b[0][0] -= x[1][2]/dy/dy;	
	     b[nx-3][0] -= x[nx-2][2]/dy/dy;
 	     break;
           }
	   break;

}
}

#endif


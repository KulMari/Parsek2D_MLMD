/***************************************************************************
BcParticles.h  -  Library to manage boundary conditions for particles

developers: Stefano Markidis, Giovanni Lapenta
***************************************************************************/
#include <stdlib.h>

#ifndef BcParticles_H
#define BcParticles_H

/** set the boundary condition  for particle in 2D
<ul>
<li>bc = 1 perfect mirror </li>
<li>bc = 2 riemission     </li>
</ul>*/
inline void BCpart(double *x,double *u,double Lx,double ut,int bcFaceXright, int bcFaceXleft){
	if (*x > Lx){
		switch(bcFaceXright){
			case 1:   // perfect mirror
				*x = 2*Lx - *x;
				*u = - *u;
				break;
			case 2:   // riemmission
				double harvest, prob,theta;
				*x = 2*Lx - *x;
				// u
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*u = - fabs(ut*prob*cos(theta));
            	
				
				
	    }
	} else if (*x < 0) {
		switch(bcFaceXleft){
			case 1:   // perfect mirror
				*x = -*x;
				*u = -*u;
				break;
			case 2:   // riemmission
				*x = -*x;
				double harvest, prob,theta;
				// u
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*u = fabs(ut*prob*cos(theta));
				break;
				
		}
		
	}
	
}
/** set the boundary condition  for particle in 2D
<ul>
<li>bc = 1 perfect mirror </li>
<li>bc = 2 riemission     </li>
</ul>*/
inline void BCpart(double *x, double *u, double *v, double *w, double Lx,double ut, double vt, double wt,  int bcFaceXright, int bcFaceXleft){
  // original, no particles in the ghost cells
  if (*x > Lx){
    switch(bcFaceXright){
    case 1:   // perfect mirror
      *x = 2*Lx -*x;
      *u = - *u;
      break;
    case 2:   // riemmission
      double harvest, prob,theta;
      *x = 2*Lx -*x;
      // u
      harvest =   rand()/(double)RAND_MAX;
      prob  = sqrt(-2.0*log(1.0-.999999*harvest));
      harvest =   rand()/(double)RAND_MAX;
      theta = 2.0*M_PI*harvest;
      *u = - fabs(ut*prob*cos(theta));
      // v
      *v = vt*prob*sin(theta);
      // w
      harvest =   rand()/(double)RAND_MAX;
      prob  = sqrt(-2.0*log(1.0-.999999*harvest));
      harvest =   rand()/(double)RAND_MAX;
      theta = 2.0*M_PI*harvest;
      *w = wt*prob*cos(theta); 
      }
  } else if (*x < 0) {
    switch(bcFaceXleft){
    case 1:   // perfect mirror
      *x = -*x;
      *u = -*u;
      break;
    case 2:   // riemmission
      *x = -*x;
      double harvest, prob,theta;
      // u
      harvest =   rand()/(double)RAND_MAX;
      prob  = sqrt(-2.0*log(1.0-.999999*harvest));
      harvest =   rand()/(double)RAND_MAX;
      theta = 2.0*M_PI*harvest;
      *u = fabs(ut*prob*cos(theta));
      // v
      *v = vt*prob*sin(theta);
      // w
      harvest =   rand()/(double)RAND_MAX;
      prob  = sqrt(-2.0*log(1.0-.999999*harvest));
      harvest =   rand()/(double)RAND_MAX;
      theta = 2.0*M_PI*harvest;
      *w = wt*prob*cos(theta);   
      break;
    }
  }
}

// BCpart, 2D with particles also in the ghost cells
inline void BCpart(double *x, double *u, double *v, double *w, double Lx,double ut, double vt, double wt,  int bcFaceXright, int bcFaceXleft, double Dx, VirtualTopology *vct, bool PERIODIC){
  // particles also in Ghost Cells
  //if (*x > Lx+ Dx ) {
  if (*x > Lx ) {
    if (PERIODIC)
      {// periodic case
	//cout << "Periodic particle BC" <<endl;
	*x= *x - (Lx);
      }
    else if (*x > Lx+ Dx )
      {// non periodic case
	switch(bcFaceXright){
	case 1:   // perfect mirror
	  //cout <<"Particle BC: perfect mirror" <<endl;
	  *x = 2*(Lx +Dx) -*x;
	  *u = - *u;
	  break;
	case 2:   // riemmission
	  //cout <<"Particle BC: reimmission" <<endl;
	  double harvest, prob,theta;
	  *x =  2*(Lx +Dx) -*x;
	  // u
	  harvest =   rand()/(double)RAND_MAX;
	  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
	  harvest =   rand()/(double)RAND_MAX;
	  theta = 2.0*M_PI*harvest;
	  *u = - fabs(ut*prob*cos(theta));
	  // v
	  *v = vt*prob*sin(theta);
	  // w
	  harvest =   rand()/(double)RAND_MAX;
	  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
	  harvest =   rand()/(double)RAND_MAX;
	  theta = 2.0*M_PI*harvest;
	  *w = wt*prob*cos(theta); 
	}// end switch
      }// end non periodic case
  } //else if (*x < - Dx) {
  else if (*x < 0) {
    if (PERIODIC)
      {// periodic case
	*x= *x + (Lx);
	//cout << "Periodic particle BC" <<endl;
      }
    else if (*x < - Dx)
      {// non periodic case
	switch(bcFaceXleft){
	case 1:   // perfect mirror
	  *x = -Dx + (-*x -Dx);
	  *u = -*u;
	  break;
	case 2:   // riemmission
	  *x = -Dx + (-*x -Dx);
	  double harvest, prob,theta;
	  // u
	  harvest =   rand()/(double)RAND_MAX;
	  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
	  harvest =   rand()/(double)RAND_MAX;
	  theta = 2.0*M_PI*harvest;
	  *u = fabs(ut*prob*cos(theta));
	  // v
	  *v = vt*prob*sin(theta);
	  // w
	  harvest =   rand()/(double)RAND_MAX;
	  prob  = sqrt(-2.0*log(1.0-.999999*harvest));
	  harvest =   rand()/(double)RAND_MAX;
	  theta = 2.0*M_PI*harvest;
	  *w = wt*prob*cos(theta);   
	  break;
	}// end switch
      }// end non periodic
  }// end x<0
}

/** set the boundary condition on boundaries for particle in 3D*/
inline void BCpart(double *x, double *y, double *z, double *u, double *v, double *w, double Lx, double Ly, double Lz, double ut, double vt, double wt,  int bcFaceXright, int bcFaceXleft, int bcFaceYright, int bcFaceYleft,int bcFaceZright,int bcFaceZleft, VirtualTopology *vct){
            //if (*x > Lx && vct->getXright_neighbor()==-1){
            if (*x > Lx && vct->getXright_neighbor()==MPI_PROC_NULL){
		switch(bcFaceXright){
            case 1:   // perfect mirror
				*x = 2*Lx - *x;
				*u = - *u;
				break;
            case 2:   // riemmission
				double harvest, prob,theta;
				*x = 2*Lx - *x;
				// u
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*u = - fabs(ut*prob*cos(theta));
				// v
				*v = vt*prob*sin(theta);
				// w
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*w = wt*prob*cos(theta);   
				break;
				
				
				
		}
	}
	    //if (*x < 0 && vct->getXleft_neighbor()==-1){
	    if (*x < 0 && vct->getXleft_neighbor()==MPI_PROC_NULL){
		switch(bcFaceXleft){
            case 1:   // perfect mirror
				*x = -*x;
				*u = -*u;
				break;
            case 2:   // riemmission
				*x = -*x;
				double harvest, prob,theta;
				// u
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*u = fabs(ut*prob*cos(theta));
				// v
				*v = vt*prob*sin(theta);
				// w
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*w = wt*prob*cos(theta);   
				break;
				
				
				
		}
	}
	    //if (*y > Ly && vct->getYright_neighbor()==-1){
	    if (*y > Ly && vct->getYright_neighbor()==MPI_PROC_NULL){
		switch(bcFaceYright){
            case 1:   // perfect mirror
				*y = 2*Ly - *y;
				*v = -*v;
				break;
            case 2:   // riemmission
				*y = 2*Ly - *y;
				double harvest, prob,theta;
				// u
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*u = ut*prob*cos(theta);
				// v
				*v = -fabs(vt*prob*sin(theta));
				// w
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*w = wt*prob*cos(theta); 
				break;
				
				
				
		}
	}
	    //if (*y < 0 && vct->getYleft_neighbor()==-1){
	    if (*y < 0 && vct->getYleft_neighbor()==MPI_PROC_NULL){ 
		switch(bcFaceYleft){
            case 1:   // perfect mirror
				*y = -*y;
				*v = -*v;
				break;
            case 2:   // riemmission
				*y = -*y;
				double harvest, prob,theta;
				// u
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*u = ut*prob*cos(theta);
				// v
				*v = fabs(vt*prob*sin(theta));
				// w
				harvest =   rand()/(double)RAND_MAX;
				prob  = sqrt(-2.0*log(1.0-.999999*harvest));
				harvest =   rand()/(double)RAND_MAX;
				theta = 2.0*M_PI*harvest;
				*w = wt*prob*cos(theta);
				break;
				
				
				
		}
	}
}
#endif


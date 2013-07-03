/***************************************************************************
                          GMRES.h  - Generalized minimal residual solver
                             -------------------
    begin             : Wed Jul 14 2004
    copyright         : (C) 2004 Los Alamos National Laboratory
    developers        : Stefano Markidis, Giovanni Lapenta
    email             : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/


#ifndef GMRES_H
#define GMRES_H

#include <iostream>

#include "../mathlib/Basic.h"


typedef  void (Field::*FIELD_IMAGE)(double*, double*, Grid*, VirtualTopology*);
typedef  void (*GENERIC_IMAGE)(double*, double*, Grid*, VirtualTopology*);

using std::cout;
using std::cerr;
using std::endl;
/**
*  GMRES  - Generalized minimal residual  solver
*
* @date Fri Jun 4 2004
* @par Copyright:
* (C) 2004 Los Alamos National Laboratory
* @author Stefano Markidis, Giovanni Lapenta
* @version 1.0
*/
/** GMRES */
/** Modified Arnoldi: method to calculate the modified Arnoldi decomposition */
inline void ModifiedArnoldi(int m, double *x,int xkrylovlen, double **H, double **V,FIELD_IMAGE FunctionImage,Grid *grid,VirtualTopology *vct,Field *field){
   // allocate v, w
   double *v = new double[xkrylovlen];
   double *w = new double[xkrylovlen];
   // initialize v. The initial direction vector xkrylov defaults to the first unit vector
   eq(v,x,xkrylovlen);
   putColumn(V,v,0,xkrylovlen);
   for (int j=0; j < m; j++){
     (field->*FunctionImage)(w,v,grid,vct);
     for(int i=0;i <=j;i++){
        getColumn(v,V,i,xkrylovlen);
        H[i][j] = dotP(w,v,xkrylovlen);
        addscale(-H[i][j],w,v,xkrylovlen);
     }
     H[j+1][j] =normP(w,xkrylovlen);
     for (int s=0; s < xkrylovlen; s++)
      v[s] = w[s]/H[j+1][j];
     putColumn(V,v,j+1,xkrylovlen);
   }
}

inline void GMRES(double *xkrylov, int xkrylovlen, double *b, int maxit, double tol, int m ,FIELD_IMAGE FunctionImage,Grid *grid,VirtualTopology *vct,Field *field){
   if (m > xkrylovlen){
     if(vct->getCartesian_rank()==0)
      cerr << "In GMRES the dimension of Krylov space(m) can't be > (length of krylov vector)/(# processors)" << endl;
     return;
   }
   
   
   bool GMRESVERBOSE = true;
   int nr, counter;
   double delta,rho,tmp,initial_error;
   // allocate r, z
   double *r  = new double[xkrylovlen];
   double *im = new double[xkrylovlen];
   double *z  = new double[xkrylovlen];
   // allocate y, c, s
   double *y = new double[m+1];
   double *c = new double[m+1];
   eqValue(0.0,c,m+1);
   double *s = new double[m+1];
   eqValue(0.0,s,m+1);
   // allocate H  for storing the results from decomposition
   double **H = newArr(double,m+1,m);
   for (int ii=0; ii < m+1;ii++)
    for (int jj=0; jj < m;jj++)
      H[ii][jj] =0;
   // allocate V
   double **V = newArr(double,xkrylovlen,m+1);
   for (int ii=0; ii < xkrylovlen;ii++)
    for (int jj=0; jj < m+1;jj++)
      V[ii][jj] =0;
   if(GMRESVERBOSE && vct->getCartesian_rank()==0){
      cout << "------------------------------------" << endl;
      cout << "-             GMRES                -" << endl;
      cout << "------------------------------------" << endl;
      cout << endl;
    }
   // initialize xkrylov: initial guessing = 0;
   eqValue(0.0,xkrylov,xkrylovlen);
   // r = b  - A*x
   (field->*FunctionImage)(im,xkrylov,grid,vct);
   sub(r,b,im,xkrylovlen);
   initial_error = normP(r,xkrylovlen);
   if (vct->getCartesian_rank() ==0)
    cout << "Initial error: " << initial_error << endl;
   if (initial_error < 1E-16)
   	return;
   for(int j=0;j<maxit;j++){
      eqValue(0.0,y,m+1);
      y[0] = normP(r,xkrylovlen);
      if (y[0] > 10E8*initial_error){
        if(vct->getCartesian_rank() ==0){
          cerr << "GMRES not converging" << endl;
          cerr << "GMRES stopped" << endl;
        }
        break;
      }
      if (vct->getCartesian_rank()==0)
        cout << "GMRES Restart # " << j << " - Norm of Residual relative to initial error = " << y[0]/initial_error << endl;
      // Normalize the residual
      for (int ii=0; ii < xkrylovlen;ii++)
        r[ii] /= y[0];
      /// MODIFIED ARNOLDI
      ModifiedArnoldi(m,r,xkrylovlen,H,V,FunctionImage,grid,vct,field);
      // Givens rotations to accomplish QR factorization
      counter = 0;
      for (int i=0;i<m;i++){
          for (int k=1;k <=i; k++){
            tmp = H[k-1][i];
            H[k-1][i] = c[k-1]*H[k-1][i] + s[k-1]*H[k][i];
            H[k][i]  = -s[k-1]*tmp + c[k-1]*H[k][i];
	    cout<<"i= "<<i<<" k= "<<k<<endl;
	    cout<<"  H[k-1][i]="<<H[k-1][i]<<"  H[k][i]="<<H[k][i]<<endl;
	  }
//	  cout<<"s[0]="<<s[0]<<" c[0] = "<<c[0];
//	  cout<<"  H[1][0]="<<H[1][0]<<endl;
          delta = sqrt(H[i][i]*H[i][i] + H[i+1][i]*H[i+1][i]);
          c[i] = H[i][i]/delta;
          s[i] = H[i+1][i]/delta;
	  cout<<"s["<<i<<"] = "<<s[i]<<endl;
          H[i][i] = c[i]*H[i][i] + s[i]*H[i+1][i];

          for (int k=i+1;k<m+1;k++)
            H[k][i] = 0.0;
          y[i+1] = -s[i]*y[i];
          y[i] = c[i]*y[i];
          rho = fabs(y[i+1]);
	  if ((GMRESVERBOSE) && (vct->getCartesian_rank()==0) && (i >=1))
              cout << "GMRES Iteration # " << i  << "\t relative error:" << rho << endl;
	  cout<<"rho = "<<rho<<endl;
	  if (rho < tol*initial_error) {
//		  if (rho < tol) {
            
		  nr = i;
		  //	    cout<<"    nr = "<<nr<<"   i ="<<i<<endl;
            counter =i;
            break;
          }
          counter++;
      }
      // Backsolve to obtain coefficients
      eqValue(0.0,z,xkrylovlen);
      if (counter>=(m-1)){
        nr = m;
        z[nr-1] = y[nr-1]/H[nr-1][nr-1];
      }
      
      for(int k=nr-2;k>=0;k--){
        z[k] = y[k];
	for (int ll=k+1;ll<nr;ll++){
		z[k] -= H[k][ll]*z[ll];	//cout<<"H[k][ll]= "<<H[k][ll]<<endl;
	}

        z[k] = z[k]/H[k][k];
      }
      // Linear combination of basis vectors of the krylov space
      for (int i=0; i < nr;i++){
        getColumn(r,V,i,xkrylovlen);
        addscale(z[i],xkrylov,r,xkrylovlen);
      }
      if(rho<tol*initial_error)
        break;
      // r = b  - A*x
      (field->*FunctionImage)(im,xkrylov,grid,vct);
      sub(r,b,im,xkrylovlen);
      
   }
   (field->*FunctionImage)(im,xkrylov,grid,vct);
   sub(r,b,im,xkrylovlen);
   double final_error = normP(r,xkrylovlen);
   delete[] r;
   delete[] im;
   delete[] z;
   delete[] y;
   delete[] c;
   delete[] s;
   delArr(H,m+1);
   delArr(V,xkrylovlen);


}




#endif

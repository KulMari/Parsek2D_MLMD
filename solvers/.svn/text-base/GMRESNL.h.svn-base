 /***************************************************************************
                           GMRESNL.h  - Generalized minimal residual solver for non linear system
                              -------------------
     begin             : Wed Jul 14 2004
     copyright         : (C) 2004 Los Alamos National Laboratory
     developers        : Stefano Markidis, Giovanni Lapenta
     email             : markidis@lanl.gov, lapenta@lanl.gov
  ***************************************************************************/


#ifndef GMRESNL_H
#define GMRESNL_H

#include <iostream>

#include "../mathlib/Basic.h"
#include "../mathlib/DirDer.h"

typedef  void (Field::*FIELD_IMAGE)(double*, double*, Grid*, VirtualTopology*);
typedef  void (*GENERIC_IMAGE)(double*, double*, Grid*, VirtualTopology*);

using std::cout;
using std::cerr;
using std::endl;
 /**
 *  GMRESNL  - Generalized minimal residual  solver for non linear system
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */
 /** Non linear GMRES. it takes as argumnt a pointer to a method of Field abstract class */
 inline void GMRESNL(double *xkrylov, int xkrylovlen, int maxit, double tol, int m,FIELD_IMAGE FunctionImage,Grid *grid, VirtualTopology *vct,Field *field){
   if (m > xkrylovlen){
     if(vct->getCartesian_rank()==0)
      cerr << "In GMRES the dimension of Krylov space(m) can't be > (length of krylov vector)/(# processors)" << endl;
     return;
   }
   bool GMRESVERBOSE = true;
   int nr, counter;
   double delta,rho,tmp,initial_error;
   // allocate r, z
   double *b  = new double[xkrylovlen];
   double *f0 = new double[xkrylovlen];
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
   (field->*FunctionImage)(b,xkrylov,grid,vct);
   neg(b,xkrylovlen);
   eq(r,b,xkrylovlen);
   initial_error = normP(r,xkrylovlen);
   rho = initial_error;
   for(int j=0;j<maxit;j++){
      (field->*FunctionImage)(f0,xkrylov,grid,vct);
      eqValue(rho,y,m+1);
      y[0] = rho;
      if (y[0] > 10E8*initial_error){
        if(vct->getCartesian_rank() ==0){
          cerr << "GMRES not converging" << endl;
          cerr << "GMRES stopped" << endl;
        }
        break;
      }
      if ((GMRESVERBOSE) && (vct->getCartesian_rank()==0))
        cout << "Restart # " << j << " - Norm of Residual relative to initial error = " << y[0]/initial_error << endl;
      ///----------------------------------
      /// MODIFIED ARNOLDI
      ///----------------------------------
      double *v = new double[xkrylovlen];
      double *w = new double[xkrylovlen];
      // initialize v. The initial direction vector xkrylov defaults to the first unit vector
      eq(v,r,xkrylovlen);
      putColumn(V,v,0,xkrylovlen);
      for (int jjj=0; jjj < m; jjj++){
          DirDer(w,xkrylov,xkrylovlen,v,f0,FunctionImage,grid,vct,field);
          for(int iii=0;iii <=jjj;iii++){
              getColumn(v,V,iii,xkrylovlen);
              H[iii][jjj] = dotP(w,v,xkrylovlen);
             addscale(-H[iii][jjj],w,v,xkrylovlen);
          }
          H[jjj+1][jjj] =normP(w,xkrylovlen);
          for (int s=0; s < xkrylovlen; s++)
            v[s] = w[s]/H[jjj+1][jjj];
          putColumn(V,v,jjj+1,xkrylovlen);
      }
      delete[] v;
      delete[] w;

      // Givens rotations to accomplish QR factorization
      counter = 0;
      for (int i=0;i<m;i++){
          for (int k=1;k <=i; k++){
            tmp = H[k-1][i];
            H[k-1][i] = c[k-1]*H[k-1][i] + s[k-1]*H[k][i];
            H[k][i]  = -s[k-1]*tmp + c[k-1]*H[k][i];
          }
          delta = sqrt(H[i][i]*H[i][i] + H[i+1][i]*H[i+1][i]);
          c[i] = H[i][i]/delta;
          s[i] = H[i+1][i]/delta;

          H[i][i] = c[i]*H[i][i] + s[i]*H[i+1][i];

          for (int k=i+1;k<m+1;k++)
            H[k][i] = 0.0;
          y[i+1] = -s[i]*y[i];
          y[i] = c[i]*y[i];
          rho = fabs(y[i+1]);
          if ((GMRESVERBOSE) && (vct->getCartesian_rank()==0) && (i >=1))
              cout << "GMRES Iteration # " << i  << endl;
          if (rho < tol*initial_error) {
            nr = i;
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
        for (int ll=k+1;ll<nr;ll++)
          z[k] -= H[k][ll]*z[ll];
        z[k] = z[k]/H[k][k];
      }
      // Linear combination of basis vectors of the krylov space
      for (int i=0; i < nr;i++){
        getColumn(r,V,i,xkrylovlen);
        addscale(z[i],xkrylov,r,xkrylovlen);
      }
      if(rho<tol*initial_error)
        break;
      
   }
   // deallocate
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


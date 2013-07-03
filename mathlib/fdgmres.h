#ifndef FDGMRES_H
#define FDGMRES_H

#include <iostream>

#include "../mathlib/Basic.h"
#include "../mathlib/DirDer.h"


typedef  void (Field::*FIELD_IMAGE)(double*, double*, Grid*, VirtualTopology*);
typedef  void (*GENERIC_IMAGE)(double*, double*, Grid*, VirtualTopology*);

using std::cout;
using std::cerr;
using std::endl;
 /**
 * - Generalized minimal residual  solver for non linear system
 *
 * @date Fri Jun 4 2004
 * @par Copyright:
 * (C) 2004 Los Alamos National Laboratory
 * @author Stefano Markidis, Giovanni Lapenta
 * @version 1.0
 */
 /** Non linear GMRES. it takes as argumnt a pointer to a method of Field abstract class */
inline void fdgmres(double *x, double *xc, int n, int kmax, double errtol,int* total_iters, double* error, double* f0, FIELD_IMAGE FunctionImage,Grid *grid, VirtualTopology *vct,Field *field){
   // right side of linear for step is -f0 if the default initial iterate is used
   double rho;
   double beta;
   int k, counter;
   double normv;
   double normv2,hr,tmp,delta;
   int nr;
   // allocate r, b
   double *r  = new double[n];
   double *b  = new double[n];
   bool GMRESVERBOSE = true;
   //allocate h,v,c,s,g
   // allocate H  for storing the results from decomposition
   double **H = newArr(double,kmax+1,kmax);
   for (int ii=0; ii < kmax+1;ii++)
    for (int jj=0; jj < kmax;jj++)
      H[ii][jj] =0;
   double *c = new double[kmax+1];
   eqValue(0.0,c,kmax+1);
   double *s = new double[kmax+1];
   eqValue(0.0,s,kmax+1);
   double *g = new double[kmax+1];
   // allocate V
   double **V = newArr(double,n,kmax);
   double *v1 = new double[n];
   double *v2 = new double[n];
   double *y  = new double[n];
   // right side of linear for step is -f0 if the default initial iterate is used
   eq(b,f0,n);
   neg(b,n);
   // use zero vector as initial iterate for Newton step
   eq(r,b,n);
   rho = normP(r,n);
   eqValue(rho,g,kmax+1);
   errtol*=normP(b,n);
   total_iters=0;
   eq(v1,r,n);
   scale(r,1/rho,n);
   putColumn(V,v1,0,n);
   beta = rho;
  
   for (k=0; k < kmax; k++){
     // call directional derivative function
     getColumn(v1,V,k,n);
     DirDer(v2,xc,n,v1,f0,FunctionImage,grid,vct,field);
     putColumn(V,v2,k+1,n);
     normv = normP(v2,n);
     // Modified Gram-Schmidt
     for (int j=0;j < k;j++){
         getColumn(v1,V,j,n);
         H[j][k] = dotP(v1,v2,n);
         addscale(-H[j][k],v2,v1,n);
         putColumn(V,v2,k+1,n);
     }
     H[k+1][k] =normP(v2,n);
     normv2 = H[k+1][k];
     // reortogonalize?
     // Brown/Hindmarsh condition
     if (normv + .001*normv2 == normv){
        for (int j=0;j < k;j++){
           getColumn(v1,V,j,n);
           hr =  dotP(v1,v2,n);
           H[j][k] += hr;
           addscale(-hr,v2,v1,n);
           putColumn(V,v2,k+1,n);
        }
     }
     // watch out for happy breakdown
     if(H[k+1][k] !=0){
          scale(v2,1/H[k+1][k+1],n);
          putColumn(V,v2,k+1,n);
     }
     // form and store the information for the new Given's rotation
     // Givens rotations to accomplish QR factorization
      counter = 0;
      for (int i=0;i<kmax;i++){
          for (int j=1;j <=i; j++){
            tmp = H[j-1][i];
            H[j-1][i] = c[j-1]*H[j-1][i] + s[j-1]*H[j][i];
            H[j][i]  = -s[j-1]*tmp + c[j-1]*H[j][i];
          }
          delta = sqrt(H[i][i]*H[i][i] + H[i+1][i]*H[i+1][i]);
          c[i] = H[i][i]/delta;
          s[i] = H[i+1][i]/delta;

          H[i][i] = c[i]*H[i][i] + s[i]*H[i+1][i];

          for (int j=i+1;j<kmax+1;j++)
            H[j][i] = 0.0;
          g[i+1] = -s[i]*g[i];
          g[i] *= c[i];
          rho = fabs(g[i+1]);
          if ((GMRESVERBOSE) && (vct->getCartesian_rank()==0) && (i >=1))
              cout << "GMRES Iteration # " << i  << endl;
          if (rho < errtol) {
            nr = i;
            counter =i;
            break;
          }
          counter++;
      }
      // Backsolve to obtain coefficients
      eqValue(0.0,y,n);
      if (counter>=(kmax-1)){
        nr = kmax;
        g[nr-1] /= H[nr-1][nr-1];
      }
      for(int j=nr-2;j>=0;j--){
        y[j] = g[j];
        for (int ll=j+1;ll<nr;ll++)
          y[j] -= H[j][ll]*y[ll];
        y[j] /= H[j][j];
      }
      // Linear combination of basis vectors of the krylov space
      for (int i=0; i < nr;i++){
        getColumn(v1,V,i,n);
        addscale(y[i],x,v1,n);
      }
      error[k] = rho;
      cout << x[0]in << endl;
  }

   
   
}
#endif


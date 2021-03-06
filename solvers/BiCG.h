 
#ifndef GMRES_new_H
#define GMRES_new_H

#include <iostream>
#include <math.h> 
#include "../mathlib/Basic.h"


typedef  void (Field::*FIELD_IMAGE)(double*, double*, Grid*, VirtualTopology*);
typedef  void (*GENERIC_IMAGE)(double*, double*, Grid*, VirtualTopology*);


using std::cout;
using std::cerr;
using std::endl;

inline void GMRES(FIELD_IMAGE FunctionImage, double *xkrylov, int xkrylovlen, double *b,int m, int max_iter,double tol, Grid *grid,VirtualTopology *vct, Field *field){

inline void BiCG(FIELD_IMAGE FunctionImage, double *x, int xkrylovlen, double *b, int max_iter, double tol, Grid *grid,VirtualTopology *vct, Field *field){
 double resid, normb;
 rho_1 = new double [2];
 rho_2 = new double [2];
 alpha = new double [2];
 beta = new double [2];
 
 double *im = new double[xkrylovlen];
 double *r = new double[xkrylovlen];
 double *rtilde = new double[xkrylovlen];
 double *p = new double[xkrylovlen];
 double *ptilde = new double[xkrylovlen];
 double *q = new double[xkrylovlen];
 double *qtilde = new double[xkrylovlen];


  normb = normP(b,xkrylovlen);
// r= b - A*x
(field->*FunctionImage)(im,xkrylov,grid,vct);
   sub(r,b,im,xkrylovlen);
   initial_error = normP(r,xkrylovlen);

  eq(rtilde,r,xkrylovlen);

  if (normb == 0.0)
    normb = 1;
  resid= initial_error/normb;
  if (resid <= tol) {
    if (vct->getCartesian_rank() ==0)
    cout<<"BiCG converged without iterations: initial error < tolerance"<<endl;
    return ;
  }

  for (int i = 1; i <= max_iter; i++) {
//    z = M.solve(r);
//    ztilde = M.trans_solve(rtilde);
    rho_1[0] = dotP(vct->getCART_COMM(),r,r,xkrylovlen);
    if (rho_1[0] == 0) { 
      tol = norm(r) / normb;
      max_iter = i;
    if (vct->getCartesian_rank() ==0)
    cout<<"BiCG: breakdown occurred !!"<<endl;
    return ;
    }
    if (i == 1) {
      eq(p,r,xkrylovlen);
      eq(ptilde,rtilde,xkrylovlen);
    } else {
      beta[0] = rho_1[0] / rho_2[0];
addscale(1,beta(0),p,r,xkrylovlen);
addscale(1,beta(0),ptilde,rtilde,xkrylovlen);
    }

//     q = A * p;
(field->*FunctionImage)(q,p,grid,vct);

    q = A * p;
    qtilde = A.trans_mult(ptilde);
    alpha(0) = rho_1(0) / dot(ptilde, q);
    x += alpha(0) * p;
    r -= alpha(0) * q;
    rtilde -= alpha(0) * qtilde;

    rho_2(0) = rho_1(0);
    if ((resid = norm(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;
}











#endif

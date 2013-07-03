#ifndef GMRES_new_H
#define GMRES_new_H

#include <iostream>
#include <math.h> 
#include "../mathlib/Basic.h"


typedef  void (Field::*FIELD_IMAGE)(double*, double*, Grid*, VirtualTopology*);
typedef  void (*GENERIC_IMAGE)(double*, double*, Grid*, VirtualTopology*);

void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn);
void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn);

using std::cout;
using std::cerr;
using std::endl;



inline void GMRES(FIELD_IMAGE FunctionImage, double *xkrylov, int xkrylovlen, double *b,int m, int max_iter,double tol, Grid *grid,VirtualTopology *vct, Field *field){
if (m > xkrylovlen){
     if(vct->getCartesian_rank()==0)
      cerr << "In GMRES the dimension of Krylov space(m) can't be > (length of krylov vector)/(# processors)" << endl;
     return;
   }
  bool GMRESVERBOSE = true;
  double resid, initial_error, normb;
  int i, j , k;
   double *r  = new double[xkrylovlen];
   double *im = new double[xkrylovlen];
  double *v = new double[xkrylovlen];
  double *w = new double[xkrylovlen];
    

  double *s = new double[m+1];
  double *cs = new double[m+1];
  double *sn = new double[m+1];
  double *y = new double[m+1];

eqValue(0.0,s,m+1);
eqValue(0.0,cs,m+1);
eqValue(0.0,sn,m+1);
eqValue(0.0,y,m+1);

  
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
/*if (vct->getCartesian_rank() ==0){
 for (int i=0; i<xkrylovlen; i++)
	    cout<<"b["<<i<<"]= "<<b[i]<<endl;}
*/
// initialize xkrylov: initial guessing = 0;
//   eqValue(0.0,xkrylov,xkrylovlen);

// r = b  - A*x
   (field->*FunctionImage)(im,xkrylov,grid,vct);
/*if (vct->getCartesian_rank() ==0){
 for (int i=0; i<xkrylovlen; i++)
	    cout<<"im["<<i<<"]= "<<im[i]<<endl;}*/
   sub(r,b,im,xkrylovlen);
   initial_error = normP(r,xkrylovlen);
   normb= normP(b,xkrylovlen);
   if (normb == 0.0)
    normb = 1.0;
   if (vct->getCartesian_rank() ==0){
    cout << "Initial error: " << initial_error << endl;
    cout << "normb: " << normb << endl;}

  if ((initial_error / normb) <= tol) {
    tol = initial_error / normb;
    if (vct->getCartesian_rank() ==0)
    cout<<"GMRES converged without iterations: initial error < tolerance"<<endl;
    return; }

  int counter=1;
  while (counter <= max_iter) {
   scale(v,r,(1.0 / initial_error),xkrylovlen);
  putColumn(V,v,0,xkrylovlen); 
  eqValue(0.0,s,m+1);
    s[0] = initial_error;
// start iterations    
  for (i = 0; i < m && counter <= max_iter; i++) {
    
// w= A*V(:,i)
      getColumn(v,V,i,xkrylovlen);
      (field->*FunctionImage)(w,v,grid,vct);
      for (k = 0; k <= i; k++) {
      getColumn(v,V,k,xkrylovlen);
      H[k][i] = dotP(w,v,xkrylovlen);
      addscale(-H[k][i],w,v,xkrylovlen);
      }

     H[i+1][i] = normP(w,xkrylovlen);
     for (register int l=0; l < xkrylovlen; l++)
     v[l] = w[l]/H[i+1][i];
     putColumn(V,v,i+1,xkrylovlen);


      for (k = 0; k < i; k++)
      ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);

      GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
      ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
      ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);
      resid=fabs(s[i+1]) / normb ;
//      (field->*FunctionImage)(im,xkrylov,grid,vct);
//      sub(r,b,im,xkrylovlen); 
//      resid=normP(r,xkrylovlen)/normb;
 //     cout<<"resid="<<resid<<" counter = "<<counter<<" i ="<<i<<" tol="<<tol<<endl;
      if ( resid < tol) {
// Update
//eqValue(0.0,y,m+1);
eq(y,s,m+1);
     for (int ii = i; ii >= 0; ii--) {
       y[ii] = s[ii]/H[ii][ii];
       for (int jj = ii - 1; jj >= 0; jj--)
         y[jj] -= H[jj][ii] * y[ii];
      }

     for (int jj = 0; jj <= i; jj++){
      getColumn(v,V,jj,xkrylovlen);
      addscale(y[jj],xkrylov,v,xkrylovlen);}
      if (vct->getCartesian_rank()==0)
      cout<<"GMRES converged at restart # "<<counter<<"; iteration #"<<i<<" with error: "<<resid<<endl;
        delete[] r;
   	delete[] im;
	delete[] s;
	delete[] v;
   	delete[] cs;
   	delete[] sn;
   	delete[] w;
	delete[] y;
	delArr(H,m+1);
   	delArr(V,xkrylovlen);
      return ;
      }
    }
    if (vct->getCartesian_rank() ==0 & GMRESVERBOSE)
    cout << "Restart: " << counter <<" error: "<<resid<< endl;

// Update
//eqValue(0.0,y,m+1);
eq(y,s,m+1);
     for (int ii = i-1; ii >= 0; ii--) {
       y[ii] = s[ii]/H[ii][ii];

       for (int jj = ii - 1; jj >= 0; jj--)
         y[jj] -= H[jj][ii] * y[ii];
      }

     for (int jj = 0; jj <= i-1; jj++){
       getColumn(v,V,jj,xkrylovlen);
       addscale(y[jj],xkrylov,v,xkrylovlen);}

//    Update(xkrylov, i - 1, H, s, V, xkrylovlen, m);

// r = b  - A*x
    (field->*FunctionImage)(im,xkrylov,grid,vct);
    sub(r,b,im,xkrylovlen);
//    for (int it=0; i<xkrylovlen;i++)
//	    cout<<"r["<<i<<"] = "<<r[i]<<endl;
    initial_error= normP(r,xkrylovlen);
    resid = initial_error/normb;
    if (resid < tol) {
     if (vct->getCartesian_rank()==0)
     cout<<"GMRES converged at restart # "<<counter<<"; iteration #"<<i<<" with error: "<<resid<<endl;
        delete[] r;
   	delete[] im;
	delete[] s;
	delete[] v;
   	delete[] cs;
   	delete[] sn;
   	delete[] w;
	delete[] y;
	delArr(H,m+1);
   	delArr(V,xkrylovlen);
     return ;
    }
counter++;
  }
  
if (vct->getCartesian_rank()==0 )
cout<<"GMRES not converged !! Final error: "<<resid<<endl;

   delete[] r;
   delete[] im;
   delete[] s;
   delete[] v;
   delete[] cs;
   delete[] sn;
   delete[] w;
   delete[] y;
   delArr(H,m+1);
   delArr(V,xkrylovlen);
return;
}


inline void GeneratePlaneRotation(double &dx, double &dy, double &cs, double &sn){
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    double temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    double temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }

}

inline void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn){
  double temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}



#endif



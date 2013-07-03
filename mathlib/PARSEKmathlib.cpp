/***************************************************************************
           PARSEKmathlib.cpp  -  Mathematics Library Implementation
                             -------------------
    begin                : Mon May 24 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#include "PARSEKmathlib.h"

/***********************************/
/*  Functions Declarations        */
/*********************************/

int min_dimension(const Vector& v1, const Vector& v2){
  int min_dim = (v1.Dimension() < v2.Dimension())?v1.Dimension():v2.Dimension();
  return(min_dim);
}

double dot(const Vector& u, const Vector& v){
  double sum = 0.0;
  int min_dim = min_dimension(u,v);
  for (int i =0;i < min_dim;i++)
    sum +=u(i)*v(i);

  return(sum);
}

double dot(int N, const Vector& u, const Vector& v){
  double sum = 0.0;

  for (int i=0; i < N; i++)
    sum += u(i)*v(i);

  return(sum);
}

double dot(int N, double *a, double *b){
  double sum =0.0;

  for(int i=0;i <N;i++)
    sum += a[i]*b[i];

  return(sum);
}

void swap(double &a, double &b){
  double tmp = a;
  a = b;
  b = tmp;
}

double Sign(double x){
  double xs;
  xs = (x>=0.0)?1.0:-1.0;
  return(xs);
}
  
  



    
  

      
  
  

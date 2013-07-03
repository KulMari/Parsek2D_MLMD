/***************************************************************************
                          PARSEKmathlib.h  -  Mathematics Library
                             -------------------
    begin                : Mon May 24 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/


#ifndef PARSEKMATHLIB_H
#define PARSEKMATHLIB_H

#include <iostream>
#include <math.h>

using std::cout;
using std::cerr;
/**
  * Vector class with vector operations. Used to simplfy the development of GMRES solver
  *
  * @date Fri Jun 4 2004
  * @par Copyright:
  * (C) 2004 Los Alamos National Laboratory
  * @author Stefano Markidis, Giovanni Lapenta
  * @version 1.0
  * 
  *
  */
class Vector{
  /** overloaded + operator */
  friend Vector operator+(const Vector& v1, const Vector& v2);
  /** overloaded - operator */
  friend Vector operator-(const Vector& v1, const Vector& v2);
  /** overloaded * operator alfa*vector */
  friend Vector operator*(const double alfa, const Vector& v);
  /** overloaded * operator vector*alfa */
  friend Vector operator*(const Vector& v, const double alfa);
  /** overloaded operator / vector/a */
  friend Vector operator/(const Vector& v, const double alfa);

  private:
  /** length of vector */
  int dimension;
  /** array with values of Vector */
  double *data;

  public:
  //*******************************/
  // Costructors & Destructors   */
  //*****************************/
  Vector();
  Vector(int dim);
  Vector(const Vector& v);
  ~Vector();

  //*******************************/
  // Overloaded Operators        */
  //*****************************/
  int     operator==(const Vector& v) const;
  int     operator!=(const Vector& v) const;
  Vector& operator=(const Vector& v);
  double  operator()(const int i) const;
  double& operator()(const int i);
  Vector  operator-(const Vector& v);
 
  //******************************
  // General Methods
  //******************************
  int Dimension() const;
  void Initialize(int dim);
  void Initialize(double a);
  void Initialize(double *v);
  double Norm_l1();
  double Norm_l2();
  double Norm_linf();
  void Normalize();
  void Print() const;
  
};

/*******************************/
/* Costructors & Destructors  */
/*****************************/
/** default constructor */
inline Vector::Vector(){
 dimension = 0;
 data = NULL;
}
/** constructor */
inline Vector::Vector(int dim){
 dimension = dim;
 data = new double[dimension];

 for (int i=0;i < dimension;i++)
  data[i] = 0.0;
}
/** copy constructor */
inline Vector::Vector(const Vector& v){
 dimension = v.Dimension();
 data = new double[dimension];

 for (int i=0;i < dimension;i++)
  data[i] = v.data[i];
}
/** destructor */
inline Vector::~Vector(){
 dimension=0;
 delete[] data;
 data = NULL;
}


/***********************************/
/*  Vector Operator Overloading   */
/*********************************/
/** assignment operator */
inline Vector& Vector::operator=(const Vector &v){
  dimension = v.Dimension();
  for (int i=0;i<dimension;i++)
    data[i] = v.data[i];
  return *this;
}
/** () : returns a double */
inline double Vector::operator()(const int i) const{
  // bounds check: when all is ok delete it
  if(i>=0 && i<dimension)
    return(data[i]);

  cerr << "Vector: Invalid index "<<i<<"for Vector of dimension " << dimension << "\n";
  return(0);
}
/** () : returns a reference */
inline double& Vector::operator()(const int i){
   // bounds check: when all is ok delete it
   if(i>=0 && i<dimension)
    return(data[i]);

  cerr << "Vector: Invalid index "<<i<<"for Vector of dimension " << dimension << "\n";
  return(data[0]);

}
/** unitary operator: - */
inline Vector Vector::operator-(const Vector& v){
  Vector x(v.dimension);
  for (int i=0;i <v.dimension;i++)
    x.data[i] = - v.data[i];
  return x;
}
/** sum of vectors: + */
inline Vector operator+(const Vector& v1, const Vector& v2){
  Vector x(v1.Dimension());
  for (int i=0;i < v1.Dimension();i++)
    x.data[i] = v1(i) + v2(i);
  return x;
}
/** difference of vectors: - */
inline Vector operator-(const Vector& v1, const Vector& v2){
  Vector x(v1.Dimension());
  for (int i=0;i < v1.Dimension();i++)
    x.data[i] = v1(i) - v2(i);
  return x;
}
/** vector scaling: alfa*v */
inline Vector operator*(const double alfa, const Vector &v){
  Vector x(v.Dimension());
  for(int i = 0;i < v.Dimension();i++)
    x.data[i] = alfa*v(i);
  return x;
}
/** vector scaling: v*alfa */
Vector operator*(const Vector &v, const double alfa){
  Vector x(v.Dimension());
  for(int i = 0;i < v.Dimension();i++)
    x.data[i] = alfa*v(i);
  return x;
}
/** vector scaling: v/a */
Vector operator/(const Vector &v, const double alfa){
  Vector x(v.Dimension());
  for(int i = 0;i < v.Dimension();i++)
    x.data[i] = v(i)/alfa;
  return x;
}
//******************************
// General Methods
//******************************
/** initialize a vector */
inline void Vector::Initialize(int dim){
  if(dimension!=0)
    delete[] data;
    dimension = dim;
    data = new double[dimension];
    for(int i=0;i<dimension;i++)
      data[i]=0.0;
}
/** initialize a vector to a value a */
inline void Vector::Initialize(double a){
  for(int i=0;i <dimension;i++)
    data[i] = a;
}
/** initialize a vector to an array v */
inline void Vector::Initialize(double *v){
  for(int i=0;i < dimension;i++)
    data[i] = v[i];
}
/** get the dimension of a vector */
inline int Vector::Dimension() const{
  return(dimension);
}
/** print the vector */
inline void Vector::Print() const{
  cout << "\n";
  cout << "[ ";
  if(dimension>0)
    cout << data[0];
  for(int i=1;i < dimension;i++)
    cout << "; " << data[i];
  cout << "]\n";
}
/**  Norm1 */
inline double Vector::Norm_l1(){
 double sum = 0.0;
 for (int i=0;i <dimension;i++)
  sum +=fabs(data[i]);
 return(sum);
}
/** Norm2 */
inline double Vector::Norm_l2(){
 double sum = 0.0;
 for (int i=0; i<dimension;i++)
  sum += data[i]*data[i];
 return(sqrt(sum));
}
/** Norm inf */
double Vector::Norm_linf(){
  double maxval = 0.0, tmp;
  for (int i=0; i < dimension; i++){
    tmp = fabs(data[i]);
    maxval = (maxval > tmp) ? maxval:tmp;

  }
  return(maxval);

}
/** Normalize the vector */
inline void Vector::Normalize(){
  double tmp = 1.0/Norm_l2();
  for(int i=0;i < dimension;i++)
    data[i] = data[i]*tmp;
}

/***********************************/
/*  Functions Declarations        */
/*********************************/
int min_dimension(const Vector& v1, const Vector& v2);
double dot(const Vector& u, const Vector& v);
double dot(int N, double *a, double *b);
double dot(int N, const Vector& u, const Vector& v);
void Swap(double &a, double &b);
double Sign(double x);

#endif

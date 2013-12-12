/***************************************************************************
                  
                             -------------------
    begin                : Fri Jun 4 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/
#ifndef Alloc_H
#define Alloc_H
/** subroutines for allocation and deallocation of arrays 2D, 3D, 4D */
/**
 2 dimensional arrays
*/


// original implementation from Stefano
/** The allocator for 2D array */
template <class type>
inline type **_new_2_array(int sz1,int sz2,type *stupid)
{
 type **foo;
 //foo = new (type *)[sz1];
 foo = new type *[sz1];
 for (int i=0;i<sz1;i++) foo[i] = new type[sz2];
 return foo;
}
/** macro for allocate 2D array */
#define newArr(type,sz1,sz2) _new_2_array((sz1),(sz2),(type *) NULL)

/** deallocator for a 2D array*/
template <class type>
inline void delArr(type **foo,int sz1)
{ for (int i=0;i<sz1;i++) delete[] foo[i]; delete[] foo; }

/**
 3 dimensional arrays
*/

/** The allocator for 3D array */
template <class type>
inline type ***_new_3_array(int sz1,int sz2,int sz3,type *stupid)
{
 type ***foo;
 foo = new type **[sz1];
 for (int i=0;i<sz1;i++) foo[i] = newArr(type,sz2,sz3);
 return foo;
}
/** macro for allocate 3D array */
//#define newArr3(type,sz1,sz2,sz3) _new_3_array((sz1),(sz2),(sz3),(type *) NULL)

///** deallocator for a 3D array*/
//template <class type>
//inline void delArr3(type ***foo,int sz1,int sz2)
//{ for (int i=0;i<sz1;i++) delArr(foo[i],sz2); delete[] foo; }

/**
4 dimensional arrays
*/

/** The allocator for 4D array */
template <class type>
inline type ****_new_4_array(int sz1,int sz2,int sz3,int sz4,type *stupid)
{
 type ****foo;
 foo = new type ***[sz1];
 for (int i=0;i<sz1;i++) foo[i] = newArr3(type,sz2,sz3,sz4);
 return foo;
}
/** macro for allocate 4D array */
#define newArr4(type,sz1,sz2,sz3,sz4) \
_new_4_array((sz1),(sz2),(sz3),(sz4),(type *) NULL);

/** deallocator for a 4D array*/
template <class type>
inline void delArr4(type ****foo,int sz1,int sz2,int sz3)
{ for (int i=0;i<sz1;i++) delArr3(foo[i],sz2,sz3); delete[] foo; }
// end original implementation from Stefano

// implementation from Thomas Ponweiser, <Thomas.Ponweiser@risc-software.at>
template <typename T>
void allocArr1(T** arr, size_t n1) {
  *arr = new T[n1];
}

template <typename T>
void allocArr2(T*** arr, size_t n1, size_t n2) {
  *arr = new T*[n1];
  allocArr1(*arr, n1 * n2);
  for(size_t i = 1; i < n1; i++) {
    (*arr)[i] = **arr + i * n2;
  }
}

template <typename T>
void allocArr3(T**** arr, size_t n1, size_t n2, size_t n3) {
  *arr = new T**[n1];
  allocArr2(*arr, n1 * n2, n3);
  for(size_t i = 1; i < n1; i++) {
    (*arr)[i] = **arr + i * n2;
  }   
}

template <typename T>
void allocArr4(T***** arr, size_t n1, size_t n2, size_t n3, size_t n4) {
  *arr = new T***[n1];
  allocArr3(*arr, n1 * n2, n3, n4);
  for(size_t i = 1; i < n1; i++) {
    (*arr)[i] = **arr + i * n2;
  }
}

template <typename T>
void freeArr1(T** arr) {
  delete[] *arr;
  *arr = NULL;
}

template <typename T>
void freeArr2(T*** arr) {
  freeArr1(*arr);
  delete[] *arr;
  *arr = NULL;
}

template <typename T>
void freeArr3(T**** arr) {
  freeArr2(*arr);
  delete[] *arr;
  *arr = NULL;
}

template <typename T>
void freeArr4(T***** arr) {
  freeArr3(*arr);
  delete[] *arr;
  *arr = NULL;
}
// end implementation from Thomas Poiweiser

#endif

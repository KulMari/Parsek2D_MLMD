/***************************************************************************
                          FFT.h  -  Fast Fourier transform method (using FFTW library)
                             -------------------
    begin             : Thu Nov 09 2006
    developers        : Enrico Camporeale
 ***************************************************************************/


#ifndef FFT_H
#define FFT_H

#include <iostream>
#include <fftw3.h>
#include "../mathlib/Basic.h"
#include "../utility/TransArraySpace.h"

typedef  void (Field::*FIELD_IMAGE)(double*, double*, Grid*, VirtualTopology*);
typedef  void (*GENERIC_IMAGE)(double*, double*, Grid*, VirtualTopology*);
using std::cout;
using std::cerr;
using std::endl;

/** FFT solver
*
* @date Thu Nov 9 2006
* @author Enrico Camporeale
* @version 1.0
*/
			
inline void FFT (int nx, int ny, double* vector){

	
fftw_plan plan;
plan=fftw_plan_r2r_2d(nx, ny, vector, vector,FFTW_RODFT00, FFTW_RODFT00,0);
for (int i=0; i<(nx-2)*(ny-2);i++)
	cout<<"vector["<<i<<"] PRIMA = "<<vector[i]<<endl;


fftw_execute(plan);
fftw_destroy_plan(plan);
for (int i=0; i<(nx-2)*(ny-2);i++)
	cout<<"vector["<<i<<"] = "<<vector[i]<<endl;

}

#endif


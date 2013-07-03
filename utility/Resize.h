/***************************************************************************
                  
                             -------------------
    begin                : Fri Jun 4 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/
#ifndef RESIZE_H
#define RESIZE_H
/** utility for resizing an array of double */ 
inline void resize(double *buffer, int old_size, int new_size){
     double* new_buffer = new double[new_size];
     for (int i=0; i < new_size;i++)
     	new_buffer[i] = buffer[i];
     delete[] buffer;
     buffer = new_buffer;
}

#endif

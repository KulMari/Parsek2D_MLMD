/***************************************************************************
    ParallelIO.h  -  Library to manage Input-Output for a parallel environment
                       -------------------
    begin                : Fri Jun 4 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/

#ifndef ParallelIO_H
#define ParallelIO_H

#include "hdf5.h"
/* dataset data type */
typedef double DATATYPE;

inline void writePHIhdf5(char filename[], double time, VirtualTopology *vct, Grid *grid,EMfieldsIMP3D *field){
   char* temp;
   char* temp1;
   strcpy(temp,filename);
   strcat(temp,"_cycle_");
   sprintf(temp1,"%d",time);
   strcat(temp,temp1);
   strcat(temp,".hdf");
   cout << "Opening file: " << temp << endl;

   hid_t fid1;               // HF5 file IDs 
   hid_t acc_tpl1;           // File access templates
   hid_t sid1;               // Dataspace ID 
   hid_t file_dataspace;     // File dataspace ID 
   hid_t mem_dataspace;      // memory dataspace ID
   hid_t dataset1, dataset2, dataset3, datset4; // Dataset ID
   hsize_t dims1[3] = {(grid->getNXC()-2)*vct->getXLEN(),(grid->getNYC()-2)*vct->getYLEN(),(grid->getNZC()-2)*vct->getZLEN()};     /* dataspace dim sizes */

   DATATYPE data_array1[(grid->getNXC()-2)*vct->getXLEN()][(grid->getNYC()-2)*vct->getYLEN()][(grid->getNZC()-2)*vct->getZLEN()];  /* data buffer*/
   DATATYPE data_array2[(grid->getNXC()-2)*vct->getXLEN()][(grid->getNYC()-2)*vct->getYLEN()][(grid->getNZC()-2)*vct->getZLEN()];
   DATATYPE data_array3[(grid->getNXC()-2)*vct->getXLEN()][(grid->getNYC()-2)*vct->getYLEN()][(grid->getNZC()-2)*vct->getZLEN()];
   DATATYPE data_array4[(grid->getNXC()-2)*vct->getXLEN()][(grid->getNYC()-2)*vct->getYLEN()][(grid->getNZC()-2)*vct->getZLEN()];
   
   hssize_t start[3];             // for hyperslab setting
   hssize_t count[3], stride[3];  // for hyperslab setting

   herr_t ret;   // generic return value 
   
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Info info = MPI_INFO_NULL;
   /*--------------------------
   * START AN HDF5 FILE
   *--------------------------- */
   // setup file access template with parallel IO acces 
   acc_tpl1 = H5Pcreate(H5P_FILE_ACCESS);
   // set Parallel access with communicator 
   H5Pset_fapl_mpio(acc_tpl1,comm,info);
   // create the file collectively 
   fid1=H5Fcreate(temp,H5F_ACC_TRUNC,H5P_DEFAULT,acc_tpl1);
   // Release file-access template 
   ret=H5Pclose(acc_tpl1);

   /*-----------------------------
   * Define the dimension of the overall datasets
   * and the slabs local to the MPI process
   *--------------------------------*/
   // setup dimensionality object 
   sid1 =H5Screate_simple(3,dims1,NULL);
   // create a dataset collectively 
   dataset1 = H5Dcreate(fid1,"Central Grid Point X Coordinate",H5T_NATIVE_DOUBLE,sid1,H5P_DEFAULT);
   dataset2 = H5Dcreate(fid1,"Central Grid Point Y Coordinate",H5T_NATIVE_DOUBLE,sid1,H5P_DEFAULT);
   dataset3 = H5Dcreate(fid1,"Central Grid Point Z Coordinate",H5T_NATIVE_DOUBLE,sid1,H5P_DEFAULT);
   dataset4 = H5Dcreate(fid1,"Phi - Electrostatic Potential",H5T_NATIVE_DOUBLE,sid1,H5P_DEFAULT);
   // setup dimensions of the slab this process accesses
   count[0] =  (grid->getNXC()-2);
   count[1] =  (grid->getNYC()-2);
   count[2] =  (grid->getNZC()-2);
   start[0] =  vct->getCartesian_rank()*count[0];
   start[1] =  0;
   start[2] =  0;
   stride[0] =1;
   stride[1] =1;
   stride[2] =1;
   /* put grid in the data_array */
   DATATYPE *dataptr1 = &data_array1[0][0][0];
   DATATYPE *dataptr2 = &data_array2[0][0][0];
   DATATYPE *dataptr3 = &data_array3[0][0][0];
   DATATYPE *dataptr4 = &data_array4[0][0][0];
   for (int i=1;i < (grid->getNXC()-2);i++)
    for (int j=1;j < (grid->getNYC()-2);j++)
      for (int k=1;k < (grid->getNZC()-2);k++){
        *dataptr1++ = grid->getXC(i,j,k);
        *dataptr2++ = grid->getYC(i,j,k);
        *dataptr3++ = grid->getZC(i,j,k);
        *dataptr4++ = field->getPHI(i,j,k);
      }
   // create a file dataspace inpendently 
   file_dataspace = H5Dget_space(dataset1);
   H5Sselect_hyperslab(file_dataspace,H5S_SELECT_SET,start,stride,count,NULL);
   // create a memory dataspace independtly 
   mem_dataspace = H5Screate_simple(3,count,NULL);
   // write data indendently:XC
   H5Dwrite(dataset1,H5T_NATIVE_DOUBLE,mem_dataspace,file_dataspace,H5P_DEFAULT,data_array1);
   // write data indendently:YC
   H5Dwrite(dataset2,H5T_NATIVE_DOUBLE,mem_dataspace,file_dataspace,H5P_DEFAULT,data_array2);
   // write data indendently:ZC
   H5Dwrite(dataset3,H5T_NATIVE_DOUBLE,mem_dataspace,file_dataspace,H5P_DEFAULT,data_array3);
   // write data indendently:PHI
   H5Dwrite(dataset4,H5T_NATIVE_DOUBLE,mem_dataspace,file_dataspace,H5P_DEFAULT,data_array4);
   // release dataspace ID */
   H5Sclose(file_dataspace);
   // close dataset collectively
   H5Dclose(dataset1);
   H5Dclose(dataset2);
   H5Dclose(dataset3);
   H5Dclose(dataset4);
   // release all IDs created
   H5Sclose(sid1);
   // close the file collectively
   H5Fclose(fid1);
   
   
    
   
}
#endif

/********************************************************************************************
ParsekES.cpp  - A main file for running a parallel Particle-in-Cell with Electrostatic field and relativistic particles
			        -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/




// MPI
#include "mpi.h"
#include "mpidata/MPIdata.h"
// topology
#include "processtopology/VirtualTopology.h"
#include "processtopology/VCtopology.h"
// input
#include "inputoutput/CollectiveIO.h"
#include "inputoutput/Collective.h"
// grid
#include "grids/Grid.h"
#include "grids/Grid1DCU.h"
// fields
#include "fields/Field.h"
#include "fields/ESfield1D.h"
// particle
#include "particles/Particles.h"
#include "particles/Particles1Dcomm.h"
#include "particles/Particles1D.h"

//  output
#include "PSKOutput2D/PSKOutput.h"
#include "PSKOutput2D/PSKhdf5adaptor.h"
#include "inputoutput/Restart.h"
#include "inputoutput/SerialIO.h"
// performance
#include "performances/Timing.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::cerr;
using std::endl;




int main (int argc, char **argv) {
 // initialize MPI environment
 int nprocs, myrank;
 MPIdata *mpi = new MPIdata(&argc,&argv);
 nprocs = mpi->nprocs; // nprocs = number of processors
 myrank = mpi->rank; // myrank = rank of tha process
 
 Timing *my_clock = new Timing(myrank); // performance and timing object
 
 Collective *col = new Collective(argc,argv); // Every proc loads the parameters of simulation from class Collective
 bool verbose = col->getVerbose(); 
 string SaveDirName = col->getSaveDirName();
 string RestartDirName = col->getRestartDirName();
 const int restart = col->getRestart_status();
 const int ns = col->getNs(); // get the number of particle species involved in simulation
 const int first_cycle = col->getLast_cycle()+1; // get the last cycle from the restart

 VCtopology *vct = new VCtopology(); // initialize the virtual cartesian topology 
 if (nprocs != vct->getNprocs()){
    if (myrank == 0 ) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << " matrix: Change XLEN,YLEN in method VCTopology.init()" << endl;
      mpi->finalize_mpi();
      return(1);
    }
  }
  // We create a new communicator with a 3D virtual Cartesian topology
  vct->setup_vctopology(MPI_COMM_WORLD);
  // Print the initial setting about INPUT values: mpi values, processor topology values, collective values 
  if (myrank==0){
	  mpi->Print();
	  vct->Print();
	  col->Print();
  }
  Grid1DCU* grid = new Grid1DCU(col,vct);  // Create the local 1D grid
  ESfield1D* ESf= new ESfield1D(col, grid);  // Create Electrostatic field
  ESf->init(vct,grid); // initialize ES field with density constant values
  Particles1D *part = new Particles1D[ns]; // particles
  for (int i=0; i < ns; i++)
     part[i].allocate(i,col,vct,grid);
  
  
  // Initial Condition for PARTICLES if you are not starting from RESTART
  for (int i=0; i < ns; i++)
   part[i].two_beams(grid,ESf,vct); // electrons - 2 beams counterstreaming
 

  PSK::OutputManager < PSK::OutputAdaptor > output_mgr; // Create an Output Manager
  myOutputAgent< PSK::HDF5OutputAdaptor > hdf5_agent; // Create an Output Agent for HDF5 output
  hdf5_agent.set_simulation_pointers(ESf, grid, vct, mpi, col);
  for (int i=0 ; i<ns ; ++i)
    hdf5_agent.set_simulation_pointers_part(&part[i]);
  output_mgr.push_back( &hdf5_agent ); // Add the HDF5 output agent to the Output Manager's list
  if (myrank ==0 & restart<2){
          hdf5_agent.open(SaveDirName+"/settings.hdf");
          output_mgr.output("collective + total_topology + proc_topology",0);
          hdf5_agent.close();
          hdf5_agent.open(RestartDirName+"/settings.hdf");
          output_mgr.output("collective + total_topology + proc_topology",0);
          hdf5_agent.close();
  }
  stringstream ss;
  ss << myrank;
  if (restart==0){ // new simulation from input file
      hdf5_agent.open(SaveDirName+"/proc"+ss.str()+".hdf");
      output_mgr.output("proc_topology ",0);
  } else{ // restart append the results to the previous simulation 
     hdf5_agent.open_append(SaveDirName+"/proc"+ss.str()+".hdf");
	 output_mgr.output("proc_topology ",0);
  }
  MPI_Barrier( MPI_COMM_WORLD ) ;
 
  //*******************************************//
  //****     Start the  Simulation!         ***//
  //*******************************************//
 
  for (int cycle = first_cycle; cycle < (col->getNcycles()+first_cycle); cycle++){
	 if (myrank==0 && verbose){
		  cout << "***********************" << endl;
		  cout << "*   cycle = " << cycle + 1 <<"        *" << endl;
		  cout << "***********************" << endl;
          }		  
        my_clock->start_interpP2G(); // for profiling
		ESf->setZeroDensities();   
        for (int i=0; i < ns; i++) 
		  part[i].interpP2G(ESf,grid,vct); // interpolation Particle to grid
        ESf->sumOverSpecies(vct);
		MPI_Barrier(MPI_COMM_WORLD);
        ESf->interpDensitiesN2C(vct,grid);
		my_clock->stop_interpP2G();    // for profiling
		MPI_Barrier(MPI_COMM_WORLD);
 	
         // Poisson Solver
		 my_clock->start_field(); // profiling
	     ESf->calculateField(grid,vct); // calculate Lap(PHI) = 4 pi rho
	     my_clock->stop_field(); // profiling
   
	    // push particles
        my_clock->start_mover(); // profiling
	    for (int i=0; i < ns; i++)   // push each species: BACKGROUND ION: move only electrons !!!!
           part[i].mover_relativistic(grid,vct,ESf);
		   //part[i].mover_explicit(grid,vct,ESf);
	    my_clock->stop_mover(); // profiling
	    MPI_Barrier(MPI_COMM_WORLD);
		// Output to proc file
		// OUTPUT to large file, called proc**
	  if (cycle%(col->getFieldOutputCycle())==0 || cycle==first_cycle){
		  //output_mgr.output("k_energy + E_energy + B_energy",cycle);
		  output_mgr.output("Eall + rhos",cycle);
	  } 
	  if (cycle%(col->getParticlesOutputCycle())==0 && col->getParticlesOutputCycle()!=1)
         output_mgr.output("position + velocity + q ",cycle, 1);
			
		 // write to ASCII file
		 //if (cycle%50==0){
		 //  writePHIascii1D("PHI2p",myrank, cycle, grid,ESf); 
	       //writeParticles1D("elec",myrank,cycle,&part[0]);
		//   writeRHOascii1D("rho",myrank,cycle,grid,ESf);
         //}      
		  if (cycle%(col->getRestartOutputCycle())==0 && cycle != first_cycle)
	       writeRESTART_ES(RestartDirName,myrank,cycle,ns,mpi,vct,col,grid,ESf,part,0); // without ,0 add to restart file
 
  }
  
  hdf5_agent.close(); //close HDF5
  
  writeRESTART_ES(RestartDirName,myrank,(col->getNcycles()+first_cycle)-1,ns,mpi,vct,col,grid,ESf,part,0);  // write last restart
  
  my_clock->stopTiming(); // stop profiler
  
  mpi->finalize_mpi();  // close MPI
  
  return(0);
  
}





















/********************************************************************************************
ParsekINJ.cpp  - A main file for injecting particles
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
#include "grids/Grid2DCU.h"
// fields
#include "fields/Field.h"
#include "fields/EMfields.h"
// particle
#include "particles/Particles.h"
#include "particles/Particles2Dcomm.h"
#include "particles/Particles2D.h"


//  output
#include "PSKOutput2D/PSKOutput.h"
#include "PSKOutput2D/PSKhdf5adaptor.h"
#include "inputoutput/Restart.h"
// performance
#include "performances/Timing.h"
// wave
//#include "perturbation/Planewave.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
using std::cerr;
using std::endl;




int main (int argc, char **argv) {
 // initialize MPI environment
 int nprocs, myrank, mem_avail;
 MPIdata *mpi = new MPIdata(&argc,&argv);
 nprocs = mpi->nprocs; // nprocs = number of processors
 myrank = mpi->rank;   // myrank = rank of tha process (ID)
 Timing *my_clock = new Timing(myrank); // performance and timing object
 Collective *col = new Collective(argc,argv); // Every proc loads the parameters of simulation from class Collective
 bool verbose = col->getVerbose();
 string SaveDirName = col->getSaveDirName();
 string RestartDirName = col->getRestartDirName();
 const int restart = col->getRestart_status();
 const int ns = col->getNs(); // get the number of particle species involved in simulation
 const int first_cycle = col->getLast_cycle()+1; // get the last cycle from the restart
 // initialize the virtual cartesian topology 
 VCtopology *vct = new VCtopology();
 if (nprocs != vct->getNprocs()){ // Check if we can map the processes into a matrix ordering defined in Collective.cpp
    if (myrank == 0 ) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << " matrix: Change XLEN,YLEN in method VCTopology.init()" << endl;
      mpi->finalize_mpi();
      return(1);
    }
  }
  vct->setup_vctopology(MPI_COMM_WORLD); // 2D Cartesian Topology
  // Print the initial settings from the INPUT file 
  if (myrank==0){
	  mpi->Print();
	  vct->Print();
	  col->Print();
  }
  MPI_Barrier( MPI_COMM_WORLD ) ;
  Grid2DCU* grid = new Grid2DCU(col,vct); // Create the local grid
  EMfields *EMf = new EMfields(col, grid); // Create Electromagnetic Fields Object
  // Initial Condition for FIELD if you are not starting from RESTART
  EMf->initUniform(vct,grid); // initialize with zero magnetic field and constant density
  //EMf->addDipole(0.1,10.0,10.0,grid); // add a dipole with B0 = 2.0 in position 10.0, 10.0
  EMf->addIMF(-0.01,80, grid); // add a IMF magnetic field
  // Allocation of particles
  Particles2D *part = new Particles2D[ns];
  for (int i=0; i < ns; i++)
    part[i].allocate(i,col,vct,grid);
  // Initial Condition for PARTICLES if you are not starting from RESTART
  if (restart==0){
     //Planewave *wave = new Planewave(col, EMf, grid, vct);
     //wave->Wave_Rotated(part); // Single Plane Wave
     for (int i=0; i < ns; i++)
	   part[i].maxwellian(grid,EMf,vct);  // all the species have Maxwellian distribution in the velocity
  }
  // Initialize the output (simulation results and restart file)
  PSK::OutputManager< PSK::OutputAdaptor > output_mgr;  // Create an Output Manager
  myOutputAgent< PSK::HDF5OutputAdaptor > hdf5_agent; // Create an Output Agent for HDF5 output
  hdf5_agent.set_simulation_pointers(EMf, grid, vct, mpi, col);
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
  }
  else{ // restart append the results to the previous simulation 
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
          // delete particles inside the sphere of radius 2 at center (5.0,5.0)
	  //for (int i=0; i < ns; i++)
	    //part[i].deleteParticlesInsideSphere(2.0,10.0,10.0); // radius 2 in (10,10)
	  // interpolation
	  my_clock->start_interpP2G();  // for profiling 
	  EMf->setZeroDensities(); // set to zero the densities
	  for (int i=0; i < ns; i++) 
		  part[i].interpP2G(EMf,grid,vct); // interpolate Particles to Grid(Nodes)
	  //EMf->adjustNonPeriodicDensities(vct); // adjust boundaries if not periodic BC
	  EMf->sumOverSpecies();               // sum all over the species
	  MPI_Barrier(MPI_COMM_WORLD);
	  EMf->interpDensitiesN2C(vct,grid);   // calculate densities on centers from nodes
	  EMf->calculateHatFunctions(grid,vct); // calculate the hat quantities for the implicit method
	  my_clock->stop_interpP2G(); // for profiling
	  MPI_Barrier(MPI_COMM_WORLD);
 	  /// OUTPUT to large file, called proc**
	  if (cycle%(col->getFieldOutputCycle())==0 || cycle==first_cycle){
		  output_mgr.output("k_energy + E_energy + B_energy",cycle);
		  output_mgr.output("Eall + Ball + rhos + Jsall",cycle);
	  } 
	  if (cycle%(col->getParticlesOutputCycle())==0 && col->getParticlesOutputCycle()!=1)
         output_mgr.output("position + velocity + q ",cycle, 1);
	  // solve Maxwell equations
	  my_clock->start_field(); // for profiling
	  EMf->calculateField(grid,vct); // calculate the EM fields
	  my_clock->stop_field();  // for profiling
	  // push the particles
	  my_clock->start_mover();  // for profiling
	  for (int i=0; i < ns; i++) // move each species
	      // dipole with B0 =.1 in position (10,10) earth has radius 2.0
		  mem_avail = part[i].moverDipole(0.1,10.0,10.0,2.0,grid,vct,EMf); // use the Predictor Corrector scheme 
	      if (mem_avail < 0){ // not enough memory space allocated for particles: stop the simulation
		     if (myrank==0){
		         cout << "*************************************************************" << endl;
				 cout << "Simulation stopped. Not enough memory allocated for particles" << endl;
		         cout << "*************************************************************" << endl;
		      }
		      cycle = col->getNcycles()+first_cycle; // exit from the time loop
	   }
	   /////////////////////////////////////
	   ///  INJECTION  /////////////////////
	   /////////////////////////////////////
	  
	        if (myrank==0){
		         cout << "Injecting Plasma from Boundaries" << endl;
			}
			// inject plasma from boundaries
			for (int i=0; i < ns; i++){ 
		       mem_avail =  part[i].inject_Xleft_Wall(0.0,0.0,0.0,col->getUth(i),col->getVth(i),col->getWth(i),grid,vct,col->getRHOinit(i),0); 
	           mem_avail += part[i].inject_Xright_Wall(0.0,0.0,0.0,col->getUth(i),col->getVth(i),col->getWth(i),grid,vct,col->getRHOinit(i),0);
			   mem_avail += part[i].inject_Yleft_Wall(0.0,0.0,0.0,col->getUth(i),col->getVth(i),col->getWth(i),grid,vct,col->getRHOinit(i),0); 
	           mem_avail += part[i].inject_Yright_Wall(0.0,0.0,0.0,col->getUth(i),col->getVth(i),col->getWth(i),grid,vct,col->getRHOinit(i),0);
			   
			}
			// inject solar wind from X left with velocity Vinj
			for (int i=0; i < ns; i++)
		       mem_avail += part[i].inject_Xleft_Wall(col->getVinj(),0.0,0.0,0.05,0.05,0.05,grid,vct,col->getRHOinit(i),0); 
	        // check if there is enough space for particles injected
			if (mem_avail < 0)
			  cycle = col->getNcycles()+first_cycle; // exit from the time loop
	  
	   // print to screen number of particles
	   //if (myrank==0){
	     //cout << "********* Particles Info ************" << endl;
		 //int npSys=0, npB=0, npD=0;
		 //for (int i=0; i < ns; i++){
		  // npSys += part[i].getNPsystem(); npB += part[i].getNPdeletedBoundary(); npD += part[i].getNPdeletedDipole();
		 // }
		  //cout << "Np Total = " << npSys << "\t N eliminated at Boundary" << npB << "\t N eliminated at Dipole" << npD << endl;
	      //cout << "************************************" << endl;
	   //}
	   my_clock->stop_mover(); // for profiling
	   MPI_Barrier(MPI_COMM_WORLD);
           // Output save a file for the RESTART
	   if (cycle%(col->getRestartOutputCycle())==0 && cycle != first_cycle)
	     writeRESTART(RestartDirName,myrank,cycle,ns,mpi,vct,col,grid,EMf,part); // without ,0 add to restart file
	   MPI_Barrier(MPI_COMM_WORLD);
 
  }  // end of the cycle
  
  // Output: close the proc** file and write a restart file
  hdf5_agent.close();
  if (mem_avail==0) // write the restart only if the simulation finished succesfully
   writeRESTART(RestartDirName,myrank,(col->getNcycles()+first_cycle)-1,ns,mpi,vct,col,grid,EMf,part,0);  
  // stop profiling
  my_clock->stopTiming();
  // close MPI
  mpi->finalize_mpi();
  
  return(0);
  
}



















































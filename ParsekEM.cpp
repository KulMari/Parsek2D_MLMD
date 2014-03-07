/********************************************************************************************
ParsekEM.cpp  - A main file for running a parallel Particle-in-Cell with Electromagnetic field
			        -------------------
developers: Stefano Markidis, Enrico Camporeale, Giovanni Lapenta, David Burgess
********************************************************************************************/


// MPI
#include "mpi.h"
#include "mpidata/MPIdata.h"
#include "hdf5.h"
// topology
#include "processtopology/VirtualTopology.h"
#include "processtopology/VCtopology.h"
#include "processtopology/VCtopologyparticles.h"
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
#include "inputoutput/SerialIO.h"
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


//#include "epik_user.h"   //for scalasca selective instrumentation- cannot find it!

int main (int argc, char **argv) {

  bool TEST=1; //this takes a lot od reduces, put to 0 normally
  bool TEST_B=0; // to eliminate the barriers; when debugging, the barriers mark the beginning/ end of each phase
  //  bool SCALASCA_SELECTIVE=0;
  
  bool RandomizeParticlePos= true;

  int proj=1;
  int interp=1;
  int solveFields=1;
  int P2Gops=1;
  int PRASendRcvOps=1; 
  int PRACollectionMethod=1;   /* =0, coarse particle to be used for repopulation collected one by one during communicate ops; WORKING, not fatser; use 1
				  =1, coarse particles collected all together after mover; USE THIS WITH SUBCYCLING */
  int RefLevelAdj=0; /* options:                                                                      
		      --0: same adjust for coarse and refined level (multiply)-- chec with double periodic reconnection                      
		      --1: interp of OS particles for the refined level; more expensive, use it for influxes */

 // initialize MPI environment
 int nprocs, myrank, mem_avail;
 MPIdata *mpi = new MPIdata(&argc,&argv);
 nprocs = mpi->nprocs; // nprocs = number of processors
 myrank = mpi->rank;   // myrank = rank of the process (ID)

 Collective *col = new Collective(argc,argv); // Every proc loads the parameters of simulation from class Collective
 bool SubCycling= col->getSubCycling();
 //cout << "SubCycling is " << SubCycling <<endl;
 int TimeRatio= col->getTimeRatio();
 //cout << "TimeRatio is " << TimeRatio <<endl;
 bool verbose = col->getVerbose();
 string SaveDirName = col->getSaveDirName();
 string RestartDirName = col->getRestartDirName();
 int level; //Level of the grid the process is living on
 const int restart = col->getRestart_status();
 const int ns = col->getNs(); // get the number of particle species involved in simulation
 const int ngrids = col->getNgrids(); // get the number of particle species involved in simulation
 const int first_cycle = col->getLast_cycle()+1; // get the last cycle from the restart
 //double dt = col->getDt();
 // initialize the virtual cartesian topology
 VCtopology *vct = new VCtopology();
 VCtopologyparticles *vctparticles = new VCtopologyparticles();
 vct->setXLEN(col->getXLEN());
 vct->setYLEN(col->getYLEN());
 vct->setNgrids(ngrids);
 vct->setDivisions();
 vct->setNprocs();
 vct->setPeriodicity(col->getBcEMfaceXleft(),col->getBcEMfaceXright(),col->getBcEMfaceYleft(),col->getBcEMfaceYright());
 vctparticles->setXLEN(col->getXLEN());
 vctparticles->setYLEN(col->getYLEN());
 vctparticles->setNgrids(ngrids);
 vctparticles->setDivisions();
 vctparticles->setNprocs();
 vctparticles->setPeriodicity(col->getBcEMfaceXleft(),col->getBcEMfaceXright(),col->getBcEMfaceYleft(),col->getBcEMfaceYright());
 vctparticles->setRefLevelAdj(RefLevelAdj);
 vct->setRefLevelAdj(RefLevelAdj);

 if (nprocs/ngrids != vct->getNprocs()){ // Check if we can map the processes into a matrix ordering defined in Collective.cpp
    if (myrank == 0 ) {
      cerr << "Error: " << nprocs << " processes cant be mapped into " << vct->getXLEN() << "x" << vct->getYLEN() << " matrix: Change XLEN,YLEN in method VCTopology.init()" << endl;
      mpi->finalize_mpi();
      return(1);
    }
  }
  vct->setup_vctopology(MPI_COMM_WORLD); // field 2D Cartesian Topology
  vctparticles->setup_vctopology(MPI_COMM_WORLD); // particles 2D Cartesian Topology
  // initialize the central cell index
  const int nx0 = col->getNxc() / vct->getXLEN();// get the number of cells in x for each processor
  const int ny0 = col->getNyc() / vct->getYLEN();// get the number of cells in y for each processor
  // Print the initial settings from the INPUT file
  if (myrank==0){
	  mpi->Print();
	  vct->Print();
	  col->Print();
  }
  

/*** Check grid and particles are on the same level ***/
int coord[3], coord_particles[3];
MPI_Cart_coords(vct->getCART_COMM_TOTAL(), myrank, 3, coord);
MPI_Cart_coords(vctparticles->getCART_COMM_TOTAL(), myrank, 3, coord_particles);
if (coord[0] != coord_particles[0]) {
    if (myrank == 0 ) {
      cerr << "Error: " << myrank << " proces is not on the same level for grid and particles" << endl;
      mpi->finalize_mpi();
      return(1);
    }
 }
  

/******************************************************/
  level = coord[0];
  Grid2DCU* grid = new Grid2DCU(col,vct,level); // Create the grids

  EMfields *EMf = new EMfields(col, grid, vct); // Create Electromagnetic Fields Object
  
  // Initial Condition for FIELD if you are not starting from RESTART
  //EMf->initUniform(vct,grid); // initialize with constant values
  EMf->initDoubleHarris(vct,grid); // initialize with constant values 
  //EMf->initLightwave(vct,grid); // initialize with a dipole
  
  
  if (vct->getNgrids()>1 && interp)
    {
      int initBCres;
      initBCres=EMf->initWeightBC(vct,grid,col);// Initialize the weight for interpolation of BC between grids
      if (initBCres<0)
	{cout <<"problems with initWeightBC, exiting...\n";
	  mpi->Abort();
	  return (-1);
	}
    }

  if (TEST_B) 
    {
      MPI_Barrier(vct->getCART_COMM_TOTAL()) ;
      if (vct->getCartesian_rank_COMMTOTAL()==0)
	{
	  cout << "initWeightBC ended well\n";
	}
    }// end TEST
  
  if (vct->getNgrids()>1 && proj)
    {
      int initProjres;
      initProjres=EMf->initWeightProj(vct,grid,col);// Initialize the weight for projection of refined fields between grids
      if (initProjres<0)
        {cout <<"problems with initWeightProj, exiting...\n";
          mpi->Abort();
          return (-1);
        }
    }

  if (TEST_B)
    {
      MPI_Barrier(vct->getCART_COMM_TOTAL()) ;
      if (vct->getCartesian_rank_COMMTOTAL()==0)
	{
	  cout << "initWeightProj ended well\n";
	}
    }  // end TEST
  // Allocation of particles
  Particles2D *part = new Particles2D[ns];
  
  if (TEST_B)
    {
      MPI_Barrier(vct->getCART_COMM_TOTAL());
      if (vct->getCartesian_rank_COMMTOTAL()==0)
	{
	  cout << "I am starting allocating particles\n";
	}
    }// end TEST


  // end for debugging
  for (int i=0; i < ns; i++)
    {
      // init operations for particle repopulation
      part[i].setPRACollectionMethod(PRACollectionMethod);
      part[i].allocate(i,col,vctparticles,grid);
      //cout << "R"<< myrank <<"Before initPRAVariables, is " <<i <<endl;
      int PRASucc=part[i].initPRAVariables(i, col, vctparticles, grid, EMf);
      //int PRASucc=1;
      //cout << "R"<< myrank <<"After initPRAVariables, is " <<i <<endl;
      if (PRASucc<0)
	{
	  cout << "*******************************************************************************************" << endl;
	  cout << "PRA parameters need to be rechecked (most likely, overlapping PRA areas), stopping the sim)" << endl;
	  cout << "*******************************************************************************************" << endl;
	  mpi->Abort();
	  return (-1);
	}
      part[i].checkAfterInitPRAVariables(i, col, vctparticles, grid);
    }
  

  // Initial Condition for PARTICLES if you are not starting from RESTART
  if (restart==0)
  {
    for (int i=0; i < ns; i++)
      {
	//part[i].maxwellian(grid,EMf,vctparticles);  // all the species have Maxwellian distribution in the velocity
	part[i].DoubleHarris(grid,EMf,vctparticles);
	//int out;
	//out =part[i].maxwellian_sameParticleInit(grid,EMf,vctparticles);
	//if (out<0)
	//  {
	//    cout << "Error with maxwellian_sameParticleInit"<<endl;
	//    mpi->Abort();
	//    return -1;
	//    }
	//part[i].maxwellian_ball(grid,EMf,vctparticles);  // all the species have Maxwellian distribution in the velocity
      } 
  }// end restart
  // Initialize the output (simulation results and restart file)
  PSK::OutputManager< PSK::OutputAdaptor > output_mgr;  // Create an Output Manager
  myOutputAgent< PSK::HDF5OutputAdaptor > hdf5_agent; // Create an Output Agent for HDF5 output
  hdf5_agent.set_simulation_pointers(EMf, grid, vct, mpi, col);
  for (int i=0 ; i<ns ; ++i)
      hdf5_agent.set_simulation_pointers_part(&part[i]);
  output_mgr.push_back( &hdf5_agent ); // Add the HDF5 output agent to the Output Manager's list
  if (myrank ==0 & restart<2){
	  hdf5_agent.open(SaveDirName+"/settings.hdf");  // write in proc dir
	  output_mgr.output("collective + total_topology + proc_topology",0);
	  hdf5_agent.close();
	  hdf5_agent.open(RestartDirName+"/settings.hdf"); // write in restart dir
	  output_mgr.output("collective + total_topology + proc_topology",0);
	  hdf5_agent.close();
  }
  stringstream num_proc;
  num_proc << myrank;
  if (restart==0){ // new simulation from input file
    hdf5_agent.open(SaveDirName+"/proc"+num_proc.str()+".hdf");
    output_mgr.output("proc_topology ",0);
	hdf5_agent.close();
    hdf5_agent.open(SaveDirName+"/part"+num_proc.str()+".hdf");
    output_mgr.output("proc_topology ",0);
	hdf5_agent.close();
  }
  else{ // restart append the results to the previous simulation
    hdf5_agent.open_append(SaveDirName+"/proc"+num_proc.str()+".hdf");
	output_mgr.output("proc_topology ",0);
	hdf5_agent.close();
    hdf5_agent.open_append(SaveDirName+"/part"+num_proc.str()+".hdf");
	output_mgr.output("proc_topology ",0);
	hdf5_agent.close();
  }
 
  //  Conserved Quantities

  double Eenergy, Benergy, TOTKenergy, TOTmomentum;
  double *Kenergy;
  double *momentum;
  string cq;
  Kenergy = new double[ns];
  momentum = new double[ns];
  if (TEST)
    {
      stringstream levelstr;
      levelstr << level;
      cq = SaveDirName + "/ConservedQuantities_"+ levelstr.str() +".txt";
      if (vct->getCartesian_rank() == 0) { // this is the rank on the level
	ofstream my_file(cq.c_str());
	my_file.close();
      }
    }
  if (TEST_B)
    {
      cout<<"R"<<myrank << " Init ops done"<<endl;
      MPI_Barrier(vct->getCART_COMM()) ;
      if (! (vct->getCartesian_rank_COMMTOTAL()%(vct->getXLEN()*vct->getYLEN()))   )
	cout << "Level " << level << ": Init ops done\n";
  
    }//end TEST
  //*******************************************//<<
  //****     Start the  Simulation!         ***//
  //*******************************************//

  // time stamping
  Timing *my_clock = new Timing(myrank);
  bool CoarseOp;
  // end time stamping

  bool TESTSUB= false;


  for (int cycle = first_cycle; cycle < (col->getNcycles()+first_cycle); cycle++)
  {
    if (SubCycling)
      {
	CoarseOp= !(cycle % TimeRatio); //sub
	if (level == ngrids-1 and vct->getCartesian_rank() ==0)
	  {
	    cout << "Cycle " << cycle << ", CoarseOp " << CoarseOp <<endl;
	  }
      }
    else
      {
	CoarseOp= true;
      }
    //if (myrank==0 && verbose)
    if (level == ngrids-1 and vct->getCartesian_rank() ==0)
    {
      cout << "***********************" << endl;
      cout << "*   cycle = " << cycle + 1 <<"        *" << endl;
      cout << "***********************" << endl;
    }

    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;
      if (level==0 and vct->getCartesian_rank() ==0)
	cout << "level 0, barrier before P2Gops\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
        cout << "level 1, barrier before P2Gops\n";
      }

    if (P2Gops and (level or (!level and CoarseOp)))
    {
     
      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, Inside P2G\n";

      if (level==1 and vct->getCartesian_rank() ==0)
        cout << "level 1, Inside P2G\n";

 
      EMf->setZeroDensities(); // set to zero the densities
      for (int i=0; i < ns; i++)
	{
	  part[i].interpP2G(EMf,grid,vct); // interpolate Particles to Grid(Nodes)
	  if (cycle >0 && level >0 && RefLevelAdj==1)
	    {
	      int out= part[i].interpP2G_OS(EMf,grid,vct); 
	      if (out<0)
		{
		  cout << "Error with interpP2G_OS"<<endl;
		  mpi->Abort();
		  return -1;
		}
	    }
	}
      
      // adjustNonPeriodicDensities taken out of sumOverSpecies
      EMf->adjustNonPeriodicDensities(vct, cycle);

      EMf->sumOverSpecies(vct);               // correct boundaries if necessary and sum the charge densities  all over the species
      MPI_Barrier(vct->getCART_COMM());
      
      EMf->interpDensitiesN2C(vct,grid);   // calculate densities on centers from nodes

      EMf->calculateHatFunctions(grid,vct); // calculate the hat quantities for the implicit method
    }// end P2Gops      

    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;
      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, barrier after P2Gops\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
	cout << "level 1, barrier after P2Gops\n";

      }

    MPI_Barrier(vct->getCART_COMM()); //leave this here!!!
	  
    if  (CoarseOp ) // sub; be careful with timing the output
      {
	// OUTPUT to large file, called proc**
	if (cycle%(col->getFieldOutputCycle())==0 || cycle==first_cycle)
	  {
	    hdf5_agent.open_append(SaveDirName+"/proc"+num_proc.str()+".hdf");
	    output_mgr.output("k_energy + E_energy + B_energy + pressure",cycle);
	    output_mgr.output("Eall + Ball + rhos + rho + phi + Jsall",cycle);
	    hdf5_agent.close();
	    
	  }
	if (cycle%(col->getParticlesOutputCycle())==0 && col->getParticlesOutputCycle()!=1)
	  {
	    hdf5_agent.open_append(SaveDirName+"/part"+num_proc.str()+".hdf");
	    output_mgr.output("position + velocity + q +ID",cycle, 1);
	    hdf5_agent.close();
	  }// end output
      }//end sub

    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;

      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, barrier before interpolating BC\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
	cout << "level 1, barrier before interpolating BC\n";

      }

    if (interp) // sub; to be done at ALL CYCLES, also when SubCycling
    {
      if (level > 0)
      {
	//cout << "receive BC" << endl;
	//fflush(stdout);
	if (vct->getCartesian_rank() ==0)
	  cout <<"Level " << level << " Cycle " << cycle << " I am receiving BC\n";
	EMf->receiveBC(grid,vct, col, cycle % TimeRatio +1, TimeRatio);
      }
    }//end interp
    // solve Maxwell equations
    //if (solveFields)

    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;

      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, between interpolation and solver\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
	cout << "level 1, between interpolation and solver\n";

      }

    if (solveFields and  (level or (!level and CoarseOp) )) // sub
    {
      EMf->calculateField(grid,vct); // calculate the EM fields
    }

    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;
      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, after solver\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
	cout << "level 1, barrier after solver\n";
      
      }

    //If needed: receive or send  boundary conditions to other levels
    //if (interp)
    if (interp and CoarseOp)  //end
    {
      if (level < ngrids - 1 && ngrids >1)
      {
	if (vct->getCartesian_rank() ==0)
          cout <<"Level " << level << " Cycle " << cycle << "I am sending BC\n";
	EMf->sendBC(grid,vct);
      }
    }//end interp


    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;

      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, barrier after sending BC\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
	cout << "level 1, barrier after sending BC\n";

      }

    //If needed: receive or send refined fields from/to finer/coarser levels 
    if (proj and (!SubCycling or (SubCycling and (cycle % (TimeRatio)== TimeRatio-1)  )  )  ) //be careful here!!!
    {
      if (level > 0 )
      {
	if (vct->getCartesian_rank() ==0)
          cout <<"Level " << level << " Cycle " << cycle << "I am sending projection\n";
	EMf->sendProjection(grid,vct);
      }
    }
    
    // to avoid deadlock, if Subcycling collect particles here
    int mem_avail1;
    if (SubCycling and level==0 and CoarseOp)
      {
	if (vct->getCartesian_rank() ==0)
          cout <<"Level " << level << " Cycle " << cycle << "I am collecting and sending particles to the refined grid, BEFORE RECEIVING THE PROJ (this are the old particles)\n";
	for (int i=0; i < ns; i++) // species                                                                       
	  {
            //cout << "Collection Method: " << part[i].getPRACollectionMethod()<< endl;                             
            if (part[i].getPRACollectionMethod()==1)
              {
                int resCollMethod;
                resCollMethod= part[i].CollectivePRARepopulationAdd(vctparticles, grid);
                if (resCollMethod<0)
                  {
                    cout << "CollectivePRARepopulationAdd failed, check buffers\n";
                    mpi->Abort();
                    return -1;
                  }
	      }
	    mem_avail1=2;
	    mem_avail1=part[i].PRASend(grid, vctparticles);
	    if (mem_avail1 <0)
	      { // just one of the 2 is tested, according to the level                                              
		fflush(stdout);
		cout << "**************************************************************************************" << endl;
		cout <<"R" << myrank <<"L" << level << endl;
		cout << "Simulation stopped. Not enough memory allocated for particles after PRA send or receive" << endl;
		cout << "***************************************************************************************" << endl;
		//cycle = col->getNcycles()+1;                                                                      
		mpi->Abort();
		return (-1);
	      }
	   
	  }//end bracket on species          
      }// end SubCYling



    if (proj and CoarseOp)
      {
      if (level < ngrids - 1)
      {
	if (!SubCycling or (SubCycling and cycle+ TimeRatio < col->getNcycles()+first_cycle )) // otherwise stuck because the coarse grid is still waiting for a proj the refined grid won't do because it's over already

	  {	
	    if (vct->getCartesian_rank() ==0)
	      cout <<"Level " << level << " Cycle " << cycle << "I am receiving projection\n";
	    EMf->receiveProjection(col,grid,vct);
	  }// end if Sub
      }
    }// end projection
    //if ((cycle%(col->getFieldOutputCycle())==0 || cycle==first_cycle) and TEST){
    //EMf->outputghost(vct, col, cycle);
    //	}
    // Here we made the assumption that level n can be updated by non updated level n+1.
      
    if (TESTSUB==true)
      {
	MPI_Barrier(vct->getCART_COMM()) ;
	if (level==0 and vct->getCartesian_rank() ==0)
	  cout << "level 0, barrier after proj\n";
	
	if (level==1 and vct->getCartesian_rank() ==0)
	  cout << "level 1, barrier after proj\n";
      }
    if  (level or (!level and CoarseOp) )// sub
      {
	for (int i=0; i < ns; i++) // move each species
	  {
	    // initialize the buffers for communicatin of PRA particles from the coarser to the finer grids
	    part[i].initPRABuffers (grid, vctparticles);
	    
	    mem_avail = part[i].mover_PC(grid,vctparticles,EMf); // use the Predictor Corrector scheme
	    if (mem_avail < 0)
	      { // not enough memory space allocated for particles: stop the simulation		
		cout << "*************************************************************" << endl;
		cout <<"R" << myrank <<"L" << level << "Simulation stopped. Not enough memory allocated for particles" << endl;
		cout << "*************************************************************" << endl;
		mpi->Abort();
		return (-1);	
		//cycle = col->getNcycles()+1; // exit from the time loop
	      }// end bracket memavail<0
	  }// end species
      } // end sub

    if (!SubCycling and level==0 and CoarseOp) // this is done here only if I am not subcyclign
      {
	for (int i=0; i < ns; i++) // species
	{	    
	    //cout << "Collection Method: " << part[i].getPRACollectionMethod()<< endl;
	    if (part[i].getPRACollectionMethod()==1)
	      {
		int resCollMethod;
		resCollMethod= part[i].CollectivePRARepopulationAdd(vctparticles, grid);
		if (resCollMethod<0)
		  {
		    cout << "CollectivePRARepopulationAdd failed, check buffers\n";
		    mpi->Abort();
		    return -1;
		  }
	      }
	  }//end bracket on species
      } // end sub

    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;
      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, barrier before PRA S/R\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
	cout << "level 1, barrier before PRA S/R\n";
      }
    
    mem_avail1=2; 
    int mem_avail2=2;
    // communicate PRA particles for all species
    //if(PRASendRcvOps)
    //if(PRASendRcvOps and (level or (!level and CoarseOp)     ) )
    if(PRASendRcvOps)
    {
      for (int i=0; i < ns; i++) 
      {
	if (CoarseOp)
	  {
	    if (level < ngrids - 1 and !SubCycling) // this is done here only if I am not SubCyclign
	      {
		mem_avail1=part[i].PRASend(grid, vctparticles);
	      }
	    
	    if (level>0)
	      {
		mem_avail2= part[i].PRAReceive(grid,vctparticles,EMf);
	      }
	    
	    if (mem_avail1 < 0 || mem_avail2<0)
	      { // just one of the 2 is tested, according to the level	 
		fflush(stdout);
		cout << "**************************************************************************************" << endl;
		cout <<"R" << myrank <<"L" << level << endl;
		cout << "Simulation stopped. Not enough memory allocated for particles after PRA send or receive" << endl;
		cout << "***************************************************************************************" << endl;
		//cycle = col->getNcycles()+1; 
		mpi->Abort();
		return (-1);
	      }
	  } // end CoarseOP
	if (level and !CoarseOp) // this has to be executed on the refined grid only when PRAReceive is NOT executed
	  {
	    if (vct->getCartesian_rank()==0)
	      {
		cout << "Cycle " <<cycle << ": I am executing SubCyclingParticles\n";
	      }
	    int RPRes= part[i].SubCyclingParticles(grid, CoarseOp,vctparticles);
	    if (RPRes==-1)
	      {
		cout <<"Simulation stopped, problems with SubCylingParticles\n";
		mpi->Abort();
		return (-1);
	      }
	    /*if (RandomizeParticlePos)
	      part[i].RandomizePositionPRAParticles(i, vctparticles, grid);*/
	  }

	if (RandomizeParticlePos and level)
	  {
	    part[i].RandomizePositionPRAParticles(i, vct, grid);
	    }
      }// end species
    }// end of PRASendRcvOps
    
    if (TESTSUB==true)
      {
      MPI_Barrier(vct->getCART_COMM()) ;
      if (level==0 and vct->getCartesian_rank() ==0)
        cout << "level 0, barrier after PRA S/R\n";
      
      if (level==1 and vct->getCartesian_rank() ==0)
	cout << "level 1, barrier after PRA S/R\n";
      }

    // Output save a file for the RESTART
      if  (CoarseOp ) // sub; be careful with timing the output
	{
	  if (cycle%(col->getRestartOutputCycle())==0 && cycle != first_cycle)
	    writeRESTART(RestartDirName,myrank,cycle,ns,mpi,vct,col,grid,EMf,part,0); // without ,0 add to restart file	   
	}

      //if (TEST)// and !(cycle%20 ))
      //if (TEST and (level or (!level and CoarseOp) )) // sub; be careful with timing the output
      if (TEST and CoarseOp and (cycle%(col->getFieldOutputCycle())==0 or cycle==first_cycle))
      {
	Eenergy= EMf->getEenergy(vct);
	Benergy= EMf->getBenergy(vct);
	TOTKenergy=0.0;
	TOTmomentum=0.0;
	for (int is=0; is<ns; is++)
	  {
	    Kenergy[is]=part[is].getKenergy(vctparticles);
	    TOTKenergy+=Kenergy[is];
	    momentum[is]=part[is].getP(vctparticles);
	    TOTmomentum+=momentum[is];
	    }

	if (vct->getCartesian_rank() == 0) {
	  ofstream my_file(cq.c_str(), fstream::app);
	  
	  my_file << cycle << "\t" << "\t" << (Eenergy + Benergy + TOTKenergy) << "\t" << TOTmomentum << "\t" << Eenergy << "\t" << Benergy << "\t" << TOTKenergy <<"\t";
	  for (int is=0; is<ns; is++)
	    {
	      my_file <<Kenergy[is] <<"\t" <<momentum[is] <<"\t" ;
	      }
	  my_file <<endl;
	  my_file.close();
	 
	}
	
      }// end TEST


  }  // end of the cycle
  if (mem_avail==0) // write the restart only if the simulation finished succesfully
    writeRESTART(RestartDirName,myrank,(col->getNcycles()+first_cycle)-1,ns,mpi,vct,col,grid,EMf,part,0);

    // time stamping
    my_clock->stopTiming();
  // end time stamping

   // close MPI
  mpi->finalize_mpi();

  delete[] Kenergy;
  delete[] momentum;
  
  return(0);

}


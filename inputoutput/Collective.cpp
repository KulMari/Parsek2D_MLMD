/***************************************************************************
                          Collective.h  -  description
                             -------------------
    begin                : Wed Jun 2 2004
    copyright            : (C) 2004 Los Alamos National Laboratory
    developers           : Stefano Markidis, Giovanni Lapenta
    email                : markidis@lanl.gov, lapenta@lanl.gov
 ***************************************************************************/



/**
*  Collective properties. Input physical parameters for the simulation.
*
* @date July 2013
* @par Copyright:
* (C) KULeuven
* @author Maria Elena Innocenti, Pierre Henry, Stefano Markidis, Giovanni Lapenta
* @version 1.0
*/
#include <iostream>
#include <string.h>
#include <stdlib.h>

#include <math.h>
#include "Collective.h"
// use hdf5 for the restart file
#include "hdf5.h"

//#include "../ConfigFile/src/ConfigFile.h"
//#include "../ConfigFile/src/input_array.h"

using std::cout;
using std::endl;
/** Read the input file from text file and put the data in a collective wrapper:
    if it's a restart read from input file basic sim data and load particles and EM field
    from restart file */
void Collective::ReadInput(string inputfile){
  using namespace std;
  int test_verbose;

  //Loading the input file 
  ConfigFile config(inputfile);
  
   // the following variables are ALWAYS taken from inputfile, even if restarting
   {
	dt = config.read<double>( "dt" );
	ncycles = config.read<int>( "ncycles" );
	th = config.read<double>( "th" );
	config.readInto(Smooth, "Smooth" ) ;

	Nvolte = config.read<int>( "Nvolte" ); // MP kv_file.get_data( "Nvolte", Nvolte );
	XLEN = config.read<int>( "XLEN" ); // MP kv_file.get_data( "XLEN", XLEN );
	YLEN = config.read<int>( "YLEN" ); // MP kv_file.get_data( "YLEN", YLEN );
	ngrids = config.read<int>( "ngrids" ); // MP kv_file.get_data( "ngrids", ngrids );
	ratio = config.read<int>( "ratio" ); // MP kv_file.get_data( "ratio", ratio );
	
	SaveDirName = config.read<string>( "SaveDirName" );
	RestartDirName = config.read<string>( "RestartDirName" );
	ns = config.read<int>( "ns" );
	NpMaxNpRatio = config.read<double>( "NpMaxNpRatio" );
	// GEM Challenge
	B0x = config.read<double>( "B0x" );
	B0y = config.read<double>( "B0y" );
	B0z = config.read<double>( "B0z" );
	delta = config.read<double>( "delta" );

	rhoINIT = new double[ns];
	array_double rhoINIT0 = config.read<array_double>( "rhoINIT" );
	rhoINIT[0]=rhoINIT0.a;
	if (ns > 1)
	  rhoINIT[1]=rhoINIT0.b;
	if (ns > 2)
	  rhoINIT[2]=rhoINIT0.c;
	if (ns > 3)
	  rhoINIT[3]=rhoINIT0.d;
	if (ns > 4)
	  rhoINIT[4]=rhoINIT0.e;
	if (ns > 5)
	  rhoINIT[5]=rhoINIT0.f;
	// take the tolerance of the solvers
	CGtol = config.read<double>( "CGtol" );
	GMREStol = config.read<double>( "GMREStol" );
	NiterMover = config.read<int>( "NiterMover" );
	
	// take the injection of the particless
	Vinj = config.read<double>( "Vinj" );

	// take the output cycles
	FieldOutputCycle = config.read<int>( "FieldOutputCycle" );
	ParticlesOutputCycle = config.read<int>( "ParticlesOutputCycle" );
	RestartOutputCycle = config.read<int>( "RestartOutputCycle" );

	//AMR variables, ME       	
	PRA_Xleft = config.read<int>( "PRA_Xleft" ); 
	PRA_Xright = config.read<int>( "PRA_Xright" ); 
	PRA_Yleft = config.read<int>( "PRA_Yleft" ); 
	PRA_Yright = config.read<int>( "PRA_Yright" ); 
	
    }

    if (RESTART1){    // you are restarting
      RestartDirName = config.read<string>( "RestartDirName" );
       ReadRestart(RestartDirName);
    }
    else
      { // this is not done if restarting
       restart_status=0;
       last_cycle=-1;
       c = config.read<double>( "c" );
       Lx = config.read<double>( "Lx" );
       Ly = config.read<double>( "Ly" );
       nxc = config.read<int>( "nxc" );
       nyc = config.read<int>( "nyc" );

       npcelx = new int[ns];
       npcely = new int[ns];
       qom = new double[ns];
       uth = new double[ns];
       vth = new double[ns];
       wth = new double[ns];
       u0 = new double[ns];
       v0 = new double[ns];
       w0 = new double[ns];

       array_int npcelx0 = config.read<array_int>( "npcelx" );
       array_int npcely0 = config.read<array_int>( "npcely" );
       array_double qom0 = config.read<array_double>( "qom" );
       array_double uth0 = config.read<array_double>( "uth" );
       array_double vth0 = config.read<array_double>( "vth" );
       array_double wth0 = config.read<array_double>( "wth" );
       array_double u00 = config.read<array_double>( "u0" );
       array_double v00 = config.read<array_double>( "v0" );
       array_double w00 = config.read<array_double>( "w0" );

       npcelx[0]=npcelx0.a;
       npcely[0]=npcely0.a;
       qom[0]=qom0.a;
       uth[0]=uth0.a;
       vth[0]=vth0.a;
       wth[0]=wth0.a;
       u0[0]=u00.a;
       v0[0]=v00.a;
       w0[0]=w00.a;

       if (ns > 1){
	 npcelx[1]=npcelx0.b;
	 npcely[1]=npcely0.b;
	 qom[1]=qom0.b;
	 uth[1]=uth0.b;
	 vth[1]=vth0.b;
	 wth[1]=wth0.b;
	 u0[1]=u00.b;
	 v0[1]=v00.b;
	 w0[1]=w00.b;
       }

       if (ns > 2){
	 npcelx[2]=npcelx0.c;
	 npcely[2]=npcely0.c;
	 qom[2]=qom0.c;
	 uth[2]=uth0.c;
	 vth[2]=vth0.c;
	 wth[2]=wth0.c;
	 u0[2]=u00.c;
	 v0[2]=v00.c;
	 w0[2]=w00.c;
       }

       if (ns > 3){
	 npcelx[3]=npcelx0.d;
	 npcely[3]=npcely0.d;
	 qom[3]=qom0.d;
	 uth[3]=uth0.d;
	 vth[3]=vth0.d;
	 wth[3]=wth0.d;
	 u0[3]=u00.d;
	 v0[3]=v00.d;
	 w0[3]=w00.d;
       }

       if (ns > 4){
	 npcelx[4]=npcelx0.e;
	 npcely[4]=npcely0.e;
	 qom[4]=qom0.e;
	 uth[4]=uth0.e;
	 vth[4]=vth0.e;
	 wth[4]=wth0.e;
	 u0[4]=u00.e;
	 v0[4]=v00.e;
	 w0[4]=w00.e;
       }

       if (ns > 5){
	 npcelx[5]=npcelx0.f;
	 npcely[5]=npcely0.f;
	 qom[5]=qom0.f;
	 uth[5]=uth0.f;
	 vth[5]=vth0.f;
	 wth[5]=wth0.f;
	 u0[5]=u00.f;
	 v0[5]=v00.f;
	 w0[1]=w00.f;
       }

       verbose = config.read<bool>( "verbose" );
       
       // PHI Electrostatic Potential  
       bcPHIfaceXright = config.read<int>( "bcPHIfaceXright" );
       bcPHIfaceXleft  = config.read<int>( "bcPHIfaceXleft" );
       bcPHIfaceYright = config.read<int>( "bcPHIfaceYright" );
       bcPHIfaceYleft  = config.read<int>( "bcPHIfaceYleft" );

       // EM field boundary condition 
       bcEMfaceXright = config.read<int>( "bcEMfaceXright" );
       bcEMfaceXleft =  config.read<int>( "bcEMfaceXleft" );
       bcEMfaceYright = config.read<int>( "bcEMfaceYright" );
       bcEMfaceYleft =  config.read<int>( "bcEMfaceYleft" );

       // Particles Boundary condition                                                  
       bcPfaceXright = config.read<int>( "bcPfaceXright" );
       bcPfaceXleft =  config.read<int>( "bcPfaceXleft" );
       bcPfaceYright = config.read<int>( "bcPfaceYright" );
       bcPfaceYleft =  config.read<int>( "bcPfaceYleft" );
      } // end this is not done if restarting

    TrackParticleID =  new bool[ns];
    array_bool TrackParticleID0 = config.read<array_bool>( "TrackParticleID" );
    TrackParticleID[0]=TrackParticleID0.a;
    if (ns > 1)
      TrackParticleID[1]=TrackParticleID0.b;
    if (ns > 2)
      TrackParticleID[2]=TrackParticleID0.c;
    if (ns > 3)
      TrackParticleID[3]=TrackParticleID0.d;
    if (ns > 4)
      TrackParticleID[4]=TrackParticleID0.e;
    if (ns > 5)
      TrackParticleID[5]=TrackParticleID0.f;
}

/**
*  Read the collective information from the RESTART file in HDF5 format
*
*  There are three restart status:
*  restart_status = 0 --->  new inputfile
*  restart_status = 1 --->  RESTART and restart and result directories does not coincide
*  restart_status = 2 --->  RESTART and restart and result directories coincide
*
*/
int Collective::ReadRestart(string inputfile){

  /*  working on the restart case right now...
  cout << "Restart has never been properly adapted to the MLMD case... Exiting..." << endl;
  cerr << "Restart has never been properly adapted to the MLMD case... Exiting..." << endl;
  return -1;
  */
  
  restart_status=1;
  // hdf stuff
  hid_t    file_id;
  hid_t    dataset_id;
  herr_t   status;
  /*
   * Open the  setting file for the restart.
   */
  file_id = H5Fopen((inputfile+"/settings.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0){
    cout << "couldn't open file: " << inputfile << endl;	\
    return -1;}
  
  // read c
  dataset_id = H5Dopen1(file_id, "/collective/c");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&c);
  status = H5Dclose(dataset_id);
  
  // read Lx
  dataset_id = H5Dopen1(file_id, "/collective/Lx");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
  status = H5Dclose(dataset_id);
  // read Ly
  dataset_id = H5Dopen1(file_id, "/collective/Ly");
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
  status = H5Dclose(dataset_id);
  // read nxc
  dataset_id = H5Dopen1(file_id, "/collective/Nxc");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nxc);
  status = H5Dclose(dataset_id);
  // read nyc
  dataset_id = H5Dopen1(file_id, "/collective/Nyc");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&nyc);
  status = H5Dclose(dataset_id);
  // read ns
  dataset_id = H5Dopen1(file_id, "/collective/Ns");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ns);
  status = H5Dclose(dataset_id);
  // read ngrids
  dataset_id = H5Dopen1(file_id, "/collective/Ngrids");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&ngrids);
  status = H5Dclose(dataset_id);
  
  /** Boundary condition information */
  // read EMfaceXleft
  dataset_id = H5Dopen1(file_id, "/collective/bc/EMfaceXleft");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcEMfaceXleft);
  status = H5Dclose(dataset_id);
  // read EMfaceXright
  dataset_id = H5Dopen1(file_id, "/collective/bc/EMfaceXright");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcEMfaceXright);
  status = H5Dclose(dataset_id);
  // read EMfaceYleft
  dataset_id = H5Dopen1(file_id, "/collective/bc/EMfaceYleft");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcEMfaceYleft);
  status = H5Dclose(dataset_id);
  // read EMfaceYright
  dataset_id = H5Dopen1(file_id, "/collective/bc/EMfaceYright");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcEMfaceYright);
  status = H5Dclose(dataset_id);
  
  // read PHIfaceXleft
  dataset_id = H5Dopen1(file_id, "/collective/bc/PHIfaceXleft");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPHIfaceXleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceXright
  dataset_id = H5Dopen1(file_id, "/collective/bc/PHIfaceXright");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPHIfaceXright);
  status = H5Dclose(dataset_id);
  // read PHIfaceYleft
  dataset_id = H5Dopen1(file_id, "/collective/bc/PHIfaceYleft");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPHIfaceYleft);
  status = H5Dclose(dataset_id);
  // read PHIfaceYright
  dataset_id = H5Dopen1(file_id, "/collective/bc/PHIfaceYright");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPHIfaceYright);
  status = H5Dclose(dataset_id);
  
  // read PfaceXleft
  dataset_id = H5Dopen1(file_id, "/collective/bc/PfaceXleft");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPfaceXleft);
  status = H5Dclose(dataset_id);
  // read PfaceXright
  dataset_id = H5Dopen1(file_id, "/collective/bc/PfaceXright");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPfaceXright);
  status = H5Dclose(dataset_id);
  // read PfaceYleft
  dataset_id = H5Dopen1(file_id, "/collective/bc/PfaceYleft");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPfaceYleft);
  status = H5Dclose(dataset_id);
  // read PfaceYright
  dataset_id = H5Dopen1(file_id, "/collective/bc/PfaceYright");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&bcPfaceYright);
  status = H5Dclose(dataset_id);
  
  
  // allocate fields depending on species
  npcelx = new int[ns];
  npcely = new int[ns];
  qom = new double[ns];
  uth = new double[ns];
  vth = new double[ns];
  wth = new double[ns];
  u0 = new double[ns];
  v0 = new double[ns];
  w0 = new double[ns];
  
  // read data  from species0, species 1, species2,...
  string* name_species = new string[ns];
  
  stringstream *ss = new stringstream[ns];
  
  for (int i=0;i <ns;i++){
    ss[i] << i;
    name_species[i] = "/collective/species_"+ ss[i].str() +"/";
    
  }
  
  // npcelx for different species
  for(int i=0;i < ns;i++){
    dataset_id = H5Dopen1(file_id, (name_species[i]+"Npcelx").c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&npcelx[i]);
    status = H5Dclose(dataset_id);
  }
  
  // npcely for different species
  for(int i=0;i < ns;i++){
    dataset_id = H5Dopen1(file_id, (name_species[i]+"Npcely").c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&npcely[i]);
    status = H5Dclose(dataset_id);
  }
  
  // qom for different species
  for(int i=0;i < ns;i++){
    dataset_id = H5Dopen1(file_id, (name_species[i]+"qom").c_str());
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&qom[i]);
    status = H5Dclose(dataset_id);
  }
  /** not needed for restart **/
  for (int i=0; i<ns; i++)
    uth[i]=0.0;
  for (int i=0; i<ns; i++)
    vth[i]=0.0;
  for (int i=0; i<ns; i++)
    wth[i]=0.0;
  for (int i=0; i<ns; i++)
    u0[i]=0.0;
  for (int i=0; i<ns; i++)
    v0[i]=0.0;
  for (int i=0; i<ns; i++)
    w0[i]=0.0;
  // verbose on
  verbose = 1;
  
  if (RestartDirName == SaveDirName){
    restart_status=2;
    // read dt
    dataset_id = H5Dopen1(file_id, "/collective/Dt");
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&dt);
    status = H5Dclose(dataset_id);
    // read th
    dataset_id = H5Dopen1(file_id, "/collective/Th");
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&th);
    status = H5Dclose(dataset_id);
    // read Smooth
    dataset_id = H5Dopen1(file_id, "/collective/Smooth");
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Smooth);
    status = H5Dclose(dataset_id);
    // read Nvolte
    dataset_id = H5Dopen1(file_id, "/collective/Nvolte");
    status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Nvolte);
    status = H5Dclose(dataset_id);
  }
  
  status = H5Fclose(file_id);
  delete[] name_species;
  
  // read last cycle (not from settings, but from restart0.hdf)
  
  file_id = H5Fopen((inputfile+"/restart0.hdf").c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id < 0){
    cout << "couldn't open file: " << inputfile << endl;	\
    return -1;}
  
  dataset_id = H5Dopen1(file_id, "/last_cycle");
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,&last_cycle);
  status = H5Dclose(dataset_id);
  status = H5Fclose(file_id);

  return 0;
  
}

Collective::Collective(int argc, char** argv) {

  if( argc < 2 ){
   inputfile = "inputfile";
   RESTART1=false;}
  else if ( argc < 3 ){
    inputfile = argv[1];
    RESTART1=false;}
  else{
      if( strcmp(argv[1],"restart")==0 ){
        inputfile = argv[2];
        RESTART1=true;}
       else if ( strcmp(argv[2],"restart")==0  ){
         inputfile = argv[1];
         RESTART1=true;}
   else{
    cout<<"Error: syntax error in mpirun arguments. Did you mean to 'restart' ?"<<endl;
    return;}
  }

    ReadInput(inputfile);

    /** fourpi = 4 greek pi */
    fourpi = 16.0*atan(1.0);

    /** dx = space step - X direction */
    dx = Lx/ (double) nxc;
    /** dy = space step - Y direction */
    dy = Ly/ (double) nyc;

    /** npcelx = number of particles per cell - Direction X
    <ul>
    	<li> npcelx[0] = for species 0</li>
    	<li> npcelx[1] = for species 1</li>
    </ul>

    */
    /** npcel = number of particles per cell */
    /** np = number of particles of different species */
    /** npMax = maximum number of particles of different species */

    npcel = new int[ns];
    np = new int[ns];
    npMax = new int[ns];

    for (int i=0; i<ns; i++){
      npcel[i] = npcelx[i]*npcely[i];
      np[i] = npcel[i]*nxc*nyc;
      npMax[i] = (int) (NpMaxNpRatio*np[i]/ (XLEN*YLEN) );  // without this, risks to go out of boundaries for int
      //npMax[i] = (int) (NpMaxNpRatio*np[i]);
      //Debug6
      /*cout <<"Debug6, species " << i <<endl;
      cout <<"npcelx[i] "<<npcelx[i] <<" npcely[i] " << npcely[i] << " npcel[i] " << npcel[i];
      cout <<"nxc " << nxc << " nyc " << nyc << " np[i] " << np[i] << endl;
      cout <<"NpMaxNpRatio " << NpMaxNpRatio << " npMax[i] "<< npMax[i] <<endl;*/
     }



}
/** destructor */
Collective::~Collective(){
  delete[] np;
  delete[] npcel;
  delete[] npcelx;
  delete[] npcely;
  delete[] npMax;
  delete[] qom;

  delete[] uth;
  delete[] vth;
  delete[] wth;

  delete[] u0;
  delete[] v0;
  delete[] w0;

  delete[] rhoINIT;

}
/** Print Simulation Parameters */
void Collective::Print(){
  cout << endl;
  cout << "Simulation Parameters" << endl;
  cout << "---------------------" << endl;
  cout << "Number of species    = " << ns << endl;
  for (int i=0; i < ns;i++)
    cout << "Number of superparticles of species " << i << " = " << np[i] << "\t (MAX = " << npMax[i] << ")" << endl;

  cout << "x-Length                 = " << Lx      << endl;
  cout << "y-Length                 = " << Ly      << endl;
  cout << "Number of grids          = " << ngrids  << endl;
  cout << "Number of cells (x)      = " << nxc     << endl;
  cout << "Number of cells (y)      = " << nyc     << endl;
  cout << "Number of processors in X = " << XLEN << endl;
  cout << "Number of processors in Y = " << YLEN << endl;
  cout << "Number of grids = " << ngrids << endl;
  cout << "Ratio between the grids = " << ratio << endl;
  cout << "Time step                = " << dt      << endl;
  cout << "Number of cycles         = " << ncycles << endl;
  cout << "Results saved in: "<< SaveDirName <<endl;
  cout << "---------------------" << endl;
  cout << "Check Simulation Constraints" << endl;
  cout << "---------------------" << endl;
  cout << "Accuracy Constraint:  " << endl;
  for (int i=0; i < ns;i++){
    cout << "u_th < dx/dt species " << i << "....." ;
    if (uth[i] < (dx/dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;

    cout << "v_th < dy/dt species " << i << "......" ;
    if (vth[i] < (dy/dt))
      cout << "OK" << endl;
    else
      cout << "NOT SATISFIED. STOP THE SIMULATION." << endl;
  }
  cout << endl;
  cout << "Finite Grid Stability Constraint:  ";
  cout << endl;
  for (int is=0; is < ns;is++){
   if (uth[is]*dt/dx > .1)
     cout << "OK u_th*dt/dx (species " << is << ") = "  << uth[is]*dt/dx << " > .1" << endl;
   else
     cout << "WARNING.  u_th*dt/dx (species "<< is <<") = "  << uth[is]*dt/dx << " < .1" << endl;

   if (vth[is]*dt/dy > .1)
     cout << "OK v_th*dt/dy (species " << is <<") = " << vth[is]*dt/dy << " > .1" << endl;
   else
     cout << "WARNING. v_th*dt/dy (species "<< is << ") = "<< vth[is]*dt/dy << " < .1" << endl;

  }
  cout << "PRA dimensions: " <<endl;
  cout << "PRA_Xleft: " << PRA_Xleft << ", PRA_Xright: "<< PRA_Xright << "PRA_Yleft: "<< PRA_Yleft << ", PRA_Yright: "<< PRA_Yright << endl;

  
}
/** get the physical space dimensions            */
int Collective::getDim(){
 return(dim);
}
/** get Lx */
double Collective::getLx(){
 return(Lx);
}
/** get Ly */
double Collective::getLy(){
 return(Ly);
}
/** get nxc */
int Collective::getNxc(){
 return(nxc);
}
/** get nyx */
int Collective::getNyc(){
  return(nyc);
}
/** get dx */
double Collective::getDx(){
  return(dx);
}
/** get dy */
double Collective::getDy(){
  return(dy);
}
/** get the light speed */
inline double Collective::getC(){
 return(c);
}
/** get the time step */
double Collective::getDt(){
 return(dt);
}
/** get the decentering parameter */
double Collective::getTh(){
 return(th);
}
/** get the smoothing parameter */
double Collective::getSmooth(){
	return(Smooth);
}
/** get the smoothing times */
int Collective::getNvolte(){
	return(Nvolte);
}
/** get the number of time cycles */
int Collective::getNcycles(){
 return(ncycles);
}
/** get the number of grids */
int Collective::getNgrids(){
 return(ngrids);
}
/** get the number of processors - x direction */
int Collective::getXLEN(){
 return(XLEN);
}
/** get the number of processors - x direction */
int Collective::getYLEN(){
 return(YLEN);
}
/** get the size ratio between grids */
int Collective::getRatio(){
 return(ratio);
}
/** get the number of species */
int Collective::getNs(){
 return(ns);
}
/** get the number of particles per cell for species nspecies */
int Collective::getNpcel(int nspecies){
 return(npcel[nspecies]);
}
/** get the number of particles per cell for species nspecies - direction X */
int Collective::getNpcelx(int nspecies){
 return(npcelx[nspecies]);
}
/** get the number of particles per cell for species nspecies - direction Y */
int Collective::getNpcely(int nspecies){
 return(npcely[nspecies]);
}

/** get the number of particles for different species */
int Collective::getNp(int nspecies){
 return(np[nspecies]);
}
/** get maximum number of particles for different species */
int Collective::getNpMax(int nspecies){
 return(npMax[nspecies]);
}
double Collective::getNpMaxNpRatio(){
  return(NpMaxNpRatio);
}
/** get charge to mass ratio for different species */
double Collective::getQOM(int nspecies){
 return(qom[nspecies]);
}
/** get the background density for GEM challenge */
double Collective::getRHOinit(int nspecies){
 return(rhoINIT[nspecies]);
}
/** get thermal velocity  - Direction X     */
double Collective::getUth(int nspecies){
 return(uth[nspecies]);
}
/** get thermal velocity  - Direction Y     */
double Collective::getVth(int nspecies){
 return(vth[nspecies]);
}
/** get thermal velocity  - Direction Z     */
double Collective::getWth(int nspecies){
return(wth[nspecies]);
}
/** get beam velocity - Direction X        */
double Collective::getU0(int nspecies){
 return(u0[nspecies]);
}
/** get beam velocity - Direction Y        */
double Collective::getV0(int nspecies){
 return(v0[nspecies]);
}
/** get beam velocity - Direction Z        */
double Collective::getW0(int nspecies){
return(w0[nspecies]);
}
/** get Boundary Condition Particles: FaceXright */
int Collective::getBcPfaceXright(){
 return(bcPfaceXright);
}
/** get Boundary Condition Particles: FaceXleft */
int Collective::getBcPfaceXleft(){
 return(bcPfaceXleft);
}
/** get Boundary Condition Particles: FaceYright */
int Collective::getBcPfaceYright(){
 return(bcPfaceYright);
}
/** get Boundary Condition Particles: FaceYleft */
int Collective::getBcPfaceYleft(){
 return(bcPfaceYleft);
}
/** get Boundary Condition Electrostatic Potential: FaceXright */
int Collective::getBcPHIfaceXright(){
 return(bcPHIfaceXright);
}
/** get Boundary Condition Electrostatic Potential:FaceXleft */
int Collective::getBcPHIfaceXleft(){
 return(bcPHIfaceXleft);
}
/** get Boundary Condition Electrostatic Potential:FaceYright */
int Collective::getBcPHIfaceYright(){
 return(bcPHIfaceYright);
}
/** get Boundary Condition Electrostatic Potential:FaceYleft */
int Collective::getBcPHIfaceYleft(){
 return(bcPHIfaceYleft);
}
/** get Boundary Condition EM Field: FaceXright */
int Collective::getBcEMfaceXright(){
 return(bcEMfaceXright);
}
/** get Boundary Condition EM Field: FaceXleft */
int Collective::getBcEMfaceXleft(){
 return(bcEMfaceXleft);
}
/** get Boundary Condition EM Field: FaceYright */
int Collective::getBcEMfaceYright(){
 return(bcEMfaceYright);
}
/** get Boundary Condition EM Field: FaceYleft */
int Collective::getBcEMfaceYleft(){
 return(bcEMfaceYleft);
}

/** Get GEM Challenge parameters */
double Collective::getDelta(){
  return(delta);
}
 double Collective::getB0x(){
  return(B0x);
}
 double Collective::getB0y(){
  return(B0y);
}
 double Collective::getB0z(){
  return(B0z);
}
/** get the boolean value for verbose results */
 bool Collective::getVerbose(){
  return(verbose);
}
/** get the boolean value for TrackParticleID */
 bool Collective::getTrackParticleID(int nspecies){
  return(TrackParticleID[nspecies]);
   }
 int Collective::getRestart_status(){
  return(restart_status);
 }
 /** get SaveDirName  */
 string Collective::getSaveDirName(){
 return(SaveDirName);
 }
 /** get RestartDirName  */
 string Collective::getRestartDirName(){
 return(RestartDirName);
 }
 /** get inputfile  */
 string Collective::getinputfile(){
 return(inputfile);
  }
/** get last_cycle  */
int Collective::getLast_cycle(){
 return(last_cycle);
}
/** get the velocity of injection of the plasma from the wall */
double Collective::getVinj(){
 return(Vinj);
}
/** get the converging tolerance for CG solver */
double Collective::getCGtol(){
 return(CGtol);
}
/** get the converging tolerance for GMRES solver */
double Collective::getGMREStol(){
 return(GMREStol);
}
/** get the numbers of iteration for the PC mover */
int Collective::getNiterMover(){
 return(NiterMover);
}
/** output of fields */
int Collective::getFieldOutputCycle(){
 return(FieldOutputCycle);
}
/** output of fields */
int Collective::getParticlesOutputCycle(){
  return(ParticlesOutputCycle);
}
/** output of fields */
int Collective::getRestartOutputCycle(){
  return(RestartOutputCycle);
}

/**AMR variables, ME*/
/**get number of cells, ghost cell INCLUDED, for particle repopulation; x left*/
int Collective::GetPRA_Xleft(){
  return(PRA_Xleft);
}
/**get number of cells, ghost cell INCLUDED, for particle repopulation; x right*/
int Collective::GetPRA_Xright(){
  return(PRA_Xright);
}
/**get number of cells, ghost cell INCLUDED, for particle repopulation; y left*/
int Collective::GetPRA_Yleft(){
  return(PRA_Yleft);
}
/**get number of cells, ghost cell INCLUDED, for particle repopulation; y right*/
int Collective::GetPRA_Yright(){
  return(PRA_Yright);
}

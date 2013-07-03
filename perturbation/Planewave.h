/***************************************************************************
                          Planewave.h  -  Add plane waves perturbation to fields and particles
                             -------------------
    begin             : Mon May 14 2004
    developers        : Enrico Camporeale, David Burgess
 ***************************************************************************/

#ifndef Planewave_H
#define Planewave_H


#include <iostream>
#include <fstream>
#include <complex>

#include <math.h>
#include <mpi.h>
#include "../kvf/src/kvfDataSource.h"


using std::cout;
using std::cerr;
using std::endl;
using std::cin;
using namespace std;



/**
*   Add plane waves perturbation to fields and particles 
*   
*/

class Planewave {
  public:
    /** constructor: initialize physical parameters with values */
    Planewave(CollectiveIO *col, Field* EMf, Grid* grid, VCtopology* vct);
    /** Resize vectors depending on number of waves*/
    void resize();
    /** Rotate XY plane in counterclockwise direction and apply a Single Wave perturbation (imposed on a maxwellian)*/
    void SingleWave_Rotated(Particles2D* part);
    void Wave_Rotated(Particles2D* part);
 
  private:
  string inputfile;
  string ampl_file;
  double initial_ampl; // this is defined as (delta_B)/B0
  int ns,Nwaves;
  vector<double>  kpar, theta, Ex_mod,Ey_mod,Ez_mod,Bx_mod,By_mod,Bz_mod, n0_mod, n1_mod ;   
  vector<double>  Ex_phase,Ey_phase,Ez_phase,Bx_phase,By_phase,Bz_phase, n0_phase,n1_phase;  
  vector<double>  jx0_phase,jy0_phase,jz0_phase,jx1_phase,jy1_phase,jz1_phase,jx0_mod,jy0_mod,jz0_mod,jx1_mod,jy1_mod,jz1_mod;
  double scaling_factor[3], B0;
  EMfieldsIMP3D* _field;
  Grid* _grid;
  VCtopology* _vct;


};


inline Planewave::Planewave(CollectiveIO *col, Field* EMf, Grid* grid, VCtopology* vct){
  using namespace KVF;
  _field=EMf;
  _grid=grid; 
  _vct=vct;
  
   // read from inputfile
   inputfile=col->getinputfile(); 
   // number of species
   ns =col->getNs();

   try{
      cout << "About to open file: "<< inputfile << " in Planewave"<<endl;
      KVFDataSource kv_file( "/home/ec/parsek2D/"+inputfile );
      if( kv_file.num_parse_errors() != 0 ){
	cout << "THERE WERE ERRORS READING THE FILE!! CHECK Planewave.h and " << "/home/ec/parsek2D/" << inputfile << endl;
	kv_file.diag_print_errors() ;
	exit(1);}

      kv_file.get_data( "wave_amplitudes_file", ampl_file );
      kv_file.get_data( "initial_amplitude", initial_ampl );
      kv_file.get_data( "Number_of_waves", Nwaves );

      vector<double> d_vvals;
      kv_file.get_data( "Amplitude_scale", d_vvals );
      for (register int i=0; i<3; i++)
		scaling_factor[i]=d_vvals[i];
}
catch( KVFException& e )
{
	cout << "Planewave:: Caught exception " << endl;
	e.diag_cout();
}

resize();
ifstream infile;
infile.open (ampl_file.c_str());
  if (infile.is_open()){
    for (register int i=0; i<Nwaves; i++){
      infile >> kpar[i] >> theta[i] >> Ex_mod[i] >> Ex_phase[i] >> Ey_mod[i] >> Ey_phase[i] >> Ez_mod[i] >> Ez_phase[i] >> Bx_mod[i]>> Bx_phase[i] >> By_mod[i] >> By_phase[i] >> Bz_mod[i] >> Bz_phase[i] >> n0_mod[i] >> n0_phase[i] >> n1_mod[i] >> n1_phase[i]>> jx0_mod[i] >> jx0_phase[i] >> jy0_mod[i] >> jy0_phase[i] >> jz0_mod[i] >> jz0_phase[i] >> jx1_mod[i] >> jx1_phase[i] >> jy1_mod[i] >> jy1_phase[i] >> jz1_mod[i] >> jz1_phase[i] ;}
  }
   else{cout<<"ERROR: wave_amplitude_file "<< ampl_file.c_str()<< " not open"<<endl;}
    infile.close();

// Rescale amplitudes
    for (register int i=0; i<Nwaves; i++){
	Ex_mod[i] *=scaling_factor[0];
	Ey_mod[i] *=scaling_factor[0];
	Ez_mod[i] *=scaling_factor[0];
	Bx_mod[i] *=scaling_factor[0];
	By_mod[i] *=scaling_factor[0];
	Bz_mod[i] *=scaling_factor[0];
	n0_mod[i] *=scaling_factor[1];	n1_mod[i] *=scaling_factor[1];
 	jx0_mod[i]*=scaling_factor[2];	jx1_mod[i]*=scaling_factor[2];	
 	jy0_mod[i]*=scaling_factor[2];	jy1_mod[i]*=scaling_factor[2];	
 	jz0_mod[i]*=scaling_factor[2];	jz1_mod[i]*=scaling_factor[2];	
	kpar[i] *= scaling_factor[1]/scaling_factor[0];}

}
/** resize */
inline void Planewave::resize(){
	kpar.resize( Nwaves );
	theta.resize( Nwaves ); 
	Ex_mod.resize ( Nwaves );
	Ey_mod.resize ( Nwaves );
	Ez_mod.resize ( Nwaves );
	Bx_mod.resize ( Nwaves );
	By_mod.resize ( Nwaves );
	Bz_mod.resize ( Nwaves );
	n0_mod.resize ( Nwaves );
	n1_mod.resize ( Nwaves );
	jx0_mod.resize ( Nwaves );   jx1_mod.resize ( Nwaves );
	jy0_mod.resize ( Nwaves );   jy1_mod.resize ( Nwaves );
	jz0_mod.resize ( Nwaves );   jz1_mod.resize ( Nwaves );
		
	Ex_phase.resize ( Nwaves );
	Ey_phase.resize ( Nwaves );
	Ez_phase.resize ( Nwaves );
	Bx_phase.resize ( Nwaves );
	By_phase.resize ( Nwaves );
	Bz_phase.resize ( Nwaves );
	n0_phase.resize ( Nwaves );
	n1_phase.resize ( Nwaves );
	jx0_phase.resize ( Nwaves ); jx1_phase.resize ( Nwaves );
	jy0_phase.resize ( Nwaves ); jy1_phase.resize ( Nwaves );
	jz0_phase.resize ( Nwaves ); jz1_phase.resize ( Nwaves );
};

/** ROTATED WAVE */
inline void Planewave::Wave_Rotated(Particles2D* part) {
  complex<double> Ex_rot,Ey_rot, Bx_rot, By_rot, jx0_rot, jy0_rot, jx1_rot, jy1_rot;
  complex<double> I(0.0,1.0);
  double theta_rad= theta[0]/180*M_PI;
  double kperp;

  B0= sqrt(_field->getBx(1,1,0)*_field->getBx(1,1,0)+_field->getBy(1,1,0)*_field->getBy(1,1,0)+_field->getBz(1,1,0)*_field->getBz(1,1,0));

for (register int i=0; i<Nwaves; i++){

{// rotation

Ex_rot= Ex_mod[i]*exp(I*Ex_phase[i])*cos(theta_rad)+Ey_mod[i]*exp(I*Ey_phase[i])*sin(theta_rad);
Ey_rot=-Ex_mod[i]*exp(I*Ex_phase[i])*sin(theta_rad)+Ey_mod[i]*exp(I*Ey_phase[i])*cos(theta_rad);

Bx_rot= Bx_mod[i]*exp(I*Bx_phase[i])*cos(theta_rad)+By_mod[i]*exp(I*By_phase[i])*sin(theta_rad);
By_rot=-Bx_mod[i]*exp(I*Bx_phase[i])*sin(theta_rad)+By_mod[i]*exp(I*By_phase[i])*cos(theta_rad);

jx0_rot= jx0_mod[i]*exp(I*jx0_phase[i])*cos(theta_rad)+jy0_mod[i]*exp(I*jy0_phase[i])*sin(theta_rad);
jy0_rot=-jx0_mod[i]*exp(I*jx0_phase[i])*sin(theta_rad)+jy0_mod[i]*exp(I*jy0_phase[i])*cos(theta_rad);

jx1_rot= jx1_mod[i]*exp(I*jx1_phase[i])*cos(theta_rad)+jy1_mod[i]*exp(I*jy1_phase[i])*sin(theta_rad);
jy1_rot=-jx1_mod[i]*exp(I*jx1_phase[i])*sin(theta_rad)+jy1_mod[i]*exp(I*jy1_phase[i])*cos(theta_rad);

Ex_mod[i]=abs(Ex_rot);Ey_mod[i]=abs(Ey_rot);
Bx_mod[i]=abs(Bx_rot);By_mod[i]=abs(By_rot);
jx0_mod[i]=abs(jx0_rot);jy0_mod[i]=abs(jy0_rot);
jx1_mod[i]=abs(jx1_rot);jy1_mod[i]=abs(jy1_rot);

Ex_phase[i]=atan2(imag(Ex_rot),real(Ex_rot));
Ey_phase[i]=atan2(imag(Ey_rot),real(Ey_rot));
Bx_phase[i]=atan2(imag(Bx_rot),real(Bx_rot));
By_phase[i]=atan2(imag(By_rot),real(By_rot));
jx0_phase[i]=atan2(imag(jx0_rot),real(jx0_rot));jx1_phase[i]=atan2(imag(jx1_rot),real(jx1_rot));
jy0_phase[i]=atan2(imag(jy0_rot),real(jy0_rot));jy1_phase[i]=atan2(imag(jy1_rot),real(jy1_rot));
kpar[i]=kpar[i]/cos(theta_rad);kperp=0.0;
}//


// perturbation on fields
 _field->AddPerturbationRho(initial_ampl,kpar[i], kperp, Bx_mod[i], By_mod[i], Bz_mod[i], -n0_mod[i], n0_phase[i], n1_mod[i], n1_phase[i], B0, _grid);
 _field->AddPerturbation(initial_ampl,kpar[i], kperp, Ex_mod[i], Ex_phase[i], Ey_mod[i], Ey_phase[i], Ez_mod[i], Ez_phase[i], Bx_mod[i], Bx_phase[i], By_mod[i], By_phase[i], Bz_mod[i], Bz_phase[i], B0,  _grid);

}

// Perturbation on particles
  part[0].maxwellian(_grid,_field,_vct);  // maxwellian
  part[1].maxwellian(_grid,_field,_vct);  // maxwellian

 for (register int i=0; i<Nwaves; i++){
  part[0].AddPerturbationJ(initial_ampl,kpar[i],kperp, Bx_mod[i], By_mod[i], Bz_mod[i], jx0_mod[i], jx0_phase[i], jy0_mod[i], jy0_phase[i], jz0_mod[i], jz0_phase[i], B0, _grid);
  part[1].AddPerturbationJ(initial_ampl,kpar[i],kperp, Bx_mod[i], By_mod[i], Bz_mod[i], jx1_mod[i], jx1_phase[i], jy1_mod[i], jy1_phase[i], jz1_mod[i], jz1_phase[i], B0, _grid);

}
}


#endif

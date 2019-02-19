#ifdef COOLING_GRACKLE



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../io.h"

#include "cool_grackle.h"
#include "../universal_constants.h"
#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

Cool_GK::Cool_GK( void ){}


void Initialize_Grackle( Cool_GK &Cool, struct parameters P,  Grav3D &Grav, Cosmology &Cosmo ){

  chprintf( "\nUsing GRACKLE for cooling...\n");


  grackle_verbose = 1;
  #ifdef MPI_CHOLLA
  // Enable output
  if (procID != 0 ) grackle_verbose = 0;
  #endif


  Cool.tiny_number = 1.e-20;
  Cool.gamma = P.gamma;

  Cool.dens_conv = Cosmo.rho_0_gas;
  Cool.momentum_conv = Cosmo.rho_0_gas * Cosmo.v_0_gas;
  Cool.energy_conv =   Cosmo.v_0_gas * Cosmo.v_0_gas ;

  Real Msun_CGS = 1.98847e33; //Msun in gr
  Real kpc_CGS = 3.086e21;  //kpc in cm
  Real km_CGS = 1e5; //km in cm

  Cool.dens_to_CGS = Cool.dens_conv * Msun_CGS / kpc_CGS / kpc_CGS / kpc_CGS * Cosmo.cosmo_h * Cosmo.cosmo_h;
  Cool.vel_to_CGS = km_CGS;
  Cool.energy_to_CGS =  km_CGS * km_CGS;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  Cool.units.comoving_coordinates = 1; // 1 if cosmological sim, 0 if not
  Cool.units.a_units = 1.0 ; // units for the expansion factor
  Cool.units.a_value = Cosmo.current_a / Cool.units.a_units;
  Cool.units.density_units = Cool.dens_to_CGS  / Cosmo.current_a / Cosmo.current_a / Cosmo.current_a ;
  Cool.units.length_units = kpc_CGS / Cosmo.cosmo_h * Cosmo.current_a;
  Cool.units.time_units = KPC / Cosmo.cosmo_h ;
  Cool.units.velocity_units = Cool.units.length_units / Cosmo.current_a / Cool.units.time_units; // since u = a * dx/dt

  // Second, create a chemistry object for parameters.  This needs to be a pointer.
  Cool.data = new chemistry_data;
  if (set_default_chemistry_parameters(Cool.data) == 0) {
    chprintf( "GRACKLE: Error in set_default_chemistry_parameters.\n");
    exit(-1) ;
  }
  // Set parameter values for chemistry.
  // Access the parameter storage with the struct you've created
  // or with the grackle_data pointer declared in grackle.h (see further below).
  Cool.data->use_grackle = 1;            // chemistry on
  Cool.data->with_radiative_cooling = 1; // G.Cooling on
  Cool.data->primordial_chemistry = 1;   // molecular network with H, He
  // #ifdef GRACKLE_METAL_COOLING
  // Cool.data->metal_cooling = 1;          // metal cooling on
  // #else
  Cool.data->metal_cooling = 1;          // metal cooling off
  // #endif
  Cool.data->UVbackground = 1;           // UV background on
  Cool.data->grackle_data_file = "src/cooling/CloudyData_UVB=HM2012.h5"; // data file
  Cool.data->use_specific_heating_rate = 0;
  Cool.data->use_volumetric_heating_rate = 0;
  Cool.data->cmb_temperature_floor = 1;
  Cool.data->omp_nthreads = 10;

  if ( Cool.data->UVbackground == 1) chprintf( "GRACKLE: Loading UV Background File: %s\n", Cool.data->grackle_data_file );

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&Cool.units) == 0) {
    chprintf( "GRACKLE: Error in initialize_chemistry_data.\n");
    exit(-1) ;
  }

  if ( Cool.data->UVbackground == 1){
    Cool.a_UVB_on = 1 / (Cool.data->UVbackground_redshift_on + 1 );
    chprintf( "GRACKLE: UVB on: %f \n", Cool.a_UVB_on  );
  }

}

#endif

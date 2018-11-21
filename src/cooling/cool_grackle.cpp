#ifdef COOLING_GRACKLE



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../io.h"

#include "cool_grackle.h"

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

Cool_GK::Cool_GK( void ){}

#define mh     1.67262171e-24
#define kboltz 1.3806504e-16

void Initialize_Grackle( Cool_GK &Cool, struct parameters P,  Grav3D &Grav, Cosmology &Cosmo ){

  chprintf( "\nUsing GRACKLE for cooling...\n");

  // Cool.mass_h = 1.67262171e-24;
  // Cool.k_boltz = 1.3806504e-16;
  //
  //
  // grackle_verbose = 1;
  //
  // #ifdef MPI_CHOLLA
  // if (procID != 0 ) grackle_verbose = 0;
  // #endif
  //
  // Cool.dens_conv = Cosmo.rho_0_gas;
  // Cool.momentum_conv = Cosmo.rho_0_gas * Cosmo.v_0_gas;
  // Cool.energy_conv = Cosmo.rho_0_gas * Cosmo.v_0_gas * Cosmo.v_0_gas;
  //
  //
  // Real Msun_CGS = 1.98847e33; //Msun in gr
  // Real kpc_CGS = 3.086e21;  //kpc in cm
  // Real km_CGS = 1e5; //km in cm
  //
  // Cool.dens_to_CGS = Msun_CGS / kpc_CGS / kpc_CGS / kpc_CGS * Cosmo.cosmo_h * Cosmo.cosmo_h;
  // Cool.vel_to_CGS = km_CGS;
  // Cool.energy_to_CGS = km_CGS * km_CGS;
  //
  // Cool.units.comoving_coordinates = 1; // 1 if cosmological sim, 0 if not
  // Cool.units.density_units = Cool.dens_conv * Cool.dens_to_CGS;
  // Cool.units.length_units = kpc_CGS;
  // Cool.units.time_units = 1.0e12;
  // Cool.units.velocity_units = km_CGS;
  // Cool.units.a_units = 1.0; // units for the expansion factor
  // // Set expansion factor to 1 for non-cosmological simulation.
  // Cool.units.a_value = Cosmo.current_a;
  //
  // Cool.temperature_units = Cool.mass_h * pow(Cool.units.a_units * Cool.units.length_units / Cool.units.time_units, 2) / Cool.k_boltz;
  //
  // Cool.data = new chemistry_data;
  //
  //
  // if (set_default_chemistry_parameters( Cool.data ) == 0) {
  //   fprintf(stderr, " GRACKLE: Error in set_default_chemistry_parameters.\n");
  //   return exit(-1);
  // }
  //
  // // Set parameter values for chemistry.
  // // Access the parameter storage with the struct you've created
  // // or with the grackle_data pointer declared in grackle.h (see further below).
  // Cool.data->use_grackle = 1;            // chemistry on
  // Cool.data->with_radiative_cooling = 1; // cooling on
  // Cool.data->primordial_chemistry = 3;   // molecular network with H, He, D
  // Cool.data->metal_cooling = 1;          // metal cooling on
  // Cool.data->UVbackground = 1;           // UV background on
  // Cool.data->grackle_data_file = "src/cooling/CloudyData_UVB=HM2012.h5"; // data file
  //
  // // Finally, initialize the chemistry object.
  // if (initialize_chemistry_data( &Cool.units ) == 0) {
  //   fprintf(stderr, " GRACKLE: Error in initialize_chemistry_data.\n");
  //   return exit(-1);
  // }
  //
  // Cool.tiny_number = 1.e-20;


  double initial_redshift = 0.;

  grackle_verbose = 0;
  #ifdef MPI_CHOLLA
  // Enable output
  if (procID != 0 ) grackle_verbose = 0;
  #endif

  // First, set up the units system.
  // These are conversions from code units to cgs.
  // code_units Cool.units;
  Cool.units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  Cool.units.density_units = 1.67e-24;
  Cool.units.length_units = 3.08567e+24;  //cm per Mpc
  Cool.units.time_units = 31557600000000.0;  //sec per Myr
  Cool.units.a_units = 1.0; // units for the expansion factor
  Cool.units.velocity_units = Cool.units.length_units / Cool.units.time_units;
  // Set expansion factor to 1 for non-cosmological simulation.
  Cool.units.a_value = 1. / (1. + initial_redshift) / Cool.units.a_units;

  // Second, create a chemistry object for parameters.  This needs to be a pointer.
  // chemistry_data *my_grackle_data;
  Cool.data = (chemistry_data *) malloc(sizeof(chemistry_data));
  if (set_default_chemistry_parameters(Cool.data) == 0) {
    chprintf( "GRACKLE: Error in set_default_chemistry_parameters.\n");
    exit(-1) ;
  }
  // Set parameter values for chemistry.
  // Access the parameter storage with the struct you've created
  // or with the grackle_data pointer declared in grackle.h (see further below).
  Cool.data->use_grackle = 1;            // chemistry on
  Cool.data->with_radiative_cooling = 0; // G.Cooling on
  Cool.data->primordial_chemistry = 1;   // molecular network with H, He, D
  Cool.data->metal_cooling = 0;          // metal cooling on
  Cool.data->UVbackground = 0;           // UV background on
  Cool.data->grackle_data_file = "src/cooling/CloudyData_UVB=HM2012.h5"; // data file

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&Cool.units) == 0) {
    chprintf( "GRACKLE: Error in initialize_chemistry_data.\n");
    exit(-1) ;
  }


}

#endif

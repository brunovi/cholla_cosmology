#ifdef COOLING_GRACKLE



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../io.h"

#include "grackle_functions.h"

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

void Initialize_Grackle_Fields( Grid3D &G ){

  // chprintf( "GRACKLE Local Fields.\n");
  //
  // G.Cool.fields.grid_rank = 3;
  // G.Cool.fields.grid_dimension = new int[3];
  // G.Cool.fields.grid_start = new int[3];
  // G.Cool.fields.grid_end = new int[3];
  //
  // G.Cool.fields.grid_dimension[0] = G.Grav.nz_local;
  // G.Cool.fields.grid_dimension[1] = G.Grav.ny_local;
  // G.Cool.fields.grid_dimension[2] = G.Grav.nx_local;
  //
  // for ( int i=0; i<3; i++ ){
  //   G.Cool.fields.grid_start[i] = G.H.n_ghost;
  // }
  //
  // G.Cool.fields.grid_dimension[0] = 0;
  // G.Cool.fields.grid_dimension[1] = 0;
  // G.Cool.fields.grid_dimension[2] = 0;
  //
  //
  // G.Cool.fields.grid_start[0] = 0;
  // G.Cool.fields.grid_start[1] = 0;
  // G.Cool.fields.grid_start[2] = 0;
  //
  // G.Cool.fields.grid_end[0] = G.Grav.nz_local ;
  // G.Cool.fields.grid_end[1] = G.Grav.ny_local ;
  // G.Cool.fields.grid_end[2] = G.Grav.nx_local ;
  // //
  // // G.Cool.fields.grid_end[0] = G.Grav.nz_local + G.H.n_ghost;
  // // G.Cool.fields.grid_end[1] = G.Grav.ny_local + G.H.n_ghost;
  // // G.Cool.fields.grid_end[2] = G.Grav.nx_local + G.H.n_ghost;
  //
  // G.Cool.field_size = ( G.Grav.nx_local + 2*G.H.n_ghost ) * ( G.Grav.ny_local + 2*G.H.n_ghost ) * ( G.Grav.nz_local + 2*G.H.n_ghost );
  // chprintf( " nx:%d   ny:%d   nz:%d \n", G.Cool.fields.grid_dimension[2], G.Cool.fields.grid_dimension[1], G.Cool.fields.grid_dimension[0] );
  // chprintf( " n_ghost:%d\n", G.H.n_ghost);
  // chprintf( " field_size: %d\n", G.Cool.field_size );
  //
  // G.Cool.fields.density         = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.internal_energy = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.x_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.y_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.z_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //
  // // for primordial_chemistry >= 1
  // if( G.Cool.data->primordial_chemistry >= 1){
  //   chprintf( " Allocating memory for: HI, HII, HeI, HeII, HeIII, e   densities\n");
  //   G.Cool.fields.HI_density      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.HII_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.HeI_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.HeII_density    = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.HeIII_density   = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.e_density       = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // }
  //
  // // for primordial_chemistry >= 2
  // if( G.Cool.data->primordial_chemistry >= 2){
  //   chprintf( " Allocating memory for: HM, H2I, H2II   densities\n");
  //   G.Cool.fields.HM_density      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.H2I_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.H2II_density    = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // }
  //
  //
  // // for primordial_chemistry >= 3
  // if( G.Cool.data->primordial_chemistry == 3){
  //   chprintf( " Allocating memory for: DI, DII, HDI   densities\n");
  //   G.Cool.fields.DI_density      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.DII_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //   G.Cool.fields.HDI_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // }
  //
  // // for metal_cooling = 1
  // if( G.Cool.data->metal_cooling == 1){
  //   chprintf( " Allocating memory for: metal density\n");
  //   G.Cool.fields.metal_density   = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // }
  //
  // // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  // G.Cool.fields.volumetric_heating_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // // specific heating rate (provide in units [egs s^-1 g^-1]
  // G.Cool.fields.specific_heating_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //
  // // radiative transfer ionization / dissociation rate fields (provide in units [1/s])
  // G.Cool.fields.RT_HI_ionization_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.RT_HeI_ionization_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.RT_HeII_ionization_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.RT_H2_dissociation_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // // radiative transfer heating rate field (provide in units [erg s^-1 cm^-3])
  // G.Cool.fields.RT_heating_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //
  //
  // chprintf( "GRACKLE initialized.\n\n");

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int nx = G.Grav.nx_local;
  int ny = G.Grav.ny_local;
  int nz = G.Grav.nz_local;
  G.Cool.field_size = nx * ny * nz;
  G.Cool.fields.grid_rank = 3;
  G.Cool.fields.grid_dimension = new int[3];
  G.Cool.fields.grid_start = new int[3];
  G.Cool.fields.grid_end = new int[3];
  for (int i = 0;i < 3;i++) {
    G.Cool.fields.grid_dimension[i] = nx; // the active dimension not including ghost zones.
    G.Cool.fields.grid_start[i] = 0;
    G.Cool.fields.grid_end[i] =  nx - 1;
  }
  // G.Cool.fields.grid_dimension[0] = G.Cool.field_size;
  // G.Cool.fields.grid_end[0] = G.Cool.field_size - 1;
  G.Cool.fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  G.Cool.fields.density         = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.internal_energy = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.x_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.y_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.z_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));

  // for primordial_chemistry >= 1
  if( G.Cool.data->primordial_chemistry >= 1){
    chprintf( " Allocating memory for: HI, HII, HeI, HeII, HeIII, e   densities\n");
    G.Cool.fields.HI_density      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.HII_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.HeI_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.HeII_density    = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.HeIII_density   = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.e_density       = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  }

  // for primordial_chemistry >= 2
  if( G.Cool.data->primordial_chemistry >= 2){
    chprintf( " Allocating memory for: HM, H2I, H2II   densities\n");
    G.Cool.fields.HM_density      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.H2I_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.H2II_density    = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  }
  // for primordial_chemistry >= 3
  if( G.Cool.data->primordial_chemistry == 3){
    chprintf( " Allocating memory for: DI, DII, HDI   densities\n");
    G.Cool.fields.DI_density      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.DII_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
    G.Cool.fields.HDI_density     = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  }
  // for metal_cooling = 1
  // if( G.Cool.data->metal_cooling == 1){
    chprintf( " Allocating memory for: metal density\n");
    G.Cool.fields.metal_density   = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // }
  // // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  // G.Cool.fields.volumetric_heating_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // // specific heating rate (provide in units [egs s^-1 g^-1]
  // G.Cool.fields.specific_heating_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  //
  // // radiative transfer ionization / dissociation rate fields (provided in units of [1/s])
  // G.Cool.fields.RT_HI_ionization_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.RT_HeI_ionization_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.RT_HeII_ionization_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // G.Cool.fields.RT_H2_dissociation_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // // radiative transfer heating rate (provide in units [erg s^-1 cm^-3])
  // G.Cool.fields.RT_heating_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));

};

#define mh     1.67262171e-24
#define kboltz 1.3806504e-16
void Copy_Fields_to_Grackle( Grid3D &G ){



  gr_float tiny_number = 1.e-20;
  // set temperature units
  double temperature_units = mh * pow(G.Cool.units.a_units *
                                      G.Cool.units.length_units /
                                      G.Cool.units.time_units, 2) / kboltz;


  int i;
  for (i = 0;i < G.Cool.field_size;i++) {
    G.Cool.fields.density[i] = 1.0;
    G.Cool.fields.x_velocity[i] = 0.0;
    G.Cool.fields.y_velocity[i] = 0.0;
    G.Cool.fields.z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    Real gamma = 1.666667;
    G.Cool.fields.internal_energy[i] = 1000./ ( gamma -1 ) / temperature_units;
    chprintf( "H fraction: %f\n", grackle_data->HydrogenFractionByMass ); 
    if( G.Cool.data->primordial_chemistry >= 1){
      G.Cool.fields.HI_density[i] = grackle_data->HydrogenFractionByMass * G.Cool.fields.density[i];
      // G.Cool.fields.HI_density[i] = G.Cool.fields.density[i];
      G.Cool.fields.HII_density[i] = tiny_number * G.Cool.fields.density[i];
      G.Cool.fields.HeI_density[i] = (1.0 - grackle_data->HydrogenFractionByMass) * G.Cool.fields.density[i];
      G.Cool.fields.HeI_density[i] = tiny_number * G.Cool.fields.density[i];
      G.Cool.fields.HeII_density[i] = tiny_number * G.Cool.fields.density[i];
      G.Cool.fields.HeIII_density[i] = tiny_number * G.Cool.fields.density[i];
      G.Cool.fields.e_density[i] = tiny_number * G.Cool.fields.density[i];
    }
    if( G.Cool.data->primordial_chemistry >= 2){
      G.Cool.fields.HM_density[i] = tiny_number * G.Cool.fields.density[i];
      G.Cool.fields.H2I_density[i] = tiny_number * G.Cool.fields.density[i];
      G.Cool.fields.H2II_density[i] = tiny_number * G.Cool.fields.density[i];
    }
    if( G.Cool.data->primordial_chemistry == 3){
      G.Cool.fields.DI_density[i] = 2.0 * 3.4e-5 * G.Cool.fields.density[i];
      G.Cool.fields.DII_density[i] = tiny_number * G.Cool.fields.density[i];
      G.Cool.fields.HDI_density[i] = tiny_number * G.Cool.fields.density[i];
    }
    // // solar metallicity
    // if( G.Cool.data->metal_cooling == 1){
    //   G.Cool.fields.metal_density[i] = grackle_data->SolarMetalFractionByMass * G.Cool.fields.density[i];
    // }

    // G.Cool.fields.volumetric_heating_rate[i] = 0.0;
    // G.Cool.fields.specific_heating_rate[i] = 0.0;
    //
    // G.Cool.fields.RT_HI_ionization_rate[i] = 0.0;
    // G.Cool.fields.RT_HeI_ionization_rate[i] = 0.0;
    // G.Cool.fields.RT_HeII_ionization_rate[i] = 0.0;
    // G.Cool.fields.RT_H2_dissociation_rate[i] = 0.0;
    // G.Cool.fields.RT_heating_rate[i] = 0.0;
  }
}
void Do_Cooling_Step( Real dt, Grid3D &G ){

  Copy_Fields_to_Grackle( G );
  //
  // Real dt_1 = 3.15e7 * 1e6 / G.Cool.units.time_units;
  // dt_1 = 0;
  //
  // if (solve_chemistry(&G.Cool.units, &G.Cool.fields, dt_1) == 0) {
  //   chprintf("GRACKLE: Error in solve_chemistry.\n");
  //   return ;
  // }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  // Evolving the chemistry.
  // some timestep
  double dt_1 = 3.15e7 * 1e6 / G.Cool.units.time_units;
  // double dt_1 = 3.15e3 * 1e6 / G.Cool.units.time_units;

  // if (solve_chemistry(&G.Cool.units, &G.Cool.fields, dt_1) == 0) {
  //   chprintf( "GRACKLE: Error in solve_chemistry.\n");
  //   return ;
  // }

  // Calculate cooling time.
  gr_float *cooling_time;
  cooling_time = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  if (calculate_cooling_time(&G.Cool.units, &G.Cool.fields, cooling_time) == 0) {
    chprintf( "GRACKLE: Error in calculate_cooling_time.\n");
    exit(-1) ;
  }

  chprintf( "Cooling time = %le s.\n", cooling_time[0] *  G.Cool.units.time_units);

  // Calculate temperature.
  gr_float *temperature;
  temperature = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  if (calculate_temperature(&G.Cool.units, &G.Cool.fields,  temperature) == 0) {
    chprintf( "GRACKLE: Error in calculate_temperature.\n");
    return ;
  }

  chprintf("Temperature = %le K.\n", temperature[0]);

  // Calculate pressure.
  gr_float *pressure;
  pressure = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  if (calculate_pressure(&G.Cool.units, &G.Cool.fields, pressure) == 0) {
    chprintf( "GRACKLE: Error in calculate_pressure.\n");
    return ;
  }

  chprintf("Pressure = %le.\n", pressure[0]);

  // Calculate gamma.
  gr_float *gamma;
  gamma = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  if (calculate_gamma(&G.Cool.units, &G.Cool.fields,  gamma) == 0) {
    chprintf("GRACKLE: Error in calculate_gamma.\n");
    return ;
  }

  chprintf( "gamma = %le.\n", gamma[0]);

}

void Clear_Data_Grackle( Grid3D &G ){
  free( G.Cool.fields.density );
  free( G.Cool.fields.internal_energy );
  free( G.Cool.fields.x_velocity );
  free( G.Cool.fields.y_velocity );
  free( G.Cool.fields.z_velocity );
  if( G.Cool.data->primordial_chemistry >= 1){
    free(G.Cool.fields.HI_density);
    free(G.Cool.fields.HII_density);
    free(G.Cool.fields.HeI_density);
    free(G.Cool.fields.HeII_density);
    free(G.Cool.fields.HeIII_density);
    free(G.Cool.fields.e_density);
  }
  if( G.Cool.data->primordial_chemistry >= 2){
    free(G.Cool.fields.HM_density);
    free(G.Cool.fields.H2I_density);
    free(G.Cool.fields.H2II_density);
  }
  if( G.Cool.data->primordial_chemistry == 3){
    free(G.Cool.fields.DI_density);
    free(G.Cool.fields.DII_density);
    free(G.Cool.fields.HDI_density);
  }
  if( G.Cool.data->metal_cooling == 1){
    free( G.Cool.fields.metal_density );
  }
}

#endif

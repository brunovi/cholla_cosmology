#ifdef COOLING_GRACKLE



#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "../io.h"

#include "grackle_functions.h"
#include "../universal_constants.h"

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif

void Initialize_Grackle_Fields( Grid3D &G ){

  chprintf( "GRACKLE Local Fields.\n");


  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int n_cells = G.H.nx * G.H.ny * G.H.nz;
  int nx = G.Grav.nx_local;
  int ny = G.Grav.ny_local;
  int nz = G.Grav.nz_local;
  G.Cool.field_size = n_cells;
  G.Cool.fields.grid_rank = 3;
  G.Cool.fields.grid_dimension = new int[3];
  G.Cool.fields.grid_start = new int[3];
  G.Cool.fields.grid_end = new int[3];
  G.Cool.fields.grid_dimension[0] = G.H.nx; // the active dimension not including ghost zones.
  G.Cool.fields.grid_dimension[1] = G.H.ny; // the active dimension not including ghost zones.
  G.Cool.fields.grid_dimension[2] = G.H.nz; // the active dimension not including ghost zones.
    for (int i = 0;i < 3;i++) {
    G.Cool.fields.grid_start[i] = G.H.n_ghost;
  }
  G.Cool.fields.grid_end[0] =  G.H.nx - G.H.n_ghost ;
  G.Cool.fields.grid_end[1] =  G.H.ny - G.H.n_ghost ;
  G.Cool.fields.grid_end[2] =  G.H.nz - G.H.n_ghost ;
  G.Cool.fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  // G.Cool.fields.density         = G.C.density;
  G.Cool.fields.density         = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.x_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.y_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.z_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.internal_energy = (Real *) malloc(G.Cool.field_size * sizeof(Real));

  chprintf( " Allocating memory for: HI, HII, HeI, HeII, HeIII, e   densities\n");
  G.Cool.fields.HI_density      = &G.C.scalar[ 0*n_cells ];
  G.Cool.fields.HII_density     = &G.C.scalar[ 1*n_cells ];;
  G.Cool.fields.HeI_density     = &G.C.scalar[ 2*n_cells ];;
  G.Cool.fields.HeII_density    = &G.C.scalar[ 3*n_cells ];;
  G.Cool.fields.HeIII_density   = &G.C.scalar[ 4*n_cells ];;
  G.Cool.fields.e_density       = &G.C.scalar[ 5*n_cells ];;

  chprintf( " Allocating memory for: metal density\n");
  G.Cool.fields.metal_density   = &G.C.scalar[ 6*n_cells ];;

  G.Cool.temperature = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.cooling_time = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  Set_Initial_Fields_Grackle( G );

  chprintf( "GRACKLE initialized.\n\n");
};


void Set_Initial_Fields_Grackle( Grid3D &G ){

  // set temperature units
  // double temperature_units = MASS_HYDROGEN * pow(G.Cool.units.a_units * G.Cool.units.length_units / G.Cool.units.time_units, 2) / K_BOLTZ;

  chprintf( "H fraction: %f\n", grackle_data->HydrogenFractionByMass );
  chprintf( "Solar Metal fraction: %f\n", grackle_data->SolarMetalFractionByMass );
  for (int i = 0;i < G.Cool.field_size;i++) {
    G.Cool.fields.x_velocity[i] = 0.0;
    G.Cool.fields.y_velocity[i] = 0.0;
    G.Cool.fields.z_velocity[i] = 0.0;
    G.Cool.fields.density[i] = G.C.density[i] * G.Cool.dens_conv ;
    G.Cool.fields.internal_energy[i] = G.C.GasEnergy[i]  / G.Cool.fields.density[i] * G.Cool.energy_conv * G.Cool.dens_conv / G.Cosmo.current_a / G.Cosmo.current_a ;

    // G.Cool.fields.HI_density[i] = grackle_data->HydrogenFractionByMass * G.Cool.fields.density[i];
    G.Cool.fields.HI_density[i] = G.Cool.fields.density[i];
    G.Cool.fields.HII_density[i] = G.Cool.tiny_number * G.Cool.fields.density[i];
    // G.Cool.fields.HeI_density[i] = (1.0 - grackle_data->HydrogenFractionByMass) * G.Cool.fields.density[i];
    G.Cool.fields.HeI_density[i] = G.Cool.tiny_number * G.Cool.fields.density[i];
    G.Cool.fields.HeII_density[i] = G.Cool.tiny_number * G.Cool.fields.density[i];
    G.Cool.fields.HeIII_density[i] = G.Cool.tiny_number * G.Cool.fields.density[i];
    G.Cool.fields.e_density[i] = G.Cool.tiny_number * G.Cool.fields.density[i];
    G.Cool.fields.metal_density[i] = G.Cool.tiny_number  * G.Cool.fields.density[i];
    // G.Cool.fields.metal_density[i] = G.Cool.data->SolarMetalFractionByMass  * G.Cool.fields.density[i];
  }

}

void Copy_Fields_to_Grackle( Grid3D &G ){
  // set temperature units
  double temperature_units = MASS_HYDROGEN * pow(G.Cool.units.a_units * G.Cool.units.length_units / G.Cool.units.time_units, 2) / K_BOLTZ;

  chprintf( "H fraction: %f\n", grackle_data->HydrogenFractionByMass );
  for (int i = 0;i < G.Cool.field_size;i++) {
    G.Cool.fields.density[i] = G.C.density[i] * G.Cool.dens_conv ;
    G.Cool.fields.internal_energy[i] = G.C.GasEnergy[i]  / G.Cool.fields.density[i] * G.Cool.energy_conv * G.Cool.dens_conv / G.Cosmo.current_a / G.Cosmo.current_a ;
    // initilize internal energy (here 1000 K for no reason)
    // G.Cool.fields.internal_energy[i] = 1000./ ( G.Cool.gamma -1 ) / temperature_units;
  }
}

void Do_Cooling_Step( Real dt, Grid3D &G ){

  Copy_Fields_to_Grackle( G );

  // Calculate temperature.
  if (calculate_temperature(&G.Cool.units, &G.Cool.fields,  G.Cool.temperature) == 0) {
    chprintf( "GRACKLE: Error in calculate_temperature.\n");
    return ;
  }
  chprintf("Temperature = %le K.\n", G.Cool.temperature[128]);

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/
  // Real dt_cool = G.Cosmo.dt_secs;
  // dt_cool = 1e-4;
  // chprintf( "dt_cool: %f s\n", dt_cool );
  // if (solve_chemistry(&G.Cool.units, &G.Cool.fields, dt_cool ) == 0) {
  //   chprintf( "GRACKLE: Error in solve_chemistry.\n");
  //   return ;
  // }
  //
  // Calculate cooling time.
  if (calculate_cooling_time(&G.Cool.units, &G.Cool.fields, G.Cool.cooling_time) == 0) {
    chprintf( "GRACKLE: Error in calculate_cooling_time.\n");
    exit(-1) ;
  }
  chprintf( "Cooling time = %le s.\n", G.Cool.cooling_time[0] *  G.Cool.units.time_units);

  //
  // // Calculate pressure.
  // gr_float *pressure;
  // pressure = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // if (calculate_pressure(&G.Cool.units, &G.Cool.fields, pressure) == 0) {
  //   chprintf( "GRACKLE: Error in calculate_pressure.\n");
  //   return ;
  // }
  //
  // chprintf("Pressure = %le.\n", pressure[0]);
  //
  // // Calculate gamma.
  // gr_float *gamma;
  // gamma = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  // if (calculate_gamma(&G.Cool.units, &G.Cool.fields,  gamma) == 0) {
  //   chprintf("GRACKLE: Error in calculate_gamma.\n");
  //   return ;
  // }
  //
  // chprintf( "gamma = %le.\n", gamma[0]);

}

void Clear_Data_Grackle( Grid3D &G ){
  free( G.Cool.fields.density );
  free( G.Cool.fields.internal_energy );
  free( G.Cool.fields.x_velocity );
  free( G.Cool.fields.y_velocity );
  free( G.Cool.fields.z_velocity );
  free( G.Cool.temperature );
  free( G.Cool.cooling_time );


}

#endif

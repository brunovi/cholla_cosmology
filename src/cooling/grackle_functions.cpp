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
  G.Cool.fields.grid_start[0] = G.H.n_ghost;
  G.Cool.fields.grid_start[1] = G.H.n_ghost;
  G.Cool.fields.grid_start[2] = G.H.n_ghost;
  G.Cool.fields.grid_end[0] =  G.H.nx - G.H.n_ghost - 1 ;
  G.Cool.fields.grid_end[1] =  G.H.ny - G.H.n_ghost - 1 ;
  G.Cool.fields.grid_end[2] =  G.H.nz - G.H.n_ghost - 1 ;
  // G.Cool.fields.grid_start[0] = 0;
  // G.Cool.fields.grid_start[1] = 0;
  // G.Cool.fields.grid_start[2] = 0;
  // G.Cool.fields.grid_end[0] =  G.H.nx - 1 ;
  // G.Cool.fields.grid_end[1] =  G.H.ny - 1 ;
  // G.Cool.fields.grid_end[2] =  G.H.nz - 1 ;

  G.Cool.fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  G.Cool.fields.density         = G.C.density;
  // G.Cool.fields.density         = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.x_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.y_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.z_velocity      = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.fields.internal_energy = (Real *) malloc(G.Cool.field_size * sizeof(Real));

  chprintf( " Allocating memory for: HI, HII, HeI, HeII, HeIII, e   densities\n");
  G.Cool.fields.HI_density      = &G.C.scalar[ 0*n_cells ];
  G.Cool.fields.HII_density     = &G.C.scalar[ 1*n_cells ];
  G.Cool.fields.HeI_density     = &G.C.scalar[ 2*n_cells ];
  G.Cool.fields.HeII_density    = &G.C.scalar[ 3*n_cells ];
  G.Cool.fields.HeIII_density   = &G.C.scalar[ 4*n_cells ];
  G.Cool.fields.e_density       = &G.C.scalar[ 5*n_cells ];

  // if ( G.Cool.data->metal_cooling == 1){
    chprintf( " Allocating memory for: metal density\n");
    G.Cool.fields.metal_density   = &G.C.scalar[ 6*n_cells ];
  // }

  #ifdef OUTPUT_TEMPERATURE
  G.Cool.temperature = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  #endif


  #ifdef OUTPUT_COOLING_RATE
  G.Cool.cooling_time = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  G.Cool.cooling_rate = (Real *) malloc(G.Cool.field_size * sizeof(Real));
  #endif

  Set_Initial_Fields_Grackle( G );

  chprintf( "GRACKLE initialized.\n\n");
};


void Set_Initial_Fields_Grackle( Grid3D &G ){


  chprintf( "H fraction: %f\n", grackle_data->HydrogenFractionByMass );
  chprintf( "Solar Metal fraction: %f\n", grackle_data->SolarMetalFractionByMass );

  int nx_g, ny_g, nz_g, nx, ny, nz, nGHST;
  nx_g = G.H.nx;
  ny_g = G.H.ny;
  nz_g = G.H.nz;
  nx = G.H.nx_real;
  ny = G.H.ny_real;
  nz = G.H.nz_real;
  nGHST = G.H.n_ghost;
  int i, j, k, i_g, j_g, k_g, id;
  int counter = 0;
  for (k=0; k<nz_g; k++) {
    for (j=0; j<ny_g; j++) {
      for (i=0; i<nx_g; i++) {
        id = i + j*nx_g + k*nx_g*ny_g;
        // G.Cool.fields.density[id] = G.C.density[id];
        G.Cool.fields.x_velocity[id] = 0.0;
        G.Cool.fields.y_velocity[id] = 0.0;
        G.Cool.fields.z_velocity[id] = 0.0;

        G.Cool.fields.internal_energy[id] = G.C.GasEnergy[id]  / G.C.density[id] * G.Cool.energy_conv / G.Cosmo.current_a / G.Cosmo.current_a ;
        // G.Cool.fields.metal_density[id] = 0  * G.Cool.fields.density[id];

      }
    }
  }

  #ifdef OUTPUT_TEMPERATURE
  if (calculate_temperature(&G.Cool.units, &G.Cool.fields,  G.Cool.temperature) == 0) {
    chprintf( "GRACKLE: Error in calculate_temperature.\n");
    return ;
  }
  chprintf("Temperature = %le K.\n", G.Cool.temperature[0]);
  #endif
}

void Copy_Fields_to_Grackle( Grid3D &G ){

  for (int id = 0;id < G.Cool.field_size;id++) {
    // G.Cool.fields.density[id] = G.C.density[id];
    G.Cool.fields.internal_energy[id] = G.C.GasEnergy[id]  / G.C.density[id] * G.Cool.energy_conv / G.Cosmo.current_a / G.Cosmo.current_a ;
  }
}

void Do_Cooling_Step( Grid3D &G ){

  Real kpc_CGS = 3.086e21;  //kpc in cm

  // Upfate the units conversion
  G.Cool.units.a_value = G.Cosmo.current_a / G.Cool.units.a_units;
  G.Cool.units.density_units = G.Cool.dens_to_CGS  / G.Cosmo.current_a / G.Cosmo.current_a / G.Cosmo.current_a ;
  G.Cool.units.length_units = kpc_CGS / G.Cosmo.cosmo_h * G.Cosmo.current_a;

  Copy_Fields_to_Grackle( G );

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/
  Real dt_cool = G.Cosmo.dt_secs;
  chprintf( " dt_cool: %e s\n", dt_cool );
  if (solve_chemistry(&G.Cool.units, &G.Cool.fields, dt_cool / G.Cool.units.time_units ) == 0) {
    chprintf( "GRACKLE: Error in solve_chemistry.\n");
    return ;
  }

  #ifdef OUTPUT_TEMPERATURE
  if (calculate_temperature(&G.Cool.units, &G.Cool.fields,  G.Cool.temperature) == 0) {
    chprintf( "GRACKLE: Error in calculate_temperature.\n");
    return ;
  }
  // for (int i=0; i<G.Cool.field_size; i++)  chprintf(" Temperature = %le K.\n", G.Cool.temperature[i]);
  #endif




  Update_Internal_Energy(G);

  //
  // // Calculate cooling time.
  // if (calculate_cooling_time(&G.Cool.units, &G.Cool.fields, G.Cool.cooling_time) == 0) {
  //   chprintf( "GRACKLE: Error in calculate_cooling_time.\n");
  //   exit(-1) ;
  // }
  // chprintf( "Cooling time = %le s.\n", G.Cool.cooling_time[5034] * G.Cool.units.time_units );
  //
  //
  // for (int i = 0;i < G.Cool.field_size;i++) {
  //   // G.Cool.cooling_rate[i] = cooling_units * G.Cool.fields.internal_energy[i] / G.Cool.cooling_time[i] / G.Cool.fields.density[i]   ;
  //   // G.Cool.cooling_rate[i] = G.Cool.cooling_time[i] ;
  //   G.Cool.cooling_rate[i] = G.Cool.cooling_time[i]   ;
  //   // G.Cool.cooling_rate[i] = G.Cool.temperature[i]   ;
  //   // chprintf( "Cooling time = %le s.\n", G.Cool.cooling_time[i] * G.Cool.units.time_units );
  // }

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


void Update_Internal_Energy( Grid3D &G ){

  int nx_g, ny_g, nz_g, nx, ny, nz, nGHST;
  nx_g = G.H.nx;
  ny_g = G.H.ny;
  nz_g = G.H.nz;
  nx = G.H.nx_real;
  ny = G.H.ny_real;
  nz = G.H.nz_real;
  nGHST = G.H.n_ghost;

  Real dens, ge_0, ge_1, delta_ge, vx, vy, vz, E, Ek;
  int k, j, i, id;
  // Real avrg_0, avrg_1;
  // avrg_0 = 0;
  // avrg_1 = 0;
  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      for (i=0; i<nx; i++) {
        id = (i+nGHST) + (j+nGHST)*nx_g + (k+nGHST)*nx_g*ny_g;
        dens = G.C.density[id];
        vx = G.C.momentum_x[id] / dens;
        vy = G.C.momentum_y[id] / dens;
        vz = G.C.momentum_z[id] / dens;
        Ek = 0.5 * dens * ( vx*vx + vy*vy + vz*vz );
        E = G.C.Energy[id];
        if ( E < Ek ) {
          std::cout << "##Negative Energy before Cooling: " <<  E   << " " << dens << std::endl;
          G.C.Energy[id] = Ek + G.C.GasEnergy[id];
        }
        ge_0 = G.C.GasEnergy[id];
        if ( ge_0 < 0 ) std::cout << "##Negative GasEnergy before Cooling: " <<  ge_0   << " " << dens << std::endl;
        if ( fabs(( E - (Ek + ge_0) )  / E ) > 1e-5 ) std::cout << "##Energy Error before cooling: " << E << "  " << Ek + ge_0 << std::endl;
        ge_1 = G.Cool.fields.internal_energy[id] * dens / G.Cool.energy_conv  * G.Cosmo.current_a * G.Cosmo.current_a;
        // avrg_0 += ge_0;
        delta_ge = ge_1 - ge_0;
        G.C.GasEnergy[id] += delta_ge ;
        G.C.Energy[id] += delta_ge ;
        ge_1 = G.C.GasEnergy[id];
        E = G.C.Energy[id];
        // avrg_1 += ge_1;
        if (ge_1 <= 0 ) std::cout << "##Negative gasEnergy after Cooling" <<  G.C.GasEnergy[id] << std::endl;
        if ( fabs(( E - (Ek + ge_1) )  / E ) > 1e-5 ) std::cout << "##Energy Error after cooling: " << E << "  " << Ek + ge_1 << std::endl;

        // if (E < 0 ) {
        //   std::cout << "##Negative Energy after Cooling: " <<  G.C.Energy[id] << "  " << Ek + G.C.GasEnergy[id]  << std::endl;
        //   G.C.Energy[id] = Ek + G.C.GasEnergy[id];
        // }
      }
    }
  }
  // avrg_0 /= (nx_g*ny_g*nz_g);
  // avrg_1 /= (nx_g*ny_g*nz_g);
  // std::cout << delta_ge_min << "   " << delta_ge_max << std::endl;
  // std::cout << avrg_0 << "   " << avrg_1 << std::endl;
}



void Clear_Data_Grackle( Grid3D &G ){
  // free( G.Cool.fields.density );
  free( G.Cool.fields.x_velocity );
  free( G.Cool.fields.y_velocity );
  free( G.Cool.fields.z_velocity );
  free( G.Cool.fields.internal_energy );


  #ifdef OUTPUT_COOLING_RATE
  free( G.Cool.cooling_rate );
  free( G.Cool.cooling_time );
  #endif

  #ifdef OUTPUT_TEMPERATURE
  free( G.Cool.temperature );
  #endif

}

#endif

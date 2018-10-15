#ifdef GRAVITY

#include "gravity_functions.h"

#ifdef COSMOLOGY
#include"../particles/particles_dynamics.h"
#endif



void Copy_Hydro_Density_to_Gravity( Grid3D &G ){
  int i, j, k, id, id_grav;
  // Copy the density array from hydro conserved to gravity density array
  for (k=0; k<G.Grav.nz_local; k++) {
    for (j=0; j<G.Grav.ny_local; j++) {
      for (i=0; i<G.Grav.nx_local; i++) {
        id = (i+G.H.n_ghost) + (j+G.H.n_ghost)*G.H.nx + (k+G.H.n_ghost)*G.H.nx*G.H.ny;
        id_grav = (i) + (j)*G.Grav.nx_local + (k)*G.Grav.nx_local*G.Grav.ny_local;
        #ifdef ONLY_PM
        G.Grav.F.density_h[id_grav] = 0;
        #else
        #ifdef COSMOLOGY
        G.Grav.F.density_h[id_grav] = G.C.density[id]*G.Cosmo.rho_0_gas;
        // G.Grav.F.density_h[id_grav] = G.C.density[id];
        #else
        G.Grav.F.density_h[id_grav] = G.C.density[id] ;
        #endif //COSMOLOGY
        #endif //ONLY_PM
      }
    }
  }
}

Real Get_Density_Average( Grid3D &G ){
  // int nGHST = G.Parts.G.n_ghost_particles_grid;
  // int nx_g = Parts.G.nx_local + 2*nGHST;
  // int ny_g = Parts.G.ny_local + 2*nGHST;
  // int nz_g = Parts.G.nz_local + 2*nGHST;
  int nx = G.Grav.nx_local;
  int ny = G.Grav.ny_local;
  int nz = G.Grav.nz_local;
  int k, j, i, id;
  Real dens_avrg=0;
  for( k=0; k<nz; k++){
    for( j=0; j<ny; j++){
      for( i=0; i<nx; i++){
        // id = (i+nGHST) + (j+nGHST)*nx_g + (k+nGHST)*nx_g*ny_g;
        id = (i) + (j)*nx + (k)*nx*ny;
        dens_avrg += G.Grav.F.density_h[id];
      }
    }
  }
  dens_avrg /= (nx*ny*nz);
  return dens_avrg;
}


#ifdef POTENTIAL_CUFFT
void Compute_Gravitational_Potential( Grid3D &G, Potential_CUFFT_3D &p_solver, Real *time_pot, Real *time_pDens, Real *time_pDens_trans, struct parameters P ){
  Real time_potential, start, stop, time_setup;

  Copy_Hydro_Density_to_Gravity( G );

  start = get_time();
  #ifdef PARTICLES
  Get_Particles_Density_CIC( G, P, time_pDens, time_pDens_trans );
  Copy_Particles_Density_to_Gravity( G );
  #endif
  stop = get_time();
  time_setup = (stop - start) * 1000.0;

  time_potential = p_solver.Get_Potential( G.Grav );
  // Extrapolate_Grav_Potential( G.Grav );
  // Copy_Potential_To_Hydro_Grid( G );

  *time_pot = time_potential;
  *time_pDens = time_setup;
}
#endif


#ifdef POTENTIAL_FFTW
void Compute_Gravitational_Potential( Grid3D &G, Potential_FFTW_3D &p_solver, Real *time_pot, Real *time_pDens, Real *time_pDens_trans, struct parameters P ){
  Real time_potential;

  Copy_Hydro_Density_to_Gravity( G );

  // start = get_time();
  #ifdef PARTICLES
  Get_Particles_Density_CIC( G, P, time_pDens, time_pDens_trans );
  Copy_Particles_Density_to_Gravity( G );
  #endif
  // stop = get_time();
  // time_setup = (stop - start) * 1000.0;

  time_potential = p_solver.Get_Potential( G.Grav );

  *time_pot = time_potential;
}
#endif


#ifdef POTENTIAL_PFFT
void Compute_Gravitational_Potential( Grid3D &G, Potential_PFFT_3D &p_solver, Real *time_pot, Real *time_pDens, Real *time_pDens_trans, struct parameters P ){
  Real time_potential;

  Copy_Hydro_Density_to_Gravity( G );

  #ifdef PARTICLES
  Get_Particles_Density_CIC( G, P, time_pDens, time_pDens_trans );
  Copy_Particles_Density_to_Gravity( G );
  #endif

  #ifdef COSMOLOGY
  Real dens_avrg;
  dens_avrg = Get_Density_Average( G );
  #ifdef MPI_CHOLLA
  dens_avrg = ReduceRealAvg( dens_avrg );
  #endif
  // G.Cosmo.dens_avrg = dens_avrg;
  G.Grav.dens_avrg = dens_avrg;
  // std::cout << "Density Averge: " << dens_avrg << std::endl;
  chprintf( "Densitty Average:  %f\n", dens_avrg);
  #endif

  time_potential = p_solver.Get_Potential( G.Grav );

  *time_pot = time_potential;
}
#endif





void Extrapolate_Grav_Potential( Grid3D &G ){
  int n_ghost_pot = N_GHOST_POTENTIAL;
  int nx_pot = G.Grav.nx_local + 2*n_ghost_pot;
  int ny_pot = G.Grav.ny_local + 2*n_ghost_pot;
  int nz_pot = G.Grav.nz_local + 2*n_ghost_pot;
  int n_ghost_grid = G.H.n_ghost;
  int nx_grid = G.Grav.nx_local + 2*n_ghost_grid;
  int ny_grid = G.Grav.ny_local + 2*n_ghost_grid;
  int nz_grid = G.Grav.nz_local + 2*n_ghost_grid;
  int nGHST = n_ghost_grid - n_ghost_pot;
  Real pot_now, pot_prev, pot_extrp;
  int k, j, i, id_pot, id_grid;
  for ( k=0; k<nz_pot; k++ ){
    for ( j=0; j<ny_pot; j++ ){
      for ( i=0; i<nx_pot; i++ ){
        id_pot = i + j*nx_pot + k*nx_pot*ny_pot;
        id_grid = (i+nGHST) + (j+nGHST)*nx_grid + (k+nGHST)*nx_grid*ny_grid;
        pot_now = G.C.Grav_potential[id_grid];
        if ( G.Grav.INITIAL ){
          pot_prev = pot_now;
          pot_extrp = pot_now;
          G.Grav.INITIAL = false;
        } else{
          pot_prev = G.Grav.F.potential_1_h[id_pot];
          pot_extrp = pot_now + 0.5 * G.Grav.dt_now * ( pot_now - pot_prev ) / G.Grav.dt_prev;
          // pot_extrp = pot_now;
        }

        #ifdef COSMOLOGY
        pot_extrp *= G.Cosmo.current_a * G.Cosmo.current_a / G.Cosmo.phi_0_gas;
        // pot_extrp *= G.Cosmo.current_a * G.Cosmo.current_a ;
        #endif

        G.C.Grav_potential[id_grid] = pot_extrp;
        G.Grav.F.potential_1_h[id_pot] = pot_now;
      }
    }
  }
}




void Copy_Potential_To_Hydro_Grid( Grid3D &G ){
  int n_ghost_pot = N_GHOST_POTENTIAL;
  int nx_pot = G.Grav.nx_local + 2*n_ghost_pot;
  int ny_pot = G.Grav.ny_local + 2*n_ghost_pot;
  int nz_pot = G.Grav.nz_local + 2*n_ghost_pot;
  int n_ghost_grid = G.H.n_ghost;
  int nx_grid = G.Grav.nx_local + 2*n_ghost_grid;
  int ny_grid = G.Grav.ny_local + 2*n_ghost_grid;
  int nz_grid = G.Grav.nz_local + 2*n_ghost_grid;
  int k, j, i, id_pot, id_grid;
  for ( k=0; k<G.Grav.nz_local; k++ ){
    for ( j=0; j<G.Grav.ny_local; j++ ){
      for ( i=0; i<G.Grav.nx_local; i++ ){
        id_pot = (i+n_ghost_pot) + (j+n_ghost_pot)*nx_pot + (k+n_ghost_pot)*nx_pot*ny_pot;
        id_grid = (i+n_ghost_grid) + (j+n_ghost_grid)*nx_grid + (k+n_ghost_grid)*nx_grid*ny_grid;
        G.C.Grav_potential[id_grid] = G.Grav.F.potential_h[id_pot];
      }
    }
  }
}

void Copy_Potential_From_Hydro_Grid( Grid3D &G ){
  int n_ghost_pot = N_GHOST_POTENTIAL;
  int nx_pot = G.Grav.nx_local + 2*n_ghost_pot;
  int ny_pot = G.Grav.ny_local + 2*n_ghost_pot;
  int nz_pot = G.Grav.nz_local + 2*n_ghost_pot;
  int n_ghost_grid = G.H.n_ghost;
  int nx_grid = G.Grav.nx_local + 2*n_ghost_grid;
  int ny_grid = G.Grav.ny_local + 2*n_ghost_grid;
  int nz_grid = G.Grav.nz_local + 2*n_ghost_grid;
  int nGHST = n_ghost_grid - n_ghost_pot;
  int k, j, i, id_pot, id_grid;
  for ( k=0; k<nz_pot; k++ ){
    for ( j=0; j<ny_pot; j++ ){
      for ( i=0; i<nx_pot; i++ ){
        id_pot = i + j*nx_pot + k*nx_pot*ny_pot;
        id_grid = (i+nGHST) + (j+nGHST)*nx_grid + (k+nGHST)*nx_grid*ny_grid;
        G.Grav.F.potential_h[id_pot] = G.C.Grav_potential[id_grid];
      }
    }
  }
}

void Set_dt( Grid3D &G, bool &output_now ){

  #ifdef COSMOLOGY
  Real delta_a_part;
  Real dt_courant, dt_gas, da_courant;
  #ifdef PECULIAR_VEL
  Real delta_t_part = Get_Particles_dt_cosmo( G );
  // chprintf( " dt: %f  \n", delta_t_part );
  delta_a_part = G.Cosmo.Get_da_from_dt( delta_t_part );
  #else
  delta_a_part = Get_Particles_da_cosmo( G );
  dt_courant = G.H.dt;
  da_courant = G.Cosmo.Get_da_courant( dt_courant);
  chprintf( "Delta_a parts: %f:   Delta_a gas: %f\n", delta_a_part, da_courant);
  da_courant = std::min(delta_a_part, da_courant);
  #endif
  // Real delta_a_part = G.Cosmo.max_delta_a;
  if( da_courant > G.Cosmo.max_delta_a) da_courant = G.Cosmo.max_delta_a;


  G.Cosmo.delta_a = da_courant;
  if ( (G.Cosmo.current_a + G.Cosmo.delta_a) >  G.Cosmo.next_output ){
    G.Cosmo.delta_a = G.Cosmo.next_output - G.Cosmo.current_a;
    output_now = true;
    // chprintf( " ################################## \n");
  }

  Real dt = G.Cosmo.Get_dt_from_da( G.Cosmo.delta_a );
  Real da_2 = G.Cosmo.Get_da_from_dt( dt/2 );
  G.Cosmo.delta_a_2 = da_2;
  G.Particles.dt = dt;

  dt_gas = G.Cosmo.Get_Cosmology_dt( G.Cosmo.delta_a );
  G.H.dt = dt_gas;
  // Real dt_courant = dt_gas*G.Cosmo.t_0_gas;
  // G.H.dt = dt;

  chprintf( "Current_a: %f    delta_a: %f   dt: %f    da_courant: %f \n", G.Cosmo.current_a, G.Cosmo.delta_a, dt_gas, da_courant  );
  chprintf( " dt: %f  \n", dt );

  #endif


  if ( G.Grav.INITIAL ){
    G.Grav.dt_prev = G.H.dt;
    G.Grav.dt_now = G.H.dt;
  }else{
    G.Grav.dt_prev = G.Grav.dt_now;
    G.Grav.dt_now = G.H.dt;
  }


  // #ifdef PARTICLES
  // dt_particles = Get_Particles_dt( G.Particles );
  // dt_min = std::min( G.H.dt, dt_particles );
  // G.H.dt = dt_min;
  // G.Particles.dt = dt_min;
  // // G.Particles.dt = 0;
  // #endif






}


#endif

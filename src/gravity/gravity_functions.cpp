#ifdef GRAVITY

#include "gravity_functions.h"



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
        G.Grav.F.density_h[id_grav] = 4 * M_PI * G.C.density[id] ;
        #endif
      }
    }
  }
}


#ifdef POTENTIAL_CUFFT
void Compute_Gravitational_Potential( Grid3D &G, Potential_CUFFT_3D &p_solver, Real *time_pot, Real *time_pDens, Real *time_pDens_trans, struct parameters P ){
  Real time_potential, start, stop, time_setup;

  Copy_Hydro_Density_to_Gravity( G );

  start = get_time();
  #ifdef PARTICLES
  Get_Particles_Density_CIC( G, P );
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
        }
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


#endif
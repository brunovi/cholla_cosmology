#ifdef GRAVITY

#include "gravity_functions.h"
#include "../universal_constants.h"

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

  G.Grav.dens_avrg = dens_avrg;
  chprintf( " Density Average:  %f\n", dens_avrg);
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
        pot_now = G.C.Grav_potential[id_grid]  ;
        if ( G.Grav.INITIAL ){
          pot_prev = pot_now;
          pot_extrp = pot_now *  G.Cosmo.current_a * G.Cosmo.current_a ;
        } else{
          pot_prev = G.Grav.F.potential_1_h[id_pot] ;
          #ifdef GRAVITY_CORRECTOR
          pot_extrp = pot_now * G.Cosmo.current_a * G.Cosmo.current_a;
          #else
          // pot_extrp = pot_now * G.Cosmo.current_a * G.Cosmo.current_a + 0.5 * G.Grav.dt_now * ( pot_now * G.Cosmo.current_a * G.Cosmo.current_a - pot_prev * (G.Cosmo.current_a - G.Cosmo.delta_a ) * (G.Cosmo.current_a - G.Cosmo.delta_a ) ) / G.Grav.dt_prev;
          pot_extrp = pot_now * G.Cosmo.current_a * G.Cosmo.current_a + 0.5 * G.Grav.dt_now * ( pot_now * G.Cosmo.current_a * G.Cosmo.current_a - pot_prev * (G.Cosmo.current_a ) * (G.Cosmo.current_a ) ) / G.Grav.dt_prev;
          #endif
        }

        #ifdef COSMOLOGY
        // pot_extrp *= 1 / G.Cosmo.phi_0_gas * G.Cosmo.current_a * G.Cosmo.current_a;
        pot_extrp *= 1 / G.Cosmo.phi_0_gas ;
        #endif

        G.C.Grav_potential[id_grid] = pot_extrp;
        G.Grav.F.potential_1_h[id_pot] = pot_now;
      }
    }
  }
  G.Grav.INITIAL = false;
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
  Real pot;
  int k, j, i, id_pot, id_grid;
  for ( k=0; k<G.Grav.nz_local; k++ ){
    for ( j=0; j<G.Grav.ny_local; j++ ){
      for ( i=0; i<G.Grav.nx_local; i++ ){
        id_pot = (i+n_ghost_pot) + (j+n_ghost_pot)*nx_pot + (k+n_ghost_pot)*nx_pot*ny_pot;
        id_grid = (i+n_ghost_grid) + (j+n_ghost_grid)*nx_grid + (k+n_ghost_grid)*nx_grid*ny_grid;
        pot = G.Grav.F.potential_h[id_pot];

        // #ifdef COSMOLOGY
        // pot *= G.Cosmo.current_a * G.Cosmo.current_a / G.Cosmo.phi_0_gas;
        // #endif

        G.C.Grav_potential[id_grid] = pot;
      }
    }
  }
}

#ifdef GRAVITY_CORRECTOR
void Get_Gavity_Corrector( Grid3D &G, int g_start, int g_end ){

  int nx_grav, ny_grav, nz_grav;
  // nGHST_grav = G.Particles.G.n_ghost_particles_grid;
  nx_grav = G.Grav.nx_local;
  ny_grav = G.Grav.ny_local;
  nz_grav = G.Grav.nz_local;

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.Grav.nx_local + 2*nGHST_grid;
  ny_grid = G.Grav.ny_local + 2*nGHST_grid;
  nz_grid = G.Grav.nz_local + 2*nGHST_grid;

  int nGHST = nGHST_grid ;

  Real dx, dy, dz;
  dx = G.Grav.dx;
  dy = G.Grav.dy;
  dz = G.Grav.dz;

  Real delta_max, delta;
  delta_max = 100000;
  Real pot_factor = 1. / G.Cosmo.phi_0_gas * (G.Cosmo.current_a - G.Cosmo.delta_a) * (G.Cosmo.current_a - G.Cosmo.delta_a) ;
  Real phi_l, phi_r;
  int k, j, i, id_l, id_r, id;
  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i-1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r = (i+1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l = G.C.Grav_potential[id_l] * pot_factor;
        phi_r = G.C.Grav_potential[id_r] * pot_factor;
        G.Grav.F.gravity_x_h[id] = -0.5 * ( phi_r - phi_l ) / dx;
        phi_l = G.C.Grav_potential_0[id_l] ;
        phi_r = G.C.Grav_potential_0[id_r] ;
        G.Grav.F.gravity_x_h_prev[id] = -0.5 * ( phi_r - phi_l ) / dx;
        delta = fabs( (G.Grav.F.gravity_x_h_prev[id] - G.Grav.F.gravity_x_h[id]) / G.Grav.F.gravity_x_h_prev[id] );
        if ( delta > delta_max ) std::cout << "### Grav X: " << delta<< " delta_max: " << delta_max << std::endl;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i + nGHST) + (j-1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r = (i + nGHST) + (j+1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l = G.C.Grav_potential[id_l] * pot_factor;
        phi_r = G.C.Grav_potential[id_r] * pot_factor;
        G.Grav.F.gravity_y_h[id] = -0.5 * ( phi_r - phi_l ) / dy;
        phi_l = G.C.Grav_potential_0[id_l] ;
        phi_r = G.C.Grav_potential_0[id_r] ;
        G.Grav.F.gravity_y_h_prev[id] = -0.5 * ( phi_r - phi_l ) / dy;
        delta = fabs( (G.Grav.F.gravity_y_h_prev[id] - G.Grav.F.gravity_y_h[id]) / G.Grav.F.gravity_y_h_prev[id] );
        if ( delta > delta_max ) std::cout << "### Grav Y: " << delta<< " delta_max: " << delta_max << std::endl;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i + nGHST) + (j + nGHST)*nx_grid + (k-1 + nGHST)*ny_grid*nz_grid;
        id_r = (i + nGHST) + (j + nGHST)*nx_grid + (k+1 + nGHST)*ny_grid*nz_grid;
        phi_l = G.C.Grav_potential[id_l] * pot_factor;
        phi_r = G.C.Grav_potential[id_r] * pot_factor;
        G.Grav.F.gravity_z_h[id] = -0.5 * ( phi_r - phi_l ) / dz;
        phi_l = G.C.Grav_potential_0[id_l];
        phi_r = G.C.Grav_potential_0[id_r];
        G.Grav.F.gravity_z_h_prev[id] = -0.5 * ( phi_r - phi_l ) / dz;
        delta = fabs( (G.Grav.F.gravity_z_h_prev[id] - G.Grav.F.gravity_z_h[id]) / G.Grav.F.gravity_z_h_prev[id] );
        if ( delta > delta_max ) std::cout << "### Grav Z: " << delta<< " delta_max: " << delta_max << std::endl;
      }
    }
  }
}

void Add_Gravity_Corrector( Grid3D &G, int g_start, int g_end ){

  int nx_grav, ny_grav, nz_grav;
  // nGHST_grav = G.Particles.G.n_ghost_particles_grid;
  nx_grav = G.Grav.nx_local;
  ny_grav = G.Grav.ny_local;
  nz_grav = G.Grav.nz_local;

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.Grav.nx_local + 2*nGHST_grid;
  ny_grid = G.Grav.ny_local + 2*nGHST_grid;
  nz_grid = G.Grav.nz_local + 2*nGHST_grid;

  int nGHST = nGHST_grid ;

  Real d, vx, vy, vz, gx, gy, gz;
  Real d_0, vx_0, vy_0, vz_0, gx_0, gy_0, gz_0;
  Real current_a_prev;
  Real E, u_floor, delta_u;

  current_a_prev = G.Cosmo.current_a - G.Cosmo.delta_a;
  int k, j, i, id_grav, id_grid;
  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id_grav = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_grid = (i + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;

        d_0 = G.C.density_0[id_grid];
        vx_0 = G.C.momentum_x_0[id_grid] / d_0;
        vy_0 = G.C.momentum_y_0[id_grid] / d_0;
        vz_0 = G.C.momentum_z_0[id_grid] / d_0;
        gx_0 = G.Grav.F.gravity_x_h_prev[id_grav];
        gy_0 = G.Grav.F.gravity_y_h_prev[id_grav];
        gz_0 = G.Grav.F.gravity_z_h_prev[id_grav];

        d = G.C.density[id_grid];
        vx = G.C.momentum_x[id_grid] / d;
        vy = G.C.momentum_y[id_grid] / d;
        vz = G.C.momentum_z[id_grid] / d;
        gx = G.Grav.F.gravity_x_h[id_grav];
        gy = G.Grav.F.gravity_y_h[id_grav];
        gz = G.Grav.F.gravity_z_h[id_grav];

        G.C.momentum_x[id_grid] -= 0.5 * G.H.dt * d_0 * gx_0;
        G.C.momentum_y[id_grid] -= 0.5 * G.H.dt * d_0 * gy_0;
        G.C.momentum_z[id_grid] -= 0.5 * G.H.dt * d_0 * gz_0;
        G.C.Energy[id_grid] -= 0.5 * G.H.dt * d_0 * ( vx_0*gx_0 + vy_0*gy_0 + vz_0*gz_0 );

        G.C.momentum_x[id_grid] += 0.5 * G.H.dt * d * gx;
        G.C.momentum_y[id_grid] += 0.5 * G.H.dt * d * gy;
        G.C.momentum_z[id_grid] += 0.5 * G.H.dt * d * gz;
        G.C.Energy[id_grid] += 0.5 * G.H.dt * ( d*vx*gx + d*vy*gy + d*vz*gz );

        // G.C.momentum_x[id_grid] +=  G.H.dt * d_0 * gx_0;
        // G.C.momentum_y[id_grid] +=  G.H.dt * d_0 * gy_0;
        // G.C.momentum_z[id_grid] +=  G.H.dt * d_0 * gz_0;
        // G.C.Energy[id_grid] +=  G.H.dt * d_0 * ( vx_0*gx_0 + vy_0*gy_0 + vz_0*gz_0 );

        // G.C.momentum_x[id_grid] +=  G.H.dt * d * gx;
        // G.C.momentum_y[id_grid] +=  G.H.dt * d * gy;
        // G.C.momentum_z[id_grid] +=  G.H.dt * d * gz;
        // G.C.Energy[id_grid] +=  G.H.dt * d * ( vx*gx + vy*gy + vz*gz );



        // u_floor = 0;
        // E = G.C.Energy[id_grid];
        // if ( E < u_floor ){
        //   delta_u += u_floor - E;
        //   printf("###Thread Energy change  %f -> %f \n", E, u_floor );
        //   G.C.GasEnergy[id_grid] += delta_u;
        //   G.C.Energy[id_grid] += delta_u;
        // }
      }
    }
  }

}


void Apply_Gavity_Corrector( Grid3D &G ){

  Get_Gavity_Corrector( G, 0, G.Grav.nz_local );

  Add_Gravity_Corrector( G, 0, G.Grav.nz_local );




}
#endif



void Set_dt( Grid3D &G, bool &output_now, int n_step ){

  #ifdef COSMOLOGY
  Real delta_a_part;
  Real dt_courant, dt_gas, da_courant;


  #ifdef ONLY_PM
  delta_a_part = Get_Particles_da_cosmo( G );
  chprintf( " Delta_a particles: %f\n", delta_a_part );
  da_courant = delta_a_part;
  #else
  delta_a_part = Get_Particles_da_cosmo( G );
  dt_courant = G.H.dt;
  da_courant = G.Cosmo.Get_da_courant( dt_courant);
  chprintf( " Delta_a particles: %f   Delta_a gas: %f\n", delta_a_part, da_courant);
  da_courant = std::min(delta_a_part, da_courant);
  #endif


  // Real delta_a_part = G.Cosmo.max_delta_a;
  if( da_courant > G.Cosmo.max_delta_a){
    da_courant = G.Cosmo.max_delta_a;
    chprintf( " Seting max delta_a: %f\n", da_courant );
  }

  // Real da_min = delta_a_part / 5000;
  // if ( da_courant < da_min ){
  //   da_courant = da_min;
  //   chprintf( " Seting min delta_a: %f\n", da_courant );
  // }


  G.Cosmo.delta_a = da_courant;
  if ( (G.Cosmo.current_a + G.Cosmo.delta_a) >  G.Cosmo.next_output ){
    G.Cosmo.delta_a = G.Cosmo.next_output - G.Cosmo.current_a;
    output_now = true;
  }
  // if ( n_step%1 == 0)  output_now = true;

  Real dt = G.Cosmo.Get_dt_from_da( G.Cosmo.delta_a );
  G.Cosmo.dt_secs = dt * G.Cosmo.time_conversion;
  G.Cosmo.t_secs += G.Cosmo.dt_secs;


  G.Particles.dt = dt;

  dt_gas = G.Cosmo.Get_Cosmology_dt( G.Cosmo.delta_a );
  G.H.dt = dt_gas;

  chprintf( " Current_a: %f    delta_a: %f     da_courant: %f  dt:  %f\n", G.Cosmo.current_a, G.Cosmo.delta_a, da_courant, dt  );
  chprintf( " t_physical: %f Myr   dt_physical: %f Myr\n", G.Cosmo.t_secs/MYR, G.Cosmo.dt_secs/MYR );


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

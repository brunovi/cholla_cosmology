#ifdef GRAVITY

#include "../io.h"
#include "gravity_functions.h"
#include "../universal_constants.h"

#ifdef PARTICLES
#include"../particles/particles_dynamics.h"
#endif

#ifdef PARALLEL_OMP
#include"../parallel_omp.h"
#endif

#ifdef GRAVITY_CORRECTOR
#include "../correction_functions.h"
#endif

#ifdef MPI_CHOLLA
void Transfer_Potential_Boundaries_MPI( Grid3D &G, struct parameters P){


  G.Grav.TRANSFER_POTENTIAL_BOUNDARIES = true;
  chprintf( " Transfering Gravitational Potential Boundaries\n");
  G.Set_Boundary_Conditions(P);
  G.Grav.TRANSFER_POTENTIAL_BOUNDARIES = false;

}
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

Real Get_Density_Average_function( Grid3D &G, int g_start, int g_end ){
  // int nGHST = G.Parts.G.n_ghost_particles_grid;
  // int nx_g = Parts.G.nx_local + 2*nGHST;
  // int ny_g = Parts.G.ny_local + 2*nGHST;
  // int nz_g = Parts.G.nz_local + 2*nGHST;
  int nx = G.Grav.nx_local;
  int ny = G.Grav.ny_local;
  int nz = G.Grav.nz_local;
  int k, j, i, id;
  Real dens_avrg=0;
  for( k=g_start; k<g_end; k++){
    for( j=0; j<ny; j++){
      for( i=0; i<nx; i++){
        // id = (i+nGHST) + (j+nGHST)*nx_g + (k+nGHST)*nx_g*ny_g;
        id = (i) + (j)*nx + (k)*nx*ny;
        dens_avrg += G.Grav.F.density_h[id];
      }
    }
  }
  // dens_avrg /= (nx*ny*nz);
  return dens_avrg;
}

Real Get_Density_Average( Grid3D &G ){

  Real dens_sum;

  #ifndef PARALLEL_OMP
  dens_sum = Get_Density_Average_function( G, 0, G.Grav.nz_local );
  #else
  dens_sum = 0;
  Real dens_sum_all[N_OMP_THREADS];
  #pragma omp parallel num_threads( N_OMP_THREADS )
  {
    int omp_id, n_omp_procs;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local,  &g_start, &g_end  );
    dens_sum_all[omp_id] = Get_Density_Average_function( G, g_start, g_end );

  }
  for ( int i=0; i<N_OMP_THREADS; i++ ){
    dens_sum += dens_sum_all[i];
  }
  #endif

  return dens_sum /  ( G.Grav.nx_local * G.Grav.ny_local * G.Grav.nz_local);
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
void Compute_Gravitational_Potential( Grid3D &G, struct parameters P ){
  Real time_potential;

  #ifndef ONLY_PM
  Copy_Hydro_Density_to_Gravity( G );
  #endif

  #ifdef PARTICLES
  Get_Particles_Density_CIC( G, P );
  Copy_Particles_Density_to_Gravity( G );
  #endif

  Real dens_avrg;
  dens_avrg = Get_Density_Average( G );
  #ifdef MPI_CHOLLA
  dens_avrg = ReduceRealAvg( dens_avrg );
  #endif

  G.Grav.dens_avrg = dens_avrg;
  chprintf( " Density Average:  %f\n", dens_avrg);

  #ifdef CPU_TIME
  G.Timer.Start_Timer();
  #endif
  G.p_solver.Get_Potential( G.Grav );
  #ifndef GRAVITY_CPU
  Copy_Potential_To_Hydro_Grid( G );
  #endif
  #ifdef CPU_TIME
  G.Timer.End_and_Record_Time( 3 );
  #endif

}
#endif





void Extrapolate_Grav_Potential_function( Grid3D &G, int g_start, int g_end ){
  int n_ghost_pot = N_GHOST_POTENTIAL;
  int nx_pot = G.Grav.nx_local + 2*n_ghost_pot;
  int ny_pot = G.Grav.ny_local + 2*n_ghost_pot;
  int nz_pot = G.Grav.nz_local + 2*n_ghost_pot;

  Real *potential;
  int n_ghost_grid, nx_grid, ny_grid, nz_grid;
  #ifdef GRAVITY_CPU
  potential = G.Grav.F.potential_h;
  n_ghost_grid = N_GHOST_POTENTIAL;
  #else
  potential = G.C.Grav_potential;
  n_ghost_grid = G.H.n_ghost;
  #endif

  nx_grid = G.Grav.nx_local + 2*n_ghost_grid;
  ny_grid = G.Grav.ny_local + 2*n_ghost_grid;
  nz_grid = G.Grav.nz_local + 2*n_ghost_grid;

  int nGHST = n_ghost_grid - n_ghost_pot;
  Real pot_now, pot_prev, pot_extrp;
  int k, j, i, id_pot, id_grid;
  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_pot; j++ ){
      for ( i=0; i<nx_pot; i++ ){
        id_pot = i + j*nx_pot + k*nx_pot*ny_pot;
        id_grid = (i+nGHST) + (j+nGHST)*nx_grid + (k+nGHST)*nx_grid*ny_grid;
        pot_now = potential[id_grid]  ;
        if ( G.Grav.INITIAL ){
          pot_prev = pot_now;
          pot_extrp = pot_now;
        } else{
          pot_prev = G.Grav.F.potential_1_h[id_pot] ;
          #ifdef GRAVITY_CORRECTOR
          pot_extrp = pot_now;
          #else //GRAVITY_CORRECTOR
          pot_extrp = pot_now  + 0.5 * G.Grav.dt_now * ( pot_now - pot_prev  ) / G.Grav.dt_prev;
          #endif //GRAVITY_CORRECTOR
        }

        #ifdef COSMOLOGY
        pot_extrp *= G.Cosmo.current_a * G.Cosmo.current_a / G.Cosmo.phi_0_gas ;
        #endif

        potential[id_grid] = pot_extrp;
        G.Grav.F.potential_1_h[id_pot] = pot_now;
      }
    }
  }
}

void Extrapolate_Grav_Potential( Grid3D &G ){

  #ifndef PARALLEL_OMP

  Extrapolate_Grav_Potential_function( G, 0, G.Grav.nz_local + 2*N_GHOST_POTENTIAL );

  #else

  #pragma omp parallel num_threads( N_OMP_THREADS )
  {
    int omp_id, n_omp_procs;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local + 2*N_GHOST_POTENTIAL,  &g_start, &g_end  );

    Extrapolate_Grav_Potential_function( G, g_start, g_end );
  }
  #endif


  G.Grav.INITIAL = false;

}

#ifndef GRAVITY_CPU
void Copy_Potential_To_Hydro_Grid_function( Grid3D &G, int g_start, int g_end ){
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
  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<G.Grav.ny_local; j++ ){
      for ( i=0; i<G.Grav.nx_local; i++ ){
        id_pot  = (i+n_ghost_pot)  + (j+n_ghost_pot)*nx_pot   + (k+n_ghost_pot)*nx_pot*ny_pot;
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

void Copy_Potential_To_Hydro_Grid( Grid3D &G ){


  #ifndef PARALLEL_OMP

  Copy_Potential_To_Hydro_Grid_function( G, 0, G.Grav.nz_local );

  #else

  #pragma omp parallel num_threads( N_OMP_THREADS )
  {
    int omp_id, n_omp_procs;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local,  &g_start, &g_end  );

    Copy_Potential_To_Hydro_Grid_function( G, g_start, g_end );
  }
  #endif
}
#endif


// void Get_Potential_Difference( Grid3D &G ){
//
//   int nx_grav, ny_grav, nz_grav, nGHST_grav;
//   nGHST_grav = N_GHOST_POTENTIAL;
//   nx_grav = G.Grav.nx_local + 2*nGHST_grav;
//   ny_grav = G.Grav.ny_local + 2*nGHST_grav;
//   nz_grav = G.Grav.nz_local + 2*nGHST_grav;
//
//   int nx_grid, ny_grid, nz_grid, nGHST_grid;
//   nGHST_grid = G.H.n_ghost;
//   nx_grid = G.Grav.nx_local + 2*nGHST_grid;
//   ny_grid = G.Grav.ny_local + 2*nGHST_grid;
//   nz_grid = G.Grav.nz_local + 2*nGHST_grid;
//
//   Real pot_grid, pot_grav, error;
//
//   int nGHST = nGHST_grid - nGHST_grav;
//
//   int i, j, k, id_grid, id_grav;
//   for( k=0; k<nz_grav; k++){
//     for( j=0; j<ny_grav; j++){
//       for( i=0; i<nx_grav; i++){
//         // if ( i < nGHST_grav || i >= nx_grav - 2*nGHST_grav - 1) continue;
//         // if ( j < nGHST_grav || j >= ny_grav - 2*nGHST_grav - 1) continue;
//         // if ( k < nGHST_grav || k >= nz_grav - 2*nGHST_grav - 1) continue;
//
//         id_grav = (i) + (j)*nx_grav + (k)*nx_grav*ny_grav;
//         id_grid = (i+nGHST) + (j+nGHST)*nx_grid + (k+nGHST)*nx_grid*ny_grid;
//         pot_grid = G.C.Grav_potential[id_grid];
//         pot_grav = G.Grav.F.potential_h[id_grav];
//         // if( pot_grid == 0 ) std::cout << i << " "<< j << " "<< k << " "<< std::endl;
//         error = fabs( ( pot_grid - pot_grav ) / pot_grid );
//         if ( error>0.01) chprintf( " %f %f %f %d %d %d \n", pot_grav, pot_grid, error, i,j,k);
//       }
//     }
//   }
// }

#ifdef GRAVITY_CPU
void Get_Gavity_Field_Function( Grid3D &G, int g_start, int g_end ){

  int nx_grav, ny_grav, nz_grav;
  // nGHST_grav = G.Particles.G.n_ghost_particles_grid;
  nx_grav = G.Grav.nx_local;
  ny_grav = G.Grav.ny_local;
  nz_grav = G.Grav.nz_local;

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  Real *potential;

  #ifdef GRAVITY_CPU
  potential = G.Grav.F.potential_h;
  nGHST_grid = N_GHOST_POTENTIAL;
  #else
  potential = G.C.Grav_potential;
  nGHST_grid = G.H.n_ghost;
  #endif

  nx_grid = G.Grav.nx_local + 2*nGHST_grid;
  ny_grid = G.Grav.ny_local + 2*nGHST_grid;
  nz_grid = G.Grav.nz_local + 2*nGHST_grid;

  int nGHST = nGHST_grid ;

  Real dx, dy, dz;
  dx = G.Grav.dx;
  dy = G.Grav.dy;
  dz = G.Grav.dz;

  Real pot_factor;
  pot_factor = 1;

  Real phi_l, phi_r;
  int k, j, i, id_l, id_r, id;
  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i-1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r = (i+1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l = potential[id_l] * pot_factor;
        phi_r = potential[id_r] * pot_factor;
        G.Grav.F.gravity_x_h[id] = -0.5 * ( phi_r - phi_l ) / dx;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i + nGHST) + (j-1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r = (i + nGHST) + (j+1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l = potential[id_l] * pot_factor;
        phi_r = potential[id_r] * pot_factor;
        G.Grav.F.gravity_y_h[id] = -0.5 * ( phi_r - phi_l ) / dy;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i + nGHST) + (j + nGHST)*nx_grid + (k-1 + nGHST)*ny_grid*nz_grid;
        id_r = (i + nGHST) + (j + nGHST)*nx_grid + (k+1 + nGHST)*ny_grid*nz_grid;
        phi_l = potential[id_l] * pot_factor;
        phi_r = potential[id_r] * pot_factor;
        G.Grav.F.gravity_z_h[id] = -0.5 * ( phi_r - phi_l ) / dz;
      }
    }
  }
}

void Add_Gavity_To_Hydro_Function( Grid3D &G, int g_start, int g_end ){

  int nx_grav, ny_grav, nz_grav;
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
  Real d_0, vx_0, vy_0, vz_0;

  #ifdef GRAVITY_DELTA_EK
  Real Ek_0, Ek_1;
  #endif

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

        d = G.C.density[id_grid];
        vx = G.C.momentum_x[id_grid] / d;
        vy = G.C.momentum_y[id_grid] / d;
        vz = G.C.momentum_z[id_grid] / d;

        gx = G.Grav.F.gravity_x_h[id_grav];
        gy = G.Grav.F.gravity_y_h[id_grav];
        gz = G.Grav.F.gravity_z_h[id_grav];

        #ifdef GRAVITY_DELTA_EK
        Ek_0 = 0.5 * d * ( vx*vx + vy*vy + vz*vz );
        #endif


        G.C.momentum_x[id_grid] += 0.5 *  G.H.dt * gx * ( d_0 + d );
        G.C.momentum_y[id_grid] += 0.5 *  G.H.dt * gy * ( d_0 + d );
        G.C.momentum_z[id_grid] += 0.5 *  G.H.dt * gz * ( d_0 + d );
        #ifndef GRAVITY_DELTA_EK
        G.C.Energy[id_grid] += 0.5 * G.H.dt * ( gx*(d_0*vx_0 + d*vx) + gy*(d_0*vy_0 + d*vy) + gz*(d_0*vz_0 + d*vz)   );
        #else
        vx = G.C.momentum_x[id_grid] / d;
        vy = G.C.momentum_y[id_grid] / d;
        vz = G.C.momentum_z[id_grid] / d;
        Ek_1 = 0.5 * d * ( vx*vx + vy*vy + vz*vz );
        G.C.Energy[id_grid] += Ek_1 - Ek_0;
        #endif
      }
    }
  }

}

void Add_Gavity_To_Hydro( Grid3D &G ){

  #ifndef PARALLEL_OMP
  Get_Gavity_Field_Function( G, 0, G.Grav.nz_local );
  Add_Gavity_To_Hydro_Function( G, 0, G.Grav.nz_local );
  #else
  #pragma omp parallel num_threads( N_OMP_THREADS )
  {
    int omp_id, n_omp_procs;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local,  &g_start, &g_end  );

    Get_Gavity_Field_Function( G, g_start, g_end );
    #pragma omp barrier
    Add_Gavity_To_Hydro_Function( G, g_start, g_end );
  }
  #endif
}




#endif







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
        // delta = fabs( (G.Grav.F.gravity_x_h_prev[id] - G.Grav.F.gravity_x_h[id]) / G.Grav.F.gravity_x_h_prev[id] );
        // if ( delta > delta_max ) std::cout << "### Grav X: " << delta<< " delta_max: " << delta_max << std::endl;
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
        // delta = fabs( (G.Grav.F.gravity_y_h_prev[id] - G.Grav.F.gravity_y_h[id]) / G.Grav.F.gravity_y_h_prev[id] );
        // if ( delta > delta_max ) std::cout << "### Grav Y: " << delta<< " delta_max: " << delta_max << std::endl;
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
        // delta = fabs( (G.Grav.F.gravity_z_h_prev[id] - G.Grav.F.gravity_z_h[id]) / G.Grav.F.gravity_z_h_prev[id] );
        // if ( delta > delta_max ) std::cout << "### Grav Z: " << delta<< " delta_max: " << delta_max << std::endl;
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

        // // Gravity first ordef
        // G.C.momentum_x[id_grid] +=  G.H.dt * d_0 * gx_0;
        // G.C.momentum_y[id_grid] +=  G.H.dt * d_0 * gy_0;
        // G.C.momentum_z[id_grid] +=  G.H.dt * d_0 * gz_0;
        // G.C.Energy[id_grid] +=  G.H.dt * d_0 * ( vx_0*gx_0 + vy_0*gy_0 + vz_0*gz_0 );

        d = G.C.density[id_grid];
        vx = G.C.momentum_x[id_grid] / d;
        vy = G.C.momentum_y[id_grid] / d;
        vz = G.C.momentum_z[id_grid] / d;
        gx = G.Grav.F.gravity_x_h[id_grav];
        gy = G.Grav.F.gravity_y_h[id_grav];
        gz = G.Grav.F.gravity_z_h[id_grav];


        G.C.momentum_x[id_grid] += 0.5 *  G.H.dt * d * gx;
        G.C.momentum_y[id_grid] += 0.5 *  G.H.dt * d * gy;
        G.C.momentum_z[id_grid] += 0.5 *  G.H.dt * d * gz;
        G.C.Energy[id_grid] += 0.5 * G.H.dt * d * ( vx*gx + vy*gy + vz*gz );

        G.C.momentum_x[id_grid] += 0.5 * G.H.dt * d_0 * gx_0;
        G.C.momentum_y[id_grid] += 0.5 * G.H.dt * d_0 * gy_0;
        G.C.momentum_z[id_grid] += 0.5 * G.H.dt * d_0 * gz_0;
        G.C.Energy[id_grid] += 0.5 *  G.H.dt * d_0 * ( vx_0*gx_0 + vy_0*gy_0 + vz_0*gz_0 );

        // G.C.momentum_x[id_grid] -= 0.5 * G.H.dt * d_0 * gx_0;
        // G.C.momentum_y[id_grid] -= 0.5 * G.H.dt * d_0 * gy_0;
        // G.C.momentum_z[id_grid] -= 0.5 * G.H.dt * d_0 * gz_0;
        // G.C.Energy[id_grid] -= 0.5 * G.H.dt * d_0 * ( vx_0*gx_0 + vy_0*gy_0 + vz_0*gz_0 );

        // G.C.momentum_x[id_grid] += 0.5 * G.H.dt * d * gx;
        // G.C.momentum_y[id_grid] += 0.5 * G.H.dt * d * gy;
        // G.C.momentum_z[id_grid] += 0.5 * G.H.dt * d * gz;
        // G.C.Energy[id_grid] += 0.5 * G.H.dt * ( d*vx*gx + d*vy*gy + d*vz*gz );





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


void Apply_Gavity_Corrector( Grid3D &G, struct parameters P ){

  Transfer_Potential_Boundaries_MPI( G, P);

  #ifndef PARALLEL_OMP
  Get_Gavity_Corrector( G, 0, G.Grav.nz_local );
  Add_Gravity_Corrector( G, 0, G.Grav.nz_local );
  #else

  #pragma omp parallel num_threads( N_OMP_THREADS )
  {
    int omp_id, n_omp_procs;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    // #pragma omp barrier

    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local,  &g_start, &g_end  );

    Get_Gavity_Corrector( G, g_start, g_end );
    #pragma omp barrier
    Add_Gravity_Corrector( G, g_start, g_end );
  }
  #endif

  #ifdef DE
  Sync_Energies_3D_Host( G );
  #endif




}
#endif



void Set_dt( Grid3D &G, bool &output_now, int n_step ){

  #ifdef COSMOLOGY
  Real da_particles;
  Real dt_hydro, da_hydro, dt_courant, da_courant;

  chprintf( " Current_z: %f \n", G.Cosmo.current_z );

  #ifdef ONLY_PM
  da_particles = Get_Particles_da_cosmo( G );
  chprintf( " Delta_a particles: %f\n", da_particles );
  da_courant = da_particles;
  #else
  da_particles = Get_Particles_da_cosmo( G );
  dt_hydro = G.H.dt;
  da_hydro = G.Cosmo.Get_da_courant( dt_hydro);
  chprintf( " Delta_a particles: %f   Delta_a gas: %f\n", da_particles, da_hydro);
  da_courant = std::min(da_particles, da_hydro);
  #endif

  if( da_courant > G.Cosmo.max_delta_a){
    da_courant = G.Cosmo.max_delta_a;
    chprintf( " Seting max delta_a: %f\n", da_courant );
  }

  Real da_min = da_particles / 10;
  if ( da_courant < da_min ){
    da_courant = da_min;
    chprintf( " Seting min delta_a: %f\n", da_courant );
  }

  #ifdef COOLING_GRACKLE
  if ( fabs(G.Cosmo.current_a + da_courant - G.Cool.a_UVB_on) < 0.001 ){
    da_courant /= 10;
    chprintf( " Starting UVB. Limiting delta_a:  %f \n", da_courant);
  }
  #endif


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

  dt_hydro = G.Cosmo.Get_Cosmology_dt( G.Cosmo.delta_a );
  G.H.dt = dt_hydro;

  chprintf( " Current_a: %f    delta_a: %f     da_courant: %f  dt:  %f\n", G.Cosmo.current_a, G.Cosmo.delta_a, da_courant, dt  );
  chprintf( " t_physical: %f Myr   dt_physical: %f Myr\n", G.Cosmo.t_secs/MYR, G.Cosmo.dt_secs/MYR );

  #else //COSMOLOGY
  Real dt_hydro, dt_courant;
  dt_hydro = G.H.dt;
  #ifndef PARTICLES
  chprintf( " Delta_t hydro: %f\n", dt_hydro);
  #endif

  #ifdef PARTICLES
  Real dt_particles;
  dt_particles = Get_Particles_dt( G.Particles );
  chprintf( " Delta_t hydro: %f    Delta_t particles: %f \n", dt_hydro, dt_particles );
  dt_courant = std::min( dt_hydro, dt_particles );
  chprintf( " Delta_t courant: %f\n", dt_courant);
  G.H.dt = dt_courant;
  G.Particles.dt = dt_courant;

  #endif //PARTICLES
  #endif //COSMOLOGY



  // Set times for potential extrapolation
  if ( G.Grav.INITIAL ){
    G.Grav.dt_prev = G.H.dt;
    G.Grav.dt_now = G.H.dt;
  }else{
    G.Grav.dt_prev = G.Grav.dt_now;
    G.Grav.dt_now = G.H.dt;
  }


  // #ifdef PARTICLES

  // #endif

}


#endif

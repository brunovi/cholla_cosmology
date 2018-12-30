#ifdef PARTICLES

#include"particles_dynamics.h"
#include"../io.h"

#ifdef MPI_CHOLLA
#include "../mpi_routines.h"
#endif



void Get_Particles_Acceleration( Grid3D &G, part_int_t p_start, part_int_t p_end, int g_start, int g_end ){
  Get_Gravity_Field( G, g_start, g_end );
  // Get_Gravity_Field_order4( G, g_start, g_end );
  Get_Gravity_CIC( G.Particles, p_start, p_end );
}


void Advance_Particles_step1( Particles_3D &Particles, part_int_t p_start, part_int_t p_end ){

  part_int_t pID;
  Real dt = Particles.dt;
  // Advance velocities by half a step
  for ( pID=p_start; pID<p_end; pID++ ){
    Particles.vel_x[pID] += 0.5 * dt * Particles.grav_x[pID];
    Particles.vel_y[pID] += 0.5 * dt * Particles.grav_y[pID];
    Particles.vel_z[pID] += 0.5 * dt * Particles.grav_z[pID];
  }

  //Advance Posiotions using advanced velocities
  for ( pID=p_start; pID<p_end; pID++ ){
    Particles.pos_x[pID] += dt * Particles.vel_x[pID];
    Particles.pos_y[pID] += dt * Particles.vel_y[pID];
    Particles.pos_z[pID] += dt * Particles.vel_z[pID];
  }
}

void Advance_Particles_step2( Particles_3D &Particles, part_int_t p_start, part_int_t p_end ){

  part_int_t pID;
  Real dt = Particles.dt;
  // Advance velocities by half a step
  for ( pID=p_start; pID<p_end; pID++ ){
    Particles.vel_x[pID] += 0.5 * dt * Particles.grav_x[pID];
    Particles.vel_y[pID] += 0.5 * dt * Particles.grav_y[pID];
    Particles.vel_z[pID] += 0.5 * dt * Particles.grav_z[pID];
  }
}
#ifdef COSMOLOGY
void Advance_Particles_step1_cosmo_LeapFrog( Particles_3D &Particles, Cosmology &Cosmo, part_int_t p_start, part_int_t p_end ){

  part_int_t pIndx;
  Real da = Cosmo.delta_a;
  Real current_a = Cosmo.current_a;

  Real scale_factor = Scale_Function( current_a, Cosmo.Omega_M, Cosmo.Omega_L, Cosmo.Omega_K ) / Cosmo.H0 * Cosmo.cosmo_h;
  Real scale_factor_1 = Scale_Function( current_a + 0.5*da, Cosmo.Omega_M, Cosmo.Omega_L, Cosmo.Omega_K  ) / Cosmo.H0 * Cosmo.cosmo_h;
  Real a2_inv = 1./( ( current_a + 0.5*da )*( current_a + 0.5*da ));
  // Advance velocities by half a step
  Real pos_x, vel_x, grav_x;
  Real pos_y, vel_y, grav_y;
  Real pos_z, vel_z, grav_z;
  for ( pIndx=p_start; pIndx<p_end; pIndx++ ){
    pos_x = Particles.pos_x[pIndx];
    pos_y = Particles.pos_y[pIndx];
    pos_z = Particles.pos_z[pIndx];
    vel_x = Particles.vel_x[pIndx];
    vel_y = Particles.vel_y[pIndx];
    vel_z = Particles.vel_z[pIndx];
    grav_x = Particles.grav_x[pIndx];
    grav_y = Particles.grav_y[pIndx];
    grav_z = Particles.grav_z[pIndx];

    vel_x += scale_factor * da * grav_x;
    vel_y += scale_factor * da * grav_y;
    vel_z += scale_factor * da * grav_z;

    pos_x += a2_inv * scale_factor_1 * da * vel_x;
    pos_y += a2_inv * scale_factor_1 * da * vel_y;
    pos_z += a2_inv * scale_factor_1 * da * vel_z;

    Particles.pos_x[pIndx] = pos_x;
    Particles.pos_y[pIndx] = pos_y;
    Particles.pos_z[pIndx] = pos_z;

    Particles.vel_x[pIndx] = vel_x;
    Particles.vel_y[pIndx] = vel_y;
    Particles.vel_z[pIndx] = vel_z;
  }
}

Real da_dt( Real a, Real Omega_M, Real Omega_L){
  Real a2 = a * a ;
  return sqrt( Omega_M/a + a2*Omega_L);
}


void Advance_Particles_step1_cosmo( Particles_3D &Particles, Cosmology &Cosmo, part_int_t p_start, part_int_t p_end ){

  part_int_t pIndx;
  Real da = Cosmo.delta_a;
  Real current_a = Cosmo.current_a;

  Real scale_factor = Scale_Function( current_a, Cosmo.Omega_M, Cosmo.Omega_L, Cosmo.Omega_K ) / Cosmo.H0 * Cosmo.cosmo_h;
  Real scale_factor_1 = Scale_Function( current_a + 0.5*da, Cosmo.Omega_M, Cosmo.Omega_L, Cosmo.Omega_K  ) / Cosmo.H0 * Cosmo.cosmo_h;
  Real a2_inv = 1./( ( current_a + 0.5*da )*( current_a + 0.5*da ));
  // Advance velocities by half a step
  Real pos_x, vel_x, grav_x;
  Real pos_y, vel_y, grav_y;
  Real pos_z, vel_z, grav_z;
  for ( pIndx=p_start; pIndx<p_end; pIndx++ ){
    pos_x = Particles.pos_x[pIndx];
    pos_y = Particles.pos_y[pIndx];
    pos_z = Particles.pos_z[pIndx];
    vel_x = Particles.vel_x[pIndx];
    vel_y = Particles.vel_y[pIndx];
    vel_z = Particles.vel_z[pIndx];
    grav_x = Particles.grav_x[pIndx];
    grav_y = Particles.grav_y[pIndx];
    grav_z = Particles.grav_z[pIndx];

    vel_x += scale_factor * 0.5 * da * grav_x;
    vel_y += scale_factor * 0.5 * da * grav_y;
    vel_z += scale_factor * 0.5 * da * grav_z;

    pos_x += a2_inv * scale_factor_1 * da * vel_x;
    pos_y += a2_inv * scale_factor_1 * da * vel_y;
    pos_z += a2_inv * scale_factor_1 * da * vel_z;

    Particles.pos_x[pIndx] = pos_x;
    Particles.pos_y[pIndx] = pos_y;
    Particles.pos_z[pIndx] = pos_z;

    Particles.vel_x[pIndx] = vel_x;
    Particles.vel_y[pIndx] = vel_y;
    Particles.vel_z[pIndx] = vel_z;
  }
}

void Advance_Particles_step2_cosmo( Particles_3D &Particles, Cosmology &Cosmo, part_int_t p_start, part_int_t p_end ){

  part_int_t pIndx;
  Real da = Cosmo.delta_a;
  Real current_a = Cosmo.current_a;

  Real scale_factor = Scale_Function( current_a , Cosmo.Omega_M, Cosmo.Omega_L, Cosmo.Omega_K ) / Cosmo.H0 * Cosmo.cosmo_h;
  // Advance velocities by half a step
  Real grav_x;
  Real grav_y;
  Real grav_z;
  for ( pIndx=p_start; pIndx<p_end; pIndx++ ){
    grav_x = Particles.grav_x[pIndx];
    grav_y = Particles.grav_y[pIndx];
    grav_z = Particles.grav_z[pIndx];
    Particles.vel_x[pIndx] += scale_factor * 0.5 * da * grav_x;
    Particles.vel_y[pIndx] += scale_factor * 0.5 * da * grav_y;
    Particles.vel_z[pIndx] += scale_factor * 0.5 * da * grav_z;

  }

}
#endif

Real Get_Particles_dt( Particles_3D &Particles ){
  part_int_t pID;
  Real dt, dt_min, vel;
  dt_min = 1e100;

  for ( pID=0; pID<Particles.n_local; pID++ ){
    vel = fabs(Particles.vel_x[pID]);
    if ( vel > 0){
      dt = Particles.G.dx / vel;
      dt_min = std::min( dt_min, dt);
    }
    vel = fabs(Particles.vel_y[pID]);
    if ( vel > 0){
      dt = Particles.G.dy / vel;
      dt_min = std::min( dt_min, dt);
    }
    vel = fabs(Particles.vel_z[pID]);
    if ( vel > 0){
      dt = Particles.G.dz / vel;
      dt_min = std::min( dt_min, dt);
    }
  }

  #ifdef MPI_CHOLLA
  Real dt_min_global = ReduceRealMin(dt_min);
  dt_min = dt_min_global;
  #endif

  return Particles.C_cfl * dt_min;
}

#ifdef COSMOLOGY
Real Get_Particles_da_cosmo( Grid3D &G ){
  part_int_t pID;
  Real da, da_min, vel;
  da_min = 1e100;
  Real scale_factor = Scale_Function( G.Cosmo.current_a , G.Cosmo.Omega_M, G.Cosmo.Omega_L, G.Cosmo.Omega_K  ) / G.Cosmo.H0 * G.Cosmo.cosmo_h;
  Real a2 = ( G.Cosmo.current_a )*( G.Cosmo.current_a  );
  Real vel_factor = a2 / scale_factor;


  for ( pID=0; pID<G.Particles.n_local; pID++ ){
    vel = fabs(G.Particles.vel_x[pID]);
    if ( vel > 0){
      da = G.Particles.G.dx * vel_factor / vel;
      da_min = std::min( da_min, da);
    }
    vel = fabs(G.Particles.vel_y[pID]);
    if ( vel > 0){
      da = G.Particles.G.dy * vel_factor / vel;
      da_min = std::min( da_min, da);
    }
    vel = fabs(G.Particles.vel_z[pID]);
    if ( vel > 0){
      da = G.Particles.G.dz * vel_factor / vel;
      da_min = std::min( da_min, da);
    }
  }

  #ifdef MPI_CHOLLA
  Real da_min_global = ReduceRealMin(da_min);
  da_min = da_min_global;
  #endif
  return G.Particles.C_cfl * da_min;
}

Real Get_Particles_dt_cosmo( Grid3D &G ){
  part_int_t pID;
  Real da, da_min, vel;
  da_min = 1e100;
  Real scale_factor = Scale_Function( G.Cosmo.current_a , G.Cosmo.Omega_M, G.Cosmo.Omega_L, G.Cosmo.Omega_K  ) / G.Cosmo.H0 * G.Cosmo.cosmo_h;
  Real a2 = ( G.Cosmo.current_a )*( G.Cosmo.current_a  );
  Real vel_factor = a2 / scale_factor;


  for ( pID=0; pID<G.Particles.n_local; pID++ ){
    vel = fabs(G.Particles.vel_x[pID]);
    if ( vel > 0){
      da = G.Particles.G.dx * vel_factor / vel;
      da_min = std::min( da_min, da);
    }
    vel = fabs(G.Particles.vel_y[pID]);
    if ( vel > 0){
      da = G.Particles.G.dy * vel_factor / vel;
      da_min = std::min( da_min, da);
    }
    vel = fabs(G.Particles.vel_z[pID]);
    if ( vel > 0){
      da = G.Particles.G.dz * vel_factor / vel;
      da_min = std::min( da_min, da);
    }
  }

  Real dt_cosmo = G.Cosmo.Get_Cosmology_dt( da_min );

  //
  // #ifdef MPI_CHOLLA
  // Real da_min_global = ReduceRealMin(da_min);
  // da_min = da_min_global;
  // #endif
  // return G.Particles.C_cfl * da_min;
  return G.Particles.C_cfl * dt_cosmo;
}
#endif

void Update_Particles( Grid3D &G, int step ){


  #ifdef REVERT_STEP
  if ( step == 1) Copy_Particles_Vectors( G );
  #endif

  #ifdef CPU_TIME
  G.Timer.Start_Timer();
  #endif

  #ifndef GRAVITY_OMP
  if ( step == 1 ){
    #ifndef COSMOLOGY
    Advance_Particles_step1( G.Particles, 0, G.Particles.n_local );
    #else
    Advance_Particles_step1_cosmo( G.Particles, G.Cosmo, 0, G.Particles.n_local );
    #endif
  }
  else if ( step == 2 ){
    Get_Particles_Acceleration( G, 0, G.Particles.n_local, 0, G.Particles.G.nz_local + 2*G.Particles.G.n_ghost_particles_grid  );
    #ifndef COSMOLOGY
    Advance_Particles_step2( G.Particles, 0, G.Particles.n_local);
    #else
    Advance_Particles_step2_cosmo( G.Particles, G.Cosmo, 0, G.Particles.n_local);
    #endif
  }
  #endif

  #ifdef GRAVITY_OMP

  #pragma omp parallel num_threads( N_OMP_GRAVITY_THREADS )
  {
    int omp_id, n_omp_procs;
    part_int_t p_start, p_end;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();

    Get_OMP_Indxs( G.Particles.n_local, N_OMP_GRAVITY_THREADS, omp_id, G.Particles.G.nz_local + 2*G.Particles.G.n_ghost_particles_grid, &p_start, &p_end, &g_start, &g_end );

    if ( step == 1 ){
      #ifndef COSMOLOGY
      Advance_Particles_step1( G.Particles, p_start, p_end );
      #else
      Advance_Particles_step1_cosmo( G.Particles, G.Cosmo, p_start, p_end );
      #endif
    }
    else if ( step == 2 ){
      Get_Particles_Acceleration( G, p_start, p_end, g_start, g_end );
      #pragma omp barrier
      #ifndef COSMOLOGY
      Advance_Particles_step2( G.Particles, p_start, p_end );
      #else
      Advance_Particles_step2_cosmo( G.Particles, G.Cosmo, p_start, p_end );
      #endif
    }

  }
  #endif

  #ifdef CPU_TIME
  if ( step == 1) G.Timer.End_and_Record_Time(6);
  if ( step == 2) G.Timer.End_and_Record_Time(7);
  #endif

  part_int_t n_total;
  if ( step == 2 ){
    // MPI_Barrier(world);
    n_total = ReducePartIntSum( G.Particles.n_local );
    chprintf( " Total Particles: %ld\n", n_total );
    if ( n_total != G.Particles.n_total_0) {
      chprintf( " WARNING:  Lost %ld particles\n", G.Particles.n_total_0 - n_total );
      // break;
    }
  }

}

#ifdef REVERT_STEP
void Copy_Particles_Vectors( Grid3D &G ){

  #ifdef PARTICLE_IDS
  G.Particles.partIDs_0.clear();
  G.Particles.partIDs_0 = G.Particles.partIDs;
  #endif
  #ifndef SINGLE_PARTICLE_MASS
  G.Particles.mass_0.clear();
  G.Particles.mass_0 = G.Particles.mass;
  #endif
  G.Particles.pos_x_0.clear();
  G.Particles.pos_y_0.clear();
  G.Particles.pos_z_0.clear();
  G.Particles.vel_x_0.clear();
  G.Particles.vel_y_0.clear();
  G.Particles.vel_z_0.clear();
  G.Particles.grav_x_0.clear();
  G.Particles.grav_y_0.clear();
  G.Particles.grav_z_0.clear();

  G.Particles.pos_x_0 = G.Particles.pos_x;
  G.Particles.pos_y_0 = G.Particles.pos_y;
  G.Particles.pos_z_0 = G.Particles.pos_z;
  G.Particles.vel_x_0 = G.Particles.vel_x;
  G.Particles.vel_y_0 = G.Particles.vel_y;
  G.Particles.vel_z_0 = G.Particles.vel_z;
  G.Particles.grav_x_0 = G.Particles.grav_x;
  G.Particles.grav_y_0 = G.Particles.grav_y;
  G.Particles.grav_z_0 = G.Particles.grav_z;



}
#endif


#endif

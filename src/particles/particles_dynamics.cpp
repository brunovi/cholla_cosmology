#ifdef PARTICLES

#include"particles_dynamics.h"
#include"../io.h"

// void Advance_Particles_LeapFrog( Particles_3D &Particles ){
//
//
// }

void Get_Particles_Acceleration( Grid3D &G ){
  Get_Gravity_Field( G );
  Get_Gravity_CIC( G.Particles );
}


void Advance_Particles_step1( Particles_3D &Particles ){

  part_int_t pID;
  Real dt = Particles.dt;
  // Advance velocities by half a step
  for ( pID=0; pID<Particles.n_local; pID++ ){
    Particles.vel_x[pID] += 0.5 * dt * Particles.grav_x[pID];
    Particles.vel_y[pID] += 0.5 * dt * Particles.grav_y[pID];
    Particles.vel_z[pID] += 0.5 * dt * Particles.grav_z[pID];
  }

  //Advance Posiotions using advanced velocities
  for ( pID=0; pID<Particles.n_local; pID++ ){
    Particles.pos_x[pID] += dt * Particles.vel_x[pID];
    Particles.pos_y[pID] += dt * Particles.vel_y[pID];
    Particles.pos_z[pID] += dt * Particles.vel_z[pID];
  }
}

void Advance_Particles_step2( Particles_3D &Particles ){

  part_int_t pID;
  Real dt = Particles.dt;
  // Advance velocities by half a step
  for ( pID=0; pID<Particles.n_local; pID++ ){
    Particles.vel_x[pID] += 0.5 * dt * Particles.grav_x[pID];
    Particles.vel_y[pID] += 0.5 * dt * Particles.grav_y[pID];
    Particles.vel_z[pID] += 0.5 * dt * Particles.grav_z[pID];
  }
}

float Get_Particles_dt( Particles_3D &Particles ){
  part_int_t pID;
  Real dt, dt_min, vel;
  dt_min = 1e100;

  for ( pID=0; pID<Particles.n_local; pID++ ){
    vel = abs(Particles.vel_x[pID]);
    if ( vel > 0){
      dt = Particles.G.dx / vel;
      dt_min = std::min( dt_min, dt);
    }
    vel = abs(Particles.vel_y[pID]);
    if ( vel > 0){
      dt = Particles.G.dy / vel;
      dt_min = std::min( dt_min, dt);
    }
    vel = abs(Particles.vel_z[pID]);
    if ( vel > 0){
      dt = Particles.G.dz / vel;
      dt_min = std::min( dt_min, dt);
    }
  }
  return C_cfl * dt_min;
}


float Update_Particles( Grid3D &G, int step ){

  if ( step == 1 ){
    Advance_Particles_step1( G.Particles );
    return 0;
  }
  if ( step == 2 ){
    Get_Particles_Acceleration( G );
    Advance_Particles_step2( G.Particles );
    return 0;
  }
}

#endif

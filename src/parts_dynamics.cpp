#ifdef PARTICLES

#include"parts_dynamics.h"




void Part3D::Advance_Particle_Step_LeapFrog( partID_t pIndx, Real dt, bool only_vel ){
  if ( only_vel){
    vel_x[pIndx] += dt*grav_x[pIndx];
    vel_y[pIndx] += dt*grav_y[pIndx];
    vel_z[pIndx] += dt*grav_z[pIndx];
  }
  else{
    // Update velocities using current forces
    vel_x[pIndx] += dt*grav_x[pIndx];
    vel_y[pIndx] += dt*grav_y[pIndx];
    vel_z[pIndx] += dt*grav_z[pIndx];
    // Update positions using updated velocities
    pos_x[pIndx] += dt*vel_x[pIndx];
    pos_y[pIndx] += dt*vel_y[pIndx];
    pos_z[pIndx] += dt*vel_z[pIndx];
  }

  // chprintf("%f ", dt);
}

#ifdef COSMOLOGY

void Part3D::Advance_Particle_Step_LeapFrog_cosmo( partID_t pIndx, bool only_vel ){
  // Real vel_proper, pos_proper;
  Real vel_comov, pos_comov, grav, vel_peculiar;

  Real delta_a_local = delta_a;
  Real grav_factor = 1.0;

  // if ( only_vel ) delta_a_local = -0.5 * delta_a;

  // Real scale_function   = scale_funct( current_a );
  // Real scale_function_1 = scale_funct( current_a + 0.5*delta_a_local );
  Real scale_function   = scale_funct( current_a )/H0*cosmo_h;
  Real scale_function_1 = scale_funct( current_a + 0.5*delta_a_local )/H0*cosmo_h;
  Real a2_inv = 1. / ( (current_a+0.5*delta_a_local) * (current_a+0.5*delta_a_local) );


  //X axis
  // pos_comov = pos_x[pIndx]/r_0 ;
  pos_comov = pos_x[pIndx] ;
  vel_peculiar = vel_x[pIndx] ;    //From GADGET guide
  vel_comov = vel_peculiar;
  grav = grav_x[pIndx] * grav_factor;
  //update velocity
  vel_comov += scale_function*grav*delta_a_local;
  //update position
  pos_comov += a2_inv*scale_function_1*vel_comov*delta_a_local ;
  // pos_x[pIndx] = pos_comov*r_0 ;
  //
  //
  //
  //grav was NAN and vel_comov was NAN (because of grav)
  /*
  if(pIndx==0)
  {
    printf("Adv Part procID %d scale fnc %e %e a2_inv %e grav %e delta_a_local %e vel_comov %e pos_x %e\n",procID,scale_function,scale_function_1,a2_inv,grav,delta_a_local,vel_comov,pos_x[pIndx]);
    fflush(stdout);
  }
 */

  pos_x[pIndx] = pos_comov ;
  vel_x[pIndx] = vel_comov ;



  //Y axis
  // pos_comov = pos_y[pIndx]/r_0 ;
  pos_comov = pos_y[pIndx];
  vel_peculiar = vel_y[pIndx] ;
  vel_comov = vel_peculiar;
  grav = grav_y[pIndx] * grav_factor;
  //update velocity
  vel_comov += scale_function*grav*delta_a_local;
  //update position
  pos_comov += a2_inv*scale_function_1*vel_comov*delta_a_local ;
  // pos_y[pIndx] = pos_comov*r_0 ;
  pos_y[pIndx] = pos_comov ;
  vel_y[pIndx] = vel_comov;

  //Z axis
  // pos_comov = pos_z[pIndx]/r_0 ;
  pos_comov = pos_z[pIndx] ;
  vel_peculiar = vel_z[pIndx] ;   //From GADGET guide
  vel_comov = vel_peculiar;
  grav = grav_z[pIndx] * grav_factor;
  //update velocity
  vel_comov += scale_function*grav*delta_a_local;
  //update position
  pos_comov += a2_inv*scale_function_1*vel_comov*delta_a_local  ;
  // pos_z[pIndx] = pos_comov*r_0 ;
  pos_z[pIndx] = pos_comov ;
  vel_z[pIndx] = vel_comov;
}
#endif


void Part3D::Advance_Step_LeapFrog( Real dt, bool only_vel, partID_t p_indx_start, partID_t p_indx_end  ){

  partID_t pIndx;

  for ( pIndx=p_indx_start; pIndx<p_indx_end; pIndx++ ){

    #ifdef COSMOLOGY
    Advance_Particle_Step_LeapFrog_cosmo( pIndx, only_vel );
    #else
    Advance_Particle_Step_LeapFrog( pIndx, dt,  only_vel );
    #endif
  }


}




#endif //PARTICLES

#ifdef PARTICLES

#include"particles_dynamics.h"
#include"../io.h"

float Advance_Particles( Grid3D &G ){

  Get_Gavity_Field( G );
  Get_Gravity_CIC( G.Particles );

  return 0;
}

#endif

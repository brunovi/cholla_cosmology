#ifdef COSMOLOGY

#include"cosmology_units.h"

void Change_DM_Frame_System( Grid3D &G, bool forward ){
  part_int_t pIndx;
  Real vel_factor;
  if (forward ) vel_factor = G.Cosmo.current_a ;
  else vel_factor =  1./G.Cosmo.current_a;
  for ( pIndx=0; pIndx<G.Particles.n_local; pIndx++ ){
    G.Particles.vel_x[pIndx] *= vel_factor;
    G.Particles.vel_y[pIndx] *= vel_factor;
    G.Particles.vel_z[pIndx] *= vel_factor;
  }
}

void Change_Cosmological_Frame_Sytem( Grid3D &G, bool forward ){

  if (forward) chprintf( " Converting to Cosmological Comoving System\n");
  else chprintf( " Converting to Cosmological Physical System\n");
  Change_DM_Frame_System( G, forward );
}


#endif

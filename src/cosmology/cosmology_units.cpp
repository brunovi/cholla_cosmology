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

void Change_GAS_Frame_System( Grid3D &G, bool forward ){
  Real dens_factor, momentum_factor, energy_factor;
  if ( forward ){
    dens_factor = 1 / G.Cosmo.rho_0_gas;
    momentum_factor = 1 / G.Cosmo.rho_0_gas / G.Cosmo.v_0_gas;
    energy_factor = 1 / G.Cosmo.rho_0_gas / G.Cosmo.v_0_gas / G.Cosmo.v_0_gas;
  }
  else{
    dens_factor = G.Cosmo.rho_0_gas;
    momentum_factor =  G.Cosmo.rho_0_gas * G.Cosmo.v_0_gas;
    energy_factor =  G.Cosmo.rho_0_gas * G.Cosmo.v_0_gas * G.Cosmo.v_0_gas;
  }
  int k, j, i, id;
  for (k=0; k<G.H.nz; k++) {
    for (j=0; j<G.H.ny; j++) {
      for (i=0; i<G.H.nx; i++) {
        id = i + j*G.H.nx + k*G.H.nx*G.H.ny;
        G.C.density[id] = G.C.density[id] * dens_factor ;
        G.C.momentum_x[id] = G.C.momentum_x[id] *  momentum_factor ;
        G.C.momentum_y[id] = G.C.momentum_y[id] *  momentum_factor ;
        G.C.momentum_z[id] = G.C.momentum_z[id] *  momentum_factor ;
        G.C.Energy[id] = G.C.Energy[id] * energy_factor ;

        #ifdef DE
        G.C.GasEnergy[id] = G.C.GasEnergy[id]  * energy_factor ;
        #endif
      }
    }
  }
}

void Change_Cosmological_Frame_Sytem( Grid3D &G, bool forward ){

  if (forward) chprintf( " Converting to Cosmological Comoving System\n");
  else chprintf( " Converting to Cosmological Physical System\n");
  Change_DM_Frame_System( G, forward );
  Change_GAS_Frame_System( G, forward );
}


#endif

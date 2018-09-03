#ifdef COSMOLOGY

#include"cosmology.h"
#include "../io.h"
#include "io_cosmology.h"

Cosmology::Cosmology( void ){}

void Initialize_Cosmology( Cosmology &Cosmo, struct parameters P, Particles_3D &Particles ){

  chprintf( "Cosmological Simulation\n");


  Cosmo.H0 = 67.74;                //[km/s / Mpc]
  Cosmo.cosmo_h = Cosmo.H0/100;
  Cosmo.H0 /= 1000;               //[km/s / kpc]
  Cosmo.Omega_M = 0.3089;
  Cosmo.Omega_L = 0.6911;
  Cosmo.Omega_K = 0.0;

  Cosmo.current_z = Particles.current_z;
  Cosmo.current_a = Particles.current_a;

  Cosmo.cosmo_G = 4.29448e-06 ;     //kpc km^2 s^-2 Msun^-1

  Cosmo.max_delta_a = 1e-4;

  Cosmo.r_0_dm   = P.xlen/P.nx;
  Cosmo.t_0_dm   = 1. / Cosmo.H0;
  Cosmo.v_0_dm   = Cosmo.r_0_dm / Cosmo.t_0_dm / Cosmo.cosmo_h;
  Cosmo.rho_0_dm = 3*Cosmo.H0*Cosmo.H0 / ( 8*M_PI*Cosmo.cosmo_G ) * Cosmo.Omega_M /Cosmo.cosmo_h/Cosmo.cosmo_h;

  chprintf( " H0: %f\n", Cosmo.H0 * 1000 );
  chprintf( " Omega_M: %f\n", Cosmo.Omega_M );
  chprintf( " Omega_L: %f\n", Cosmo.Omega_L );
  chprintf( " Omega_K: %f\n", Cosmo.Omega_K );
  chprintf( " Current_a: %f\n", Cosmo.current_a );
  chprintf( " Current_z: %f\n", Cosmo.current_z );
  chprintf( " r0: %f\n", Cosmo.r_0_dm );

  Load_Scale_Outputs( P, Cosmo );

  // chprintf( "\n");
}


#endif

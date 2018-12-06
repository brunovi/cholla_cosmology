#ifdef COSMOLOGY

#include"cosmology.h"
#include "../io.h"
#include "io_cosmology.h"
#include "../universal_constants.h"

Cosmology::Cosmology( void ){}

void Initialize_Cosmology( Cosmology &Cosmo, struct parameters P, Particles_3D &Particles, Grav3D &Grav ){

  chprintf( "Cosmological Simulation\n");


  Cosmo.H0 = 67.74;                //[km/s / Mpc]
  Cosmo.cosmo_h = Cosmo.H0/100;
  Cosmo.H0 /= 1000;               //[km/s / kpc]
  Cosmo.Omega_M = 0.3089;
  Cosmo.Omega_L = 0.6911;
  Cosmo.Omega_K = 0.0;

  Cosmo.current_z = Particles.current_z;
  Cosmo.current_a = Particles.current_a;

  Grav.current_a = Cosmo.current_a;

  // Cosmo.cosmo_G = 4.29448e-06 ;     //kpc km^2 s^-2 Msun^-1
  // Cosmo.cosmo_G = 4.3019425e-06;
  Cosmo.cosmo_G = 4.300927161e-06;
  // 4.3019425e-9

  Grav.Gconst = Cosmo.cosmo_G;

  Cosmo.max_delta_a = 0.0001;
  Cosmo.delta_a = Cosmo.max_delta_a;


  Cosmo.time_conversion = KPC;
  Cosmo.t_secs = 0;

  Cosmo.r_0_dm   = P.xlen/P.nx;
  Cosmo.t_0_dm   = 1. / Cosmo.H0;
  Cosmo.v_0_dm   = Cosmo.r_0_dm / Cosmo.t_0_dm / Cosmo.cosmo_h;
  Cosmo.rho_0_dm = 3*Cosmo.H0*Cosmo.H0 / ( 8*M_PI*Cosmo.cosmo_G ) * Cosmo.Omega_M /Cosmo.cosmo_h/Cosmo.cosmo_h;
  // Cosmo.dens_avrg = 0;

  Cosmo.r_0_gas = 1.0;
  Cosmo.rho_0_gas = 3*Cosmo.H0*Cosmo.H0 / ( 8*M_PI*Cosmo.cosmo_G ) * Cosmo.Omega_M /Cosmo.cosmo_h/Cosmo.cosmo_h;
  Cosmo.t_0_gas = 1/Cosmo.H0*Cosmo.cosmo_h;
  Cosmo.v_0_gas = Cosmo.r_0_gas / Cosmo.t_0_gas;
  Cosmo.phi_0_gas = Cosmo.v_0_gas * Cosmo.v_0_gas;
  Cosmo.p_0_gas = Cosmo.rho_0_gas * Cosmo.v_0_gas * Cosmo.v_0_gas;
  Cosmo.e_0_gas = Cosmo.v_0_gas * Cosmo.v_0_gas;

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

// Real Cosmology::Get_da_from_dt( Real dt ){
//   Real a2, a3, d_a;
//   a2 = current_a * current_a;
//   a3 = a2 * current_a;
//   d_a = H0 * current_a * sqrt( Omega_M/a3 + Omega_K/a2 + Omega_L);
//   return d_a * dt;
// }

Real Cosmology::Get_da_from_dt( Real dt ){
  Real a2 = current_a * current_a;
  Real a_dot = sqrt( Omega_M/current_a + a2*Omega_L) * H0 ;
  return a_dot * dt;
}


Real Cosmology::Get_dt_from_da( Real da ){
  Real a2 = current_a * current_a;
  Real a_dot = sqrt( Omega_M/current_a + a2*Omega_L) * H0 ;
  return da / a_dot;
}

Real Scale_Function( Real a, Real Omega_M, Real Omega_L, Real Omega_K ){
  Real a3 = a * a * a;
  Real factor = ( Omega_M + a*Omega_K + a3*Omega_L ) / a;
  return 1./sqrt(factor);
}

Real Cosmology::Get_Cosmology_dt( Real da ){
  Real a2, a3, d_a, dt;
  a2 = current_a * current_a;
  a3 = a2 * current_a;
  d_a = current_a * sqrt( Omega_M/a3 + Omega_K/a2 + Omega_L);
  dt = da / d_a / current_a / current_a;
  return dt;
}

Real Cosmology::Get_da_courant( Real dt ){
  Real a2, a3, da_dt, da;
  a2 = current_a * current_a;
  a3 = a2 * current_a;
  da_dt = current_a * sqrt( Omega_M/a3 + Omega_K/a2 + Omega_L);
  da = dt * da_dt * current_a * current_a;
  return da;
}

#endif

#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include<stdio.h>
#include <cmath>
#include "io.h"
#include"global.h"
#include"universal_constants.h"
#include"grid3D.h"
#include"part3D.h"
#include "io_cosmo.h"

#include <vector>

// #define H0 68.7
// #define Omega_M 0.3
// #define Omega_L 0.7
// #define Omega_K 0.0


extern Real H0;
extern Real Omega_M;
extern Real Omega_L;
extern Real Omega_K;

extern Real cosmo_G;
extern Real cosmo_h;
extern Real current_z;
extern Real current_a;
extern Real max_delta_a;
extern Real delta_a;
extern Real r_0;
extern Real t_0;
extern Real v_0;
extern Real rho_0;
extern Real phi_0;

extern Real phi_0_gas;
extern Real rho_0_gas;
extern Real r_0_gas;
extern Real v_0_gas;
extern Real t_0_gas;
extern Real p_0_gas;
extern Real e_0_gas;

extern Real delta_a_gas;
extern Real current_a_gas;
extern Real a0;
// extern Real L_mpc;
// extern Real particle_DM_mass;
// extern Real Lambda;

extern ch_vector_t scale_outputs;
extern int n_outputs;

void normalize_dm_units( Grid3D &G, Part3D &Parts, bool inverse );
void normalize_gas_units( Grid3D &G, bool inverse );
Real get_dt_from_da( Real da);
Real get_da_from_dt( Real dt);
Real scale_funct( Real a );
void density_proper_to_comoving( Part3D &Parts );
void density_comoving_to_proper( Part3D &Parts );

void Gas_density_comoving_to_proper( Grid3D &G );
void Gas_density_proper_to_comoving( Grid3D &G );
// void potential_proper_to_cosmological( Part3D &Parts );

void potential_comoving_to_proper( Grid3D &G,  Part3D &Part );

void Initialize_Cosmology( Grid3D &G, struct parameters P );
// void Make_Dimentionless( Part3D &Parts );

#endif //COSMOLOGY_H

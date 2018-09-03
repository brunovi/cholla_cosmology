#ifdef GRAVITY

#ifndef GRAV_FUNC_H
#define GRAV_FUNC_H

#include"../grid3D.h"
#include"../global.h"
// #include"mpi_pfft.h"
#include "poisson_solver_3D.h"

#ifdef POTENTIAL_CUFFT
#include "potential_CUFFT_3D.h"
#endif

#ifdef POTENTIAL_FFTW
#include "potential_FFTW_3D.h"
#endif

#ifdef PARTICLES
#include "../particles/density_CIC.h"
#endif

void Copy_Hydro_Density_to_Gravity( Grid3D &G );

#ifdef POTENTIAL_CUFFT
void Compute_Gravitational_Potential( Grid3D &G,  Potential_CUFFT_3D &p_solver, Real *time_pot, Real *time_set);
#endif

#ifdef POTENTIAL_FFTW
void Compute_Gravitational_Potential( Grid3D &G,  Potential_FFTW_3D &p_solver, Real *time_pot, Real *time_set);
#endif

void Copy_Potential_To_Hydro_Grid( Grid3D &G );

void Copy_Potential_From_Hydro_Grid( Grid3D &G );

void Extrapolate_Grav_Potential( Grid3D &G  );

#endif //SELF_GRAV_FUNC_H
#endif //SELF_GRAVITY

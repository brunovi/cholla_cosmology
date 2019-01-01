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

#ifdef POTENTIAL_PFFT
#include "potential_PFFT_3D.h"
#endif

#ifdef PARTICLES
#include "../particles/density_CIC.h"
#endif

void Copy_Hydro_Density_to_Gravity( Grid3D &G );

#ifdef POTENTIAL_CUFFT
void Compute_Gravitational_Potential( Grid3D &G,  Potential_CUFFT_3D &p_solver, Real *time_pot, Real *time_pDens, Real *time_pDens_trans, struct parameters P);
#endif

#ifdef POTENTIAL_FFTW
void Compute_Gravitational_Potential( Grid3D &G,  Potential_FFTW_3D &p_solver, Real *time_pot, Real *time_pDens, Real *time_pDens_trans, struct parameters P);
#endif

#ifdef POTENTIAL_PFFT
void Compute_Gravitational_Potential( Grid3D &G,  Potential_PFFT_3D &p_solver, struct parameters P);
#endif

Real Get_Density_Average( Grid3D &G );

void Copy_Potential_To_Hydro_Grid( Grid3D &G );

// void Copy_Potential_From_Hydro_Grid( Grid3D &G );

#ifdef GRAVITY_CPU
void Add_Gavity_To_Hydro( Grid3D &G );
#endif


#ifdef GRAVITY_CORRECTOR
void Get_Gavity_Corrector( Grid3D &G, int g_start, int g_end);
void Apply_Gavity_Corrector( Grid3D &G, struct parameters P );
#endif

#ifdef MPI_CHOLLA
void Transfer_Potential_Boundaries_MPI( Grid3D &G, struct parameters P);
#endif

void Extrapolate_Grav_Potential( Grid3D &G  );

void Set_dt( Grid3D &G, bool &output_now, int n_step );

#endif //SELF_GRAV_FUNC_H
#endif //SELF_GRAVITY

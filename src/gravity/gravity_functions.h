#ifdef GRAVITY

#ifndef GRAV_FUNC_H
#define GRAV_FUNC_H

#include"../grid3D.h"
#include"../global.h"
// #include"mpi_pfft.h"
#include "poisson_solver_3D.h"
#include "potential_CUFFT_3D.h"

void Copy_Hydro_Density_to_Gravity( Grid3D &G );

void Compute_Gravitational_Potential( Grid3D &G,  Potential_CUFFT_3D &p_solver);

#endif //SELF_GRAV_FUNC_H
#endif //SELF_GRAVITY

/*! \file CTU_3D_cuda.h
 *  \brief Declarations for the cuda version of the 3D CTU algorithm. */

#ifdef CUDA

#ifndef CTU_3D_CUDA_H
#define CTU_3D_CUDA_H

#include"global.h"

Real CTU_Algorithm_3D_CUDA(Real *host_conserved0, Real *host_conserved1, int nx, int ny, int nz, int x_off, int y_off, int z_off, int n_ghost, Real dx, Real dy, Real dz, Real xbound, Real ybound, Real zbound, Real dt, int n_fields, Real dens_floor, Real temp_floor);


#endif //CTU_3D_CUDA_H
#endif //CUDA

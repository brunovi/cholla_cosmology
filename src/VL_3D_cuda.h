/*! \file VL_3D_cuda.h
 *  \brief Declarations for the cuda version of the 3D VL algorithm. */

#ifdef CUDA

#ifndef VL_3D_CUDA_H
#define VL_3D_CUDA_H

#include"global.h"

void Free_GPU_Memory_VL_3D( void );

Real VL_Algorithm_3D_CUDA(Real *host_conserved0, Real *host_conserved1, int nx, int ny, int nz, int x_off, int y_off, int z_off, int n_ghost, Real dx, Real dy, Real dz, Real xbound, Real ybound, Real zbound, Real dt, int n_fields, Real dens_floor, Real temp_floor );



#endif //VL_3D_CUDA_H
#endif //CUDA

#ifdef PARTICLES
#ifdef PARTICLES_CUDA

#include"density_CIC_CUDA.h"
#include"../io.h"
#include<cuda.h>
#include"../global_cuda.h"


void Get_Density_CIC_CUDA( Particles_3D &Parts ){

  part_int_t n_local = Parts.n_local;

  Real *mass_d;
  Real *pos_x_d;
  Real *pos_y_d;
  Real *pos_z_d;
  Real *density_d;

  // allocate memory on the GPU
  CudaSafeCall( cudaMalloc((void**)&mass_d,  n_local*sizeof(Real)) );
  CudaSafeCall( cudaMalloc((void**)&pos_x_d, n_local*sizeof(Real)) );
  CudaSafeCall( cudaMalloc((void**)&pos_y_d, n_local*sizeof(Real)) );
  CudaSafeCall( cudaMalloc((void**)&pos_z_d, n_local*sizeof(Real)) );
  CudaSafeCall( cudaMalloc((void**)&density_d, Parts.G.n_cells*sizeof(Real)) );
  //
  //
  //
  // free the GPU memory
  cudaFree(mass_d);
  cudaFree(pos_x_d);
  cudaFree(pos_y_d);
  cudaFree(pos_z_d);
  cudaFree(density_d);






}


#endif
#endif

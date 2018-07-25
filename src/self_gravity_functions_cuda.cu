#ifdef CUDA
#ifdef SELF_GRAVITY

#include<math.h>
#include<cuda.h>
#include"global_cuda.h"

//void CopyField_Host_To_Device( Real *field_h, Real *field_d, int n_cells ){
void CopyField_Host_To_Device( Real *field_h, Real *field_d, long n_cells ){
  CudaSafeCall( cudaMemcpy( field_d, field_h , n_cells*sizeof(Real), cudaMemcpyHostToDevice) );
}

#endif //SELF_GRAVITY
#endif

#ifdef GRAVITY
#ifdef POTENTIAL_CUFFT

#include "potential_CUFFT_3D.h"


Potential_CUFFT_3D::Potential_CUFFT_3D( void ){}

void Potential_CUFFT_3D::Initialize( Grav3D Grav){

  Lbox_x = Grav.Lbox_x;
  Lbox_y = Grav.Lbox_y;
  Lbox_z = Grav.Lbox_z;

  nx_total = Grav.nx_total;
  ny_total = Grav.ny_total;
  nz_total = Grav.nz_total;

  nx_local = Grav.nx_local;
  ny_local = Grav.ny_local;
  nz_local = Grav.nz_local;

  dx = Grav.dx;
  dy = Grav.dy;
  dz = Grav.dz;

  n_cells_local = nx_local*ny_local*nz_local;
  n_cells_total = nx_total*ny_total*nz_total;
  chprintf( " Using Poisson Solver: CUFFT\n");
  chprintf( "  CUFFT: L[ %f %f %f ] N[ %d %d %d ] dx[ %f %f %f ]\n", Lbox_x, Lbox_y, Lbox_z, nx_local, ny_local, nz_local, dx, dy, dz );

  AllocateMemory_CPU();

  chprintf( "  CUFFT: Creating FFT plan...\n");
  cufftPlan3d( &plan_cufft_fwd,  nz_local, ny_local,  nx_local, CUFFT_D2Z);
  cufftPlan3d( &plan_cufft_bwd,  nz_local, ny_local,  nx_local, CUFFT_Z2D);


}

void Potential_CUFFT_3D::AllocateMemory_CPU( void ){
  F.output_h = (Real *) malloc(n_cells_local*sizeof(Real));

}

void Potential_CUFFT_3D::AllocateMemory_GPU( void ){

  cudaMalloc( (void**)&F.transform_d, n_cells_local*sizeof(Complex_cufft));
  cudaMalloc( (void**)&F.input_d, n_cells_local*sizeof(Real_cufft));
  cudaMalloc( (void**)&F.output_d, n_cells_local*sizeof(Real_cufft));

}

void Potential_CUFFT_3D::FreeMemory_GPU( void ){
  cudaFree( F.input_d );
  cudaFree( F.output_d );
  cudaFree( F.transform_d );
}

void Potential_CUFFT_3D::Copy_Input( Grav3D &Grav ){
  cudaMemcpy( F.input_d, Grav.F.density_h, n_cells_local*sizeof(Real_cufft), cudaMemcpyHostToDevice );
}

void Potential_CUFFT_3D::Copy_Output( Grav3D &Grav ){

  cudaMemcpy( F.output_h, F.output_d, n_cells_local*sizeof(Real_cufft), cudaMemcpyDeviceToHost );
  // cudaMemcpy( F.output_h, F.input_d, n_cells_local*sizeof(Real_cufft), cudaMemcpyDeviceToHost );

  int id, id_pot;
  int i, k, j;
  for (k=0; k<nz_local; k++) {
    for (j=0; j<ny_local; j++) {
      for (i=0; i<nx_local; i++) {
        id = i + j*nx_local + k*nx_local*ny_local;
        id_pot = (i+N_GHOST_POTENTIAL) + (j+N_GHOST_POTENTIAL)*(nx_local+2*N_GHOST_POTENTIAL) + (k+N_GHOST_POTENTIAL)*(nx_local+2*N_GHOST_POTENTIAL)*(ny_local+2*N_GHOST_POTENTIAL);
        Grav.F.potential_h[id_pot] = F.output_h[id] / n_cells_local;
        // chprintf( "%f\n", Grav.F.potential_h[id]);
      }
    }
  }
}

void Potential_CUFFT_3D::Get_Potential( Grav3D &Grav ){

  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start);

  AllocateMemory_GPU();
  Copy_Input( Grav );

  cufftExecD2Z( plan_cufft_fwd, F.input_d, F.transform_d );
  cufftExecZ2D( plan_cufft_bwd, F.transform_d, F.output_d );
  Copy_Output( Grav );

  FreeMemory_GPU();

  cudaEventRecord(stop);
  cudaEventSynchronize(stop);
  float milliseconds = 0;
  cudaEventElapsedTime(&milliseconds, start, stop);
  chprintf( " CUFFT: Potential Time = %f   msecs\n", milliseconds);
}




#endif //POTENTIAL_CUFFT
#endif //GRAVITY

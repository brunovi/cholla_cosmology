#ifdef GRAVITY
#ifdef POTENTIAL_PFFT

#include "potential_PFFT_3D.h"
#include<iostream>

Potential_PFFT_3D::Potential_PFFT_3D( void ){}

void Potential_PFFT_3D::Initialize( Grav3D Grav){

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

  index_0 = -1;

  n_cells_local = nx_local*ny_local*nz_local;
  n_cells_total = nx_total*ny_total*nz_total;
  chprintf( " Using Poisson Solver: PFFT\n");
  chprintf( "  Index 0: %d\n", index_0);

  nproc_pfft = nproc;

  chprintf(" Initializing PFFT:  %d processes\n", nproc_pfft );
  pfft_init();

  chprintf( "  PFFT: L[ %f %f %f ] N_local[ %d %d %d ] dx[ %f %f %f ]\n", Lbox_x, Lbox_y, Lbox_z, nx_local, ny_local, nz_local, dx, dy, dz );

  nprocs_grid_pfft[0] = nproc_z;
  nprocs_grid_pfft[1] = nproc_y;
  nprocs_grid_pfft[2] = nproc_x;

  n_pfft[0] = Grav.nz_total;
  n_pfft[1] = Grav.ny_total;
  n_pfft[2] = Grav.nx_total;

  // /* Create 3-dimensional process grid of size np_pfft[0] x np_pfft[1] x np_pfft[2] */
  pfft_create_procmesh(3, MPI_COMM_WORLD, nprocs_grid_pfft, &comm_pfft);
  MPI_Comm_rank( comm_pfft, &procID_pfft );
  MPI_Cart_coords( comm_pfft, procID_pfft, 3, pcoords_pfft);

  /* Get parameters of data distribution */
  alloc_local = pfft_local_size_dft_r2c_3d(
      n_pfft, comm_pfft, PFFT_TRANSPOSED_OUT ,
      local_ni_pfft, local_i_start_pfft, local_ntrans_pfft, local_trans_start_pfft);

  chprintf("  PFFT process: nz:%d ny:%d nx:%d \n", nprocs_grid_pfft[0], nprocs_grid_pfft[1], nprocs_grid_pfft[2]);
  chprintf("  PFFT cells:   nz:%d ny:%d nx:%d \n", n_pfft[0], n_pfft[1], n_pfft[2]);
  chprintf("  PFFT local:   nz:%d ny:%d nx:%d \n", local_ni_pfft[0], local_ni_pfft[1], local_ni_pfft[2] );

  AllocateMemory_CPU();

  chprintf( "  PFFT: Creating FFT plan...\n");
  // /* Plan parallel  FFTs */
  plan_fwd = pfft_plan_dft_r2c_3d(n_pfft, F.input, F.transform, comm_pfft,
      PFFT_FORWARD, PFFT_TRANSPOSED_OUT | PFFT_MEASURE);

  plan_bwd = pfft_plan_dft_c2r_3d(n_pfft, F.transform, F.output, comm_pfft,
      PFFT_BACKWARD, PFFT_TRANSPOSED_IN | PFFT_MEASURE);
  chprintf( "  PFFT: Computing K for Gravity Green Funtion\n");
  Get_K_for_Green_function();

}

void Potential_PFFT_3D::AllocateMemory_CPU( void ){
  /* Allocate memory */
  F.input  = pfft_alloc_real(2*alloc_local);
  F.transform = pfft_alloc_complex(alloc_local);
  F.output = pfft_alloc_real(2*alloc_local);
  F.G = (Real *) malloc(n_cells_local*sizeof(Real));

}

void Potential_PFFT_3D::Reset( void ){
  free( F.input );
  free( F.output );
  free( F.transform );
  free( F.G );
}

void Potential_PFFT_3D::Copy_Input( Grav3D &Grav ){
  int i, k, j, id;
  for (k=0; k<nz_local; k++) {
    for (j=0; j<ny_local; j++) {
      for (i=0; i<nx_local; i++) {
        id = i + j*nx_local + k*nx_local*ny_local;
        F.input[id] = Grav.F.density_h[id] ;
      }
    }
  }
}

void Potential_PFFT_3D::Copy_Output( Grav3D &Grav ){
  int id, id_pot;
  int i, k, j;
  for (k=0; k<nz_local; k++) {
    for (j=0; j<ny_local; j++) {
      for (i=0; i<nx_local; i++) {
        id = i + j*nx_local + k*nx_local*ny_local;
        id_pot = (i+N_GHOST_POTENTIAL) + (j+N_GHOST_POTENTIAL)*(nx_local+2*N_GHOST_POTENTIAL) + (k+N_GHOST_POTENTIAL)*(nx_local+2*N_GHOST_POTENTIAL)*(ny_local+2*N_GHOST_POTENTIAL);
        Grav.F.potential_h[id_pot] = F.output[id] / n_cells_total;
      }
    }
  }
}

void Potential_PFFT_3D::Get_K_for_Green_function( void){
  int m = 0;
  int k_0, k_1, k_2;
  Real k_x, k_y, k_z, G_x, G_y, G_z, G;

  double ksqrd;
  for(ptrdiff_t k1=local_trans_start_pfft[1]; k1<local_trans_start_pfft[1]+local_ntrans_pfft[1]; k1++){
    for(ptrdiff_t k2=local_trans_start_pfft[2]; k2<local_trans_start_pfft[2]+local_ntrans_pfft[2]; k2++){
      for(ptrdiff_t k0=local_trans_start_pfft[0]; k0<local_trans_start_pfft[0]+local_ntrans_pfft[0]; k0++){
        //  printf("out[%td, %td, %td] = %.2f + I * %.2f\n", k0, k1, k2, transf[m][0],transf[m][1]);
        k_0 = k0;
        k_1 = k1;
        k_2 = k2;
        k_x = 2*M_PI*k_2/nx_total;
        k_z = 2*M_PI*k_0/nz_total;
        k_y = 2*M_PI*k_1/nx_total;
        G_x = sin( k_x/2  );
        G_y = sin( k_y/2 );
        G_z = sin( k_z/2 );
        G = -1 / ( G_x*G_x + G_y*G_y + G_z*G_z ) * dx * dx /4 ;
        if ( k_0==0 && k_1==0 && k_2 == 0 ){
          index_0 = m ;
          G = 1;
          std::cout << " ###########  K=0   index: " << index_0 << std::endl;
        }
        F.G[m] = G;
        m += 1;
      }
    }
  }
}



void Potential_PFFT_3D::Apply_K2_Funtion( void ){
  Real kx, ky, kz, k2;
  int id;
  for (int k=0; k<nz_local; k++){
    for (int j=0; j<ny_local; j++){
      for ( int i=0; i<nx_local; i++){
        id = i + j*nx_local + k*nx_local*ny_local;
        kz = k;
        ky = j;
        kx = i;
        if ( kz > nz_local/2) kz -= nz_local;
        if ( ky > ny_local/2) ky -= ny_local;
        if ( kx > nx_local/2) kx -= nx_local;
        k2  = 4 * M_PI * M_PI * ( kx*kx + ky*ky + kz*kz );
        if ( id == 0 ) k2 = 1;
        F.transform[id][0] *= -1/k2;
        F.transform[id][1] *= -1/k2;
      }
    }
  }
  if (index_0 > -1){
    F.transform[index_0][0] = 0;
    F.transform[index_0][1] = 0;
  }
}

void Potential_PFFT_3D::Apply_G_Funtion( void ){
  Real G_val;
  for ( int i=0; i<n_cells_local; i++ ){
    G_val = F.G[i];
    F.transform[i][0] *= G_val;
    F.transform[i][1] *= G_val;
  }
  if (index_0 > -1){
    F.transform[index_0][0] = 0;
    F.transform[index_0][1] = 0;
  }
}


Real Potential_PFFT_3D::Get_Potential( Grav3D &Grav ){

  double start = get_time();

  Copy_Input( Grav );

  pfft_execute( plan_fwd );
  Apply_G_Funtion();
  pfft_execute( plan_bwd );
  Copy_Output( Grav );

  double stop = get_time();
  double milliseconds = (stop - start) * 1000.0;
  chprintf( " PFFT: Potential Time = %f   msecs\n", milliseconds);
  return milliseconds;
}

#endif //POTENTIAL_CUFFT
#endif //GRAVITY

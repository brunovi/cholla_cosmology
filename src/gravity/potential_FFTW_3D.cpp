#ifdef GRAVITY
#ifdef POTENTIAL_FFTW

#include "potential_FFTW_3D.h"

Potential_FFTW_3D::Potential_FFTW_3D( void ){}

void Potential_FFTW_3D::Initialize( Grav3D Grav){

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
  chprintf( " Using Poisson Solver: FFTW\n");
  chprintf( "  FFTW: L[ %f %f %f ] N[ %d %d %d ] dx[ %f %f %f ]\n", Lbox_x, Lbox_y, Lbox_z, nx_local, ny_local, nz_local, dx, dy, dz );

  AllocateMemory_CPU();

  chprintf( "  FFTW: Creating FFT plan...\n");
  fftw_plan_fwd = fftw_plan_dft_r2c_3d( nz_local, ny_local, nx_local, F.input, F.transform, FFTW_ESTIMATE );
  fftw_plan_bwd = fftw_plan_dft_c2r_3d( nz_local, ny_local, nx_local, F.transform, F.output, FFTW_ESTIMATE );

  chprintf( "  CUFFT: Computing K for Gravity Green Funtion\n");
  Get_K_for_Green_function();

}

void Potential_FFTW_3D::AllocateMemory_CPU( void ){
  F.input = (Real_fftw *) malloc(n_cells_local*sizeof(Real_fftw));
  F.output = (Real_fftw *) malloc(n_cells_local*sizeof(Real_fftw));
  F.transform = (Complex_fftw *) malloc(n_cells_local*sizeof(Complex_fftw));
  F.G = (Real *) malloc(n_cells_local*sizeof(Real));

}

void Potential_FFTW_3D::Reset( void ){
  free( F.input );
  free( F.output );
  free( F.transform );
  free( F.G );
}

void Potential_FFTW_3D::Copy_Input( Grav3D &Grav ){
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

void Potential_FFTW_3D::Copy_Output( Grav3D &Grav ){
  int id, id_pot;
  int i, k, j;
  for (k=0; k<nz_local; k++) {
    for (j=0; j<ny_local; j++) {
      for (i=0; i<nx_local; i++) {
        id = i + j*nx_local + k*nx_local*ny_local;
        id_pot = (i+N_GHOST_POTENTIAL) + (j+N_GHOST_POTENTIAL)*(nx_local+2*N_GHOST_POTENTIAL) + (k+N_GHOST_POTENTIAL)*(nx_local+2*N_GHOST_POTENTIAL)*(ny_local+2*N_GHOST_POTENTIAL);
        Grav.F.potential_h[id_pot] = F.output[id] / n_cells_local;
      }
    }
  }
}

void Potential_FFTW_3D::Get_K_for_Green_function( void){
  Real kx, ky, kz, Gx, Gy, Gz, G;
  int id;
  for (int k=0; k<nz_local; k++){
    kz =  2*M_PI*k/nz_local;
    Gz = sin( kz/2 );
    for (int j=0; j<ny_local; j++){
      ky =  2*M_PI*j/ny_local;
      Gy = sin( ky/2 );
      for ( int i=0; i<nx_local; i++){
        id = i + j*nx_local + k*nx_local*ny_local;
        kx =  2*M_PI*i/nx_local;
        Gx = sin( kx/2 );
        G = -1 / ( Gx*Gx + Gy*Gy + Gz*Gz ) * dx * dx / 4 ;
        if ( id == 0 ) G = 1;
        F.G[id] = G;
      }
    }
  }
}


void Potential_FFTW_3D::Apply_G_Funtion( void ){
  for ( int i=0; i<n_cells_local; i++ ){
    F.transform[i][0] *= F.G[i];
    F.transform[i][1] *= F.G[i];
  }
  F.transform[0][0] = 0;
  F.transform[0][1] = 0;
}


void Potential_FFTW_3D::Get_Potential( Grav3D &Grav ){

  double start = get_time();

  Copy_Input( Grav );

  fftw_execute_dft_r2c( fftw_plan_fwd, F.input, F.transform );
  Apply_G_Funtion();
  fftw_execute_dft_c2r( fftw_plan_bwd, F.transform, F.output );
  Copy_Output( Grav );

  double stop = get_time();
  double milliseconds = (stop - start) * 1000.0;
  chprintf( " CUFFT: Potential Time = %f   msecs\n", milliseconds);
}

#endif //POTENTIAL_CUFFT
#endif //GRAVITY

#ifdef SELF_GRAVITY

#include"self_gravity_functions.h"
#include"cosmology.h"


void Grid3D::CopyDensity_To_Gravity(void){
  //int i, j, k, id, id_grav;
  long i, j, k, id, id_grav;
  // Copy the density array to the memory buffer
  Real factor = 1.0;

  #ifdef COSMOLOGY
  factor = rho_0_gas;
  #endif

  for (k=0; k<Grav.nz_local; k++) {
    for (j=0; j<Grav.ny_local; j++) {
      for (i=0; i<Grav.nx_local; i++) {
        id = (i+H.n_ghost) + (j+H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
        id_grav = (i) + (j)*Grav.nx_local + (k)*Grav.nx_local*Grav.ny_local;
        Grav.F.density_h[id_grav] = C.density[id] *factor;

	//now all report as real
	/*
	if(k==0 && j==0 && i==0)
	{
	  printf("CopyDensity TGrav procID %d id_Grav %ld %ld density %e factor %e rho_0_gas %e\n",procID,id,id_grav,C.density[id],factor,rho_0_gas);
	}
	*/
      }
    }
  }
}


Real Grav3D::Get_Potential_PFFT( void ){
  //Copy density from conserved array to gravity density field
  // MPI_Barrier( comm_pfft );
  time_start_pfft = get_time();
  // // // if ( copy_from_device ) CopyDensity_To_Gravity_device( dev_conserved,  )
  // // // G.CopyDensity_To_Gravity_host();
  // // // MPI_Barrier( comm_pfft );
  CopyDensity_To_FFT( nx_local, ny_local, nz_local, F.density_h );
  // // MPI_Barrier( comm_pfft );
  Exec_PFFT_Forw();
  // // MPI_Barrier( comm_pfft );
  // Divide_by_K2();
  Apply_G_Funtion( Lbox_x, nx_total, ny_total, nz_total );
  // // MPI_Barrier( comm_pfft );
  Exec_PFFT_Back();
  // // MPI_Barrier( comm_pfft );
  // if (copy) Copy_Potential( nx_local, ny_local, nz_local, potential );
  Copy_Potential( nx_local, ny_local, nz_local, F.potential_h );
  // Copy_Extrapolate_Potential( nx_local, ny_local, nz_local, dt, dt_1, potential, potential_1 );
  // // Set_Potential_Boundaries( nx_local, ny_local, nz_local, potential );
  time_stop_pfft = get_time();
  return time_stop_pfft - time_start_pfft;
}


#endif //SELF_GRAVITY

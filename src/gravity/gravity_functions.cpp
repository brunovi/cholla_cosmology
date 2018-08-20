#ifdef GRAVITY

#include"gravity_functions.h"


void Copy_Hydro_Density_to_Gravity( Grid3D &G ){
  int i, j, k, id, id_grav;
  // Copy the density array from hydro conserved to gravity density array
  for (k=0; k<G.Grav.nz_local; k++) {
    for (j=0; j<G.Grav.ny_local; j++) {
      for (i=0; i<G.Grav.nx_local; i++) {
        id = (i+G.H.n_ghost) + (j+G.H.n_ghost)*G.H.nx + (k+G.H.n_ghost)*G.H.nx*G.H.ny;
        id_grav = (i) + (j)*G.Grav.nx_local + (k)*G.Grav.nx_local*G.Grav.ny_local;
        G.Grav.F.density_h[id_grav] = G.C.density[id] ;
      }
    }
  }
}


void Compute_Gravitational_Potential( Grid3D &G){

  Copy_Hydro_Density_to_Gravity( G );
}




#endif

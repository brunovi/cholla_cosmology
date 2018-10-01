#ifdef PARTICLES

#include "density_boundaries.h"
#include "../io.h"

void Transfer_Particles_Density_Boundaries( Particles_3D &Parts ){
  int nGHST, nx_g, ny_g, nz_g, nx, ny, nz;
  nGHST = Parts.G.n_ghost_particles_grid;
  nx_g = Parts.G.nx_local + 2*nGHST;
  ny_g = Parts.G.ny_local + 2*nGHST;
  nz_g = Parts.G.nz_local + 2*nGHST;
  nx = Parts.G.nx_local;
  ny = Parts.G.ny_local;
  nz = Parts.G.nz_local;

  int i, j, k, id_src, id_dst;

  for ( k=0; k<nGHST; k++ ){
    for ( j=0; j<ny_g; j++ ){
      for ( i=0; i<nx_g; i++ ){
        id_src = (i) + (j)*nx_g + (k)*nx_g*ny_g;
        id_dst = (i) + (j)*nx_g + (nz_g - 2*nGHST + k )* nx_g*ny_g;
        Parts.G.density[id_dst] += Parts.G.density[id_src];
      }
    }
  }

  for ( k=0; k<nGHST; k++ ){
    for ( j=0; j<ny_g; j++ ){
      for ( i=0; i<nx_g; i++ ){
        id_src = (i) + (j)*nx_g + (nz_g - nGHST + k)*nx_g*ny_g;
        id_dst = (i) + (j)*nx_g + (k+nGHST)*nx_g*ny_g;
        Parts.G.density[id_dst] += Parts.G.density[id_src];
      }
    }
  }
  for ( k=0; k<nz_g; k++ ){
    for ( j=0; j<nGHST; j++ ){
      for ( i=0; i<nx_g; i++ ){
        id_src = (i) + (j)*nx_g + (k)*nx_g*ny_g;
        id_dst = (i) + (ny_g - 2*nGHST + j)*nx_g + (k)* nx_g*ny_g;
        Parts.G.density[id_dst] += Parts.G.density[id_src];
      }
    }
  }
  for ( k=0; k<nz_g; k++ ){
    for ( j=0; j<nGHST; j++ ){
      for ( i=0; i<nx_g; i++ ){
        id_src = (i) + (ny_g - nGHST + j)*nx_g + (k)*nx_g*ny_g;
        id_dst = (i) + (j+nGHST)*nx_g + (k)*nx_g*ny_g;
        Parts.G.density[id_dst] += Parts.G.density[id_src];
      }
    }
  }
  for ( k=0; k<nz_g; k++ ){
    for ( j=0; j<ny_g; j++ ){
      for ( i=0; i<nGHST; i++ ){
        id_src = (i) + (j)*nx_g + (k)*nx_g*ny_g;
        id_dst = (nx_g - 2*nGHST + i) + (j)*nx_g + (k)* nx_g*ny_g;
        Parts.G.density[id_dst] += Parts.G.density[id_src];
      }
    }
  }
  for ( k=0; k<nz_g; k++ ){
    for ( j=0; j<ny_g; j++ ){
      for ( i=0; i<nGHST; i++ ){
        id_src = (nx_g - nGHST + i) + (j)*nx_g + (k)*nx_g*ny_g;
        id_dst = (i+nGHST) + (j)*nx_g + (k)*nx_g*ny_g;
        Parts.G.density[id_dst] += Parts.G.density[id_src];
      }
    }
  }
}

#ifdef MPI_CHOLLA

void Transfer_Particles_Density_Boundaries_MPI( Grid3D &G, struct parameters P){


  G.Particles.TRANSFER_DENSITY_BOUNDARIES = true;
  // chprintf( " Transfering Particles Density Boundaries\n");
  G.Set_Boundary_Conditions(P);
  G.Particles.TRANSFER_DENSITY_BOUNDARIES = false;

}
#endif
#endif

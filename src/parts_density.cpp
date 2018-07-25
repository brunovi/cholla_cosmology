#ifdef PARTICLES

#include"parts_density.h"



void Part3D::Clear_Density( void ){
  int nGHST = G.grid_ghost;
  int n_cells = (G.nx+2*nGHST) * (G.ny+2*nGHST) * (G.nz+2*nGHST);

  for( int i=0; i<n_cells; i++ ) G.density_h[i] = 0;
}



void Part3D::Get_Indexes_CIC( Real pos_x, Real pos_y, Real pos_z, int &indx_x, int &indx_y, int &indx_z ){

  // Real Lx = G.xMax - G.xMin;
  // Real Ly = G.yMax - G.yMin;
  // Real Lz = G.zMax - G.zMin;

  indx_x = (int) floor( ( pos_x - G.xMin - 0.5*G.dx )/G.dx );
  indx_y = (int) floor( ( pos_y - G.yMin - 0.5*G.dy )/G.dy );
  indx_z = (int) floor( ( pos_z - G.zMin - 0.5*G.dz )/G.dz );

  // indx_x += G.grid_ghost;
  // indx_y += G.grid_ghost;
  // indx_z += G.grid_ghost;

}

void Part3D::Get_Density_CIC( partID_t p_indx_start, partID_t p_indx_end ){
  int nGHST = G.grid_ghost;
  int nx_g = G.nx+ 2*nGHST;
  int ny_g = G.ny+ 2*nGHST;
  int nz_g = G.nz+ 2*nGHST;
  int n_cells = nx_g*ny_g*nz_g;
  partID_t pIndx;
  int indx_x, indx_y, indx_z, indx;
  Real pMass, x_pos, y_pos, z_pos;
  Real cell_center_x, cell_center_y, cell_center_z;
  Real delta_x, delta_y, delta_z;
  Real dV_inv = 1./(G.dx*G.dy*G.dz);

  // #pragma omp parallel for
  for ( pIndx=0; pIndx < n_local; pIndx++ ){
    pMass = mass[pIndx] * dV_inv;
    x_pos = pos_x[pIndx];
    y_pos = pos_y[pIndx];
    z_pos = pos_z[pIndx];
    Get_Indexes_CIC( x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );
    // if ( indx_x < -1 ) std::cout << "Negative xIndx" << std::endl;
    // if ( indx_y < -1 ) std::cout << "Negative yIndx" << std::endl;
    // if ( indx_z < -1 ) std::cout << "Negative zIndx" << std::endl;
    // if ( indx_x > nx_g-2  ) std::cout << "Excess xIndx" << std::endl;
    // if ( indx_y > ny_g-2  ) std::cout << "Excess yIndx" << std::endl;
    // if ( indx_y > nz_g-2  ) std::cout << "Excess zIndx" << std::endl;
    cell_center_x = G.xMin + indx_x*G.dx + 0.5*G.dx;
    cell_center_y = G.yMin + indx_y*G.dy + 0.5*G.dy;
    cell_center_z = G.zMin + indx_z*G.dz + 0.5*G.dz;
    delta_x = 1 - ( x_pos - cell_center_x )/G.dx;
    delta_y = 1 - ( y_pos - cell_center_y )/G.dy;
    delta_z = 1 - ( z_pos - cell_center_z )/G.dz;
    indx_x += nGHST;
    indx_y += nGHST;
    indx_z += nGHST;
    indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass  * delta_x * delta_y * delta_z;
    }

    indx = (indx_x+1) + indx_y*nx_g + indx_z*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass  * (1-delta_x) * delta_y * delta_z;
    }

    indx = indx_x + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass  * delta_x * (1-delta_y) * delta_z;
    }

    indx = indx_x + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass  * delta_x * delta_y * (1-delta_z);
    }

    indx = (indx_x+1) + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass  * (1-delta_x) * (1-delta_y) * delta_z;
    }

    indx = (indx_x+1) + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass  * (1-delta_x) * delta_y * (1-delta_z);
    }

    indx = indx_x + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass  * delta_x * (1-delta_y) * (1-delta_z);
    }

    indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells ){
      // #pragma omp atomic
      G.density_h[indx] += pMass * (1-delta_x) * (1-delta_y) * (1-delta_z);
    }
  }
}


void Part3D::Add_Density_to_Grid( Grid3D &Grid, bool clearPrevious ){

  int nGHST = G.grid_ghost;
  int nx_g = G.nx+ 2*nGHST;
  int ny_g = G.ny+ 2*nGHST;
  int nz_g = G.nz+ 2*nGHST;

  //int i, j, k, id_grid, id_part;
  long i, j, k, id_grid, id_part;
  // Copy the density array to the memory buffer
  for (k=0; k<Grid.Grav.nz_local; k++) {
    for (j=0; j<Grid.Grav.ny_local; j++) {
      for (i=0; i<Grid.Grav.nx_local; i++) {
        id_part = (i+nGHST) + (j+nGHST)*nx_g + (k+nGHST)*nx_g*ny_g;
        id_grid = (i) + (j)*Grid.Grav.nx_local + (k)*Grid.Grav.nx_local*Grid.Grav.ny_local;
        if ( clearPrevious ) Grid.Grav.F.density_h[id_grid] = 0;
        Grid.Grav.F.density_h[id_grid] += G.density_h[id_part];

      }
    }
  }
}





#endif //PARTICLES

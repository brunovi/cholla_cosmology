#ifdef PARTICLES

#include"parts_gravity.h"

#ifdef COSMOLOGY
#include"cosmology.h"
#endif


void Part3D::Copy_potential_from_Grid( Grid3D &Grid ){

  //int i, j, k, indx_grid, indx_pot;
  //int nGHST_g = 1;
  //int nGHST_pot = 2;
  //int nx_g, ny_g, nz_g, nx_p, ny_p, nz_p;
  long i, j, k, indx_grid, indx_pot;
  long nGHST_g = 1;
  long nGHST_pot = 2;
  long nx_g, ny_g, nz_g, nx_p, ny_p, nz_p;
  nx_g = G.nx+2*nGHST_g;
  ny_g = G.ny+2*nGHST_g;
  nz_g = G.nz+2*nGHST_g;
  nx_p = G.nx+2*nGHST_pot;
  ny_p = G.ny+2*nGHST_pot;
  nz_p = G.nz+2*nGHST_pot;

  for (k=0; k<G.nz; k++) {
    for (j=0; j<G.ny; j++) {
      for (i=0; i<G.nx; i++) {
        indx_grid =(i+nGHST_g) + (j+nGHST_g)*nx_g + (k+nGHST_g)*nx_g*ny_g;
        indx_pot = (i+nGHST_pot) + (j+nGHST_pot)*nx_p + (k+nGHST_pot)*nx_p*ny_p;
        G.potential[indx_pot] = Grid.Grav.F.potential_h[indx_grid];

/*
	//fixed Already nan here
	if(k==0 && j==0 && i==0)
	{
	  printf("Copy potential from G procID %d indx_grid %ld indx_pot %ld pot %e\n",procID,indx_grid,indx_pot,Grid.Grav.F.potential_h[indx_grid]);
	  fflush(stdout);
	}
*/

        // #ifdef COSMOLOGY
        // Grid.Grav.F.potential_h[indx_grid] *= current_a_gas*current_a_gas / phi_0_gas;
        // #endif
      }
    }
  }
}

void Part3D::Copy_potential_to_Grid( Grid3D &Grid ){

  //int i, j, k, indx_grid, indx_pot;
  //int nGHST_g = 1;
  //int nGHST_pot = 2;
  //int nx_g, ny_g, nz_g, nx_p, ny_p, nz_p;
  long i, j, k, indx_grid, indx_pot;
  long nGHST_g = 1;
  long nGHST_pot = 2;
  long nx_g, ny_g, nz_g, nx_p, ny_p, nz_p;
  nx_g = G.nx+2*nGHST_g;
  ny_g = G.ny+2*nGHST_g;
  nz_g = G.nz+2*nGHST_g;
  nx_p = G.nx+2*nGHST_pot;
  ny_p = G.ny+2*nGHST_pot;
  nz_p = G.nz+2*nGHST_pot;

  for (k=0; k<nz_g; k++) {
    for (j=0; j<ny_g; j++) {
      for (i=0; i<nx_g; i++) {
        indx_grid =(i) + (j)*nx_g + (k)*nx_g*ny_g;
        indx_pot = (i+1) + (j+1)*nx_p + (k+1)*nx_p*ny_p;
        Grid.Grav.F.potential_h[indx_grid] = G.potential[indx_pot];

        #ifdef COSMOLOGY
        Grid.Grav.F.potential_h[indx_grid] *= current_a_gas*current_a_gas / phi_0_gas;
        #endif
/*
	if(k==0 && j==0 && i==0)
	{
	  printf("Copy potential to G procID %d indx_grid %ld indx_pot %ld pot %e current_a_gas %e phi_0_gas %e\n",procID,indx_grid,indx_pot,Grid.Grav.F.potential_h[indx_grid], current_a_gas, phi_0_gas);
	  fflush(stdout);
	}
*/
      }
    }
  }
}

void Part3D::Get_Potential_Gradient( bool cosmological, int grid_indx_start, int grid_indx_end ){

  //int i, j, k, indx_grav, indx_pot, indx_l, indx_r;
  long i, j, k, indx_grav, indx_pot, indx_l, indx_r;
  Real pot_r, pot_l, pot_c, grav;
  //int nGHST = G.grid_ghost;
  //int nGHST_pot = G.grid_ghost_pot;
  //int nx_g, ny_g, nz_g, nx_p, ny_p, nz_p;
  long nGHST = G.grid_ghost;
  long nGHST_pot = G.grid_ghost_pot;
  long nx_g, ny_g, nz_g, nx_p, ny_p, nz_p;
  nx_g = G.nx+2*nGHST;
  ny_g = G.ny+2*nGHST;
  nz_g = G.nz+2*nGHST;
  nx_p = G.nx+2*nGHST_pot;
  ny_p = G.ny+2*nGHST_pot;
  nz_p = G.nz+2*nGHST_pot;
  Real *pot_array;

  //G.potential is nans here
  if (cosmological) pot_array = G.potential;
  else pot_array = G.potential;

  for (k=grid_indx_start; k<grid_indx_end; k++) {
  // for (k=0; k<nz_g; k++) {
    for (j=0; j<ny_g; j++) {
      for (i=0; i<nx_g; i++) {
        indx_grav = (i) + (j)*nx_g + (k)*nx_g*ny_g;
        indx_pot = (i+1) + (j+1)*nx_p + (k+1)*nx_p*ny_p;
        pot_c = pot_array[ indx_pot ];
        // X component
        indx_r = (i+2) + (j+1)*nx_p + (k+1)*nx_p*ny_p;
        indx_l = (i)   + (j+1)*nx_p + (k+1)*nx_p*ny_p;
        pot_r = pot_array[ indx_r ];
        pot_l = pot_array[ indx_l ];
        grav = -0.5*(pot_r - pot_l);
        if (!cosmological) grav /= G.dx;
        #ifdef COSMO_UNITS
        grav /= G.dx;
        #endif
        G.gravity_x[indx_grav] = grav;

/*
        if((k==grid_indx_start) && (i==0) &&( j==0 ) )
	{
	  printf("Pot Grad procID %ld n?_g %d %d %d indx_pot %ld %ld %ld pot %e %e %e\n",procID,nx_g,ny_g,nz_g,indx_pot,indx_r,indx_l,pot_c,pot_r,pot_l);
	  fflush(stdout);
	}
*/

        // Y component
        indx_r = (i+1) + (j+2)*nx_p + (k+1)*nx_p*ny_p;
        indx_l = (i+1) + (j)*nx_p + (k+1)*nx_p*ny_p;
        pot_r = pot_array[ indx_r ];
        pot_l = pot_array[ indx_l ];
        grav = -0.5*(pot_r - pot_l);
        if (!cosmological) grav /= G.dy;
        #ifdef COSMO_UNITS
        grav /= G.dy;
        #endif
        G.gravity_y[indx_grav] = grav;

        // Z component
        indx_r = (i+1) + (j+1)*nx_p + (k+2)*nx_p*ny_p;
        indx_l = (i+1) + (j+1)*nx_p + (k)*nx_p*ny_p;
        pot_r = pot_array[ indx_r ];
        pot_l = pot_array[ indx_l ];
        grav = -0.5*(pot_r - pot_l);
        if (!cosmological) grav /= G.dz;
        #ifdef COSMO_UNITS
        grav /= G.dz;
        #endif
        G.gravity_z[indx_grav] = grav;
      }
    }
  }

}



void Part3D::Get_Gravity_CIC(  partID_t p_indx_start, partID_t p_indx_end ){
  int nGHST = G.grid_ghost;
  int nx_g = G.nx+ 2*nGHST;
  int ny_g = G.ny+ 2*nGHST;
  int nz_g = G.nz+ 2*nGHST;
  //int n_cells = nx_g*ny_g*nz_g;
  long n_cells = nx_g*ny_g*nz_g;

  // #pragma omp parallel
  // {

  partID_t pIndx;
  //int indx_x, indx_y, indx_z, indx;
  int indx_x, indx_y, indx_z;
  long indx;
  Real x_pos, y_pos, z_pos;
  Real cell_center_x, cell_center_y, cell_center_z;
  Real delta_x, delta_y, delta_z;
  // Real dV_inv = 1./(G.dx*G.dy*G.dz);
  Real g_x_bl, g_x_br, g_x_bu, g_x_bru, g_x_tl, g_x_tr, g_x_tu, g_x_tru;
  Real g_y_bl, g_y_br, g_y_bu, g_y_bru, g_y_tl, g_y_tr, g_y_tu, g_y_tru;
  Real g_z_bl, g_z_br, g_z_bu, g_z_bru, g_z_tl, g_z_tr, g_z_tu, g_z_tru;
  Real g_x, g_y, g_z;

  // #pragma omp parallel for
  for ( pIndx=p_indx_start; pIndx < p_indx_end; pIndx++ ){
    // pMass = mass[pIndx];
    x_pos = pos_x[pIndx];
    y_pos = pos_y[pIndx];
    z_pos = pos_z[pIndx];
    Get_Indexes_CIC( x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );
    if ( indx_x < -1 ) std::cout << "Negative xIndx" << std::endl;
    if ( indx_y < -1 ) std::cout << "Negative yIndx" << std::endl;
    if ( indx_z < -1 ) std::cout << "Negative zIndx" << std::endl;
    if ( indx_x > nx_g-1  ) std::cout << "Excess xIndx" << std::endl;
    if ( indx_y > ny_g-1  ) std::cout << "Excess yIndx" << std::endl;
    if ( indx_y > nz_g-1  ) std::cout << "Excess zIndx" << std::endl;
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
    // if ( indx >= 0 && indx < n_cells )
    g_x_bl = G.gravity_x[indx];
    g_y_bl = G.gravity_y[indx];
    g_z_bl = G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y)*nx_g + (indx_z)*nx_g*ny_g;
    // if ( indx >= 0 && indx < n_cells )
    g_x_br = G.gravity_x[indx];
    g_y_br = G.gravity_y[indx];
    g_z_br = G.gravity_z[indx];

    indx = (indx_x) + (indx_y+1)*nx_g + (indx_z)*nx_g*ny_g;
    // if ( indx >= 0 && indx < n_cells )
    g_x_bu = G.gravity_x[indx];
    g_y_bu = G.gravity_y[indx];
    g_z_bu = G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z)*nx_g*ny_g;
    // if ( indx >= 0 && indx < n_cells )
    g_x_bru = G.gravity_x[indx];
    g_y_bru = G.gravity_y[indx];
    g_z_bru = G.gravity_z[indx];

    indx = (indx_x) + (indx_y)*nx_g + (indx_z+1)*nx_g*ny_g;
    // if ( indx >= 0 && indx < n_cells )
    g_x_tl = G.gravity_x[indx];
    g_y_tl = G.gravity_y[indx];
    g_z_tl = G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y)*nx_g + (indx_z+1)*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells )
    g_x_tr = G.gravity_x[indx];
    g_y_tr = G.gravity_y[indx];
    g_z_tr = G.gravity_z[indx];

    indx = (indx_x) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells )
    g_x_tu = G.gravity_x[indx];
    g_y_tu = G.gravity_y[indx];
    g_z_tu = G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    if ( indx >= 0 && indx < n_cells )
    g_x_tru = G.gravity_x[indx];
    g_y_tru = G.gravity_y[indx];
    g_z_tru = G.gravity_z[indx];

    g_x = g_x_bl*(delta_x)*(delta_y)*(delta_z)   + g_x_br*(1-delta_x)*(delta_y)*(delta_z) +
          g_x_bu*(delta_x)*(1-delta_y)*(delta_z) + g_x_bru*(1-delta_x)*(1-delta_y)*(delta_z) +
          g_x_tl*(delta_x)*(delta_y)*(1-delta_z) + g_x_tr*(1-delta_x)*(delta_y)*(1-delta_z) +
          g_x_tu*(delta_x)*(1-delta_y)*(1-delta_z) + g_x_tru*(1-delta_x)*(1-delta_y)*(1-delta_z);

    g_y = g_y_bl*(delta_x)*(delta_y)*(delta_z)   + g_y_br*(1-delta_x)*(delta_y)*(delta_z) +
          g_y_bu*(delta_x)*(1-delta_y)*(delta_z) + g_y_bru*(1-delta_x)*(1-delta_y)*(delta_z) +
          g_y_tl*(delta_x)*(delta_y)*(1-delta_z) + g_y_tr*(1-delta_x)*(delta_y)*(1-delta_z) +
          g_y_tu*(delta_x)*(1-delta_y)*(1-delta_z) + g_y_tru*(1-delta_x)*(1-delta_y)*(1-delta_z);

    g_z = g_z_bl*(delta_x)*(delta_y)*(delta_z)   + g_z_br*(1-delta_x)*(delta_y)*(delta_z) +
          g_z_bu*(delta_x)*(1-delta_y)*(delta_z) + g_z_bru*(1-delta_x)*(1-delta_y)*(delta_z) +
          g_z_tl*(delta_x)*(delta_y)*(1-delta_z) + g_z_tr*(1-delta_x)*(delta_y)*(1-delta_z) +
          g_z_tu*(delta_x)*(1-delta_y)*(1-delta_z) + g_z_tru*(1-delta_x)*(1-delta_y)*(1-delta_z);

    grav_x[pIndx] = g_x;
    grav_y[pIndx] = g_y;
    grav_z[pIndx] = g_z;

/*
    //no more Gravity_x are nans here
    if(pIndx==p_indx_start)
    {
      printf("Get_Grav_CIC procID %ld G.gravity_? %e %e %e grav_? %e %e %e\n",procID,g_x_tru,g_y_tru,g_z_tru,grav_x[pIndx],grav_y[pIndx],grav_z[pIndx]);
      fflush(stdout);
    }
*/
  }
}








#endif //PARTICLES

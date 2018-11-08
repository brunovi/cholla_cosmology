#ifdef PARTICLES

#include"gravity_CIC.h"
#include"../io.h"


void Get_Gravity_Field( Grid3D &G, int g_start, int g_end ){
  // int nGHST_pot = N_GHOST_POTENTIAL;
  int nx_grav, ny_grav, nz_grav, nGHST_grav;
  nGHST_grav = G.Particles.G.n_ghost_particles_grid;
  nx_grav = G.Particles.G.nx_local + 2*nGHST_grav;
  ny_grav = G.Particles.G.ny_local + 2*nGHST_grav;
  nz_grav = G.Particles.G.nz_local + 2*nGHST_grav;

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.Grav.nx_local + 2*nGHST_grid;
  ny_grid = G.Grav.ny_local + 2*nGHST_grid;
  nz_grid = G.Grav.nz_local + 2*nGHST_grid;

  int nGHST = nGHST_grid - nGHST_grav;

  Real dx, dy, dz;
  dx = G.Particles.G.dx;
  dy = G.Particles.G.dy;
  dz = G.Particles.G.dz;


  Real phi_l, phi_r;
  int k, j, i, id_l, id_r, id;
  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i-1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r = (i+1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l = G.C.Grav_potential[id_l];
        phi_r = G.C.Grav_potential[id_r];
        G.Particles.G.gravity_x[id] = -0.5 * ( phi_r - phi_l ) / dx;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i + nGHST) + (j-1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r = (i + nGHST) + (j+1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l = G.C.Grav_potential[id_l];
        phi_r = G.C.Grav_potential[id_r];
        G.Particles.G.gravity_y[id] = -0.5 * ( phi_r - phi_l ) / dy;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l = (i + nGHST) + (j + nGHST)*nx_grid + (k-1 + nGHST)*ny_grid*nz_grid;
        id_r = (i + nGHST) + (j + nGHST)*nx_grid + (k+1 + nGHST)*ny_grid*nz_grid;
        phi_l = G.C.Grav_potential[id_l];
        phi_r = G.C.Grav_potential[id_r];
        G.Particles.G.gravity_z[id] = -0.5 * ( phi_r - phi_l ) / dz;
      }
    }
  }

}



void Get_Gravity_Field_order4( Grid3D &G, int g_start, int g_end ){
  // int nGHST_pot = N_GHOST_POTENTIAL;
  int nx_grav, ny_grav, nz_grav, nGHST_grav;
  nGHST_grav = G.Particles.G.n_ghost_particles_grid;
  nx_grav = G.Particles.G.nx_local + 2*nGHST_grav;
  ny_grav = G.Particles.G.ny_local + 2*nGHST_grav;
  nz_grav = G.Particles.G.nz_local + 2*nGHST_grav;

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.Grav.nx_local + 2*nGHST_grid;
  ny_grid = G.Grav.ny_local + 2*nGHST_grid;
  nz_grid = G.Grav.nz_local + 2*nGHST_grid;

  int nGHST = nGHST_grid - nGHST_grav;

  Real dx, dy, dz;
  dx = G.Particles.G.dx;
  dy = G.Particles.G.dy;
  dz = G.Particles.G.dz;


  Real phi_l_1, phi_r_1, phi_l_2, phi_r_2;
  int k, j, i, id_l_1, id_r_1,  id_l_2, id_r_2,  id;
  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l_1 = (i-1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_l_2 = (i-2 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r_1 = (i+1 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r_2 = (i+2 + nGHST) + (j + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l_1 = G.C.Grav_potential[id_l_1];
        phi_r_1 = G.C.Grav_potential[id_r_1];
        phi_l_2 = G.C.Grav_potential[id_l_2];
        phi_r_2 = G.C.Grav_potential[id_r_2];
        // G.Particles.G.gravity_x[id] = -2/3 * ( phi_r_1 - phi_l_1 ) / dx + 1/12 * ( phi_r_2 - phi_l_2 ) / dx  ;
        G.Particles.G.gravity_x[id] = -(  8*phi_r_1  - 8*phi_l_1 - phi_r_2 + phi_l_2 ) / dx / 12  ;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l_1 = (i + nGHST) + (j-1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_l_2 = (i + nGHST) + (j-2 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r_1 = (i + nGHST) + (j+1 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        id_r_2 = (i + nGHST) + (j+2 + nGHST)*nx_grid + (k + nGHST)*ny_grid*nz_grid;
        phi_l_1 = G.C.Grav_potential[id_l_1];
        phi_r_1 = G.C.Grav_potential[id_r_1];
        phi_l_2 = G.C.Grav_potential[id_l_2];
        phi_r_2 = G.C.Grav_potential[id_r_2];
        // G.Particles.G.gravity_y[id] = -2/3 * ( phi_r_1 - phi_l_1 ) / dy + 1/12 * ( phi_r_2 - phi_l_2 ) / dy  ;
        G.Particles.G.gravity_y[id] = -(  8*phi_r_1  - 8*phi_l_1 - phi_r_2 + phi_l_2 ) / dy / 12 ;
      }
    }
  }

  for ( k=g_start; k<g_end; k++ ){
    for ( j=0; j<ny_grav; j++ ){
      for ( i=0; i<nx_grav; i++ ){
        id   = (i) + (j)*nx_grav + (k)*ny_grav*nz_grav;
        id_l_1 = (i + nGHST) + (j + nGHST)*nx_grid + (k-1 + nGHST)*ny_grid*nz_grid;
        id_l_2 = (i + nGHST) + (j + nGHST)*nx_grid + (k-2 + nGHST)*ny_grid*nz_grid;
        id_r_1 = (i + nGHST) + (j + nGHST)*nx_grid + (k+1 + nGHST)*ny_grid*nz_grid;
        id_r_2 = (i + nGHST) + (j + nGHST)*nx_grid + (k+2 + nGHST)*ny_grid*nz_grid;
        phi_l_1 = G.C.Grav_potential[id_l_1];
        phi_r_1 = G.C.Grav_potential[id_r_1];
        phi_l_2 = G.C.Grav_potential[id_l_2];
        phi_r_2 = G.C.Grav_potential[id_r_2];
        // G.Particles.G.gravity_z[id] = -2/3 * ( phi_r_1 - phi_l_1 ) / dz + 1/12 * ( phi_r_2 - phi_l_2 ) / dz  ;
        G.Particles.G.gravity_z[id] = -(  8*phi_r_1  - 8*phi_l_1 - phi_r_2 + phi_l_2 ) / dz / 12 ;
      }
    }
  }

}


void Get_Gravity_CIC( Particles_3D &Particles, part_int_t p_start, part_int_t p_end ){

  int nx_g, ny_g, nz_g, nGHST;
  nGHST = Particles.G.n_ghost_particles_grid;
  nx_g = Particles.G.nx_local + 2*nGHST;
  ny_g = Particles.G.ny_local + 2*nGHST;
  nz_g = Particles.G.nz_local + 2*nGHST;

  Real xMin, yMin, zMin, dx, dy, dz;
  xMin = Particles.G.xMin;
  yMin = Particles.G.yMin;
  zMin = Particles.G.zMin;
  dx = Particles.G.dx;
  dy = Particles.G.dy;
  dz = Particles.G.dz;

  part_int_t pIndx;
  int indx_x, indx_y, indx_z, indx;
  Real x_pos, y_pos, z_pos;
  Real cell_center_x, cell_center_y, cell_center_z;
  Real delta_x, delta_y, delta_z;
  Real g_x_bl, g_x_br, g_x_bu, g_x_bru, g_x_tl, g_x_tr, g_x_tu, g_x_tru;
  Real g_y_bl, g_y_br, g_y_bu, g_y_bru, g_y_tl, g_y_tr, g_y_tu, g_y_tru;
  Real g_z_bl, g_z_br, g_z_bu, g_z_bru, g_z_tl, g_z_tr, g_z_tu, g_z_tru;
  Real g_x, g_y, g_z;
  bool ignore;
  for ( pIndx=p_start; pIndx < p_end; pIndx++ ){
    ignore = false;
    // pMass = Particles.mass[pIndx] * dV_inv;
    x_pos = Particles.pos_x[pIndx];
    y_pos = Particles.pos_y[pIndx];
    z_pos = Particles.pos_z[pIndx];
    Get_Indexes_CIC( xMin, yMin, zMin, dx, dy, dz, x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );
    if ( indx_x < -1 ) ignore = true;
    if ( indx_y < -1 ) ignore = true;
    if ( indx_z < -1 ) ignore = true;
    if ( indx_x > nx_g-2  ) ignore = true;
    if ( indx_y > ny_g-2  ) ignore = true;
    if ( indx_y > nz_g-2  ) ignore = true;
    if ( ignore ){
      #ifdef PARTICLE_IDS
      std::cout << "ERROR GRAVITY_CIC Index    pID: " << Particles.partIDs[pIndx] << std::endl;
      #else
      std::cout << "ERROR GRAVITY_CIC Index " << std::endl;
      #endif
      std::cout << "Negative xIndx: " << x_pos << "  " << indx_x << std::endl;
      std::cout << "Negative zIndx: " << z_pos << "  " << indx_z << std::endl;
      std::cout << "Negative yIndx: " << y_pos << "  " << indx_y << std::endl;
      std::cout << "Excess xIndx: " << x_pos << "  " << indx_x << std::endl;
      std::cout << "Excess yIndx: " << y_pos << "  " << indx_y << std::endl;
      std::cout << "Excess zIndx: " << z_pos << "  " << indx_z << std::endl;
      std::cout << std::endl;
      continue;
    }
    cell_center_x = xMin + indx_x*dx + 0.5*dx;
    cell_center_y = yMin + indx_y*dy + 0.5*dy;
    cell_center_z = zMin + indx_z*dz + 0.5*dz;
    delta_x = 1 - ( x_pos - cell_center_x ) / dx;
    delta_y = 1 - ( y_pos - cell_center_y ) / dy;
    delta_z = 1 - ( z_pos - cell_center_z ) / dz;
    indx_x += nGHST;
    indx_y += nGHST;
    indx_z += nGHST;

    indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
    g_x_bl = Particles.G.gravity_x[indx];
    g_y_bl = Particles.G.gravity_y[indx];
    g_z_bl = Particles.G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y)*nx_g + (indx_z)*nx_g*ny_g;
    g_x_br = Particles.G.gravity_x[indx];
    g_y_br = Particles.G.gravity_y[indx];
    g_z_br = Particles.G.gravity_z[indx];

    indx = (indx_x) + (indx_y+1)*nx_g + (indx_z)*nx_g*ny_g;
    g_x_bu = Particles.G.gravity_x[indx];
    g_y_bu = Particles.G.gravity_y[indx];
    g_z_bu = Particles.G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z)*nx_g*ny_g;
    g_x_bru = Particles.G.gravity_x[indx];
    g_y_bru = Particles.G.gravity_y[indx];
    g_z_bru = Particles.G.gravity_z[indx];

    indx = (indx_x) + (indx_y)*nx_g + (indx_z+1)*nx_g*ny_g;
    g_x_tl = Particles.G.gravity_x[indx];
    g_y_tl = Particles.G.gravity_y[indx];
    g_z_tl = Particles.G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y)*nx_g + (indx_z+1)*nx_g*ny_g;
    g_x_tr = Particles.G.gravity_x[indx];
    g_y_tr = Particles.G.gravity_y[indx];
    g_z_tr = Particles.G.gravity_z[indx];

    indx = (indx_x) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    g_x_tu = Particles.G.gravity_x[indx];
    g_y_tu = Particles.G.gravity_y[indx];
    g_z_tu = Particles.G.gravity_z[indx];

    indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
    g_x_tru = Particles.G.gravity_x[indx];
    g_y_tru = Particles.G.gravity_y[indx];
    g_z_tru = Particles.G.gravity_z[indx];

    g_x = g_x_bl*(delta_x)*(delta_y)*(delta_z)     + g_x_br*(1-delta_x)*(delta_y)*(delta_z) +
          g_x_bu*(delta_x)*(1-delta_y)*(delta_z  ) + g_x_bru*(1-delta_x)*(1-delta_y)*(delta_z) +
          g_x_tl*(delta_x)*(delta_y)*(1-delta_z)   + g_x_tr*(1-delta_x)*(delta_y)*(1-delta_z) +
          g_x_tu*(delta_x)*(1-delta_y)*(1-delta_z) + g_x_tru*(1-delta_x)*(1-delta_y)*(1-delta_z);

    g_y = g_y_bl*(delta_x)*(delta_y)*(delta_z)     + g_y_br*(1-delta_x)*(delta_y)*(delta_z) +
          g_y_bu*(delta_x)*(1-delta_y)*(delta_z)   + g_y_bru*(1-delta_x)*(1-delta_y)*(delta_z) +
          g_y_tl*(delta_x)*(delta_y)*(1-delta_z)   + g_y_tr*(1-delta_x)*(delta_y)*(1-delta_z) +
          g_y_tu*(delta_x)*(1-delta_y)*(1-delta_z) + g_y_tru*(1-delta_x)*(1-delta_y)*(1-delta_z);

    g_z = g_z_bl*(delta_x)*(delta_y)*(delta_z)     + g_z_br*(1-delta_x)*(delta_y)*(delta_z) +
          g_z_bu*(delta_x)*(1-delta_y)*(delta_z)   + g_z_bru*(1-delta_x)*(1-delta_y)*(delta_z) +
          g_z_tl*(delta_x)*(delta_y)*(1-delta_z)   + g_z_tr*(1-delta_x)*(delta_y)*(1-delta_z) +
          g_z_tu*(delta_x)*(1-delta_y)*(1-delta_z) + g_z_tru*(1-delta_x)*(1-delta_y)*(1-delta_z);

    Particles.grav_x[pIndx] = g_x;
    Particles.grav_y[pIndx] = g_y;
    Particles.grav_z[pIndx] = g_z;
  }
}

#endif

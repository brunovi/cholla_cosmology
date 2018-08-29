#ifdef PARTICLES

#include"gravity_CIC.h"
#include"../io.h"


void Get_Gavity_Field( Grid3D &G ){
  // int nGHST_pot = N_GHOST_POTENTIAL;
  int nx_g, ny_g, nz_g, nGHST;
  nGHST = G.Particles.G.n_ghost_particles_grid;
  nx_g = G.Particles.G.nx_local + 2*nGHST;
  ny_g = G.Particles.G.ny_local + 2*nGHST;
  nz_g = G.Particles.G.nz_local + 2*nGHST;

  Real dx, dy, dz;
  dx = G.Particles.G.dx;
  dy = G.Particles.G.dy;
  dz = G.Particles.G.dz;


  Real phi_l, phi_r;
  int k, j, i, id_l, id_r, id;
  for ( k=1; k<nz_g+1; k++ ){
    for ( j=1; j<ny_g+1; j++ ){
      for ( i=1; i<nx_g+1; i++ ){
        id   = (i) + (j)*nx_g + (k)*ny_g*nz_g;
        id_l = (i-1) + (j)*nx_g + (k)*ny_g*nz_g;
        id_r = (i+1) + (j)*nx_g + (k)*ny_g*nz_g;
        phi_l = G.C.Grav_potential[id_l];
        phi_r = G.C.Grav_potential[id_r];
        G.Particles.G.gravity_x[id] = -0.5 * ( phi_r - phi_l ) / dx;
      }
    }
  }

  for ( k=1; k<nz_g+1; k++ ){
    for ( j=1; j<ny_g+1; j++ ){
      for ( i=1; i<nx_g+1; i++ ){
        id   = (i) + (j)*nx_g + (k)*ny_g*nz_g;
        id_l = (i) + (j-1)*nx_g + (k)*ny_g*nz_g;
        id_r = (i) + (j+1)*nx_g + (k)*ny_g*nz_g;
        phi_l = G.C.Grav_potential[id_l];
        phi_r = G.C.Grav_potential[id_r];
        G.Particles.G.gravity_y[id] = -0.5 * ( phi_r - phi_l ) / dy;
      }
    }
  }

  for ( k=1; k<nz_g+1; k++ ){
    for ( j=1; j<ny_g+1; j++ ){
      for ( i=1; i<nx_g+1; i++ ){
        id   = (i) + (j)*nx_g + (k)*ny_g*nz_g;
        id_l = (i) + (j)*nx_g + (k-1)*ny_g*nz_g;
        id_r = (i) + (j)*nx_g + (k+1)*ny_g*nz_g;
        phi_l = G.C.Grav_potential[id_l];
        phi_r = G.C.Grav_potential[id_r];
        G.Particles.G.gravity_z[id] = -0.5 * ( phi_r - phi_l ) / dz;
      }
    }
  }
}

void Get_Gravity_CIC( Grid3D &G ){

  // Real xMin, yMin, zMin, dx, dy, dz;
  // xMin = G.Particles.G.xMin;
  // yMin = G.Particles.G.yMin;
  // zMin = G.Particles.G.zMin;
  // dx = G.Particles.G.dx;
  // dy = G.Particles.G.dy;
  // dz = G.Particles.G.dz;
  //
  // part_int_t pIndx;
  // int indx_x, indx_y, indx_z, indx;
  // Real pMass, x_pos, y_pos, z_pos;
  //
  // Real cell_center_x, cell_center_y, cell_center_z;
  // Real delta_x, delta_y, delta_z;
  // Real dV_inv = 1./(G.Particles.G.dx*G.Particles.G.dy*G.Particles.G.dz);
  // bool ignore;
  // for ( pIndx=0; pIndx < G.Particles.n_local; pIndx++ ){
  //   ignore = false;
  //   pMass = G.Particles.mass[pIndx] * dV_inv;
  //   x_pos = G.Particles.pos_x[pIndx];
  //   y_pos = G.Particles.pos_y[pIndx];
  //   z_pos = G.Particles.pos_z[pIndx];
  //   Get_Indexes_CIC( xMin, yMin, zMin, dx, dy, dz, x_pos, y_pos, z_pos, indx_x, indx_y, indx_z );
  //   if ( indx_x < -1 ) ignore = true;
  //   if ( indx_y < -1 ) ignore = true;
  //   if ( indx_z < -1 ) ignore = true;
  //   if ( indx_x > nx_g-2  ) ignore = true;
  //   if ( indx_y > ny_g-2  ) ignore = true;
  //   if ( indx_y > nz_g-2  ) ignore = true;
  //   if ( ignore ){
  //     std::cout << "ERROR CIC Index    pID: " << G.Particles.partIDs[pIndx] << std::endl;
  //     std::cout << "Negative xIndx: " << x_pos << "  " << indx_x << std::endl;
  //     std::cout << "Negative zIndx: " << z_pos << "  " << indx_z << std::endl;
  //     std::cout << "Negative yIndx: " << y_pos << "  " << indx_y << std::endl;
  //     std::cout << "Excess xIndx: " << x_pos << "  " << indx_x << std::endl;
  //     std::cout << "Excess yIndx: " << y_pos << "  " << indx_y << std::endl;
  //     std::cout << "Excess zIndx: " << z_pos << "  " << indx_z << std::endl;
  //     std::cout << std::endl;
  //     continue;
  //   }
  //   cell_center_x = xMin + indx_x*dx + 0.5*dx;
  //   cell_center_y = yMin + indx_y*dy + 0.5*dy;
  //   cell_center_z = zMin + indx_z*dz + 0.5*dz;
  //   delta_x = 1 - ( x_pos - cell_center_x ) / dx;
  //   delta_y = 1 - ( y_pos - cell_center_y ) / dy;
  //   delta_z = 1 - ( z_pos - cell_center_z ) / dz;
  //   indx_x += nGHST;
  //   indx_y += nGHST;
  //   indx_z += nGHST;
  //
  //   indx = indx_x + indx_y*nx_g + indx_z*nx_g*ny_g;
  //
  //
  //   indx = (indx_x+1) + indx_y*nx_g + indx_z*nx_g*ny_g;
  //
  //
  //   indx = indx_x + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
  //
  //
  //   indx = indx_x + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
  //
  //
  //   indx = (indx_x+1) + (indx_y+1)*nx_g + indx_z*nx_g*ny_g;
  //
  //
  //   indx = (indx_x+1) + indx_y*nx_g + (indx_z+1)*nx_g*ny_g;
  //
  //
  //   indx = indx_x + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
  //
  //
  //   indx = (indx_x+1) + (indx_y+1)*nx_g + (indx_z+1)*nx_g*ny_g;
  //
  // }
}

#endif

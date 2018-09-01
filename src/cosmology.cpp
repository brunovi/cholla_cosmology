// #ifdef COSMOLOGY
//
// #include"cosmology.h"
//
// Real H0;
// Real Omega_M;
// Real Omega_L;
// Real Omega_K;
//
//
// Real cosmo_G;
// Real cosmo_h;
// Real current_z;
// Real current_a;
// Real delta_a;
// Real max_delta_a;
// Real r_0;
// Real t_0;
// Real v_0;
// Real rho_0;
// Real phi_0;
//
// Real a0;
// Real delta_a_gas;
// Real current_a_gas;
// Real phi_0_gas;
// Real rho_0_gas;
// Real r_0_gas;
// Real v_0_gas;
// Real t_0_gas;
// Real p_0_gas;
// Real e_0_gas;
// // Real L_mpc;
// // Real particle_DM_mass;
// Real Lambda;
//
// ch_vector_t scale_outputs;
// int n_outputs;
//
//
// // Real get_da_from_dt( Real dt){
// //   Real d_a, a3, a0_3 ;
// //   a3 = current_a*current_a*current_a;
// //   a0_3 = a0 * a0 * a0;
// //   d_a = H0 / cosmo_h  * t_0_gas * sqrt( Omega_M*a3*a0_3  + Omega_L*a3*a3 );
// //   // d_a = current_a * sqrt( Omega_M/a3  + Omega_L );
// //   return d_a * dt;
// // }
// //
// //
// // Real get_dt_from_da( Real da){
// //   Real d_a, a3, a0_3 ;
// //   a3 = current_a_gas*current_a_gas*current_a_gas;
// //   a0_3 = a0 * a0 * a0;
// //   d_a = H0 / cosmo_h  * t_0_gas * sqrt( Omega_M*a3*a0_3  + Omega_L*a3*a3 );
// //   // d_a = current_a * sqrt( Omega_M/a3  + Omega_L );
// //   return da / d_a ;
// // }
//
// Real get_da_from_dt( Real dt){
//   Real d_a, a3 ;
//   a3 = current_a*current_a*current_a;
//   d_a = H0 * current_a * sqrt( Omega_M/a3  + Omega_L );
//   // d_a = current_a * sqrt( Omega_M/a3  + Omega_L );
//   return d_a * dt;
// }
//
// Real get_dt_from_da( Real da){
//   Real d_a, a3 ;
//   a3 = current_a*current_a*current_a;
//   d_a =  current_a * sqrt( Omega_M/a3  + Omega_L );
//   // d_a = current_a * sqrt( Omega_M/a3  + Omega_L );
//   return da / d_a / current_a_gas / current_a_gas ;
// }
//
//
// Real scale_funct( Real a ){
//   Real a3 = a*a*a;
//   Real val = ( Omega_M + a*Omega_K + a3*Omega_L ) / a;
//   return 1.0/sqrt( val );
// }
//
// //
// // void density_proper_to_comoving( Part3D &Parts ){
// //   Real rho_crit = 3*H0*H0/(8*M_PI*cosmo_G);
// //   Real a3 = current_a * current_a * current_a;
// //
// //   int nGHST = Parts.G.grid_ghost;
// //   int n_cells = (Parts.G.nx+2*nGHST) * (Parts.G.ny+2*nGHST) * (Parts.G.nz+2*nGHST);
// //
// //   for( int i=0; i<n_cells; i++ ) Parts.G.density_h[i] *= a3;
// // }
//
//
// void normalize_dm_units( Grid3D &G, Part3D &Parts, bool inverse ){
//   int i, j, k, id;
//   // Real pos_norm = r_0;
//   // Real vel_norm = v_0 /( current_a * sqrt(current_a));  //From GADGET guide to get peculiar velocities
//
//   Real vel_norm = 1/( current_a * sqrt(current_a));  //From GADGET guide to get peculiar velocities
//   chprintf(" Normalize dm units \n");
//   partID_t pIndx;
//
//   if (!inverse){
//     for ( pIndx=0; pIndx<Parts.n_local; pIndx++ ){
//       // Parts.pos_x[pIndx] /= pos_norm;
//       // Parts.pos_y[pIndx] /= pos_norm;
//       // Parts.pos_z[pIndx] /= pos_norm;
//       Parts.vel_x[pIndx] /= vel_norm;
//       Parts.vel_y[pIndx] /= vel_norm;
//       Parts.vel_z[pIndx] /= vel_norm;
//     }
//   }
//
//   if (inverse){
//     for ( pIndx=0; pIndx<Parts.n_local; pIndx++ ){
//       // Parts.pos_x[pIndx] *= pos_norm;
//       // Parts.pos_y[pIndx] *= pos_norm;
//       // Parts.pos_z[pIndx] *= pos_norm;
//       Parts.vel_x[pIndx] *= vel_norm;
//       Parts.vel_y[pIndx] *= vel_norm;
//       Parts.vel_z[pIndx] *= vel_norm;
//     }
//   }
// }
//
// void normalize_gas_units( Grid3D &G, bool inverse ){
//   int i, j, k, id;
//   chprintf(" Normalize gas units %d \n", G.H.nz);
//   // Copy the density array to the memory buffer
//   // Real a3 = current_a*current_a*current_a;
//   // Real vel_factor = 15;
//   Real dens, vel_x, vel_y, vel_z, p, e;
//   if (!inverse){
//     for (k=0; k<G.H.nz; k++) {
//       for (j=0; j<G.H.ny; j++) {
//         for (i=0; i<G.H.nx; i++) {
//           id = i + j*G.H.nx + k*G.H.nx*G.H.ny;
//           // id_grav = (i) + (j)*Grav.nx_local + (k)*Grav.nx_local*Grav.ny_local;
//           // dens = G.C.density[id];
//           // vel_x = G.C.momentum_x[id] / dens *  ( current_a * sqrt( current_a ) ;
//           // vel_y = G.C.momentum_y[id] / dens *  ( current_a * sqrt( current_a ) ;
//           // vel_z = G.C.momentum_z[id] / dens *  ( current_a * sqrt( current_a ) ;
//           // vel2 = vel_x*vel_x + vel_y*vel_y + vel_z*velz;
//           // e =
//           //
//           //
//
//           G.C.density[id] = G.C.density[id] / rho_0_gas ;
//           G.C.momentum_x[id] = G.C.momentum_x[id] / v_0_gas / rho_0_gas  ;
//           G.C.momentum_y[id] = G.C.momentum_y[id] / v_0_gas / rho_0_gas ;
//           G.C.momentum_z[id] = G.C.momentum_z[id] / v_0_gas / rho_0_gas ;
//           G.C.Energy[id] = G.C.Energy[id]  / rho_0_gas / v_0_gas / v_0_gas ;
//         }
//       }
//     }
//   }
//   if ( inverse ){
//     for (k=0; k<G.H.nz; k++) {
//       for (j=0; j<G.H.ny; j++) {
//         for (i=0; i<G.H.nx; i++) {
//           id = i + j*G.H.nx + k*G.H.nx*G.H.ny;
//           // id_grav = (i) + (j)*Grav.nx_local + (k)*Grav.nx_local*Grav.ny_local;
//           G.C.density[id] = G.C.density[id] * rho_0_gas ;
//           G.C.momentum_x[id] = G.C.momentum_x[id] *  v_0_gas * rho_0_gas ;
//           G.C.momentum_y[id] = G.C.momentum_y[id] *  v_0_gas * rho_0_gas ;
//           G.C.momentum_z[id] = G.C.momentum_z[id] *  v_0_gas * rho_0_gas ;
//           G.C.Energy[id] = G.C.Energy[id] * rho_0_gas * v_0_gas * v_0_gas   ;
//         }
//       }
//     }
//   }
// }
//
//
//
// void Initialize_Cosmology( Grid3D &G, struct parameters P ){
//
//   H0 = 67.74;                //[km/s / Mpc]
//   cosmo_h = H0/100;
//   // cosmo_h = 0.7;
//   H0 /= 1000;               //[km/s / kpc]
//   Omega_M = 0.3089;
//   Omega_L = 0.6911;
//   Omega_K = 0.0;
//
//   cosmo_G = 4.29448e-06 ;     //kpc km^2 s^-2 Msun^-1
//
//   current_z = 100;
//   current_a = 1./(current_z+1);
//   // max_delta_a = 5e-3;
//   max_delta_a = 1e-4;
//
//   r_0   = G.Grav.Lbox_x/G.Grav.nx_total;
//   t_0   = 1./H0;
//   v_0   = r_0/t_0/cosmo_h;
//   rho_0 = 3*H0*H0 / ( 8*M_PI*cosmo_G ) * Omega_M /cosmo_h/cosmo_h;
//   // rho_0 =1.0;
//   phi_0 = v_0*v_0;
//
//   a0 = 1;
//   // a0 = pow( Omega_L/Omega_M, 1./3);
//   current_a_gas = current_a * a0;
//
//   r_0_gas = 1.0;
//   rho_0_gas = 3*H0*H0 / ( 8*M_PI*cosmo_G ) * Omega_M /cosmo_h/cosmo_h;
//   // t_0_gas = 2/H0*cosmo_h * sqrt(1./Omega_M/(a0*a0*a0));
//   t_0_gas = 1/H0*cosmo_h  ;
//   v_0_gas = r_0_gas / t_0_gas;
//   phi_0_gas = v_0_gas * v_0_gas;
//   p_0_gas = rho_0_gas * v_0_gas * v_0_gas;
//   e_0_gas = v_0_gas * v_0_gas;
//
//
//   // Lambda = 3 * H0*H0 * Omega_L;
//   // L_mpc = 115.0;
//
//   //NOTE: This assumes same number of particles per side as grid_cells per side
//   // particle_DM_mass =  1.32e5 * ( Omega_M * cosmo_h*cosmo_h ) * std::pow( L_mpc / G.Grav.nx_total*128, 3);
//   // particle_DM_mass = CRITICAL_DENSITY * Omega_M * L_mpc * L_mpc * L_mpc;
//
//
//   int io_result = Load_scale_outputs( P );
//
//   chprintf(" Cosmology Initialized \n");
//   // chprintf(" L_box: %0.2f Mpc\n", L_mpc);
//   chprintf(" H_0: %0.5f  km/s / kpc\n", H0);
//   chprintf(" h:   %0.4f \n", cosmo_h);
//   chprintf(" Omega_M: %0.3f \n", Omega_M);
//   chprintf(" Omega_L: %0.3f \n", Omega_L);
//   chprintf(" Omega_K: %0.3f \n", Omega_K);
//   chprintf(" rho_0: %0.3f  Msun/kpc^3\n", rho_0);
//   chprintf(" v_0: %0.3f  km/s\n", v_0);
//   chprintf(" phi_0: %0.3f  (km/s)^2\n", phi_0);
//   chprintf(" phi_0_gas: %0.5f  (km/s)^2\n", phi_0_gas);
//   chprintf(" v_0_gas: %0.5f  (km/s)^2\n", v_0_gas);
//   chprintf(" rho_0_gas: %0.5f\n", rho_0_gas);
//   chprintf(" cosmo_G : %0.5e\n", cosmo_G);
//
//
// }
//
//
//
//
//
//
// // void Gas_density_comoving_to_proper( Grid3D &G ){
// //   int i, j, k, id;
// //   // chprintf(" Density comoving to proper %d \n", G.H.nz);
// //   // Copy the density array to the memory buffer
// //   Real a3 = current_a*current_a*current_a;
// //   for (k=0; k<G.H.nz; k++) {
// //     for (j=0; j<G.H.ny; j++) {
// //       for (i=0; i<G.H.nx; i++) {
// //         id = i + j*G.H.nx + k*G.H.nx*G.H.ny;
// //         // id_grav = (i) + (j)*Grav.nx_local + (k)*Grav.nx_local*Grav.ny_local;
// //         G.C.density[id] = G.C.density[id] / a3 ;
// //         // G.C.density[id] = 0;
// //       }
// //     }
// //   }
// // }
// //
// // void Gas_density_proper_to_comoving( Grid3D &G ){
// //   int i, j, k, id;
// //   // chprintf(" Density comoving to proper %d \n", G.H.nz);
// //   // Copy the density array to the memory buffer
// //   Real a3 = current_a*current_a*current_a;
// //   for (k=0; k<G.H.nz; k++) {
// //     for (j=0; j<G.H.ny; j++) {
// //       for (i=0; i<G.H.nx; i++) {
// //         id = i + j*G.H.nx + k*G.H.nx*G.H.ny;
// //         // id_grav = (i) + (j)*Grav.nx_local + (k)*Grav.nx_local*Grav.ny_local;
// //         G.C.density[id] = G.C.density[id] * a3 ;
// //         // G.C.density[id] = 0;
// //       }
// //     }
// //   }
// // }
//
//
//
//
//
//
//
//
//
//
// // void density_comoving_to_proper( Part3D &Parts ){
// //   Real rho_crit = 3*H0*H0/(8*M_PI*cosmo_G);
// //   Real a3 = current_a * current_a * current_a;
// //
// //   int nGHST = Parts.G.grid_ghost;
// //   int n_cells = (Parts.G.nx+2*nGHST) * (Parts.G.ny+2*nGHST) * (Parts.G.nz+2*nGHST);
// //
// //   for( int i=0; i<n_cells; i++ ) Parts.G.density_h[i] /= a3;
// // }
//
//
// // void normalize_potential_units( Grid3D &G, Part3D &Part ){
// //
// //   int nx, ny, nz;
// //   nx = Part.G.nx;
// //   ny = Part.G.ny;
// //   nz = Part.G.nz;
// //
// //   int i, j, k, id_pot_grid, id_pot_part;
// //   int nGHST_pot_part, nGHST_pot_grid;
// //
// //
// //
// // }
//
//
//
// void potential_comoving_to_proper( Grid3D &G, Part3D &Part ){
//
//   int nx, ny, nz;
//   nx = Part.G.nx;
//   ny = Part.G.ny;
//   nz = Part.G.nz;
//
//   Real xMin, yMin, zMin, dx, dy, dz;
//   xMin = Part.G.xMin;
//   yMin = Part.G.yMin;
//   zMin = Part.G.zMin;
//   dx = Part.G.dx;
//   dy = Part.G.dy;
//   dz = Part.G.dz;
//
//   int i, j, k, id_pot_grid, id_pot_part;
//   int nGHST_pot_part, nGHST_pot_grid;
//   nGHST_pot_grid = 1;
//   nGHST_pot_part = 2;
//   // int nGHST_pot_part = Part.G.grid_ghost_pot;
//   // nGHST_pot = 1;
//   Real pos_z, pos_y, pos_x, r2;
//   Real pot_comoving, pot_lambda;
//   Real a2 = current_a*current_a;
//   for (k=0; k<nz+2*nGHST_pot_grid; k++) {
//     // pos_z = zMin + (k-nGHST_pot_grid)*dz + 0.5*dz;
//     // if ( pos_z < domainMin_z ) pos_z += Lbox_z;
//     // if ( pos_z > domainMax_z ) pos_z -= Lbox_z;
//     for (j=0; j<ny+2*nGHST_pot_grid; j++) {
//       // pos_y = yMin + (j-nGHST_pot_grid)*dy + 0.5*dy;
//       // if ( pos_y < domainMin_y ) pos_y += Lbox_y;
//       // if ( pos_y > domainMax_y ) pos_y -= Lbox_y;
//       for (i=0; i<nx+2*nGHST_pot_grid; i++) {
//         //       // pId = (i) + (j)*G.nx + (k)*G.nx*G.ny;
//         id_pot_part= (i+1) + (j+1)*(nx+2*nGHST_pot_part) + (k+1)*(nx+2*nGHST_pot_part)*(ny+2*nGHST_pot_part);
//         id_pot_grid= (i) + (j)*(nx+2*nGHST_pot_grid) + (k)*(nx+2*nGHST_pot_grid)*(ny+2*nGHST_pot_grid);
//         // pos_x = xMin + (i-nGHST_pot_grid)*dx + 0.5*dx;
//         // if ( pos_x < domainMin_x ) pos_x += Lbox_x;
//         // if ( pos_x > domainMax_x ) pos_x -= Lbox_x;
//   //
//         // r2 =  pos_x*pos_x + pos_y*pos_y + pos_z*pos_z ;
//         // Part.G.potential[id_pot] = r2;
//         pot_comoving = Part.G.potential[id_pot_part] * H0 * H0 * G.Grav.Lbox_x * G.Grav.Lbox_x ;
//         pot_comoving = Part.G.potential[id_pot_part] * H0 * H0 * G.Grav.Lbox_x * G.Grav.Lbox_x ;
//
//         // pot_comoving = 100.;
//         // pot_lambda = 0.5*H0*H0 * ( a2*Omega_L - 0.5*Omega_M/current_a) / a2 * r2;
//         // pot_lambda = 0.5*H0*H0 * ( a2*Omega_L - 0.5*Omega_M/current_a) / a2 ;
//         // pot_lambda = 0.;
//         // G.Grav.F.potential_h[id_pot_grid] = pot_comoving - pot_lambda;
//         G.Grav.F.potential_h[id_pot_grid] = pot_comoving;
//       }
//     }
//   }
// }
//
// // void Part3D::potential_proper_to_cosmological( void ){
// //   int nx, ny, nz;
// //   nx = G.nx;
// //   ny = G.ny;
// //   nz = G.nz;
// //   int nGHST_pot = G.grid_ghost_pot;
// //   int n_cells_pot = (nx+2*nGHST_pot) * (ny+2*nGHST_pot) * (nz+2*nGHST_pot);
// //   Real xMin, yMin, zMin, dx, dy, dz;
// //   Real domainMin_x, domainMin_y, domainMin_z;
// //   Real domainMax_x, domainMax_y, domainMax_z;
// //   Real Lbox_x, Lbox_y, Lbox_z;
// //   xMin = G.xMin;
// //   yMin = G.yMin;
// //   zMin = G.zMin;
// //   dx = G.dx;
// //   dy = G.dy;
// //   dz = G.dz;
// //   domainMin_x = G.domainMin_x;
// //   domainMin_y = G.domainMin_y;
// //   domainMin_z = G.domainMin_z;
// //   domainMax_x = G.domainMax_x;
// //   domainMax_y = G.domainMax_y;
// //   domainMax_z = G.domainMax_z;
// //   Lbox_x = domainMax_x - domainMin_x;
// //   Lbox_y = domainMax_y - domainMin_y;
// //   Lbox_z = domainMax_z - domainMin_z;
// //
// //   Real pos_x, pos_y, pos_z, r2;
// //   // Real pos_x_loc, pos_y_loc, pos_z_loc;
// //   Real pot_proper, pot_cosmo;
// //   Real a3 = current_a * current_a * current_a;
// //   int i, j, k, id_pot;
// //   for (k=0; k<nz+2*nGHST_pot; k++) {
// //     pos_z = zMin + (k-nGHST_pot)*dz + 0.5*dz;
// //     if ( pos_z < domainMin_z ) pos_z += Lbox_z;
// //     if ( pos_z > domainMax_z ) pos_z -= Lbox_z;
// //     for (j=0; j<ny+2*nGHST_pot; j++) {
// //       pos_y = yMin + (j-nGHST_pot)*dy + 0.5*dy;
// //       if ( pos_y < domainMin_y ) pos_y += Lbox_y;
// //       if ( pos_y > domainMax_y ) pos_y -= Lbox_y;
// //       for (i=0; i<nx+2*nGHST_pot; i++) {
// //         //       // pId = (i) + (j)*G.nx + (k)*G.nx*G.ny;
// //         id_pot = (i) + (j)*(nx+2*nGHST_pot) + (k)*(nx+2*nGHST_pot)*(ny+2*nGHST_pot);
// //         pos_x = xMin + (i-nGHST_pot)*dx + 0.5*dx;
// //         if ( pos_x < domainMin_x ) pos_x += Lbox_x;
// //         if ( pos_x > domainMax_x ) pos_x -= Lbox_x;
// //   //
// //         r2 =  pos_x*pos_x + pos_y*pos_y + pos_z*pos_z ;
// //         pot_proper = G.potential[id_pot];
// //         // pot_cosmo = pot_proper + 0.5*H0*H0 * ( Omega_L - 0.5/a3*Omega_M )*r2;
// //         pot_cosmo = pot_proper;
// //         // G.potential_cosmo[id_pot] = pot_cosmo / ( v_0*v_0 );
// //         G.potential_cosmo[id_pot] = pot_cosmo;
// //       }
// //     }
// //   }
// // }
//
//
//
//
// #endif

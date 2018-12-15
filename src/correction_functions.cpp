#include "correction_functions.h"




#ifdef DE
void Sync_Energies_3D_Host(Grid3D &G ){

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.H.nx;
  ny_grid = G.H.ny;
  nz_grid = G.H.nz;

  int nx_grav, ny_grav, nz_grav, id_grav;
  nx_grav = G.Grav.nx_local;
  ny_grav = G.Grav.ny_local;
  nz_grav = G.Grav.nz_local;

  int nGHST = nGHST_grid ;
  Real d, d_inv, vx, vy, vz, E, Ek, ge1, ge2, Emax;
  int k, j, i, id;

  int n_cells = 0;
  Real Ek_mean = 0;
  for ( k=0; k<nz_grid; k++ ){
    for ( j=0; j<ny_grid; j++ ){
      for ( i=0; i<nx_grid; i++ ){
        if ( (i < nGHST) || (i > (nx_grid - nGHST - 1) ) ) continue;
        if ( (j < nGHST) || (j > (ny_grid - nGHST - 1) ) ) continue;
        if ( (k < nGHST) || (k > (nz_grid - nGHST - 1) ) ) continue;

        id  = (i) + (j)*nx_grid + (k)*ny_grid*nz_grid;
        d = G.C.density[id];
        d_inv = 1/d;
        vx = G.C.momentum_x[id] * d_inv;
        vy = G.C.momentum_y[id] * d_inv;
        vz = G.C.momentum_z[id] * d_inv;
        Ek = 0.5*d*(vx*vx + vy*vy + vz*vz);
        Ek_mean += Ek;
        n_cells += 1;
      }
    }
  }
  Ek_mean /= n_cells;


  int imo, ipo, jmo, jpo, kmo, kpo;
  for ( k=0; k<nz_grid; k++ ){
    for ( j=0; j<ny_grid; j++ ){
      for ( i=0; i<nx_grid; i++ ){
        if ( (i < nGHST) || (i > (nx_grid - nGHST - 1) ) ) continue;
        if ( (j < nGHST) || (j > (ny_grid - nGHST - 1) ) ) continue;
        if ( (k < nGHST) || (k > (nz_grid - nGHST - 1) ) ) continue;

        id_grav = (i-nGHST) + (j-nGHST)*nx_grav + (k-nGHST)*ny_grav*nz_grav;
        id  = (i) + (j)*nx_grid + (k)*ny_grid*nz_grid;
        imo = (i-1) + (j)*nx_grid + (k)*ny_grid*nz_grid;
        ipo = (i+1) + (j)*nx_grid + (k)*ny_grid*nz_grid;
        jmo = (i) + (j-1)*nx_grid + (k)*ny_grid*nz_grid;
        jpo = (i) + (j+1)*nx_grid + (k)*ny_grid*nz_grid;
        kmo = (i) + (j)*nx_grid + (k-1)*ny_grid*nz_grid;
        kpo = (i) + (j)*nx_grid + (k+1)*ny_grid*nz_grid;

        d = G.C.density[id];
        d_inv = 1/d;
        vx = G.C.momentum_x[id] * d_inv;
        vy = G.C.momentum_y[id] * d_inv;
        vz = G.C.momentum_z[id] * d_inv;
        E = G.C.Energy[id];

        if (E < 0.0 || E != E) continue;
        Ek = 0.5*d*(vx*vx + vy*vy + vz*vz);
        ge1 = G.C.GasEnergy[id];
        ge2 = E - 0.5*d*(vx*vx + vy*vy + vz*vz);

        G.Grav.F.extra_field[id_grav] = 0;
        if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.001 ) {
        // if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.001 && Ek/Ek_mean > 0.4 ) {
          G.C.GasEnergy[id] = ge2;
          ge1 = ge2;
          G.Grav.F.extra_field[id_grav] = 1;
          // std::cout << ge2/E << std::endl;
        }

        //find the max nearby total energy
        Emax = E;
        Emax = std::max(G.C.Energy[imo], E);
        Emax = std::max(Emax, G.C.Energy[ipo]);
        Emax = std::max(Emax, G.C.Energy[jmo]);
        Emax = std::max(Emax, G.C.Energy[jpo]);
        Emax = std::max(Emax, G.C.Energy[kmo]);
        Emax = std::max(Emax, G.C.Energy[kpo]);

        if (ge2/Emax > 0.1 && ge2 > 0.0 && Emax > 0.0) {
          G.C.GasEnergy[id] = ge2;
          G.Grav.F.extra_field[id_grav] = 2;
        }
        // sync the total energy with the internal energy
        else {
          if (ge1 > 0.0) G.C.Energy[id] += ge1 - ge2;
          else G.C.GasEnergy[id] = ge2;
        }

        // if ( ge2 < 0 ) G.C.Energy[id] = 0.5*d*(vx*vx + vy*vy + vz*vz) + ge1;

        //InternalEnergy Floor at u=0.2
        Real dens, u, u_physical;
        // Real phi_0_gas = 0.01;                           //Unit Conversion
        dens  =  d;
        u = G.C.GasEnergy[id];
        u_physical = u  * G.Cosmo.v_0_gas * G.Cosmo.v_0_gas / G.Cosmo.current_a / G.Cosmo.current_a;  //convert to physical km^2/s^2

        //Boltazman constant
        Real K_b = 1.38064852e-23; //m2 kg s-2 K-1

        //Mass of proton
        Real M_p = 1.6726219e-27; //kg

        Real gamma = 1.6666667;

        Real temp = u_physical / dens * 1e6 * (gamma - 1) * M_p / K_b ;

        Real temp_0 = 1.0;
        Real u_new, delta_u;
        if ( temp < temp_0 ){
          temp = temp_0;
          u_new = temp * dens * 1e-6 / (gamma - 1) / M_p * K_b ;
          delta_u = u_new - u_physical;
          delta_u = delta_u / G.Cosmo.v_0_gas / G.Cosmo.v_0_gas * G.Cosmo.current_a * G.Cosmo.current_a;
          G.C.GasEnergy[id] += delta_u;
          G.C.Energy[id] += delta_u;
        }
      }
    }
  }
}
#endif

// #ifdef REVERT_STEP
void Get_Delta_Conserved( Grid3D &G ){

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.H.nx;
  ny_grid = G.H.ny;
  nz_grid = G.H.nz;
  int nGHST = nGHST_grid ;

  Real max_delta_d, max_delta_E, max_delta_GE;
  Real max_thr_d, max_thr_E, max_thr_GE;
  max_delta_d = 0.1;
  max_delta_E = 0.9;
  max_delta_GE = 0.9;

  max_thr_d = 0.1;
  max_thr_E = 0.1;
  max_thr_GE = 0.1;


  // Real d, d_inv, vx, vy, vz, E, ge1, ge2, Emax;
  Real v_1, v_0, delta_v, delta_max;
  delta_max = 0.1;
  int k, j, i, id;
  for ( k=0; k<nz_grid; k++ ){
    for ( j=0; j<ny_grid; j++ ){
      for ( i=0; i<nx_grid; i++ ){
        if ( (i < nGHST) || (i > (nx_grid - nGHST - 1) ) ) continue;
        if ( (j < nGHST) || (j > (ny_grid - nGHST - 1 ) ) ) continue;
        if ( (k < nGHST) || (k > (nz_grid - nGHST -1 ) ) ) continue;

        id  = (i) + (j)*nx_grid + (k)*ny_grid*nz_grid;
        //
        //Compare density
        delta_max = 0.1;
        v_0 = G.C.density_0[id];
        v_1 = G.C.density[id];
        delta_v = fabs( (v_1 - v_0 )/ v_0 );
        max_delta_d = std::max( max_delta_d, delta_v );
        if ( delta_v > delta_max ){
          std::cout << "### Delta density: " << delta_v << " delta_max: " << delta_max << std::endl;
        }

        // //Compare momentum_x
        // v_0 = G.C.momentum_x_0[id];
        // v_1 = G.C.momentum_x[id];
        // delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // if ( delta_v > delta_max ){
        //   std::cout << "### Delta momentum_x: " << delta_v << " delta_max: " << delta_max << " v0: " << v_0 << " v1: " << v_1 << std::endl;
        // }
        //
        // //Compare momentum_x
        // v_0 = G.C.momentum_y_0[id];
        // v_1 = G.C.momentum_y[id];
        // delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // if ( delta_v > delta_max ){
        //   std::cout << "### Delta momentum_y: " << delta_v << " delta_max: " << delta_max << " v0: " << v_0 << " v1: " << v_1 << std::endl;
        // }
        //
        // //Compare momentum_x
        // v_0 = G.C.momentum_z_0[id];
        // v_1 = G.C.momentum_z[id];
        // delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // if ( delta_v > delta_max ){
        //   std::cout << "### Delta momentum_z: " << delta_v << " delta_max: " << delta_max << " v0: " << v_0 << " v1: " << v_1 << std::endl;
        // }

        //Compare Energy
        delta_max = 0.8;
        v_0 = G.C.Energy_0[id];
        v_1 = G.C.Energy[id];
        delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // max_delta_E = std::max( max_delta_E, delta_v );
        if ( delta_v > delta_max ){
          std::cout << "### Delta Energy: " << delta_v << " delta_max: " << delta_max << std::endl;
        }
        // if ( delta_v > 0.05 ){
        //   G.C.Energy[id] = v_0 * ( 1 + 1*sgn( v_1 - v_0 )*0.05 );
        // }

        #ifdef DE
        //Compare GasEnergy
        delta_max = 0.1;
        v_0 = G.C.GasEnergy_0[id];
        v_1 = G.C.GasEnergy[id];
        delta_v = fabs( (v_1 - v_0 )/ v_0 );
        max_delta_GE = std::max( max_delta_GE, delta_v );
        if ( delta_v > delta_max ){
          std::cout << "### Delta GasEnergy: " << delta_v << " delta_max: " << delta_max << std::endl;
        }
        // if ( delta_v > 0.01 ){
        //   G.C.GasEnergy[id] = v_0 * ( 1 + 1*sgn( v_1 - v_0 )*0.01 );
        // }

        #endif

        // //Compare Grav_potential
        // v_0 = G.C.Grav_potential_0[id];
        // v_1 = G.C.Grav_potential[id] /  G.Cosmo.phi_0_gas * (G.Cosmo.current_a-G.Cosmo.delta_a) * (G.Cosmo.current_a-G.Cosmo.delta_a) ;
        // delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // if ( delta_v > delta_max ){
        //   std::cout << "### Grav Potential: " << delta_v << " delta_max: " << delta_max << " v0: " << v_0 << " v1: " << v_1 << std::endl;
        // }


      }
    }
  }
  // Real time_factor;
  // Real dt_next_d, dt_next_E, dt_next_GE, dt_next;
  // dt_next_d = 1e10;
  // dt_next_E = 1e10;
  // dt_next_GE = 1e10;
  //
  // // Comare max deltas to max theresholds
  // if ( max_delta_d > max_thr_d ){
  //   std::cout << "### Max Delta Density: " << max_delta_d << "   Max thereshold:" << max_thr_d << std::endl;
  //   G.H.send_revert_signal = true;
  //   time_factor = max_delta_d / max_thr_d;
  //   dt_next_d = G.H.dt / time_factor;
  //   std::cout << dt_next_d << std::endl;
  // }
  // if ( max_delta_E > max_thr_E ){
  //   std::cout << "### Max Delta Energy: " << max_delta_E << "   Max thereshold:" << max_thr_E << std::endl;
  //   G.H.send_revert_signal = true;
  //   time_factor = max_delta_E / max_thr_E;
  //   dt_next_E = G.H.dt / time_factor;
  //   // std::cout << dt_next_E << std::endl;
  // }
  // if ( max_delta_GE > max_thr_GE ){
  //   std::cout << "### Max Delta GasEnergy: " << max_delta_GE << "   Max thereshold:" << max_thr_GE << std::endl;
  //   G.H.send_revert_signal = true;
  //   G.H.revert_dt = G.H.dt / time_factor;
  //   time_factor = max_delta_GE / max_thr_GE;
  //   dt_next_GE = G.H.dt / time_factor;
  //   // std::cout << dt_next_GE << std::endl;
  // }
  //
  // if ( G.H.send_revert_signal ){
  //   dt_next = std::min( dt_next_E, dt_next_d);
  //   dt_next = std::min( dt_next_GE, dt_next);
  //   G.H.revert_dt = dt_next;
  //   std::cout << "### Sending Revert Step Signal      dt_now:" << G.H.dt << "   dt_next:" << dt_next << std::endl;
  // }


}
// #endif

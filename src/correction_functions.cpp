#ifdef GRAVITY_CORRECTOR

#include "correction_functions.h"


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

#endif

#include "correction_functions.h"




void Sync_Energies_3D_Host(Grid3D &G ){

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.H.nx;
  ny_grid = G.H.ny;
  nz_grid = G.H.nz;

  int nGHST = nGHST_grid ;

  Real d, d_inv, vx, vy, vz, E, ge1, ge2, Emax;
  int k, j, i, id;
  int imo, ipo, jmo, jpo, kmo, kpo;
  for ( k=0; k<nz_grid; k++ ){
    for ( j=0; j<ny_grid; j++ ){
      for ( i=0; i<nx_grid; i++ ){
        if ( (i < 1) || (i > (nx_grid - 2) ) ) continue;
        if ( (j < 1) || (j > (ny_grid - 2) ) ) continue;
        if ( (k < 1) || (k > (nz_grid - 2) ) ) continue;

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

        ge1 = G.C.GasEnergy[id];
        ge2 = E - 0.5*d*(vx*vx + vy*vy + vz*vz);
        // std::cout << ge2/E << std::endl;
        if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.001) {
          G.C.GasEnergy[id] = ge2;
          ge1 = ge2;
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
        }
        // sync the total energy with the internal energy
        else {
          if (ge1 > 0.0) G.C.Energy[id] += ge1 - ge2;
          else G.C.GasEnergy[id] = ge2;
        }

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

void Get_Delta_Conserved( Grid3D &G ){

  int nx_grid, ny_grid, nz_grid, nGHST_grid;
  nGHST_grid = G.H.n_ghost;
  nx_grid = G.H.nx;
  ny_grid = G.H.ny;
  nz_grid = G.H.nz;

  int nGHST = nGHST_grid ;

  // Real d, d_inv, vx, vy, vz, E, ge1, ge2, Emax;
  Real v_1, v_0, delta_v, delta_max;
  delta_max = 0.00001;
  int k, j, i, id;
  for ( k=0; k<nz_grid; k++ ){
    for ( j=0; j<ny_grid; j++ ){
      for ( i=0; i<nx_grid; i++ ){
        // if ( (i < nGHST) || (i > (nx_grid - nGHST - 1) ) ) continue;
        // if ( (j < nGHST) || (j > (ny_grid - nGHST - 1 ) ) ) continue;
        // if ( (k < nGHST) || (k > (nz_grid - nGHST -1 ) ) ) continue;

        id  = (i) + (j)*nx_grid + (k)*ny_grid*nz_grid;
        //
        // //Compare density
        // v_0 = G.C.density_0[id];
        // v_1 = G.C.density[id];
        // delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // if ( delta_v > delta_max ){
        //   std::cout << "### Delta density: " << delta_v << " delta_max: " << delta_max << std::endl;
        // }

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

        // //Compare Energy
        // v_0 = G.C.Energy_0[id];
        // v_1 = G.C.Energy[id];
        // delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // if ( delta_v > delta_max ){
        //   std::cout << "### Delta Energy: " << delta_v << " delta_max: " << delta_max << std::endl;
        // }
        //
        // #ifdef DE
        // //Compare GasEnergy
        // v_0 = G.C.GasEnergy_0[id];
        // v_1 = G.C.GasEnergy[id];
        // delta_v = fabs( (v_1 - v_0 )/ v_0 );
        // if ( delta_v > delta_max ){
        //   std::cout << "### Delta GasEnergy: " << delta_v << " delta_max: " << delta_max << std::endl;
        // }
        // #endif

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


}

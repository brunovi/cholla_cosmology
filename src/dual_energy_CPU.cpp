#ifdef GRAVITY_CPU
#ifdef DE

#include "dual_energy_CPU.h"

#ifdef DE_EKINETIC_LIMIT

void Get_Mean_Kinetic_Energy( Grid3D &G ){

  int nx_grid, ny_grid, nz_grid, nGHST;
  nGHST = G.H.n_ghost;
  nx_grid = G.H.nx;
  ny_grid = G.H.ny;
  nz_grid = G.H.nz;

  int n_cells = 0;
  Real Ek_mean = 0;

  // Real d, d_inv, vx, vy, vz, E, Ek, ge1, ge2, Emax;
  // int k, j, i, id;
  // for ( k=0; k<nz_grid; k++ ){
  //   for ( j=0; j<ny_grid; j++ ){
  //     for ( i=0; i<nx_grid; i++ ){
  //       if ( (i < nGHST) || (i > (nx_grid - nGHST - 1) ) ) continue;
  //       if ( (j < nGHST) || (j > (ny_grid - nGHST - 1) ) ) continue;
  //       if ( (k < nGHST) || (k > (nz_grid - nGHST - 1) ) ) continue;
  //
  //       id  = (i) + (j)*nx_grid + (k)*ny_grid*nz_grid;
  //       d = G.C.density[id];
  //       d_inv = 1/d;
  //       vx = G.C.momentum_x[id] * d_inv;
  //       vy = G.C.momentum_y[id] * d_inv;
  //       vz = G.C.momentum_z[id] * d_inv;
  //       Ek = 0.5*d*(vx*vx + vy*vy + vz*vz);
  //       Ek_mean += Ek;
  //       n_cells += 1;
  //     }
  //   }
  // }
  // Ek_mean /= n_cells;

  G.H.Ekin_mean = Ek_mean;

}


#endif

void Sync_Energies_3D_Host_Function(Grid3D &G, int g_start, int g_end ){

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
  int k_g, j_g, i_g;


  int imo, ipo, jmo, jpo, kmo, kpo;
  for ( k_g=g_start; k_g<g_end; k_g++ ){
    for ( j_g=0; j_g<ny_grav; j_g++ ){
      for ( i_g=0; i_g<nx_grav; i_g++ ){

        id_grav = (i_g) + (j_g)*nx_grav + (k_g)*ny_grav*nz_grav;

        i = i_g + nGHST;
        j = j_g + nGHST;
        k = k_g + nGHST;

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

        #ifdef EXTRA_FIELD
        G.Grav.F.extra_field[id_grav] = 0;
        #endif

        #ifdef DE_EKINETIC_LIMIT
        if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.001 && Ek/G.H.Ekin_mean > 0.4 ) {
        #else
        if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.001 ) {
        #endif //DE_EKINETIC_LIMIT

          G.C.GasEnergy[id] = ge2;
          ge1 = ge2;
          #ifdef EXTRA_FIELD
          G.Grav.F.extra_field[id_grav] = 1;
          #endif
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
          #ifdef EXTRA_FIELD
          G.Grav.F.extra_field[id_grav] = 2;
          #endif
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

void Sync_Energies_3D_Host(Grid3D &G ){

  #ifndef PARALLEL_OMP
  Sync_Energies_3D_Host_Function( G, 0, G.Grav.nz_local );
  #else
  #pragma omp parallel num_threads( N_OMP_THREADS )
  {
    int omp_id, n_omp_procs;
    part_int_t p_start, p_end;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local,  &g_start, &g_end  );

    Sync_Energies_3D_Host_Function( G, g_start, g_end );
    // #pragma omp barrier
  }
  #endif

}



#endif //DE
#endif //GRAVITY_CPU

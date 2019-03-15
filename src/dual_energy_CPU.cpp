#ifdef GRAVITY_CPU
#ifdef DE

#include "dual_energy_CPU.h"
#include "universal_constants.h"

#ifdef DE_EKINETIC_LIMIT

Real Get_Mean_Kinetic_Energy_function( Grid3D &G, int g_start, int g_end ){

  int nx_grid, ny_grid, nz_grid, nGHST;
  nGHST = G.H.n_ghost;
  nx_grid = G.H.nx;
  ny_grid = G.H.ny;
  nz_grid = G.H.nz;

  int nx_grav, ny_grav, nz_grav, id_grav;
  nx_grav = G.Grav.nx_local;
  ny_grav = G.Grav.ny_local;
  nz_grav = G.Grav.nz_local;

  Real Ek_mean = 0;
  Real d, d_inv, vx, vy, vz, E, Ek;

  int k, j, i, id;
  int k_g, j_g, i_g;

  for ( k_g=g_start; k_g<g_end; k_g++ ){
    for ( j_g=0; j_g<ny_grav; j_g++ ){
      for ( i_g=0; i_g<nx_grav; i_g++ ){

        i = i_g + nGHST;
        j = j_g + nGHST;
        k = k_g + nGHST;

        id  = (i) + (j)*nx_grid + (k)*ny_grid*nx_grid;

        d = G.C.density[id];
        d_inv = 1/d;
        vx = G.C.momentum_x[id] * d_inv;
        vy = G.C.momentum_y[id] * d_inv;
        vz = G.C.momentum_z[id] * d_inv;
        E = G.C.Energy[id];
        Ek = 0.5*d*(vx*vx + vy*vy + vz*vz);
        Ek_mean += Ek;
      }
    }
  }
  return Ek_mean;
}

void Get_Mean_Kinetic_Energy( Grid3D &G ){

  Real Ek_sum;

  #ifndef PARALLEL_OMP
  Ek_sum = Get_Mean_Kinetic_Energy_function( G, 0, G.Grav.nz_local );
  #else
  Ek_sum = 0;
  Real Ek_sum_all[N_OMP_THREADS];
  #pragma omp parallel num_threads( N_OMP_THREADS )
  {
    int omp_id, n_omp_procs;
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local,  &g_start, &g_end  );
    Ek_sum_all[omp_id] = Get_Mean_Kinetic_Energy_function( G, g_start, g_end );

  }
  for ( int i=0; i<N_OMP_THREADS; i++ ){
    Ek_sum += Ek_sum_all[i];
  }
  #endif

  G.H.Ekin_mean = Ek_sum / ( G.Grav.nx_local * G.Grav.ny_local * G.Grav.nz_local);

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

        id  = (i) + (j)*nx_grid + (k)*ny_grid*nx_grid;
        imo = (i-1) + (j)*nx_grid + (k)*ny_grid*nx_grid;
        ipo = (i+1) + (j)*nx_grid + (k)*ny_grid*nx_grid;
        jmo = (i) + (j-1)*nx_grid + (k)*ny_grid*nx_grid;
        jpo = (i) + (j+1)*nx_grid + (k)*ny_grid*nx_grid;
        kmo = (i) + (j)*nx_grid + (k-1)*ny_grid*nx_grid;
        kpo = (i) + (j)*nx_grid + (k+1)*ny_grid*nx_grid;

        d = G.C.density[id];
        d_inv = 1/d;
        vx = G.C.momentum_x[id] * d_inv;
        vy = G.C.momentum_y[id] * d_inv;
        vz = G.C.momentum_z[id] * d_inv;
        E = G.C.Energy[id];

        // if (E < 0.0 || E != E) continue; // BUG: This leads to negative Energy
        Ek = 0.5*d*(vx*vx + vy*vy + vz*vz);
        ge1 = G.C.GasEnergy[id];
        ge2 = E - 0.5*d*(vx*vx + vy*vy + vz*vz);

        // #ifdef DE_EKINETIC_LIMIT
        // if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.02 && Ek/G.H.Ekin_mean > 0.4 ) {
        // // if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.02  ) {
        // // if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.075  ) {
        // #else
        // if (ge2 > 0.0 && E > 0.0 && ge2/E > 0.001 ) {
        // #endif //DE_EKINETIC_LIMIT
        // 
        //   G.C.GasEnergy[id] = ge2;
        //   ge1 = ge2;
        // }

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
        // else {
        //   if (ge1 > 0.0) G.C.Energy[id] += ge1 - ge2;
        //   else G.C.GasEnergy[id] = ge2;
        // }
        // G.C.Energy[id] = Ek + G.C.GasEnergy[id];
        // if ( fabs(( G.C.Energy[id] - (Ek + G.C.GasEnergy[id]) )  / G.C.Energy[id] ) > 1e-5 ) std::cout << "##Energy Error: " << G.C.Energy[id] << "  " << Ek + G.C.GasEnergy[id] << std::endl;
        // if (G.C.Energy[id] < 0 ) std::cout << "##Negative Energy after E_sync: " <<  G.C.Energy[id] << "  " << G.C.GasEnergy[id] << " " << d << " " << Ek << std::endl;
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
    int g_start, g_end;

    omp_id = omp_get_thread_num();
    n_omp_procs = omp_get_num_threads();
    Get_OMP_Grid_Indxs( n_omp_procs, omp_id, G.Grav.nz_local,  &g_start, &g_end  );

    Sync_Energies_3D_Host_Function( G, g_start, g_end );
    #ifdef TEMPERATURE_FLOOR
    #pragma omp barrier
    Apply_Temperature_Floor_Host( G, g_start, g_end );
    #endif
  }
  #endif

}


#ifdef TEMPERATURE_FLOOR
void Apply_Temperature_Floor_Host( Grid3D &G, int g_start, int g_end ){

  Real temp_floor = G.H.temperature_floor;

  #ifdef COSMOLOGY
  temp_floor *= 1 / (gama - 1) / MASS_HYDROGEN * K_BOLTZ * 1e-10;
  temp_floor /=  G.Cosmo.v_0_gas * G.Cosmo.v_0_gas / G.Cosmo.current_a / G.Cosmo.current_a;
  #endif

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
  Real d, u, temp;
  Real d_inv, vx, vy, vz, Ekin, E;
  Real u_new, delta_u;
  int k, j, i, id;
  int k_g, j_g, i_g;
  for ( k_g=g_start; k_g<g_end; k_g++ ){
    for ( j_g=0; j_g<ny_grav; j_g++ ){
      for ( i_g=0; i_g<nx_grav; i_g++ ){
        i = i_g + nGHST;
        j = j_g + nGHST;
        k = k_g + nGHST;
        id  = (i) + (j)*nx_grid + (k)*ny_grid*nx_grid;

        d = G.C.density[id];
        u = G.C.GasEnergy[id];

        d_inv = 1/d;
        vx = G.C.momentum_x[id] * d_inv;
        vy = G.C.momentum_y[id] * d_inv;
        vz = G.C.momentum_z[id] * d_inv;
        Ekin = 0.5 * d * (vx*vx + vy*vy + vz*vz);
        E = G.C.Energy[id];

        if ( fabs(( E - (Ekin + u) )  / E ) > 1e-5 ) std::cout << "##Energy Error: " << E << "  " << Ekin + u << std::endl;

        temp = u / d;
        if ( temp < temp_floor ){
          temp = temp_floor;
          u_new = temp * d  ;
          delta_u = u_new - u;
          G.C.GasEnergy[id] += delta_u;
        }
        
        temp = (E - Ekin) / d;
        if ( temp < temp_floor ){
          temp = temp_floor;
          u_new = temp * d  ;
          delta_u = u_new - u;
          G.C.Energy[id] += delta_u;
        }
        
        E = G.C.Energy[id];
        u = G.C.GasEnergy[id];
        if ( fabs(( E - (Ekin + u) )  / E ) > 1e-5 ) std::cout << "##Energy Error: " << E << "  " << Ekin + u << std::endl;
        if (G.C.Energy[id] < 0 ) std::cout << "##Negative Energy after Temp_floor: " <<  G.C.Energy[id] << "  " << G.C.GasEnergy[id] << " " << d << " " << Ekin << std::endl;
      }
    }
  }
}


#endif //TEMPERATURE_FLOOR
#endif //DE
#endif //GRAVITY_CPU

// #ifdef PARTICLES
// #include <stdio.h>
// #include <stdlib.h>
// #include <sstream>
// #include <iostream>
// #include"omp_functions.h"
// #include <math.h>
//
// // using namespace std;
//
// bool complete_potential;
//
// partID_t get_nParticles_omp( partID_t n_parts_local, int n_omp_procs, int omp_proc_id ){
//   partID_t n_parts_omp, parts_reminder;
//
//   parts_reminder = n_parts_local % n_omp_procs;
//   n_parts_omp = (partID_t) int(n_parts_local) / int(n_omp_procs);
//
//   if ( omp_proc_id == n_omp_procs-1 ) n_parts_omp += parts_reminder;
//
//   // std::cout << omp_proc_id << " "<< n_parts_local  << " " << n_omp_procs << " " << n_parts_omp << " "  << std::endl;
//
//   return n_parts_omp;
// }
//
// void get_omp_indxs( partID_t *omp_pIndx_start, partID_t *omp_pIndx_end, int *omp_gridIndx_start, int *omp_gridIndx_end, partID_t n_parts_local, int n_omp_procs, int omp_proc_id, int nGrid ){
//   partID_t n_parts_omp, n_parts_omp_local;
//   int nGrid_indxs, nGrid_indxs_reminder;
//   n_parts_omp_local = get_nParticles_omp( n_parts_local, n_omp_procs, omp_proc_id );
//   n_parts_omp = (partID_t) n_parts_local / n_omp_procs;
//
//   *omp_pIndx_start = omp_proc_id * n_parts_omp;
//   *omp_pIndx_end = *omp_pIndx_start + n_parts_omp_local;
//
//   nGrid_indxs = (int) floor(int(nGrid) / int(n_omp_procs));
//   nGrid_indxs_reminder =  nGrid % n_omp_procs;
//   if ( omp_proc_id == n_omp_procs-1 ) nGrid_indxs += nGrid_indxs_reminder;
//   *omp_gridIndx_start = (int) omp_proc_id * (int) floor(int(nGrid) / int(n_omp_procs));
//   *omp_gridIndx_end = *omp_gridIndx_start + nGrid_indxs;
// }
// #endif

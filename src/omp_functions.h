#ifdef PARTICLES
#ifndef OMP_FUNCTIONS_H
#define OMP_FUNCTIONS_H


#include "global.h"

extern bool complete_potential;

partID_t get_nParticles_omp( partID_t nLocal, int n_omp_procs, int omp_proc_id );
void get_omp_indxs( partID_t *omp_pIndx_start, partID_t *omp_pIndx_end, int *omp_gridIndx_start, int *omp_gridIndx_end, partID_t n_parts_local, int n_omp_procs, int omp_proc_id, int nGrid );
#endif //OMP_FUNCTIONS_H
#endif

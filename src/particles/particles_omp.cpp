#ifdef PARTICLES
#ifdef PARTICLES_OMP

#include "particles_omp.h"


part_int_t get_nParticles_omp( part_int_t n_parts_local, int n_omp_procs, int omp_proc_id ){



}

void Get_OMP_Indxs( part_int_t n_parts_local, int n_omp_procs, int omp_proc_id, int nGrid, part_int_t *omp_pIndx_start, part_int_t *omp_pIndx_end, int *omp_gridIndx_start, int *omp_gridIndx_end){

  part_int_t n_parts_omp, parts_reminder;

  parts_reminder = n_parts_local % n_omp_procs;
  n_parts_omp = n_parts_local/n_omp_procs;

  // return n_parts_omp

}


#endif
#endif

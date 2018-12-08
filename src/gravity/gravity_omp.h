#ifdef GRAVITY_OMP

#ifndef GRAVITY_OMP_H
#define GRAVITY_OMP_H

#include<stdio.h>
#include<stdlib.h>
#include"math.h"
#include"../global.h"
#include <iostream>
#include <omp.h>

void Get_OMP_Indxs( part_int_t n_parts_local, int n_omp_procs, int omp_proc_id, int nGrid, part_int_t *omp_pIndx_start, part_int_t *omp_pIndx_end, int *omp_gridIndx_start, int *omp_gridIndx_end);

void Get_OMP_Grid_Indxs(  int n_omp_procs, int omp_proc_id, int nGrid,  int *omp_gridIndx_start, int *omp_gridIndx_end  );
#endif
#endif

#ifdef GRAVITY
#ifdef POTENTIAL_PFFT

#ifndef POTENTIAL_PFFT_3D_H
#define POTENTIAL_PFFT_3D_H

#include "../global.h"
#include "../io.h"
#include "grav3D.h"
#include "poisson_solver_3D.h"
#include <stdlib.h>
#include <cmath>
#include <time.h>

#include <pfft.h>


// #if PRECISION == 1
// typedef fftReal Real_fftw;
// typedef cufftComplex Complex_cufft;
// #endif
//
// #if PRECISION == 2
// typedef double Real_fftw;
// typedef fftw_complex Complex_fftw;
// #endif


class Potential_PFFT_3D : public Poisson_Solver_3D
{
  public:

  int procID_pfft;
  int nproc_pfft;
  MPI_Comm comm_pfft;

  int nprocs_grid_pfft[3];
  int pcoords_pfft[3];
  int poffset_pfft[3];
  ptrdiff_t n_pfft[3];
  ptrdiff_t alloc_local;
  ptrdiff_t local_ni_pfft[3], local_i_start_pfft[3];
  ptrdiff_t local_no_pfft[3], local_o_start_pfft[3];
  ptrdiff_t local_ntrans_pfft[3], local_trans_start_pfft[3];

  pfft_plan plan_fwd;
  pfft_plan plan_bwd;

  int index_0;

  Real xMin;
  Real yMin;
  Real zMin;

  struct Fields
  {

    pfft_complex *transform;
    double *input;
    double *output;
    // Complex_fftw *transform;
    double *G;


  } F;

  Potential_PFFT_3D( void );

  virtual void Initialize( Grav3D Grav );

  void AllocateMemory_CPU( void );

  void Reset( void );

  void Copy_Input( Grav3D &Grav );
  void Copy_Output( Grav3D &Grav );
  //
  virtual void Get_K_for_Green_function( void );
  void Apply_G_Funtion( void );
  void Apply_K2_Funtion( void );
  virtual Real Get_Potential( Grav3D &Grav );

  void Get_Index_Global(int i, int j, int k, int *i_global, int *j_global, int *k_global);


};




#endif //POTENTIAL_CUFFT_H
#endif //POTENTIAL_CUFFT
#endif //GRAVITY

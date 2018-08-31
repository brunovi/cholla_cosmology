#ifdef GRAVITY
#ifdef POTENTIAL_FFTW

#ifndef POTENTIAL_FFTW_3D_H
#define POTENTIAL_FFTW_3D_H

#include "../global.h"
#include "../io.h"
#include "grav3D.h"
#include "poisson_solver_3D.h"
#include <fftw3.h>
#include <stdlib.h>
#include <cmath>
#include <time.h>


// #if PRECISION == 1
// typedef fftReal Real_fftw;
// typedef cufftComplex Complex_cufft;
// #endif
//
#if PRECISION == 2
typedef double Real_fftw;
typedef fftw_complex Complex_fftw;
#endif


class Potential_FFTW_3D : public Poisson_Solver_3D
{
  public:

  fftw_plan fftw_plan_fwd;
  fftw_plan fftw_plan_bwd;

  struct Fields
  {

    Complex_fftw *input;
    Complex_fftw *output;
    Complex_fftw *transform;

    Real *G;

  } F;

  Potential_FFTW_3D( void );

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

};




#endif //POTENTIAL_CUFFT_H
#endif //POTENTIAL_CUFFT
#endif //GRAVITY

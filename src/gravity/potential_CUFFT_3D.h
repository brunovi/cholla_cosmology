#ifdef GRAVITY
#ifdef POTENTIAL_CUFFT

#ifndef POTENTIAL_CUFFT_3D_H
#define POTENTIAL_CUFFT_3D_H

#include "../global.h"
#include "../io.h"
#include "grav3D.h"
#include "poisson_solver_3D.h"
#include <cufft.h>


#if PRECISION == 1
typedef cufftReal Real_cufft;
typedef cufftComplex Complex_cufft;
#endif

#if PRECISION == 2
typedef cufftDoubleReal Real_cufft;
typedef cufftDoubleComplex Complex_cufft;
#endif


class Potential_CUFFT_3D : public Poisson_Solver_3D
{
  public:

  cufftHandle plan_cufft_fwd;
  cufftHandle plan_cufft_bwd;

  struct Fields
  {

    Real_cufft *input_d;
    Real_cufft *output_d;
    Complex_cufft *transform_d;

    Real *output_h;

  } F;

  Potential_CUFFT_3D( void );

  virtual void Initialize( Grav3D Grav );

  void AllocateMemory_CPU( void );
  void AllocateMemory_GPU( void );

  void FreeMemory_GPU( void );

  void Copy_Input( Grav3D &Grav );
  void Copy_Output( Grav3D &Grav );

  virtual void Get_Potential( Grav3D &Grav );
};




#endif //POTENTIAL_CUFFT_H
#endif //POTENTIAL_CUFFT
#endif //GRAVITY

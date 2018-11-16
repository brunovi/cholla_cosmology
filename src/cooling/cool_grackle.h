#ifdef  COSMOLOGY
#ifdef COOLING_GRACKLE

#ifndef INIT_GRACKLE_H
#define INIT_GRACKLE_H

#include"../global.h"
// extern "C" {
#include <grackle.h>

class Cool_GK
{
  public:

  // int procID_pfft;
  // int nproc_pfft;
  // MPI_Comm comm_pfft;


Cool_GK( void );
void Initialize( );

};

#endif
#endif
#endif

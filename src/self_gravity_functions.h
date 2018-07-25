#ifdef SELF_GRAVITY

#ifndef SELF_GRAV_FUNC_H
#define SELF_GRAV_FUNC_H

#include"grid3D.h"
#include"global.h"
#include"mpi_pfft.h"

void CopyField_Host_To_Device( Real *field_h, Real *field_d, int n_cells );



#endif //SELF_GRAV_FUNC_H
#endif //SELF_GRAVITY

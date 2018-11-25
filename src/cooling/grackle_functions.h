#ifdef  COSMOLOGY
#ifdef COOLING_GRACKLE

#ifndef GRACKLE_FUCTIONS_H
#define GRACKLE_FUCTIONS_H

#include"../global.h"
#include"../grid3D.h"

extern "C" {
#include <grackle.h>
}

void Initialize_Grackle_Fields( Grid3D &G );

void Clear_Data_Grackle( Grid3D &G );

void Copy_Fields_to_Grackle( Grid3D &G );

void Do_Cooling_Step( Real dt, Grid3D &G );

void Set_Initial_Fields_Grackle( Grid3D &G );

#endif
#endif
#endif

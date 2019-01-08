#ifndef DUAL_ENERGY_H
#define DUAL_ENERGY_H

#include"grid3D.h"
#include"global.h"

#ifdef PARALLEL_OMP
#include"parallel_omp.h"
#endif

#ifdef DE_EKINETIC_LIMIT
void Get_Mean_Kinetic_Energy( Grid3D &G );
#endif

void Sync_Energies_3D_Host(Grid3D &G );

#ifdef TEMPERATURE_FLOOR
void Apply_Temperature_Floor_Host( Grid3D &G, int g_start, int g_end );
#endif

#endif

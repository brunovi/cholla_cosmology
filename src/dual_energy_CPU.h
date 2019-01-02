#ifndef DUAL_ENERGY_H
#define DUAL_ENERGY_H

#include"grid3D.h"
#include"global.h"

#ifdef DE_EKINETIC_LIMIT
void Get_Mean_Kinetic_Energy( Grid3D &G );
#endif

void Sync_Energies_3D_Host(Grid3D &G );

#endif

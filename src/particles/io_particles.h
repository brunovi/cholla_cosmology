#ifndef IO_PARTICLES_H
#define IO_PARTICLES_H

#include"../global.h"
#include"../grid3D.h"
#include"../io.h"
#include"particles_3D.h"


/* Write the data */
void WriteData_Particles( Grid3D &G, struct parameters P, int nfile);

void OutputData_Particles( Grid3D &G, struct parameters P, int nfile);

void Load_Particles_Data( Particles_3D &Particles, struct parameters P);

#ifdef HDF5
void Load_Particles_Data_HDF5(hid_t file_id, int nfile, Particles_3D &Particles );
#endif
// #ifdef   MPI_CHOLLA
// void OutputDataMPI_Particles(Part3D &Particles, Grid3D G, struct parameters P, int nfile);
// #endif /*MPI_CHOLLA*/

#endif /*IO_PARTICLES_H*/

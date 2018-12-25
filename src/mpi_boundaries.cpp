#include"grid3D.h"
#include"mpi_routines.h"
#include"io.h"
#include"error_handling.h"
#include <iostream>
#ifdef MPI_CHOLLA

#ifdef PARTICLES
#include "particles/particles_boundaries.h"


void Grid3D::Finish_Particles_Transfer( void ){

  Particles.Remove_Transfered_Particles();

}
#endif


void Grid3D::Set_Boundaries_MPI(struct parameters P)
{
  int flags[6] = {0,0,0,0,0,0};

  if(Check_Custom_Boundary(&flags[0],P))
  {
    //perform custom boundaries
    Custom_Boundary(P.custom_bcnd);
  }

  switch(flag_decomp)
  {
    case SLAB_DECOMP:
      Set_Boundaries_MPI_SLAB(flags,P);
      break;
    case BLOCK_DECOMP:
      #ifdef PARTICLES
      if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
        Particles.Clear_Particles_For_Transfer();
        // Particles.Select_Particles_to_Transfer_All();
      }
      #endif
      Set_Boundaries_MPI_BLOCK(flags,P);
      break;
  }
}


void Grid3D::Set_Boundaries_MPI_SLAB(int *flags, struct parameters P)
{
  //perform boundary conditions around
  //edges of communication region on x faces

  if(flags[0]==5)
    Set_Edge_Boundaries(0,flags);
  if(flags[1]==5)
    Set_Edge_Boundaries(1,flags);

  //1) load and post comm for buffers
  Load_and_Send_MPI_Comm_Buffers(0, flags);

  //2) perform any additional boundary conditions
  //including whether the x face is non-MPI
  if(H.nx>1)
  {
    Set_Boundaries(0,flags);
    Set_Boundaries(1,flags);
  }
  if(H.ny>1)
  {
    Set_Boundaries(2,flags);
    Set_Boundaries(3,flags);
  }
  if(H.nz>1)
  {
    Set_Boundaries(4,flags);
    Set_Boundaries(5,flags);
  }

  //3) wait for sends and receives to finish
  //   and then load ghost cells from comm buffers
  if(flags[0]==5 || flags[1]==5)
    Wait_and_Unload_MPI_Comm_Buffers_SLAB(flags);
}


void Grid3D::Set_Boundaries_MPI_BLOCK(int *flags, struct parameters P)
{
  if (H.nx > 1) {

    /* Step 1 - Send MPI x-boundaries */
    if (flags[0]==5 || flags[1]==5) {
      Load_and_Send_MPI_Comm_Buffers(0, flags);
    }

    /* Step 2 - Set non-MPI x-boundaries */
    Set_Boundaries(0, flags);
    Set_Boundaries(1, flags);

    /* Step 3 - Receive MPI x-boundaries */
    if (flags[0]==5 || flags[1]==5) {
      Wait_and_Unload_MPI_Comm_Buffers_BLOCK(0, flags);
      #ifdef PARTICLES
      if (Particles.TRANSFER_PARTICLES_BOUNDARIES) Wait_and_Unload_MPI_Comm_Particles_Buffers_BLOCK(0, flags);
      #endif
    }
  }
  MPI_Barrier(world);
  if (H.ny > 1) {

    /* Step 4 - Send MPI y-boundaries */
    if (flags[2]==5 || flags[3]==5) {
      Load_and_Send_MPI_Comm_Buffers(1, flags);
    }

    /* Step 5 - Set non-MPI y-boundaries */
    Set_Boundaries(2, flags);
    Set_Boundaries(3, flags);

    /* Step 6 - Receive MPI y-boundaries */
    if (flags[2]==5 || flags[3]==5) {
      Wait_and_Unload_MPI_Comm_Buffers_BLOCK(1, flags);
      #ifdef PARTICLES
      if (Particles.TRANSFER_PARTICLES_BOUNDARIES) Wait_and_Unload_MPI_Comm_Particles_Buffers_BLOCK(1, flags);
      #endif
    }
  }
  MPI_Barrier(world);
  if (H.nz > 1) {

    /* Step 7 - Send MPI z-boundaries */
    if (flags[4]==5 || flags[5]==5) {
      Load_and_Send_MPI_Comm_Buffers(2, flags);
    }

    /* Step 8 - Set non-MPI z-boundaries */
    Set_Boundaries(4, flags);
    Set_Boundaries(5, flags);

    /* Step 9 - Receive MPI z-boundaries */
    if (flags[4]==5 || flags[5]==5) {
      Wait_and_Unload_MPI_Comm_Buffers_BLOCK(2, flags);
      #ifdef PARTICLES
      if (Particles.TRANSFER_PARTICLES_BOUNDARIES) Wait_and_Unload_MPI_Comm_Particles_Buffers_BLOCK(2, flags);
      #endif
    }
  }

  #ifdef PARTICLES
  if ( Particles.TRANSFER_PARTICLES_BOUNDARIES)  Finish_Particles_Transfer();
  #endif

}


void Grid3D::Set_Edge_Boundaries(int dir, int *flags)
{
  int iedge;
  int i, j, k;
  int imin[3] = {0,0,0};
  int imax[3] = {H.nx,H.ny,H.nz};
  Real a[3]   = {1,1,1};  //sign of momenta
  int idx;    //index of a real cell
  int gidx;   //index of a ghost cell

  int nedge = 0;

  if(H.ny>1)
    nedge = 2;
  if(H.ny>1 && H.nz>1)
    nedge = 8;

  for(iedge=0;iedge<nedge;iedge++)
  {
    //set the edge or corner extents
    Set_Edge_Boundary_Extents(dir, iedge, &imin[0], &imax[0]);

    /*set ghost cells*/
    for (i=imin[0]; i<imax[0]; i++) {
      for (j=imin[1]; j<imax[1]; j++) {
        for (k=imin[2]; k<imax[2]; k++) {

          //reset sign of momenta
          a[0] = 1.;
          a[1] = 1.;
          a[2] = 1.;

          //find the ghost cell index
          gidx = i + j*H.nx + k*H.nx*H.ny;

          //find the corresponding real cell index and momenta signs
          idx  = Set_Boundary_Mapping(i,j,k,flags,&a[0]);

          //idx will be >= 0 if the boundary mapping function has
          //not set this ghost cell by hand, for instance for analytical
          //boundary conditions
          //
          //Otherwise, the boundary mapping function will set idx<0
          //if this ghost cell has been set by hand
          if(idx>=0)
          {
            //set the ghost cell value
            C.density[gidx]    = C.density[idx];
            C.momentum_x[gidx] = C.momentum_x[idx]*a[0];
            C.momentum_y[gidx] = C.momentum_y[idx]*a[1];
            C.momentum_z[gidx] = C.momentum_z[idx]*a[2];
            C.Energy[gidx]     = C.Energy[idx];
            #ifdef DE
            C.GasEnergy[gidx]  = C.GasEnergy[idx];
            #endif
            #ifdef SCALAR
            for (int ii=0; ii<NSCALARS; ii++) {
              C.scalar[gidx + ii*H.n_cells]  = C.scalar[idx+ii*H.n_cells];
            }
            #endif
            #ifdef GRAVITY
            C.Grav_potential[gidx]  = C.Grav_potential[idx];
            #endif
          }
        }
      }
    }
  }
}


/*! \fn Set_Edge_Boundary_Extents(int dir, int edge, int *imin, int *imax)
 *  \brief Set the extents of the edge and corner ghost regions we want to
 *         initialize first, so we can then do communication */
void Grid3D::Set_Edge_Boundary_Extents(int dir, int edge, int *imin, int *imax)
{

  int i, j, k;
  int ni, nj, nk;


  //two D case
  if(H.ny>1 && H.nz==1)
  {
    if(dir==0 || dir==1)
    {
      ni = H.nx;
      i  = 0;
      nj = H.ny;
      j  = 1;
      nk = 1;
      k  = 2;
    }

    if(dir==2 || dir==3)
    {
      ni = H.ny;
      i  = 1;
      nj = H.nx;
      j  = 0;
      nk = 1;
      k  = 2;
    }

    //upper or lower face?
    if(!(dir%2))
    {
      *(imin+i) = 0;
      *(imax+i) = 2*H.n_ghost;


    }else{
      *(imin+i) = ni-2*H.n_ghost;
      *(imax+i) = ni;
    }

    if(edge==0)
    {
      //lower square
      *(imin+j) = 0;
      *(imax+j) = 2*H.n_ghost;
    }else{
      //upper square
      *(imin+j) = nj-2*H.n_ghost;
      *(imax+j) = nj;
    }
    return;
  }

  //three D case
  if(H.ny>1 && H.nz>1)
  {

    if(dir==0 || dir==1)
    {
      ni = H.nx;
      i  = 0;
      nj = H.ny;
      j  = 1;
      nk = H.nz;
      k  = 2;
    }

    if(dir==2 || dir==3)
    {
      ni = H.ny;
      i  = 1;
      nj = H.nz;
      j  = 2;
      nk = H.nx;
      k  = 0;
    }

    if(dir==4 || dir==5)
    {
      ni = H.nz;
      i  = 2;
      nj = H.nx;
      j  = 0;
      nk = H.ny;
      k  = 1;
    }


    //upper or lower face?
    if(!(dir%2))
    {
      *(imin+i) = H.n_ghost;
      *(imax+i) = 2*H.n_ghost;
    }else{
      *(imin+i) = ni-2*H.n_ghost;
      *(imax+i) = ni-H.n_ghost;
    }

    //edges and corners clockwise from lower left corner
    switch(edge)
    {
      //lower left corner
      case 0:
        *(imin+j) = 0;
        *(imax+j) = H.n_ghost;
        *(imin+k) = 0;
        *(imax+k) = H.n_ghost;
        break;

      //left edge
      case 1:
        *(imin+j) = H.n_ghost;
        *(imax+j) = nj-H.n_ghost;
        *(imin+k) = 0;
        *(imax+k) = H.n_ghost;
        break;

      //upper left corner
      case 2:
        *(imin+j) = nj-H.n_ghost;
        *(imax+j) = nj;
        *(imin+k) = 0;
        *(imax+k) = H.n_ghost;
        break;

      //upper edge
      case 3:
        *(imin+j) = nj-H.n_ghost;
        *(imax+j) = nj;
        *(imin+k) = H.n_ghost;
        *(imax+k) = nk-H.n_ghost;
        break;

      //upper right corner
      case 4:
        *(imin+j) = nj-H.n_ghost;
        *(imax+j) = nj;
        *(imin+k) = nk-H.n_ghost;
        *(imax+k) = nk;
        break;

      //right edge
      case 5:
        *(imin+j) = H.n_ghost;
        *(imax+j) = nj-H.n_ghost;
        *(imin+k) = nk-H.n_ghost;
        *(imax+k) = nk;
        break;

      //lower right corner
      case 6:
        *(imin+j) = 0;
        *(imax+j) = H.n_ghost;
        *(imin+k) = nk-H.n_ghost;
        *(imax+k) = nk;
        break;

      //lower edge
      case 7:
        *(imin+j) = 0;
        *(imax+j) = H.n_ghost;
        *(imin+k) = H.n_ghost;
        *(imax+k) = nk-H.n_ghost;
        break;
    }
  }
}

void Grid3D::Load_and_Send_MPI_Comm_Buffers(int dir, int *flags)
{

  switch(flag_decomp)
  {
    case SLAB_DECOMP:
      /*load communication buffers*/
      Load_and_Send_MPI_Comm_Buffers_SLAB(flags);
      break;
    case BLOCK_DECOMP:
      /*load communication buffers*/
      #ifdef PARTICLES
      if ( Particles.TRANSFER_PARTICLES_BOUNDARIES){
        Particles.Select_Particles_to_Transfer( dir );
      }
      #endif
      Load_and_Send_MPI_Comm_Buffers_BLOCK(dir, flags);
      break;
  }

}

void Grid3D::Load_and_Send_MPI_Comm_Buffers_SLAB(int *flags)
{

  int i, j, k, ii;
  int gidx;
  int idx;
  int ireq = 0;

  int offset = H.n_ghost*H.ny*H.nz;

  /*check left side*/
  if(flags[0]==5)
  {
    //load left x communication buffer
    for(i=0;i<H.n_ghost;i++)
    {
      for(j=0;j<H.ny;j++)
      {
        for(k=0;k<H.nz;k++)
        {
          idx  = (i+H.n_ghost) + j*H.nx      + k*H.nx*H.ny;
          gidx = i             + j*H.n_ghost + k*H.n_ghost*H.ny;

          for (ii=0; ii<H.n_fields; ii++) {
            *(send_buffer_0 + gidx + ii*offset) = C.density[ii*H.n_cells + idx];
          }
        }
      }
    }


    //post non-blocking receive left x communication buffer
    MPI_Irecv(recv_buffer_0, recv_buffer_length, MPI_CHREAL, source[0], 0, world, &recv_request[ireq]);

    //non-blocking send left x communication buffer
    MPI_Isend(send_buffer_0, send_buffer_length, MPI_CHREAL, dest[0],   1,    world, &send_request[0]);

    //remember how many recv's this proc expects
    ireq++;
  }

  /*check right side*/
  if(flags[1]==5)
  {
    //load right x communication buffer
    for(i=0;i<H.n_ghost;i++)
    {
      for(j=0;j<H.ny;j++)
      {
        for(k=0;k<H.nz;k++)
        {
          idx  = (i+H.nx-2*H.n_ghost) + j*H.nx      + k*H.nx*H.ny;
          gidx = i                    + j*H.n_ghost + k*H.n_ghost*H.ny;

          for (ii=0; ii<H.n_fields; ii++) {
            *(send_buffer_1 + gidx + ii*offset) = C.density[ii*H.n_cells + idx];
          }
        }
      }
    }


    //post non-blocking receive right x communication buffer
    MPI_Irecv(recv_buffer_1, recv_buffer_length, MPI_CHREAL, source[1], 1, world, &recv_request[ireq]);

    //non-blocking send right x communication buffer
    MPI_Isend(send_buffer_1, send_buffer_length, MPI_CHREAL, dest[1], 0, world, &send_request[1]);

    //remember how many recv's this proc expects
    ireq++;
  }

  //done!
}

#ifdef GRAVITY
int Grid3D::Load_Potential_To_Buffer( int direction, int side, Real *buffer, int buffer_start  ){


  int i, j, k, indx, indx_buff, length;
  int nGHST, nx_g, ny_g, nz_g, nx, ny, nz;
  nGHST = H.n_ghost;
  nx = Grav.nx_local;
  ny = Grav.ny_local;
  nz = Grav.nz_local;
  nx_g = Grav.nx_local + 2*nGHST;
  ny_g = Grav.ny_local + 2*nGHST;
  nz_g = Grav.nz_local + 2*nGHST;

  //Load X boundaries
  if (direction == 0){
    for ( k=nGHST; k<nz+nGHST; k++ ){
      for ( j=nGHST; j<ny+nGHST; j++ ){
        for ( i=0; i<nGHST; i++ ){
          if ( side == 0 ) indx = (i+nGHST) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (nx_g - 2*nGHST + i) + (j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = (j-nGHST) + (k-nGHST)*ny + i*ny*nz ;
          buffer[buffer_start+indx_buff] = C.Grav_potential[indx];
          // buffer[buffer_start+indx_buff] = -10;
          length = nGHST * nz * ny;
        }
      }
    }
  }

  //Load Y boundaries
  if (direction == 1){
    for ( k=nGHST; k<nz+nGHST; k++ ){
      for ( j=0; j<nGHST; j++ ){
        for ( i=nGHST; i<nx+nGHST; i++ ){
          if ( side == 0 ) indx = (i) + (j+nGHST)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (ny_g - 2*nGHST + j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = (i-nGHST) + (k-nGHST)*nx + j*nx*nz ;
          buffer[buffer_start+indx_buff] = C.Grav_potential[indx];
          // buffer[buffer_start+indx_buff] = -10;
          length = nGHST * nz * nx;
        }
      }
    }
  }

  //Load Z boundaries
  if (direction == 2){
    for ( k=0; k<nGHST; k++ ){
      for ( j=nGHST; j<ny+nGHST; j++ ){
        for ( i=nGHST; i<nx+nGHST; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k+nGHST)*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (j)*nx_g + (nz_g - 2*nGHST + k)*nx_g*ny_g;
          indx_buff = (i-nGHST) + (j-nGHST)*nx + k*nx*ny ;
          buffer[buffer_start+indx_buff] = C.Grav_potential[indx];
          // buffer[buffer_start+indx_buff] = -10;
          length = nGHST * nx * ny;
        }
      }
    }
  }
  return length;
}


void Grid3D::Unload_Potential_from_Buffer( int direction, int side, Real *buffer, int buffer_start  ){


  int i, j, k, indx, indx_buff;
  int nGHST, nx_g, ny_g, nz_g, nx, ny, nz;
  nGHST = H.n_ghost;
  nx = Grav.nx_local;
  ny = Grav.ny_local;
  nz = Grav.nz_local;
  nx_g = Grav.nx_local + 2*nGHST;
  ny_g = Grav.ny_local + 2*nGHST;
  nz_g = Grav.nz_local + 2*nGHST;

  //Load X boundaries
  if (direction == 0){
    for ( k=nGHST; k<nz+nGHST; k++ ){
      for ( j=nGHST; j<ny+nGHST; j++ ){
        for ( i=0; i<nGHST; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (nx_g - nGHST + i) + (j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = (j-nGHST) + (k-nGHST)*ny + i*ny*nz ;
          C.Grav_potential[indx] = buffer[buffer_start+indx_buff];
        }
      }
    }
  }

  //Load Y boundaries
  if (direction == 1){
    for ( k=nGHST; k<nz+nGHST; k++ ){
      for ( j=0; j<nGHST; j++ ){
        for ( i=nGHST; i<nx+nGHST; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (ny_g - nGHST + j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = (i-nGHST) + (k-nGHST)*nx + j*nx*nz ;
          C.Grav_potential[indx] = buffer[buffer_start+indx_buff];
        }
      }
    }
  }

  //Load Z boundaries
  if (direction == 2){
    for ( k=0; k<nGHST; k++ ){
      for ( j=nGHST; j<ny+nGHST; j++ ){
        for ( i=nGHST; i<nx+nGHST; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (j)*nx_g + (nz_g - nGHST + k)*nx_g*ny_g;
          indx_buff = (i-nGHST) + (j-nGHST)*nx + k*nx*ny ;
          C.Grav_potential[indx] = buffer[buffer_start+indx_buff];
        }
      }
    }
  }
}





#endif




#ifdef PARTICLES
int Grid3D::Load_Particles_Density_Boundary_to_Buffer( int direction, int side, Real *buffer  ){

  int i, j, k, indx, indx_buff, buffer_length;
  int nGHST, nx_g, ny_g, nz_g, nx, ny, nz;
  nGHST = Particles.G.n_ghost_particles_grid;
  nx = Particles.G.nx_local;
  ny = Particles.G.ny_local;
  nz = Particles.G.nz_local;
  nx_g = Particles.G.nx_local + 2*nGHST;
  ny_g = Particles.G.ny_local + 2*nGHST;
  nz_g = Particles.G.nz_local + 2*nGHST;

  //Load Z boundaries
  if (direction == 2){
    for ( k=0; k<nGHST; k++ ){
      for ( j=0; j<ny_g; j++ ){
        for ( i=0; i<nx_g; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (j)*nx_g + (nz_g - nGHST + k)*nx_g*ny_g;
          indx_buff = i + j*nx_g + k*nx_g*ny_g ;
          buffer[indx_buff] = Particles.G.density[indx];
        }
      }
    }
    buffer_length = nGHST * nx_g * ny_g;
  }

  //Load Y boundaries
  if (direction == 1){
    for ( k=0; k<nz_g; k++ ){
      for ( j=0; j<nGHST; j++ ){
        for ( i=0; i<nx_g; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (ny_g - nGHST + j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = i + k*nx_g + j*nx_g*nz_g ;
          buffer[indx_buff] = Particles.G.density[indx];
        }
      }
    }
    buffer_length = nGHST * nx_g * nz_g;
  }

  //Load X boundaries
  if (direction == 0){
    for ( k=0; k<nz_g; k++ ){
      for ( j=0; j<ny_g; j++ ){
        for ( i=0; i<nGHST; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (nx_g - nGHST + i) + (j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = j + k*ny_g + i*ny_g*nz_g ;
          buffer[indx_buff] = Particles.G.density[indx];
        }
      }
    }
    buffer_length = nGHST * ny_g * nz_g;
  }


  return buffer_length;
}

void Grid3D::Unload_Particles_Density_Boundary_From_Buffer( int direction, int side, Real *buffer  ){

  int i, j, k, indx, indx_buff, buffer_length;
  int nGHST, nx_g, ny_g, nz_g, nx, ny, nz;
  nGHST = Particles.G.n_ghost_particles_grid;
  nx = Particles.G.nx_local;
  ny = Particles.G.ny_local;
  nz = Particles.G.nz_local;
  nx_g = Particles.G.nx_local + 2*nGHST;
  ny_g = Particles.G.ny_local + 2*nGHST;
  nz_g = Particles.G.nz_local + 2*nGHST;

  // //Unload Z boundaries
  if (direction == 2){
    for ( k=0; k<nGHST; k++ ){
      for ( j=0; j<ny_g; j++ ){
        for ( i=0; i<nx_g; i++ ){
          if ( side == 0 ) indx = (i) + (j)*nx_g + (k + nGHST )*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (j)*nx_g + (nz_g - 2*nGHST + k)*nx_g*ny_g;
          indx_buff = i + j*nx_g + k*nx_g*ny_g ;
          Particles.G.density[indx] += buffer[indx_buff];
        }
      }
    }
  }

  //Unload Y boundaries
  if (direction == 1){
    for ( k=0; k<nz_g; k++ ){
      for ( j=0; j<nGHST; j++ ){
        for ( i=0; i<nx_g; i++ ){
          if ( side == 0 ) indx = (i) + (j + nGHST)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (i) + (ny_g - 2*nGHST + j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = i + k*nx_g + j*nx_g*nz_g ;
          Particles.G.density[indx] += buffer[indx_buff];
        }
      }
    }
  }

  //Unload X boundaries
  if (direction == 0){
    for ( k=0; k<nz_g; k++ ){
      for ( j=0; j<ny_g; j++ ){
        for ( i=0; i<nGHST; i++ ){
          if ( side == 0 ) indx = (i+nGHST) + (j)*nx_g + (k)*nx_g*ny_g;
          if ( side == 1 ) indx = (nx_g - 2*nGHST + i) + (j)*nx_g + (k)*nx_g*ny_g;
          indx_buff = j + k*ny_g + i*ny_g*nz_g ;
          Particles.G.density[indx] += buffer[indx_buff];
          // Particles.G.density[indx] += 1;
        }
      }
    }
  }

}

void Grid3D::Unload_Particles_from_Buffer_X0(){
  Particles.Unload_Particles_from_Buffer( 0, 0, recv_buffer_x0_particles, Particles.n_recv_x0, send_buffer_y0_particles, send_buffer_y1_particles, send_buffer_z0_particles, send_buffer_z1_particles, buffer_length_particles_y0_send , buffer_length_particles_y1_send, buffer_length_particles_z0_send, buffer_length_particles_z1_send);
}

void Grid3D::Unload_Particles_from_Buffer_X1(){
  Particles.Unload_Particles_from_Buffer( 0, 1, recv_buffer_x1_particles, Particles.n_recv_x1, send_buffer_y0_particles, send_buffer_y1_particles, send_buffer_z0_particles, send_buffer_z1_particles, buffer_length_particles_y0_send , buffer_length_particles_y1_send, buffer_length_particles_z0_send, buffer_length_particles_z1_send);
}

void Grid3D::Unload_Particles_from_Buffer_Y0(){
  Particles.Unload_Particles_from_Buffer( 1, 0, recv_buffer_y0_particles, Particles.n_recv_y0, send_buffer_y0_particles, send_buffer_y1_particles, send_buffer_z0_particles, send_buffer_z1_particles, buffer_length_particles_y0_send , buffer_length_particles_y1_send, buffer_length_particles_z0_send, buffer_length_particles_z1_send);
}

void Grid3D::Unload_Particles_from_Buffer_Y1(){
  Particles.Unload_Particles_from_Buffer( 1, 1, recv_buffer_y1_particles, Particles.n_recv_y1, send_buffer_y0_particles, send_buffer_y1_particles, send_buffer_z0_particles, send_buffer_z1_particles, buffer_length_particles_y0_send , buffer_length_particles_y1_send, buffer_length_particles_z0_send, buffer_length_particles_z1_send);
}

void Grid3D::Unload_Particles_from_Buffer_Z0(){
  Particles.Unload_Particles_from_Buffer( 2, 0, recv_buffer_z0_particles, Particles.n_recv_z0, send_buffer_y0_particles, send_buffer_y1_particles, send_buffer_z0_particles, send_buffer_z1_particles, buffer_length_particles_y0_send , buffer_length_particles_y1_send, buffer_length_particles_z0_send, buffer_length_particles_z1_send);
}

void Grid3D::Unload_Particles_from_Buffer_Z1(){
  Particles.Unload_Particles_from_Buffer( 2, 1, recv_buffer_z1_particles, Particles.n_recv_z1, send_buffer_y0_particles, send_buffer_y1_particles, send_buffer_z0_particles, send_buffer_z1_particles, buffer_length_particles_y0_send , buffer_length_particles_y1_send, buffer_length_particles_z0_send, buffer_length_particles_z1_send);
}

#endif

int Grid3D::Load_Hydro_Buffer_X0(){

  int i, j, k, ii;
  int gidx;
  int idx;
  int offset;

  // load left x communication buffer
  // 1D
  if (H.ny == 1 && H.nz == 1) {
    offset = H.n_ghost;
    for (i=0;i<H.n_ghost;i++) {
      idx = (i+H.n_ghost);
      gidx = i;
      for (ii=0; ii<H.n_fields; ii++) {
        *(send_buffer_x0 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
      }
    }
  }
  // 2D
  if (H.ny > 1 && H.nz == 1) {
    offset = H.n_ghost*(H.ny-2*H.n_ghost);
    for (i=0;i<H.n_ghost;i++) {
      for (j=0;j<H.ny-2*H.n_ghost;j++) {
        idx = (i+H.n_ghost) + (j+H.n_ghost)*H.nx;
        gidx = i + j*H.n_ghost;
        for (ii=0; ii<H.n_fields; ii++) {
          *(send_buffer_x0 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
        }
      }
    }
  }
  // 3D
  if (H.ny > 1 && H.nz > 1) {
    offset = H.n_ghost*(H.ny-2*H.n_ghost)*(H.nz-2*H.n_ghost);
    for(i=0;i<H.n_ghost;i++)
    {
      for(j=0;j<H.ny-2*H.n_ghost;j++)
      {
        for(k=0;k<H.nz-2*H.n_ghost;k++)
        {
          idx  = (i+H.n_ghost) + (j+H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
          gidx = i + j*H.n_ghost + k*H.n_ghost*(H.ny-2*H.n_ghost);
          for (ii=0; ii<H.n_fields; ii++) {
            *(send_buffer_x0 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
          }
        }
      }
    }
  }
  return x_buffer_length;
}

int Grid3D::Load_Hydro_Buffer_X1(){
  int i, j, k, ii;
  int gidx;
  int idx;
  int offset;

  // load right x communication buffer
  // 1D
  if (H.ny == 1 && H.nz == 1) {
    offset = H.n_ghost;
    for (i=0;i<H.n_ghost;i++) {
      idx = (i+H.nx-2*H.n_ghost);
      gidx = i;
      for (ii=0; ii<H.n_fields; ii++) {
        *(send_buffer_x1 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
      }
    }
  }
  // 2D
  if (H.ny > 1 && H.nz == 1) {
    offset = H.n_ghost*(H.ny-2*H.n_ghost);
    for (i=0;i<H.n_ghost;i++) {
      for (j=0;j<H.ny-2*H.n_ghost;j++) {
        idx = (i+H.nx-2*H.n_ghost) + (j+H.n_ghost)*H.nx;
        gidx = i + j*H.n_ghost;
        for (ii=0; ii<H.n_fields; ii++) {
          *(send_buffer_x1 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
        }
      }
    }
  }
  // 3D
  if (H.ny > 1 && H.nz > 1) {
    offset = H.n_ghost*(H.ny-2*H.n_ghost)*(H.nz-2*H.n_ghost);
    for(i=0;i<H.n_ghost;i++)
    {
      for(j=0;j<H.ny-2*H.n_ghost;j++)
      {
        for(k=0;k<H.nz-2*H.n_ghost;k++)
        {
          idx  = (i+H.nx-2*H.n_ghost) + (j+H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
          gidx = i + j*H.n_ghost + k*H.n_ghost*(H.ny-2*H.n_ghost);
          for (ii=0; ii<H.n_fields; ii++) {
            *(send_buffer_x1 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
          }
        }
      }
    }
  }
  return x_buffer_length;
}

int Grid3D::Load_Hydro_Buffer_Y0(){

  int i, j, k, ii;
  int gidx;
  int idx;
  int offset;

  // load left y communication buffer
  // 2D
  if (H.nz == 1) {
    offset = H.n_ghost*H.nx;
    for (i=0;i<H.nx;i++) {
      for (j=0;j<H.n_ghost;j++) {
        idx = i + (j+H.n_ghost)*H.nx;
        gidx = i + j*H.nx;
        for (ii=0; ii<H.n_fields; ii++) {
          *(send_buffer_y0 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
        }
      }
    }
  }
  // 3D
  if (H.nz > 1) {
    offset = H.n_ghost*H.nx*(H.nz-2*H.n_ghost);
    for(i=0;i<H.nx;i++)
    {
      for(j=0;j<H.n_ghost;j++)
      {
        for(k=0;k<H.nz-2*H.n_ghost;k++)
        {
          idx  = i + (j+H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
          gidx = i + j*H.nx + k*H.nx*H.n_ghost;
          for (ii=0; ii<H.n_fields; ii++) {
            *(send_buffer_y0 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
          }
        }
      }
    }
  }
  return y_buffer_length;
}

int Grid3D::Load_Hydro_Buffer_Y1(){
  int i, j, k, ii;
  int gidx;
  int idx;
  int offset;
  // load right y communication buffer
  // 2D
  if (H.nz == 1) {
    offset = H.n_ghost*H.nx;
    for (i=0;i<H.nx;i++) {
      for (j=0;j<H.n_ghost;j++) {
        idx = i + (j+H.ny-2*H.n_ghost)*H.nx;
        gidx = i + j*H.nx;
        for (ii=0; ii<H.n_fields; ii++) {
          *(send_buffer_y1 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
        }
      }
    }
  }
  // 3D
  if (H.nz > 1) {
    offset = H.n_ghost*H.nx*(H.nz-2*H.n_ghost);
    for(i=0;i<H.nx;i++)
    {
      for(j=0;j<H.n_ghost;j++)
      {
        for(k=0;k<H.nz-2*H.n_ghost;k++)
        {
          idx  = i + (j+H.ny-2*H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
          gidx = i + j*H.nx + k*H.nx*H.n_ghost;
          for (ii=0; ii<H.n_fields; ii++) {
            *(send_buffer_y1 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
          }
        }
      }
    }
  }
  return y_buffer_length;
}

int Grid3D::Load_Hydro_Buffer_Z0(){
  int i, j, k, ii;
  int gidx;
  int idx;
  int offset;

  // left z communication buffer
  // 3D
  offset = H.n_ghost*H.nx*H.ny;
  for(i=0;i<H.nx;i++)
  {
    for(j=0;j<H.ny;j++)
    {
      for(k=0;k<H.n_ghost;k++)
      {
        idx  = i + j*H.nx + (k+H.n_ghost)*H.nx*H.ny;
        gidx = i + j*H.nx + k*H.nx*H.ny;
        for (ii=0; ii<H.n_fields; ii++) {
          *(send_buffer_z0 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
        }
      }
    }
  }
  return z_buffer_length;
}

int Grid3D::Load_Hydro_Buffer_Z1(){

  int i, j, k, ii;
  int gidx;
  int idx;
  int offset;

  // load right z communication buffer
  offset = H.n_ghost*H.nx*H.ny;
  for(i=0;i<H.nx;i++)
  {
    for(j=0;j<H.ny;j++)
    {
      for(k=0;k<H.n_ghost;k++)
      {
        idx  = i + j*H.nx + (k+H.nz-2*H.n_ghost)*H.nx*H.ny;
        gidx = i + j*H.nx + k*H.nx*H.ny;
        for (ii=0; ii<H.n_fields; ii++) {
          *(send_buffer_z1 + gidx + ii*offset) = C.density[idx + ii*H.n_cells];
        }
      }
    }
  }
  return z_buffer_length;
}



void Grid3D::Load_and_Send_MPI_Comm_Buffers_BLOCK(int dir, int *flags)
{

  int ireq;
  ireq = 0;
  bool transfer_buffer;
  transfer_buffer = true;


  //Flag to omit hydro transfer when doing particles_density Transfering
  int buffer_length;


  #ifdef PARTICLES
  int ireq_n_particles, ireq_particles_transfer;
  ireq_n_particles = 0;
  ireq_particles_transfer = 0;
  int n_transfer_secondary, buffer_length_secondary;
  MPI_Status status_particles_secondary_0, status_particles_secondary_1;
  #endif


  /* x boundaries */
  if(dir == 0)
  {
    if (flags[0]==5) {

      if( H.TRANSFER_HYDRO_BOUNDARIES ) buffer_length = Load_Hydro_Buffer_X0();

      #ifdef GRAVITY
      if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ){
        buffer_length = Load_Potential_To_Buffer( 0, 0, send_buffer_x0, 0 );
      }
      #endif

      #ifdef PARTICLES
      if ( Particles.TRANSFER_DENSITY_BOUNDARIES) {
        buffer_length = Load_Particles_Density_Boundary_to_Buffer( 0, 0, send_buffer_x0  );
      }
      else if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
        transfer_buffer = false;
        MPI_Irecv(&Particles.n_recv_x0, 1, MPI_PART_INT, source[0], 0, world, &recv_request_n_particles[ireq_n_particles]);
        MPI_Isend(&Particles.n_send_x0, 1, MPI_PART_INT, dest[0],   1, world, &send_request_n_particles[0]);
        ireq_n_particles ++;
        // if ( Particles.n_send_x0 > 0 ){
        //   if ( Particles.n_send_x0 > 0 ) std::cout << " Sent X0:  " << Particles.n_send_x0 <<  "  " << procID <<  "  to  "  <<  dest[0] <<  std::endl;
          buffer_length = Particles.n_send_x0 * N_DATA_PER_PARTICLE_TRANSFER;
          Check_and_Grow_Particles_Buffer( &send_buffer_x0_particles , &buffer_length_particles_x0_send, buffer_length );
          Particles.Load_Particles_to_Buffer( 0, 0, send_buffer_x0_particles,  buffer_length_particles_x0_send );
          MPI_Isend(send_buffer_x0_particles, buffer_length, MPI_CHREAL, dest[0],   1, world, &send_request_particles_transfer[ireq_particles_transfer]);
          ireq_particles_transfer ++;
        // }
      }
      #endif

      if ( transfer_buffer ){
        //post non-blocking receive left x communication buffer
        MPI_Irecv(recv_buffer_x0, buffer_length, MPI_CHREAL, source[0], 0, world, &recv_request[ireq]);
        //non-blocking send left x communication buffer
        MPI_Isend(send_buffer_x0, buffer_length, MPI_CHREAL, dest[0],   1, world, &send_request[0]);
        //keep track of how many sends and receives are expected
        ireq++;
      }
    }

    if(flags[1]==5)
    {

      if (H.TRANSFER_HYDRO_BOUNDARIES) buffer_length = Load_Hydro_Buffer_X1();

      #ifdef GRAVITY
      if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ){
        buffer_length = Load_Potential_To_Buffer( 0, 1, send_buffer_x1, 0 );
      }
      #endif

      #ifdef PARTICLES
      if ( Particles.TRANSFER_DENSITY_BOUNDARIES) {
        buffer_length = Load_Particles_Density_Boundary_to_Buffer( 0, 1, send_buffer_x1  );
      }
      else if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
        transfer_buffer = false;
        MPI_Irecv(&Particles.n_recv_x1, 1, MPI_PART_INT, source[1], 1, world, &recv_request_n_particles[ireq_n_particles]);
        MPI_Isend(&Particles.n_send_x1, 1, MPI_PART_INT, dest[1],   0, world, &send_request_n_particles[1]);
        ireq_n_particles++;
        // if ( Particles.n_send_x1 > 0 ){
        //   std::cout << " Sent X1: " << Particles.n_send_x1 << std::endl;
          buffer_length = Particles.n_send_x1 * N_DATA_PER_PARTICLE_TRANSFER;
          Check_and_Grow_Particles_Buffer( &send_buffer_x1_particles , &buffer_length_particles_x1_send, buffer_length );
          Particles.Load_Particles_to_Buffer( 0, 1, send_buffer_x1_particles,  buffer_length_particles_x1_send  );
          MPI_Isend(send_buffer_x1_particles, buffer_length, MPI_CHREAL, dest[1],   0, world, &send_request_particles_transfer[ireq_particles_transfer]);\
          ireq_particles_transfer ++;
        // }
      }
      #endif

      if ( transfer_buffer ){
        //post non-blocking receive right x communication buffer
        MPI_Irecv(recv_buffer_x1, buffer_length, MPI_CHREAL, source[1], 1, world, &recv_request[ireq]);
        //non-blocking send right x communication buffer
        MPI_Isend(send_buffer_x1, buffer_length, MPI_CHREAL, dest[1],   0, world, &send_request[1]);
        //keep track of how many sends and receives are expected
        ireq++;
      }
    }

    // Receive the number of particles transfer
    #ifdef PARTICLES
    if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
      Wait_and_Recv_Particles_Transfer_BLOCK( dir, flags );
    }
    #endif

  }

  /* y boundaries */
  if (dir==1) {
    if(flags[2] == 5)
    {

      if ( H.TRANSFER_HYDRO_BOUNDARIES ) buffer_length = Load_Hydro_Buffer_Y0();

      #ifdef GRAVITY
      if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ){
        buffer_length = Load_Potential_To_Buffer( 1, 0, send_buffer_y0, 0 );
      }
      #endif

      #ifdef PARTICLES
      if ( Particles.TRANSFER_DENSITY_BOUNDARIES) {
       buffer_length = Load_Particles_Density_Boundary_to_Buffer( 1, 0, send_buffer_y0  );
      }
      else if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
        transfer_buffer = false;
        MPI_Isend(&Particles.n_send_y0, 1, MPI_PART_INT, dest[2],   3, world, &send_request_n_particles[0]);
        MPI_Irecv(&Particles.n_recv_y0, 1, MPI_PART_INT, source[2], 2, world, &recv_request_n_particles[ireq_n_particles]);
        ireq_n_particles++;
        // if ( Particles.n_send_y0 > 0 ){
        //   std::cout << " Sent Y0: " << Particles.n_send_y0 << std::endl;
          buffer_length = Particles.n_send_y0 * N_DATA_PER_PARTICLE_TRANSFER;
          Check_and_Grow_Particles_Buffer( &send_buffer_y0_particles , &buffer_length_particles_y0_send, buffer_length );
          Particles.Load_Particles_to_Buffer( 1, 0, send_buffer_y0_particles,  buffer_length_particles_y0_send  );
          MPI_Isend(send_buffer_y0_particles, buffer_length, MPI_CHREAL, dest[2],   3, world, &send_request_particles_transfer[ireq_particles_transfer]);
          ireq_particles_transfer ++;
        // }
      }
      #endif

      if ( transfer_buffer ){
        //post non-blocking receive left y communication buffer
        MPI_Irecv(recv_buffer_y0, buffer_length, MPI_CHREAL, source[2], 2, world, &recv_request[ireq]);
        //non-blocking send left y communication buffer
        MPI_Isend(send_buffer_y0, buffer_length, MPI_CHREAL, dest[2],   3, world, &send_request[0]);
        //keep track of how many sends and receives are expected
        ireq++;
      }
    }

    if(flags[3]==5)
    {

      if ( H.TRANSFER_HYDRO_BOUNDARIES ) buffer_length = Load_Hydro_Buffer_Y1();

      #ifdef GRAVITY
      if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ){
        buffer_length = Load_Potential_To_Buffer( 1, 1, send_buffer_y1, 0 );
      }
      #endif

      #ifdef PARTICLES
      if ( Particles.TRANSFER_DENSITY_BOUNDARIES) {
        buffer_length = Load_Particles_Density_Boundary_to_Buffer( 1, 1, send_buffer_y1  );
      }
      else if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
        transfer_buffer = false;
        MPI_Isend(&Particles.n_send_y1, 1, MPI_PART_INT, dest[3],   2, world, &send_request_n_particles[1]);
        MPI_Irecv(&Particles.n_recv_y1, 1, MPI_PART_INT, source[3], 3, world, &recv_request_n_particles[ireq_n_particles]);
        ireq_n_particles++;
        // if ( Particles.n_send_y1 > 0 ){
        //   std::cout << " Sent Y1: " << Particles.n_send_y1 << std::endl;
          buffer_length = Particles.n_send_y1 * N_DATA_PER_PARTICLE_TRANSFER;
          Check_and_Grow_Particles_Buffer( &send_buffer_y1_particles , &buffer_length_particles_y1_send, buffer_length );
          Particles.Load_Particles_to_Buffer( 1, 1, send_buffer_y1_particles,  buffer_length_particles_y1_send  );
          MPI_Isend(send_buffer_y1_particles, buffer_length, MPI_CHREAL, dest[3],   2, world, &send_request_particles_transfer[ireq_particles_transfer]);
          ireq_particles_transfer ++;
        // }
      }
      #endif

      if ( transfer_buffer ){
        //post non-blocking receive right y communication buffer
        MPI_Irecv(recv_buffer_y1, buffer_length, MPI_CHREAL, source[3], 3, world, &recv_request[ireq]);
        //non-blocking send right y communication buffer
        MPI_Isend(send_buffer_y1, buffer_length, MPI_CHREAL, dest[3],   2, world, &send_request[1]);
        //keep track of how many sends and receives are expected
        ireq++;
      }
    }

    // Receive the number of particles transfer
    #ifdef PARTICLES
    if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
      Wait_and_Recv_Particles_Transfer_BLOCK( dir, flags );
    }
    #endif
  }

  /* z boundaries */
  if (dir==2) {

    if(flags[4]==5)
    {
      if ( H.TRANSFER_HYDRO_BOUNDARIES ) buffer_length = Load_Hydro_Buffer_Z0();

      #ifdef GRAVITY
      if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ){
        buffer_length = Load_Potential_To_Buffer( 2, 0, send_buffer_z0, 0 );
      }
      #endif

      #ifdef PARTICLES
      if ( Particles.TRANSFER_DENSITY_BOUNDARIES) {
        buffer_length = Load_Particles_Density_Boundary_to_Buffer( 2, 0, send_buffer_z0  );
      }
      else if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
        transfer_buffer = false;
        MPI_Isend(&Particles.n_send_z0, 1, MPI_PART_INT, dest[4],   5, world, &send_request_n_particles[0]);
        MPI_Irecv(&Particles.n_recv_z0, 1, MPI_PART_INT, source[4], 4, world, &recv_request_n_particles[ireq_n_particles]);
        ireq_n_particles++;
        // if ( Particles.n_send_z0 > 0 ){
        //   std::cout << " Sent Z0: " << Particles.n_send_z0 << std::endl;
          buffer_length = Particles.n_send_z0 * N_DATA_PER_PARTICLE_TRANSFER;
          Check_and_Grow_Particles_Buffer( &send_buffer_z0_particles , &buffer_length_particles_z0_send, buffer_length );
          Particles.Load_Particles_to_Buffer( 2, 0, send_buffer_z0_particles,  buffer_length_particles_z0_send  );
          MPI_Isend(send_buffer_z0_particles, buffer_length, MPI_CHREAL, dest[4],   5, world, &send_request_particles_transfer[ireq_particles_transfer]);
          ireq_particles_transfer ++;
        // }
      }
      #endif

      if ( transfer_buffer ){
        //post non-blocking receive left z communication buffer
        MPI_Irecv(recv_buffer_z0, buffer_length, MPI_CHREAL, source[4], 4, world, &recv_request[ireq]);
        //non-blocking send left z communication buffer
        MPI_Isend(send_buffer_z0, buffer_length, MPI_CHREAL, dest[4],   5, world, &send_request[0]);
        //keep track of how many sends and receives are expected
        ireq++;
      }
    }

    if(flags[5]==5)
    {

      if ( H.TRANSFER_HYDRO_BOUNDARIES ) buffer_length = Load_Hydro_Buffer_Z1();

      #ifdef GRAVITY
      if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ){
        buffer_length = Load_Potential_To_Buffer( 2, 1, send_buffer_z1, 0 );
      }
      #endif

      #ifdef PARTICLES
      if ( Particles.TRANSFER_DENSITY_BOUNDARIES) {
        buffer_length = Load_Particles_Density_Boundary_to_Buffer( 2, 1, send_buffer_z1  );
      }
      else if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
        transfer_buffer = false;
        MPI_Isend(&Particles.n_send_z1, 1, MPI_CHREAL, dest[5],   4, world, &send_request_n_particles[1]);
        MPI_Irecv(&Particles.n_recv_z1, 1, MPI_CHREAL, source[5], 5, world, &recv_request_n_particles[ireq_n_particles]);
        ireq_n_particles++;
        // if ( Particles.n_send_z1 > 0 ){
        //   std::cout << " Sent Z1: " << Particles.n_send_z1 << std::endl;
          buffer_length = Particles.n_send_z1 * N_DATA_PER_PARTICLE_TRANSFER;
          Check_and_Grow_Particles_Buffer( &send_buffer_z1_particles , &buffer_length_particles_z1_send, buffer_length );
          Particles.Load_Particles_to_Buffer( 2, 1, send_buffer_z1_particles,  buffer_length_particles_z1_send  );
          MPI_Isend(send_buffer_z1_particles, buffer_length, MPI_CHREAL, dest[5],   4, world, &send_request_particles_transfer[ireq_particles_transfer]);
          ireq_particles_transfer ++;
        // }
      }
      #endif

      if ( transfer_buffer ){
        // if (H.TRANSFER_HYDRO_BOUNDARIES ) std::cout << " N Loaded Z1: " << send_buffer_z1[z_buffer_length_hydro]  << std::endl;
        //post non-blocking receive right x communication buffer
        MPI_Irecv(recv_buffer_z1, buffer_length, MPI_CHREAL, source[5], 5, world, &recv_request[ireq]);
        //non-blocking send right x communication buffer
        MPI_Isend(send_buffer_z1, buffer_length, MPI_CHREAL, dest[5],   4, world, &send_request[1]);
        //keep track of how many sends and receives are expected
        ireq++;
      }
    }
    // Receive the number of particles transfer
    #ifdef PARTICLES
    if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ){
      Wait_and_Recv_Particles_Transfer_BLOCK( dir, flags );
    }
    #endif
  }

}


void Grid3D::Wait_and_Unload_MPI_Comm_Buffers_SLAB(int *flags)
{
  int iwait;
  int index = 0;
  int wait_max=0;
  MPI_Status status;

  //find out how many recvs we need to wait for
  for(iwait=0;iwait<6;iwait++)
    if(flags[iwait] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm

  //wait for any receives to complete
  for(iwait=0;iwait<wait_max;iwait++)
  {
    //wait for recv completion
    MPI_Waitany(wait_max,recv_request,&index,&status);

    //depending on which face arrived, load the buffer into the ghost grid
    Unload_MPI_Comm_Buffers(status.MPI_TAG);
  }

}


void Grid3D::Wait_and_Unload_MPI_Comm_Buffers_BLOCK(int dir, int *flags)
{
  #ifdef PARTICLES
  if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ) return;
  #endif

  int iwait;
  int index = 0;
  int wait_max=0;
  MPI_Status status;


  //find out how many recvs we need to wait for
  if (dir==0) {
    if(flags[0] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[1] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }
  if (dir==1) {
    if(flags[2] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[3] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }
  if (dir==2) {
    if(flags[4] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[5] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }

  //wait for any receives to complete
  for(iwait=0;iwait<wait_max;iwait++)
  {
    //wait for recv completion
    MPI_Waitany(wait_max,recv_request,&index,&status);
    //if (procID==1) MPI_Get_count(&status, MPI_CHREAL, &count);
    //if (procID==1) printf("Process 1 unloading direction %d, source %d, index %d, length %d.\n", status.MPI_TAG, status.MPI_SOURCE, index, count);
    //depending on which face arrived, load the buffer into the ghost grid
    Unload_MPI_Comm_Buffers(status.MPI_TAG);
  }


}

#ifdef PARTICLES

void Grid3D::Wait_and_Recv_Particles_Transfer_BLOCK(int dir, int *flags)
{
  #ifdef PARTICLES
  if ( !Particles.TRANSFER_PARTICLES_BOUNDARIES ) return;
  #endif

  int iwait;
  int index = 0;
  int wait_max=0;
  MPI_Status status;


  //find out how many recvs we need to wait for
  if (dir==0) {
    if(flags[0] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[1] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }
  if (dir==1) {
    if(flags[2] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[3] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }
  if (dir==2) {
    if(flags[4] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[5] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }

  int ireq_particles_transfer = 0;

  //wait for any receives to complete
  for(iwait=0;iwait<wait_max;iwait++)
  {
    //wait for recv completion
    MPI_Waitany(wait_max,recv_request_n_particles,&index,&status);
    //if (procID==1) MPI_Get_count(&status, MPI_CHREAL, &count);
    //if (procID==1) printf("Process 1 unloading direction %d, source %d, index %d, length %d.\n", status.MPI_TAG, status.MPI_SOURCE, index, count);
    //depending on which face arrived, load the buffer into the ghost grid
    // Unload_MPI_Comm_Buffers(status.MPI_TAG);
    Load_N_Particles_Transfer(status.MPI_TAG, &ireq_particles_transfer);
  }

}

void Grid3D::Wait_and_Unload_MPI_Comm_Particles_Buffers_BLOCK(int dir, int *flags)
{

  int iwait;
  int index = 0;
  int wait_max=0;
  MPI_Status status;


  //find out how many recvs we need to wait for
  if (dir==0) {
    if(flags[0] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[1] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }
  if (dir==1) {
    if(flags[2] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[3] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }
  if (dir==2) {
    if(flags[4] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
    if(flags[5] == 5) //there is communication on this face
      wait_max++;   //so we'll need to wait for its comm
  }

  //wait for any receives to complete
  for(iwait=0;iwait<wait_max;iwait++)
  {
    //wait for recv completion
    MPI_Waitany(wait_max,recv_request_particles_transfer,&index,&status);
    //if (procID==1) MPI_Get_count(&status, MPI_CHREAL, &count);
    //if (procID==1) printf("Process 1 unloading direction %d, source %d, index %d, length %d.\n", status.MPI_TAG, status.MPI_SOURCE, index, count);
    //depending on which face arrived, load the buffer into the ghost grid
    // Unload_MPI_Comm_Buffers(status.MPI_TAG);
    // Load_N_Particles_Transfer(status.MPI_TAG, &ireq_particles_transfer);
    // std::cout << " Unloading Particles :  " << status.MPI_TAG <<  std::endl;
    Unload_Particles_From_Buffers_BLOCK(status.MPI_TAG);
  }

}

void Grid3D::Load_N_Particles_Transfer(int index, int *ireq_particles_transfer){

  int buffer_length;

  // std::cout << "ireq: " << *ireq_particles_transfer << std::endl;

  if ( index == 0){
    buffer_length = Particles.n_recv_x0 * N_DATA_PER_PARTICLE_TRANSFER;
    Check_and_Grow_Particles_Buffer( &recv_buffer_x0_particles , &buffer_length_particles_x0_recv, buffer_length );
    // if ( Particles.n_recv_x0 ){
      // std::cout << " Recv X0: " << Particles.n_recv_x0 << std::endl;
      MPI_Irecv(recv_buffer_x0_particles, buffer_length, MPI_CHREAL, source[0], 0, world, &recv_request_particles_transfer[*ireq_particles_transfer]);
    // }
  }

  if ( index == 1){
    buffer_length = Particles.n_recv_x1 * N_DATA_PER_PARTICLE_TRANSFER;
    Check_and_Grow_Particles_Buffer( &recv_buffer_x1_particles , &buffer_length_particles_x1_recv, buffer_length );
    // if ( Particles.n_recv_x1 ){
      // if ( Particles.n_recv_x1 > 0 ) std::cout << " Recv X1:  " << Particles.n_recv_x1 <<  "  " << procID <<  "  from "  <<  source[1] <<  std::endl;
      MPI_Irecv(recv_buffer_x1_particles, buffer_length, MPI_CHREAL, source[1], 1, world, &recv_request_particles_transfer[*ireq_particles_transfer]);
    // }
  }

  if ( index == 2){
    buffer_length = Particles.n_recv_y0 * N_DATA_PER_PARTICLE_TRANSFER;
    Check_and_Grow_Particles_Buffer( &recv_buffer_y0_particles , &buffer_length_particles_y0_recv, buffer_length );
    // if ( Particles.n_recv_y0 ){
      // std::cout << " Recv Y0: " << Particles.n_recv_y0 << std::endl;
      MPI_Irecv(recv_buffer_y0_particles, buffer_length, MPI_CHREAL, source[2], 2, world, &recv_request_particles_transfer[*ireq_particles_transfer]);
    // }
  }

  if ( index == 3){
    buffer_length = Particles.n_recv_y1 * N_DATA_PER_PARTICLE_TRANSFER;
    Check_and_Grow_Particles_Buffer( &recv_buffer_y1_particles , &buffer_length_particles_y1_recv, buffer_length );
    // if ( Particles.n_recv_y1 ){
      // std::cout << " Recv Y1: " << Particles.n_recv_y1 << std::endl;
      MPI_Irecv(recv_buffer_y1_particles, buffer_length, MPI_CHREAL, source[3], 3, world, &recv_request_particles_transfer[*ireq_particles_transfer]);
    // }
  }

  if ( index == 4){
    buffer_length = Particles.n_recv_z0 * N_DATA_PER_PARTICLE_TRANSFER;
    Check_and_Grow_Particles_Buffer( &recv_buffer_z0_particles , &buffer_length_particles_z0_recv, buffer_length );
    // if ( Particles.n_recv_z0 ){
      // std::cout << " Recv Z0: " << Particles.n_recv_z0 << std::endl;
      MPI_Irecv(recv_buffer_z0_particles, buffer_length, MPI_CHREAL, source[4], 4, world, &recv_request_particles_transfer[*ireq_particles_transfer]);
    // }
  }

  if ( index == 5){
    buffer_length = Particles.n_recv_z1 * N_DATA_PER_PARTICLE_TRANSFER;
    Check_and_Grow_Particles_Buffer( &recv_buffer_z1_particles , &buffer_length_particles_z1_recv, buffer_length );
    // if ( Particles.n_recv_z1 ){
      // std::cout << " Recv Z1: " << Particles.n_recv_z1 << std::endl;
      MPI_Irecv(recv_buffer_z1_particles, buffer_length, MPI_CHREAL, source[5], 5, world, &recv_request_particles_transfer[*ireq_particles_transfer]);
    // }
  }

  *ireq_particles_transfer += 1;


}

#endif



void Grid3D::Unload_MPI_Comm_Buffers(int index)
{
  switch(flag_decomp)
  {
    case SLAB_DECOMP:
      Unload_MPI_Comm_Buffers_SLAB(index);
      break;
    case BLOCK_DECOMP:
      Unload_MPI_Comm_Buffers_BLOCK(index);
      // #ifdef PARTICLES
      // Unload_Particles_From_Buffers_BLOCK(index);
      // #endif
      break;
  }
}


void Grid3D::Unload_MPI_Comm_Buffers_SLAB(int index)
{
  int i, j, k, ii;
  int idx;
  int gidx;
  int offset = H.n_ghost*H.ny*H.nz;

  //left face
  if(index==0)
  {
    //load left x communication buffer
    for(i=0;i<H.n_ghost;i++) {
      for(j=0;j<H.ny;j++) {
        for(k=0;k<H.nz;k++)
        {
          idx  = i + j*H.nx + k*H.nx*H.ny;
          gidx = i + j*H.n_ghost + k*H.n_ghost*H.ny;
          for (ii=0; ii<H.n_fields; ii++) {
            C.density[idx + ii*H.n_cells] = *(recv_buffer_0 + gidx + ii*offset);
          }
        }
      }
    }
  }

  //right face
  if(index==1)
  {
    //load left x communication buffer
    for(i=0;i<H.n_ghost;i++) {
      for(j=0;j<H.ny;j++) {
        for(k=0;k<H.nz;k++)
        {
          idx  = (i + H.nx - H.n_ghost) + j*H.nx + k*H.nx*H.ny;
          gidx = i                      + j*H.n_ghost + k*H.n_ghost*H.ny;
          for (ii=0; ii<H.n_fields; ii++) {
            C.density[idx + ii*H.n_cells] = *(recv_buffer_1 + gidx + ii*offset);
          }
        }
      }
    }
  }

  //done unloading
}


void Grid3D::Unload_MPI_Comm_Buffers_BLOCK(int index)
{
  int i, j, k, ii;
  int idx;
  int gidx;
  int offset;

  // bool transfer_hydro = H.TRANSFER_HYDRO_BOUNDARIES;

  #ifdef PARTICLES
  // if ( Particles.TRANSFER_DENSITY_BOUNDARIES ) H.TRANSFER_HYDRO_BOUNDARIES = false;
  // if ( Particles.TRANSFER_PARTICLES_BOUNDARIES ) H.TRANSFER_HYDRO_BOUNDARIES = false;
  #endif
  #ifdef GRAVITY
  // if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ) H.TRANSFER_HYDRO_BOUNDARIES = false;
  #endif

  if ( H.TRANSFER_HYDRO_BOUNDARIES ){
    //unload left x communication buffer
    if(index==0)
    {
      // 1D
      if (H.ny == 1 && H.nz == 1) {
        offset = H.n_ghost;
        for(i=0;i<H.n_ghost;i++) {
          idx  = i;
          gidx = i;
          for (ii=0; ii<H.n_fields; ii++) {
            C.density[idx + H.n_cells] = *(recv_buffer_x0 + gidx + ii*offset);
          }
        }
      }
      // 2D
      if (H.ny > 1 && H.nz == 1) {
        offset = H.n_ghost*(H.ny-2*H.n_ghost);
        for(i=0;i<H.n_ghost;i++) {
          for (j=0;j<H.ny-2*H.n_ghost;j++) {
            idx  = i + (j+H.n_ghost)*H.nx;
            gidx = i + j*H.n_ghost;
            for (ii=0; ii<H.n_fields; ii++) {
              C.density[idx + ii*H.n_cells] = *(recv_buffer_x0 + gidx + ii*offset);
            }
          }
        }
      }
      // 3D
      if (H.nz > 1) {
        offset = H.n_ghost*(H.ny-2*H.n_ghost)*(H.nz-2*H.n_ghost);
        for(i=0;i<H.n_ghost;i++) {
          for(j=0;j<H.ny-2*H.n_ghost;j++) {
            for(k=0;k<H.nz-2*H.n_ghost;k++) {
              idx  = i + (j+H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
              gidx = i + j*H.n_ghost + k*H.n_ghost*(H.ny-2*H.n_ghost);
              for (ii=0; ii<H.n_fields; ii++) {
                C.density[idx + ii*H.n_cells] = *(recv_buffer_x0 + gidx + ii*offset);
              }
            }
          }
        }
      }
      // std::cout << "N_recv X_0: " << recv_buffer_x0[x_buffer_length_hydro] << std::endl;
    }

    //unload right x communication buffer
    if(index==1)
    {
      // 1D
      if (H.ny == 1 && H.nz == 1) {
        offset = H.n_ghost;
        for(i=0;i<H.n_ghost;i++) {
          idx  = i+H.nx-H.n_ghost;
          gidx = i;
          for (ii=0; ii<H.n_fields; ii++) {
            C.density[idx + ii*H.n_cells] = *(recv_buffer_x1 + gidx + ii*offset);
          }
        }
      }
      // 2D
      if (H.ny > 1 && H.nz == 1) {
        offset = H.n_ghost*(H.ny-2*H.n_ghost);
        for(i=0;i<H.n_ghost;i++) {
          for (j=0;j<H.ny-2*H.n_ghost;j++) {
            idx  = i+H.nx-H.n_ghost + (j+H.n_ghost)*H.nx;
            gidx = i + j*H.n_ghost;
            for (ii=0; ii<H.n_fields; ii++) {
              C.density[idx + ii*H.n_cells] = *(recv_buffer_x1 + gidx + ii*offset);
            }
          }
        }
      }
      // 3D
      if (H.nz > 1) {
        offset = H.n_ghost*(H.ny-2*H.n_ghost)*(H.nz-2*H.n_ghost);
        for(i=0;i<H.n_ghost;i++) {
          for(j=0;j<H.ny-2*H.n_ghost;j++) {
            for(k=0;k<H.nz-2*H.n_ghost;k++) {
              idx  = i+H.nx-H.n_ghost + (j+H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
              gidx = i + j*H.n_ghost + k*H.n_ghost*(H.ny-2*H.n_ghost);
              for (ii=0; ii<H.n_fields; ii++) {
                C.density[idx + ii*H.n_cells] = *(recv_buffer_x1 + gidx + ii*offset);
              }
            }
          }
        }
      }
      // std::cout << "N_recv X_1: " << recv_buffer_x1[x_buffer_length_hydro] << std::endl;
    }


    //unload left y communication buffer
    if(index==2)
    {
      // 2D
      if (H.nz == 1) {
        offset = H.n_ghost*H.nx;
        for(i=0;i<H.nx;i++) {
          for (j=0;j<H.n_ghost;j++) {
            idx  = i + j*H.nx;
            gidx = i + j*H.nx;
            for (ii=0; ii<H.n_fields; ii++) {
              C.density[idx + ii*H.n_cells] = *(recv_buffer_y0 + gidx + ii*offset);
            }
          }
        }
      }
      // 3D
      if (H.nz > 1) {
        offset = H.n_ghost*H.nx*(H.nz-2*H.n_ghost);
        for(i=0;i<H.nx;i++) {
          for(j=0;j<H.n_ghost;j++) {
            for(k=0;k<H.nz-2*H.n_ghost;k++) {
              idx  = i + j*H.nx + (k+H.n_ghost)*H.nx*H.ny;
              gidx = i + j*H.nx + k*H.nx*H.n_ghost;
              for (ii=0; ii<H.n_fields; ii++) {
                C.density[idx + ii*H.n_cells] = *(recv_buffer_y0 + gidx + ii*offset);
              }
            }
          }
        }
      }
      // std::cout << "N_recv Y_0: " << recv_buffer_y0[y_buffer_length_hydro] << std::endl;
    }

    //unload right y communication buffer
    if(index==3)
    {
      // 2D
      if (H.nz == 1) {
        offset = H.n_ghost*H.nx;
        for(i=0;i<H.nx;i++) {
          for (j=0;j<H.n_ghost;j++) {
            idx  = i + (j+H.ny-H.n_ghost)*H.nx;
            gidx = i + j*H.nx;
            for (ii=0; ii<H.n_fields; ii++) {
              C.density[idx + ii*H.n_cells] = *(recv_buffer_y1 + gidx + ii*offset);
            }
          }
        }
      }
      // 3D
      if (H.nz > 1) {
        offset = H.n_ghost*H.nx*(H.nz-2*H.n_ghost);
        for(i=0;i<H.nx;i++) {
          for(j=0;j<H.n_ghost;j++) {
            for(k=0;k<H.nz-2*H.n_ghost;k++) {
              idx  = i + (j+H.ny-H.n_ghost)*H.nx + (k+H.n_ghost)*H.nx*H.ny;
              gidx = i + j*H.nx + k*H.nx*H.n_ghost;
              for (ii=0; ii<H.n_fields; ii++) {
                C.density[idx + ii*H.n_cells] = *(recv_buffer_y1 + gidx + ii*offset);
              }
            }
          }
        }
      }
      // std::cout << "N_recv Y_1: " << recv_buffer_y1[y_buffer_length_hydro] << std::endl;
    }

    //unload left z communication buffer
    if(index==4)
    {
      offset = H.n_ghost*H.nx*H.ny;
      for(i=0;i<H.nx;i++) {
        for(j=0;j<H.ny;j++) {
          for(k=0;k<H.n_ghost;k++) {
            idx  = i + j*H.nx + k*H.nx*H.ny;
            gidx = i + j*H.nx + k*H.nx*H.ny;
            for (ii=0; ii<H.n_fields; ii++) {
              C.density[idx + ii*H.n_cells] = *(recv_buffer_z0 + gidx + ii*offset);
            }
          }
        }
      }
      // std::cout << "N_recv Z_0: " << recv_buffer_z0[z_buffer_length_hydro] << std::endl;
    }

    //unload right z communication buffer
    if(index==5)
    {
      offset = H.n_ghost*H.nx*H.ny;
      for(i=0;i<H.nx;i++) {
        for(j=0;j<H.ny;j++) {
          for(k=0;k<H.n_ghost;k++) {
            idx  = i + j*H.nx + (k+H.nz-H.n_ghost)*H.nx*H.ny;
            gidx = i + j*H.nx + k*H.nx*H.ny;
            for (ii=0; ii<H.n_fields; ii++) {
              C.density[idx + ii*H.n_cells] = *(recv_buffer_z1 + gidx + ii*offset);
            }
          }
        }
      }
      // std::cout << "N_recv Z_0: " << recv_buffer_z1[z_buffer_length_hydro] << std::endl;
    }
  }

  #ifdef PARTICLES
  if (  Particles.TRANSFER_DENSITY_BOUNDARIES ){
    if ( index == 0 ) Unload_Particles_Density_Boundary_From_Buffer( 0, 0, recv_buffer_x0 );
    if ( index == 1 ) Unload_Particles_Density_Boundary_From_Buffer( 0, 1, recv_buffer_x1 );
    if ( index == 2 ) Unload_Particles_Density_Boundary_From_Buffer( 1, 0, recv_buffer_y0 );
    if ( index == 3 ) Unload_Particles_Density_Boundary_From_Buffer( 1, 1, recv_buffer_y1 );
    if ( index == 4 ) Unload_Particles_Density_Boundary_From_Buffer( 2, 0, recv_buffer_z0 );
    if ( index == 5 ) Unload_Particles_Density_Boundary_From_Buffer( 2, 1, recv_buffer_z1 );
  }
  #endif

  #ifdef GRAVITY
  if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ){
    if ( index == 0 ) Unload_Potential_from_Buffer( 0, 0, recv_buffer_x0, 0  );
    if ( index == 1 ) Unload_Potential_from_Buffer( 0, 1, recv_buffer_x1, 0  );
    if ( index == 2 ) Unload_Potential_from_Buffer( 1, 0, recv_buffer_y0, 0  );
    if ( index == 3 ) Unload_Potential_from_Buffer( 1, 1, recv_buffer_y1, 0  );
    if ( index == 4 ) Unload_Potential_from_Buffer( 2, 0, recv_buffer_z0, 0  );
    if ( index == 5 ) Unload_Potential_from_Buffer( 2, 1, recv_buffer_z1, 0  );
  }
  #endif

}


#ifdef PARTICLES
void Grid3D::Unload_Particles_From_Buffers_BLOCK(int index){

  if ( Particles.TRANSFER_DENSITY_BOUNDARIES ) return;
  if ( H.TRANSFER_HYDRO_BOUNDARIES ) return;
  #ifdef GRAVITY
  if ( Grav.TRANSFER_POTENTIAL_BOUNDARIES ) return;
  #endif

  if( index == 0) {
    Unload_Particles_from_Buffer_X0(  );
  }

  if( index == 1){
    Unload_Particles_from_Buffer_X1();
  }

  if( index == 2){
    Unload_Particles_from_Buffer_Y0();
  }

  if( index == 3){
    Unload_Particles_from_Buffer_Y1();
  }

  if( index == 4){
    Unload_Particles_from_Buffer_Z0();
  }

  if( index == 5){
    Unload_Particles_from_Buffer_Z1();
  }

}



#endif


#endif /*MPI_CHOLLA*/

// #ifdef PARTICLES
//
// #include"parts_boundaries.h"
// #include <unistd.h>
//
//
// int real_to_int( Real inVal ){
//   int outVal = (int) inVal;
//   if ( (inVal - outVal) > 0.1 ) outVal += 1;
//   if ( fabs(outVal - inVal) > 0.5 ) outVal -= 1;
//   return outVal;
// }
//
// Real Get_and_Remove_Real( partID_t indx, ch_vector_t &vec ){
//   Real value = vec[indx];
//   vec[indx] = vec.back();
//   vec.pop_back();
//   return value;
// }
//
// Real Get_and_Remove_partId( partID_t indx, ids_vector_t &vec ){
//   Real value = (Real) vec[indx];
//   vec[indx] = vec.back();
//   vec.pop_back();
//   return value;
// }
//
// void Remove_id( partID_t indx, ids_vector_t &vec ){
//   vec[indx] = vec.back();
//   vec.pop_back();
// }
//
// void Remove_Real( partID_t indx, ch_vector_t &vec ){
//   vec[indx] = vec.back();
//   vec.pop_back();
// }
//
// void Part3D::Select_Particles_to_Transfer( void ){
//
//   partID_t pIndx;
//
//   for ( pIndx=0; pIndx<n_local; pIndx++ ){
//     // if (procID == 0) std::cout << pIndx << std::endl;
//     if ( pos_x[pIndx] <= G.xMin ){
//       out_indxs_vec_x0.push_back( pIndx );
//       continue;
//     }
//     if ( pos_x[pIndx] >= G.xMax ){
//       out_indxs_vec_x1.push_back( pIndx );
//       continue;
//     }
//     if ( pos_y[pIndx] <= G.yMin ){
//       out_indxs_vec_y0.push_back( pIndx );
//       continue;
//     }
//     if ( pos_y[pIndx] >= G.yMax ){
//       out_indxs_vec_y1.push_back( pIndx );
//       continue;
//     }
//     if ( pos_z[pIndx] <= G.zMin ){
//       out_indxs_vec_z0.push_back( pIndx );
//       continue;
//     }
//     if ( pos_z[pIndx] >= G.zMax ){
//       out_indxs_vec_z1.push_back( pIndx );
//       continue;
//     }
//   }
//
// /*
//   if(out_indxs_vec_x0.size()>0)
//   {
//     printf("procID %d xMin %e xMax %e xpart %e\n",procID,G.xMin,G.xMax,pos_x[out_indxs_vec_x0[0]]);
//     printf("procID %d Selected : x0 %ld x1 %ld y0 %ld y1 %ld z0 %ld z1 %ld\n",procID,out_indxs_vec_x0.size(),out_indxs_vec_x1.size(),out_indxs_vec_y0.size(),out_indxs_vec_y1.size(),out_indxs_vec_z0.size(),out_indxs_vec_z1.size());
//     fflush(stdout);
//   }
// */
// }
//
// void Part3D::Load_Particles_to_Buffer_X( void ){
//   Load_Particles_to_Buffer( 0, 0, out_indxs_vec_x0 );
//   Load_Particles_to_Buffer( 0, 1, out_indxs_vec_x1 );
//   Load_Particles_to_Buffer( 1, 0, out_indxs_vec_y0 );
//   Load_Particles_to_Buffer( 1, 1, out_indxs_vec_y1 );
//   Load_Particles_to_Buffer( 2, 0, out_indxs_vec_z0 );
//   Load_Particles_to_Buffer( 2, 1, out_indxs_vec_z1 );
// }
//
// void Part3D::Load_Particles_to_Buffer_Y( void ){
//   // Load_Particles_to_Buffer( 1, 0, out_indxs_vec_y0 );
//   // Load_Particles_to_Buffer( 1, 1, out_indxs_vec_y1 );
//   Unload_Particles_from_Buffer( 0, 0);
//   Unload_Particles_from_Buffer( 0, 1);
// }
//
// void Part3D::Load_Particles_to_Buffer_Z( void ){
//   // Load_Particles_to_Buffer( 2, 0, out_indxs_vec_z0 );
//   // Load_Particles_to_Buffer( 2, 1, out_indxs_vec_z1 );
//   Unload_Particles_from_Buffer( 1, 0);
//   Unload_Particles_from_Buffer( 1, 1);
// }
//
// void Part3D::Finish_Unload_from_Transfers( void ){
//   Unload_Particles_from_Buffer( 2, 0);
//   Unload_Particles_from_Buffer( 2, 1);
//
// }
//
//
// void Add_Particle_To_Buffer( Real *buffer, Real pId, Real pMass,
//                             Real pPos_x, Real pPos_y, Real pPos_z,
//                             Real pVel_x, Real pVel_y, Real pVel_z ){
//   int n_in_buffer, offset_out;
//   n_in_buffer = real_to_int( buffer[0] );
//   offset_out = 2 + n_in_buffer*N_DATA_PER_PARTICLE;
//   buffer[offset_out + 0] = pId;
//   buffer[offset_out + 1] = pMass;
//   buffer[offset_out + 2] = pPos_x;
//   buffer[offset_out + 3] = pPos_y;
//   buffer[offset_out + 4] = pPos_z;
//   buffer[offset_out + 5] = pVel_x;
//   buffer[offset_out + 6] = pVel_y;
//   buffer[offset_out + 7] = pVel_z;
//   buffer[0] += 1;
// }
//
// void Part3D::Add_Particle_To_Vectors( Real pId, Real pMass,
//                             Real pPos_x, Real pPos_y, Real pPos_z,
//                             Real pVel_x, Real pVel_y, Real pVel_z ){
//   partIDs.push_back( (partID_t) pId );
//   mass.push_back( pMass );
//   pos_x.push_back( pPos_x );
//   pos_y.push_back( pPos_y );
//   pos_z.push_back( pPos_z );
//   vel_x.push_back( pVel_x );
//   vel_y.push_back( pVel_y );
//   vel_z.push_back( pVel_z );
//   grav_x.push_back(0);
//   grav_y.push_back(0);
//   grav_z.push_back(0);
//
//   n_local += 1;
// }
//
//
// void Part3D::Unload_Particles_from_Buffer( int direction, int side ){
//
//   int n_recv;
//   Real *recv_buffer;
//   if ( direction == 0 ){
//     if ( side == 0 ){
//       recv_buffer = recv_buffer_x0_particles;
//     }
//     if ( side == 1 ){
//       recv_buffer = recv_buffer_x1_particles;
//     }
//   }
//   if ( direction == 1 ){
//     if ( side == 0 ){
//       recv_buffer = recv_buffer_y0_particles;
//     }
//     if ( side == 1 ){
//       recv_buffer = recv_buffer_y1_particles;
//     }
//   }
//   if ( direction == 2 ){
//     if ( side == 0 ){
//       recv_buffer = recv_buffer_z0_particles;
//     }
//     if ( side == 1 ){
//       recv_buffer = recv_buffer_z1_particles;
//     }
//   }
//
//   int offset_buff, indx, offset_out, n_in_out;
//   partID_t pId;
//   Real pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z;
//   n_recv = real_to_int( recv_buffer[0] );
//
//   // fflush(stdout);
//   // MPI_Barrier(world);
//   // for(int i=0;i<nproc;i++)
//   // {
//   //   if(i==procID)
//   //   {
//   //     printf("  [procID %d] Receiving: ( %d %d ) %d \n",procID, direction, side, n_recv );
//   //   }
//   //   fflush(stdout);
//   //   usleep(1e2);
//   //   MPI_Barrier(world);
//   // }
//
//   offset_buff = 2;
//   for ( indx = 0; indx<n_recv; indx++ ){
//     pId    = recv_buffer[ offset_buff + 0 ];
//     pMass  = recv_buffer[ offset_buff + 1 ];
//     pPos_x = recv_buffer[ offset_buff + 2 ];
//     pPos_y = recv_buffer[ offset_buff + 3 ];
//     pPos_z = recv_buffer[ offset_buff + 4 ];
//     pVel_x = recv_buffer[ offset_buff + 5 ];
//     pVel_y = recv_buffer[ offset_buff + 6 ];
//     pVel_z = recv_buffer[ offset_buff + 7 ];
//     offset_buff += N_DATA_PER_PARTICLE;
//     if ( pPos_x < G.domainMin_x ) pPos_x += ( G.domainMax_x - G.domainMin_x );
//     if ( pPos_x > G.domainMax_x ) pPos_x -= ( G.domainMax_x - G.domainMin_x );
//     if ( ( pPos_x < G.xMin ) || ( pPos_x > G.xMax )  ){
//       std::cout << "ERROR Particle out of X domain" << std::endl;
//       std::cout << " posX: " << pPos_x << "velX: " << pVel_x << std::endl;
//       std::cout << " posY: " << pPos_y << "velY: " << pVel_y << std::endl;
//
//       continue;
//     }
//
//     if (  direction == 1 ){
//       if ( pPos_y < G.domainMin_y ) pPos_y += ( G.domainMax_y - G.domainMin_y );
//       if ( pPos_y > G.domainMax_y ) pPos_y -= ( G.domainMax_y - G.domainMin_y );
//     }
//     if ( (direction==1 || direction==2) && (( pPos_y < G.yMin ) || ( pPos_y > G.yMax ))  ){
//       std::cout << "ERROR Particle out of Y domain" << std::endl;
//       std::cout << " posX: " << pPos_x << "velX: " << pVel_x << std::endl;
//       std::cout << " posY: " << pPos_y << "velY: " << pVel_y << std::endl;
//       std::cout << " posZ: " << pPos_z << "velZ: " << pVel_z << std::endl;
//       continue;
//     }
//
//     if (direction  !=1 ){
//       if ( pPos_y <= G.yMin ){
//         // std::cout << "  Adding to Y0 " << std::endl;
//         Add_Particle_To_Buffer( send_buffer_y0_particles, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
//         continue;
//       }
//       if ( pPos_y >= G.yMax ){
//         // std::cout << "  Adding to Y1 " << std::endl;
//         Add_Particle_To_Buffer( send_buffer_y1_particles, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
//         continue;
//       }
//     }
//
//     if (  direction == 2 ){
//       if ( pPos_z < G.domainMin_z ) pPos_z += ( G.domainMax_z - G.domainMin_z );
//       if ( pPos_z > G.domainMax_z ) pPos_z -= ( G.domainMax_z - G.domainMin_z );
//     }
//     if ( (direction==2) && (( pPos_z < G.zMin ) || ( pPos_z > G.zMax ))  ){
//       std::cout << "ERROR Particle out of Z domain" << std::endl;
//       std::cout << " posX: " << pPos_x << "velX: " << pVel_x << std::endl;
//       std::cout << " posY: " << pPos_y << "velY: " << pVel_y << std::endl;
//       std::cout << " posZ: " << pPos_z << "velZ: " << pVel_z << std::endl;
//       continue;
//     }
//
//     if (direction  !=2 ){
//       if ( pPos_z <= G.zMin ){
//         // std::cout << "  Adding to Z0 " << std::endl;
//         Add_Particle_To_Buffer( send_buffer_z0_particles, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
//         continue;
//       }
//       if ( pPos_z >= G.zMax ){
//         // std::cout << "  Adding to Z1 " << std::endl;
//         Add_Particle_To_Buffer( send_buffer_z1_particles, pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pPos_z);
//         continue;
//       }
//     }
//     Add_Particle_To_Vectors( pId, pMass, pPos_x, pPos_y, pPos_z, pVel_x, pVel_y, pVel_z );
//   }
// }
//
//
//
//
//
//
//
// void Part3D::Load_Particles_to_Buffer( int direction, int side, ids_vector_t &out_indxs_vec ){
//
//   int n_out;
//   // int n_send;
//   Real *send_buffer;
//
//   if ( direction == 0 ){
//     if ( side == 0 ){
//       n_out = out_indxs_vec_x0.size();
//       // n_send = std::min( n_out, N_PARTS_PER_TRANSFER  );
//       send_buffer = send_buffer_x0_particles;
//     }
//     if ( side == 1 ){
//       n_out = out_indxs_vec_x1.size();
//       // n_send = std::min( n_out, N_PARTS_PER_TRANSFER  );
//       send_buffer = send_buffer_x1_particles;
//     }
//   }
//   if ( direction == 1 ){
//     if ( side == 0 ){
//       n_out = out_indxs_vec_y0.size();
//       // n_send = std::min( n_out, N_PARTS_PER_TRANSFER  );
//       send_buffer = send_buffer_y0_particles;
//     }
//     if ( side == 1 ){
//       n_out = out_indxs_vec_y1.size();
//       // n_send = std::min( n_out, N_PARTS_PER_TRANSFER  );
//       send_buffer = send_buffer_y1_particles;
//     }
//   }
//   if ( direction == 2 ){
//     if ( side == 0 ){
//       n_out = out_indxs_vec_z0.size();
//       // n_send = std::min( n_out, N_PARTS_PER_TRANSFER  );
//       send_buffer = send_buffer_z0_particles;
//     }
//     if ( side == 1 ){
//       n_out = out_indxs_vec_z1.size();
//       // n_send = std::min( n_out, N_PARTS_PER_TRANSFER  );
//       send_buffer = send_buffer_z1_particles;
//     }
//   }
//
//
//   //HERE
//   /*
//   int offset, n_in_buffer;
//   partID_t indx, pIndx;
//   n_in_buffer = real_to_int( send_buffer[0] );
//   // send_buffer[0] = n_out;
//   // send_buffer[1] = ( n_out - n_send );
//   offset = 2 + n_in_buffer*N_DATA_PER_PARTICLE;
//   for ( indx = 0; indx<n_out; indx++ ){
//     pIndx = out_indxs_vec[indx];
//     send_buffer[ offset + 0 ] = (Real) partIDs[pIndx];
//     send_buffer[ offset + 1 ] = mass[ pIndx ];
//     send_buffer[ offset + 2 ] = pos_x[ pIndx ];
//     send_buffer[ offset + 3 ] = pos_y[ pIndx ];
//     send_buffer[ offset + 4 ] = pos_z[ pIndx ];
//     send_buffer[ offset + 5 ] = vel_x[ pIndx ];
//     send_buffer[ offset + 6 ] = vel_y[ pIndx ];
//     send_buffer[ offset + 7 ] = vel_z[ pIndx ];
//     send_buffer[0] += 1;
//     offset += N_DATA_PER_PARTICLE;
//   }
//  */
//
//   // fflush(stdout);
//   // MPI_Barrier(world);
//   // for(int i=0;i<nproc;i++)
//   // {
//   //   if(i==procID)
//   //   {
//   //     printf("  [procID %d] Sending: ( %d %d ) %d from %d\n",procID, direction, side, (int)send_buffer[0], n_out );
//   //   }
//   //   fflush(stdout);
//   //   usleep(1e2);
//   //   MPI_Barrier(world);
//   // }
// }
//
//
// void Clear_Buffers_For_Transfers( void ){
//
//   int n_data = N_PARTS_PER_TRANSFER * N_DATA_PER_PARTICLE;
//   send_buffer_x0_particles[0] = 0;
//   send_buffer_x1_particles[0] = 0;
//   send_buffer_y0_particles[0] = 0;
//   send_buffer_y1_particles[0] = 0;
//   send_buffer_z0_particles[0] = 0;
//   send_buffer_z1_particles[0] = 0;
// }
//
// void Part3D::Remove_Particles( void ){
//
//   partID_t n_delete = 0;
//   n_delete += out_indxs_vec_x0.size();
//   n_delete += out_indxs_vec_x1.size();
//   n_delete += out_indxs_vec_y0.size();
//   n_delete += out_indxs_vec_y1.size();
//   n_delete += out_indxs_vec_z0.size();
//   n_delete += out_indxs_vec_z1.size();
//
//   ids_vector_t delete_indxs_vec;
//   // delete_indxs_vec.reserve( n_delete );
//   delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_x0.begin(), out_indxs_vec_x0.end() );
//   delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_x1.begin(), out_indxs_vec_x1.end() );
//   delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_y0.begin(), out_indxs_vec_y0.end() );
//   delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_y1.begin(), out_indxs_vec_y1.end() );
//   delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_z0.begin(), out_indxs_vec_z0.end() );
//   delete_indxs_vec.insert( delete_indxs_vec.end(), out_indxs_vec_z1.begin(), out_indxs_vec_z1.end() );
//
//   out_indxs_vec_x0.clear();
//   out_indxs_vec_x1.clear();
//   out_indxs_vec_y0.clear();
//   out_indxs_vec_y1.clear();
//   out_indxs_vec_z0.clear();
//   out_indxs_vec_z1.clear();
//
//   std::sort(delete_indxs_vec.begin(), delete_indxs_vec.end());
//
//   partID_t indx, pIndx;
//   for ( indx=0; indx<n_delete; indx++ ){
//     pIndx = delete_indxs_vec.back();
//     Remove_id( pIndx, partIDs );
//     Remove_Real( pIndx, mass );
//     Remove_Real( pIndx, pos_x );
//     Remove_Real( pIndx, pos_y );
//     Remove_Real( pIndx, pos_z );
//     Remove_Real( pIndx, vel_x );
//     Remove_Real( pIndx, vel_y );
//     Remove_Real( pIndx, vel_z );
//     Remove_Real( pIndx, grav_x );
//     Remove_Real( pIndx, grav_y );
//     Remove_Real( pIndx, grav_z );
//     delete_indxs_vec.pop_back();
//     n_local -= 1;
//   }
//   if ( delete_indxs_vec.size() != 0 ) std::cout << "ERROR: deleting particles from vectors" << std::endl;
// }
//
//
//
// void Part3D::Allocate_Buffers_For_Transfers( void ){
//   int send_buffer_size = N_PARTS_PER_TRANSFER*N_DATA_PER_PARTICLE + 2;
//   send_buffer_x0_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   send_buffer_x1_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   recv_buffer_x0_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   recv_buffer_x1_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   send_buffer_y0_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   send_buffer_y1_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   recv_buffer_y0_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   recv_buffer_y1_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   send_buffer_z0_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   send_buffer_z1_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   recv_buffer_z0_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   recv_buffer_z1_particles = (Real *) malloc(send_buffer_size*sizeof(Real));
//   send_buffer_length_particles_x0 = send_buffer_size;
//   send_buffer_length_particles_x1 = send_buffer_size;
//   recv_buffer_length_particles_x0 = send_buffer_size;
//   recv_buffer_length_particles_x1 = send_buffer_size;
//   send_buffer_length_particles_y0 = send_buffer_size;
//   send_buffer_length_particles_y1 = send_buffer_size;
//   recv_buffer_length_particles_y0 = send_buffer_size;
//   recv_buffer_length_particles_y1 = send_buffer_size;
//   send_buffer_length_particles_z0 = send_buffer_size;
//   send_buffer_length_particles_z1 = send_buffer_size;
//   recv_buffer_length_particles_z0 = send_buffer_size;
//   recv_buffer_length_particles_z1 = send_buffer_size;
//
// }
//
//
//
//
//
//
// #endif //PARTICLES

#ifdef COSMOLOGY

#include <iostream>
#include <fstream>
#include "io_cosmology.h"

using namespace std;

void Load_Scale_Outputs( struct parameters P, Cosmology &Cosmo ) {

  char filename_1[100];
  // create the filename to read from
  strcpy(filename_1, P.scale_outputs_file);
  chprintf( " Loading Scale_Factor Outpus: %s\n", filename_1);

  int n_outputs = 0;
  ifstream file_out ( filename_1 );
  string line;
  Real a_value;
  if (file_out.is_open()){
    while ( getline (file_out,line) ){
      a_value = atof( line.c_str() );
      Cosmo.scale_outputs.push_back( a_value );
      n_outputs += 1;
      // chprintf("%f\n", a_value);
    }
    file_out.close();
    Cosmo.n_outputs = Cosmo.scale_outputs.size();
    Cosmo.next_output_indx = 0;
    chprintf("  Loaded %d scale outputs \n", Cosmo.n_outputs);
  }
  else{
    chprintf("  Error: Unable to open cosmology outputs file\n");
    exit(1);
  }

}

void Set_Next_Scale_Output( Cosmology &Cosmo ){

  int scale_indx = Cosmo.next_output_indx;
  Real a_value = Cosmo.scale_outputs[scale_indx];
  if  ( ( scale_indx == 0 ) && ( abs(a_value - Cosmo.current_a )<1e-3 ) )scale_indx = 1;
  a_value = Cosmo.scale_outputs[scale_indx];
  while ( a_value <= Cosmo.current_a ){
    chprintf( "%f   %f\n", a_value, Cosmo.current_a);
    scale_indx += 1;
    a_value = Cosmo.scale_outputs[scale_indx];
  }
  Cosmo.next_output_indx = scale_indx;
  Cosmo.next_output = a_value;
  chprintf( " Next output scale_factor: %f\n", Cosmo.next_output);
}

#endif

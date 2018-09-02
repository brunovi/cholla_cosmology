#ifdef COSMOLOGY

#include <iostream>
#include <fstream>
#include "io_cosmology.h"

using namespace std;

int Load_Scale_Outputs( struct parameters P, Cosmology &Cosmo ) {

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
    chprintf("  Loaded %d scale outputs \n", Cosmo.n_outputs);
  }
  else{
    chprintf("  Error: Unable to open cosmology outputs file\n");
    exit(1);
  }

}


#endif

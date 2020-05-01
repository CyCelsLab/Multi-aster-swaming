//RCS: $Id: readwritedata.cc,v 2.1 2005/04/22 19:31:06 nedelec Exp $
// This programs reads and write parameter file
// it can be useful to get a standard formatting,
// or to make modifications of parameters on the command line

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include "iomessages.h"
#include "sim_param.h"

void showHelp()
{
  printf("USAGE:\n");
  printf("  readwritedata input_data_file [ parameter modifications ]\n");
  printf("\n");
  printf("     reads parameters from input_data_file,\n");
  printf("     modify them according to the command line ( eg. mtmax=4 )\n");
  printf("     and print parameters to stdout\n");
  printf("\n");
  printf("  you can add '> output_data_file' to redirect output\n");
  printf("  example:  rwdata data0001 dt=0.001 boxsize[2]=3 > data0002\n");
  printf("  bug: warning messages from sim_param are not displayed\n");
}

//-------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  MSG.shutUp();

  if (argc < 2) {
    showHelp();
    return EXIT_SUCCESS;
  }
  
  if (NO_ERROR != MP.parseFile(argv[1])) {
    printf("Could not open file [%s] for input\n", argv[1]);
    return EXIT_FAILURE;
  }

  MP.compatibilityOperations();
  
  //parse the command line for modifications on parameter values: 
  for( int ii = 2; ii < argc ; ++ii )
    MP.parseLine( argv[ii] );
  
  //print parameters to standard output:
  MP.printFile(stdout);

  return EXIT_SUCCESS;
}

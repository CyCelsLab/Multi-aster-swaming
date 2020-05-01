//RCS: $Id: datagen6.cc,v 2.4 2005/04/12 07:24:27 nedelec Exp $
//====================================datagen.cc=============================
// Vale & Gohta simulations of Kinetochore-Fiber interactions with aster, 
// A complex of 2 dynein, and EB1 - MAP. We do a dilution of dynein concentration
// F. Nedelec, April 2005

#include <cstdio>
#include "main.h"
#include "iomessages.h"
#include "random.h"
#include "sim_param.h"
#include "smath.h"

//===========================================================================

///write the parameters in MP onto a new file data????
int writeParameters( int nb )
{
  char filename[128];
  snprintf(filename, sizeof(filename), "data%04d", nb++);
  if ( NO_ERROR == MP.printFile(filename, PARAM_NOT_DEFAULT ) ) {
    printf("generated %s\n", filename);
    return 1;
  } else {
    printf("error writing %s\n", filename);
    return 0;
  }
}

//===========================================================================

void generateParameterFiles(int first, int dilution)
{
  int cxmax0 = MP.cxmax[0];
  int cxmax1 = MP.cxmax[1];
  
  real dilutionf = 1;
  
  for( int ii = 0; ii < dilution; ++ii ) {
    
    MP.cxmax[ 0 ] = cxmax0;                 //NCD-Eb1
    MP.cxmax[ 1 ] = cxmax1 / dilutionf;     //dynein
    writeParameters( first++ );
    
    //write the same parameters without NCD
    MP.cxmax[ 0 ] = 0;
    writeParameters( first++ );
    
    dilutionf = 2 * dilutionf;
  }
  
}

//===========================================================================

int main(int argc, char* argv[])
{
  
  //read the template file into the parameter list MP:
  if (argc > 1) 
    MP.parseFile(argv[1]);
  else
    MSG.error("datagen6", "you must provide the name of a template file");
  
  MP.clearToDefaults("datafile");
  
  //complain if no value have been set:
  if ( MP.nbValuesRead() == 0 )
    MSG.error("datagen6", "wrong data.in template file, or no value read from it");
  
  //read the first number to use:
  int first = 0;
  if ( argc > 2 )
    sscanf(argv[2], "%i", &first);
  
  //re-seed the random number generator:
  RNG.seedTimer();
  
  //we force a line with randseed=0 in output:
  MP.randseed = 0;
  MP.showThisManyValues("randseed", 1);
  
  //generate files with a range of numbers:
  generateParameterFiles(first, 10);
  
  return 1;
}

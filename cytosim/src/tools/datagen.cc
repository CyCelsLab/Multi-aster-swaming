//RCS: $Id: datagen.cc,v 2.9 2005/03/03 10:49:12 nedelec Exp $
//====================================datagen.cc=============================
// example for a datagen program
// Usage: datagen data.gen 0 1000 
// produces 1000 data files numbered from 0000 to 0999
// with data.gen as default input

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
//===========================================================================

void generateParameterFiles( int first, int last )
{
  int cnt = first;
  
  while( cnt < last ) {
        
    for(int mo = 0 ; mo < 2 ; ++mo )	{
      
      if ( mo == 1 )
        MP.haspeed[ mo ]         = RNG.choose( -0.5, -1.0 );
      else
        MP.haspeed[ mo ]         = RNG.choose( 0.5, 1.0 );
      
      MP.haattachrate[ mo ]    = RNG.choose( 1, 2, 4, 8, 16 );
      MP.hadetachrate[ mo ]    = RNG.choose( 0.01, 0.02, 0.08, 0.32, 1.28 );
    }
    
    MP.mtinit = 35;
    writeParameters( cnt++ );

    MP.mtinit = 25;
    writeParameters( cnt++ );
}
  
}


//===========================================================================

int main(int argc, char* argv[])
{
  
  //read the template file into the parameter list MP:
  if (argc > 1) 
    MP.parseFile(argv[1]);
  else
    MP.parseFile("data00.gen");

  //we forget the name of the datagen file:
  MP.clearToDefaults("datafile");
  
  //we force a line with randseed=0 in output:
  MP.randseed = 0;
  MP.showThisManyValues("randseed", 1);
  
  //complain if no value have been set:
  if ( MP.nbValuesRead() == 0 )
    MSG.error("datagen.cc", "wrong data.in template file, or no value read from it");
  
  //re-seed the random number generator:
  RNG.seedTimer();
  
  //read the first number to use:
  int first = 0;
  if ( argc > 2 )
    sscanf(argv[2], "%i", &first);
  int last = first+1000;
  if ( argc > 3 )
    sscanf(argv[3], "%i", &last);

  //generate files with a range of numbers:
  generateParameterFiles(first, last);
  
  return 1;
}

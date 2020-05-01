//RCS: $Id: datagen2.cc,v 2.11 2005/04/22 20:00:19 nedelec Exp $
//====================================datagen.cc=============================
// Vale & Gohta simulations of Kinetochore-Fiber interactions with aster, 
// A complex of plus-end binding with minus-end motor
// F. Nedelec, Feb 2005

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

void generateParameterFiles( int first, int last )
{
  int cnt = first;
  
  while( cnt < last ) {
    
start:
    MP.asmtmax[0]  = RNG.int_range( 30, 70 );
    MP.km          = RNG.choose( 40, 80, 160 );
    
    for(int mo = 0 ; mo < 3 ; ++mo )	{
      MP.hamodel[ mo ]         = 0;
      MP.hadetachmodel[ mo ]   = 2;
      MP.haforce[ mo ]         = RNG.choose( 1, 2, 3, 4, 5 );
      MP.haattachdist[ mo ]    = MP.haforce[mo] / MP.km;

      MP.haattachrate[ mo ]    = RNG.choose( 1, 2, 4, 8, 16 );
      MP.hadetachrate[ mo ]    = RNG.choose( 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56 );
      
      //the EB1 should have an unbinding rate faster than 0.05,
      //in order for the commet to be less than 1 um.
      if ( MP.haendattachdist[mo]  )
        MP.hadetachrate[ mo ]  = RNG.choose( 0.16, 0.32 );

      if ( mo == 2 )
        MP.haenddetachrate[mo] = RNG.choose( 2, 4, 8, 16, 32, 64, 128 );

      //that is the condition on the efficiency
      if ( 2*sqrt( DIM * MP.kT * MP.km ) > MP.haforce[mo] ) 
        goto start;
    }
    
    //write the parameters with both NCD and Dynein
    MP.cxmax[ 0 ] = RNG.int_range( 100, 32000 );  //NCD-Eb1
    MP.cxmax[ 1 ] = RNG.int_range( 0, 12000 );     //dynein
    writeParameters( cnt++ );
    
    //write the same parameters without NCD
    MP.cxmax[ 0 ] = 0;
    writeParameters( cnt++ );

    //write the same with less dynein
    MP.cxmax[1] /= 2;
    writeParameters( cnt++ );

    //write the same with less dynein
    MP.cxmax[1] /= 2;
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
    MP.parseFile("data02.gen");
  
  MP.clearToDefaults("datafile");

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
  
  //we force a line with randseed=0 in output:
  MP.randseed = 0;
  MP.showThisManyValues("randseed", 1);
  
  //generate files with a range of numbers:
  generateParameterFiles(first, last);
  
  return 1;
}

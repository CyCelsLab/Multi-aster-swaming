//RCS: $Id: datagen8.cc,v 1.4 2005/03/01 14:01:35 nedelec Exp $
//==================================datagen8.cc=============================
// Gradient of MT catastrophe & Gradient of Grafted-Motors density
// PredocPractical, Oct 2004


#include <stdio.h>
#include "main.h"
#include "iomessages.h"
#include "random.h"
#include "sim_param.h"
#include "smath.h"

//===========================================================================
//===========================================================================
//===========================================================================
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

                   

int main(int argc, char* argv[])
{
  if (argc > 1) 
    MP.parseFile(argv[1]);
  else
    MP.parseFile("data08.gen");

  if ( MP.nbValuesRead() == 0 )
    MSG.error("datagen8.cc", "wrong template file, or no value read from it");

  MP.showThisManyValues("randseed", 1);
  MP.randseed = 0;
  RNG.seedTimer();

  int start = 1;
  if ( argc >= 3 )
    sscanf(argv[2], "%i", &start );
  int filecnt = start;
  start += 2000;

  
  while ( filecnt <= start ) {
    
    start:

    real left_value  = RNG.real_uniform_range( 0.05, 0.4 );
    real right_value = RNG.real_uniform_range( 0.05, 0.4 );
    
    MP.ghmax[0]       = RNG.int_range( 0, 20000 );
    
    for(int mo = 0 ; mo < 1 ; ++mo ) {
        
      MP.haspeed[ mo ]         = -0.01 * RNG.pint_inc(100);
      MP.haattachdist[ mo ]    = RNG.choose( 0.005, 0.01, 0.02, 0.04 );
      MP.haforce[ mo ]         = RNG.choose( 1, 2, 3, 4, 5);
      MP.haattachrate[ mo ]    = RNG.choose( 1, 2, 4, 8, 16, 32 );
      MP.hadetachrate[ mo ]    = RNG.choose( 0.01, 0.02, 0.08, 0.32, 1.28, 2.56 );
      
      if ( 2*sqrt( DIM * MP.kT * MP.km ) > MP.haforce[mo] )  {
	      //printf("%f %f %f\n", MP.km, 2*sqrt(DIM*MP.kT*MP.km), 
	      // MP.haforce[mo]); 
	      goto start; 
	    }
    }
    
    MP.ghinit[0]         = 20; //gradient

    MP.mtcatagradient[0] = left_value;
    MP.mtcatagradient[1] = right_value;
    writeParameters( filecnt++ );
    
    MP.mtcatagradient[1] = left_value;
    MP.mtcatagradient[0] = right_value;
    writeParameters( filecnt++ );
    
    MP.mtcatagradient[0] = right_value;
    MP.mtcatagradient[1] = right_value;
    writeParameters( filecnt++ );
    
    MP.mtcatagradient[1] = left_value;
    MP.mtcatagradient[0] = left_value;
    writeParameters( filecnt++ );

    MP.mtcatagradient[0] = left_value;
    MP.mtcatagradient[1] = right_value;
    MP.ghinit[0]   = 0;  //no gradient
    writeParameters( filecnt++ );

    
  }
  
  return 1;
}

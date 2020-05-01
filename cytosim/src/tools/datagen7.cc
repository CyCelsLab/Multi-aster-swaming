//RCS: $Id: datagen7.cc,v 2.3 2005/03/01 14:01:31 nedelec Exp $
//==================================datagen7.cc=============================
// Two centrosome-asters with one hetero-complex with cxlength > 0
// F. Nedelec, Sept 16 2004


#include <stdio.h>
#include "main.h"
#include "iomessages.h"
#include "random.h"
#include "sim_param.h"
#include "smath.h"

//===========================================================================
//===========================================================================
//===========================================================================

int main(int argc, char* argv[])
{
  if (argc > 1) 
    MP.parseFile(argv[1]);
  else
    MP.parseFile("data07.gen");

  if ( MP.nbValuesRead() == 0 )
    MSG.error("datagen7.cc", "wrong template file, or no value read from it");

  MP.showThisManyValues("randseed", 1);
  MP.randseed = 0;
  RNG.seedTimer();

  int start = 1;
  if ( argc >= 3 )
    sscanf(argv[2], "%i", &start );
  int filecnt = start;
  start += 5000;

  char filename[128];
  while ( filecnt <= start ) {
    
    start:

    MP.km          = RNG.choose( 10, 20, 40, 80 );
    MP.asmtmax[0]  = RNG.int_range( 12, 80 );
    MP.cxmax[0]    = RNG.int_range( 1000, 30000 );
    MP.cxlength[0] = RNG.choose( 0.025, 0.05, 0.1, 0.2 );
    
    for(int mo = 0 ; mo < 2 ; ++mo ) {
        
      MP.haspeed[ mo ]         = 0.01 * RNG.sint_inc(100);
      MP.haattachdist[ mo ]    = MP.cxlength[0] + 0.005;
      MP.haforce[ mo ]         = RNG.choose( 1, 2, 3, 4, 5);
      MP.haattachrate[ mo ]    = RNG.choose( 1, 2, 4, 8, 16, 32 );
      MP.hadetachrate[ mo ]    = RNG.choose( 0.01, 0.02, 0.08, 0.32, 1.28 );
      
      if ( 2*sqrt( DIM * MP.kT * MP.km ) > MP.haforce[mo] )  {
	      //printf("%f %f %f\n", MP.km, 2*sqrt(DIM*MP.kT*MP.km), 
	      // MP.haforce[mo]); 
	      goto start; 
	    }
    }

    snprintf(filename, sizeof(filename), "data%04d", filecnt++);
    MP.printFile(filename, PARAM_NOT_DEFAULT );
    printf("generated %s\n", filename);
  }
  
  return 1;
}

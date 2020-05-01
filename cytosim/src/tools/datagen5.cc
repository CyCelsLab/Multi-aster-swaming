//RCS: $Id: datagen5.cc,v 2.3 2005/03/01 14:01:23 nedelec Exp $
//============================   datagen5.cc   =============================
// two asters + two soluble homo-complexes with treadmilling flux inwards


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
    MP.parseFile("data05.gen");

  if ( MP.nbValuesRead() == 0 )
    MSG.error("datagen5.cc", "wrong template file, or no value read from it");

  MP.showThisManyValues("randseed", 1);
  MP.randseed = 0;
  RNG.seedTimer();

  int start = 1;
  if ( argc >= 3 ) sscanf(argv[2], "%i", &start );
  int filecnt = start;
  start += 1000;

  char filename[128];
  while ( filecnt <= start ) 
    {
    start:

      MP.km          = 10 * ( 1 << RNG.pint_inc(4) );
      MP.asmtmax[0]  = 20 + RNG.pint_inc( 70 );

      for(int mo = 0 ; mo < 2 ; ++mo )
	{
	  MP.hamodel[ mo ]         = 3;
	  MP.haspeed[ mo ]         = 0.01 * RNG.sint_inc(100);

	  MP.haattachdist[ mo ]    = 0.005 * ( 1 << RNG.pint_inc(4) );
	  MP.haforce[ mo    ]      = MP.km * MP.haattachdist[ mo ];
	  
	  MP.haattachrate[ mo ]    = 2 + 50 * RNG.preal();
	  MP.haattachrate[ mo ]    = floor( MP.haattachrate[mo] );
	  
	  MP.hadetachrate[ mo ]    = 0.01 + (exp(RNG.preal())-1)/(exp(1)-1);
	  MP.hadetachrate[ mo ]    = floor( 100 * MP.hadetachrate[mo] )/100.0;

	  MP.haenddetachrate[ mo ] = 999;

	  if ( 2*sqrt( DIM * MP.kT * MP.km ) > MP.haforce[mo] ) 
	    {
	      //printf("%f %f %f\n", MP.km, 2*sqrt(DIM*MP.kT*MP.km), 
	      // MP.haforce[mo]); 
	      goto start; 
	    }
	}

      if ( MP.haspeed[0] * MP.haspeed[1] > 0 ) 
	goto start;

//complexes:

      MP.cxmax[ 0 ] = RNG.pint_inc( 15000 );
      MP.cxmax[ 1 ] = RNG.pint_inc( 15000 );
      
      snprintf(filename, sizeof(filename), "data%04d", filecnt++);
      MP.printFile(filename, PARAM_NOT_DEFAULT);
      printf("generated %s\n", filename);
    }

  return 1;
}

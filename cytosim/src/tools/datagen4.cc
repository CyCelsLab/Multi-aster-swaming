//RCS: $Id: datagen4.cc,v 2.4 2005/03/01 14:01:20 nedelec Exp $
//==============================  datagen4.cc  =============================
// Sept 2003: Complex EndBindingProtein-motor + a soluble homo complex

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
  RNG.seedTimer();

  if (argc > 1) 
    MP.parseFile(argv[1]);
  else
    MP.parseFile("data04.gen");

  if ( MP.nbValuesRead() == 0 )
    MSG.error("datagen4.cc", "wrong template file, or no value read from it");

  int mo, filecnt = 1, start = 1;

  if ( argc >= 3 ) sscanf(argv[2], "%i", &start );
  filecnt += start;

  MP.randseed = 0;

  char filename[128];

  start += 2000;
  while ( filecnt <= start ) {
    
start:

    MP.km          = 10 * ( 1 << RNG.pint_inc(3) );
    MP.asmtmax[0]  = 20 + RNG.pint_inc( 70 );

    for( mo = 0 ; mo < 2 ; ++mo )	{
      
      MP.haspeed[ mo ]         = 0.01 * RNG.sint_inc(100);
      
      MP.haattachdist[ mo ]    = 0.01 * ( 1+RNG.pint_inc(4) );
      MP.haforce[ mo ]         = MP.km * MP.haattachdist[ mo ];
      
      MP.haattachrate[ mo ]    = 25 * RNG.preal();
      MP.haattachrate[ mo ]    = floor( MP.haattachrate[mo] );
      
      MP.hadetachrate[ mo ]    = 0.01 + (exp(RNG.preal())-1)/(exp(1)-1);
      MP.hadetachrate[ mo ]    = floor( 100 * MP.hadetachrate[mo] )/100.0;
      
      if ( 2*sqrt( DIM * MP.kT * MP.km ) > MP.haforce[mo] ) {
	      //printf("%f %f %f\n", MP.km, 2*sqrt(DIM*MP.kT*MP.km), MP.haforce[mo]); 
	      goto start; 
      }
    }

    MP.haendattachdist[2]  = 0.5 * ( 1 + RNG.pint_exc(4) );

    //if ( MP.haspeed[0] * MP.haspeed[1] > 0 ) 	goto start;
    //if ( MP.haspeed[0] > 0 )                 goto start;
    
    //complexes:

    MP.cxmax[ 0 ] = RNG.pint_inc( 15000 );
    MP.cxmax[ 1 ] = RNG.pint_inc( 15000 );
      
    snprintf(filename, sizeof(filename), "data%04d", filecnt++);
    printf("generated %s\n", filename);
    MP.printFile(filename, PARAM_NOT_DEFAULT);
  }

  return 1;
}

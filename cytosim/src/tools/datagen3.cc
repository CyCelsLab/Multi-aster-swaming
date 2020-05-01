//RCS: $Id: datagen3.cc,v 2.10 2005/03/03 10:49:23 nedelec Exp $
//====================================datagen.cc=============================
// Vale & Gohta simulations of Kinetochore-Fiber interactions with aster, 
// F. Nedelec, Jan 2005

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
    MP.km = RNG.choose( 10, 20, 40, 80, 160 );
    
    for(int mo = 0 ; mo < 3 ; ++mo )	{
      MP.hamodel[ mo ]         = 0;
      MP.hadetachmodel[ mo ]   = 2;
      MP.haforce[ mo ]         = RNG.choose( 1, 2, 3, 4, 5 );
      MP.haattachdist[ mo ]    = MP.haforce[mo] / MP.km;

      if ( mo == 0 )
        MP.haspeed[ mo ]       = RNG.choose( -0.1, -0.2, -0.4, -0.8 );
      MP.haattachrate[ mo ]    = RNG.choose( 1, 2, 4, 8, 16 );
      MP.hadetachrate[ mo ]    = RNG.choose( 0.01, 0.02, 0.08, 0.32, 1.28, 2.56, 5.12 );
      
      //MP.haenddetachrate[ mo ] = RNG.choose( 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 50, 100 );
      if ( 2*sqrt( DIM * MP.kT * MP.km ) > MP.haforce[mo] ) 
        goto start;
    }
    
    MP.cxmax[ 0 ] = RNG.int_range( 100, 54000 );
    MP.cxmax[ 1 ] = RNG.int_range( 0, 6000 );
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
    MP.parseFile("data03.gen");
  
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

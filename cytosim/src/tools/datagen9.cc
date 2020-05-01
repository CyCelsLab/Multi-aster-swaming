//RCS: $Id: datagen9.cc,v 1.9 2005/04/08 12:04:34 foethke Exp $
//====================================datagen8.cc=============================
// different dynamics for the microtubules in S.pombe
// Usage: datagen9 data.gen 0 1000 
// produces 1000 data files numbered from 0000 to 0999
// with data.gen as default input

#include <stdio.h>
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
    MP.parseFile("data_pom1.gen");

  if ( MP.nbValuesRead() == 0 )
    MSG.error("datagen9.cc", "wrong template data file, or no value read from it");
  
  RNG.seedTimer();
  

  //we forget the name of the datagen file:
  MP.clearToDefaults("datafile");
  
  //we force a line with randseed=0 in output:
  MP.randseed = 0;
  MP.showThisManyValues("randseed", 1);
  //MP.showThisManyValues("mtdynforce", 1);

  int start = 1;
  if ( argc >= 3 ) sscanf(argv[2], "%i", &start );
  int filecnt = start;
  if ( argc >= 4 ) {
    sscanf(argv[3], "%i", &start );
    start += filecnt;
  } else
    start += 10;

  char filename[128];
  /* Values taken from "Dynamics of interphase microtubules in Schizosaccharomyces pombe",
     Douglas R. Drummond and Robert A. Cross, Curr. Biol. 10:766-775, 2000 */
  real growth      = 0.05;
  real shrinkage   = -0.0745;
  real catastrophe = 0.017;
  real rescue      = 0.031;
  //real stallForce  = 2.0;
  while ( filecnt < start ) 
    {
    //start:
	MP.initsize[2]       = RNG.real_uniform_range(-7.5, 7.5);
	MP.mtdynspeed[0]     = 0.;
    MP.mtdyntrans[0]     = 0.;
    MP.mtdynspeed[1]     = 0.;
    MP.mtdyntrans[1]     = 0.;
    MP.mtdynspeed[2]     = 0.02 + 0.1*RNG.preal();
    MP.mtdyntrans[2]     = 0.017;
    MP.mtdynspeed[3]     = -0.0745;
    MP.mtdyntrans[3]     = 0.;
    MP.mtcatgrowthrel[0] = 0.;
    MP.mtcatgrowthrel[1] = 1000 + 11000*RNG.preal();

	writeParameters( filecnt++ );
    }

  return 1;
}


//RCS: $Id: createstart.cc,v 2.0 2004/08/16 16:45:20 nedelec Exp $
//------------------------create state--------------------------------
#include "sim_param.h"
#include "iomessages.h"
#include "sim.h"
#include "exceptions.h"


void createConfiguration()
{
  unsigned long ouf = 0;
  
  real len;
  Grafted * gh;
  Microtub * mt;
  Vecteur place, box, dir;
  
  MSG.shutUp();
  if (NO_ERROR != MP.parseFile(MP.datafile)) 
    MSG.error("createstart.cc","Cannot read file %s", MP.datafile);
  MP.massageParameters();
  
  sim.resetCounters();
  sim.setSpace();
  
  MSG.setVerboseLevel(4);
  
  int grafted_type = 0;
  for (int jj=0; jj < 300; ++jj)
    {
      do {
        place = Grafted::space->randomPlaceInVolume();
        
        if ( ++ouf > MAX_TRIALS_BEFORE_STUCK )
          throw StuckException("reached MAX_TRIALS_BEFORE_STUCK while placing grafteds");

      } while ( !( (place.XX > -5) && (place.XX < 5)));
      gh=new Grafted(grafted_type, place); 
      sim.link(gh);
    }


  for (int jj=0; jj < 400; ++jj) {
    place = Microtub::space->randomPlaceInVolume();
    dir   = Vecteur::randSphere();
    len   = Microtub::initLength();
    mt = new Microtub(place,dir,len, MT_ORIGIN);
    sim.link(mt);
  }

  sim.printShortDescription();
  sim.writeState("start.in");
}

int main( int argc, char * argv[] )
{
  try {
    createConfiguration();
  } catch( StuckException e ) {
    MSG("createConfiguration %s\n", e.getMessage() );
  }
      
  return 0;
}

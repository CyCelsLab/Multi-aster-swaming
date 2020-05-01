//RCS: $Id: reducedata.cc,v 2.0 2004/08/16 16:46:00 nedelec Exp $
//=======================formatdata.cc=============================
// This read and write the parameter file, hiding some parameters,
// which have the default value.
// This was used to get a standard format for web pages.

#include <stdio.h>
#include "iomessages.h"
#include "sim_param.h"


void reduce()
{
  MP.hide("version");
  MP.hide("dimension");
  MP.hide("record");
  MP.hide("dt");
  MP.hide("nbiter");
  MP.hide("kT");
  MP.hide("visc");
  MP.hide("boxshape");
  MP.hide("mtmax");
  MP.showThisManyValues("asmtmax");
  //MP.hide("mtrodlength");
  MP.hide("mtrigid");
  MP.hide("mtminlength");
  //MP.hide("mtinitlength");
  MP.hide("mtmonlength");
  MP.hide("mtmaxlength");
  MP.hide("ghmax");
  MP.hide("assize");
  MP.hide("initcode");
  MP.hide("randseed");
  MP.hide("mtconfine");

  //Dirty: we shift dynamic instability values of plus-end at beginning:

  if (( MP.mtdyntrans[0] == 0 ) && ( MP.mtdyntrans[1] == 0 ))  {
    for(int i=0; i<2; ++i) {
      MP.mtdyntrans[i] = MP.mtdyntrans[i+2];
      MP.mtdynspeed[i] = MP.mtdynspeed[i+2];
    }
    MP.showThisManyValues("mtdyntrans", 2);
    MP.showThisManyValues("mtdynspeed", 2);
  }
}

int main(int argc, char* argv[])
{
  if (argc == 1) {
    if (NO_ERROR==MP.parseFile(DATA_OUT)) { 
      MP.compatibilityOperations();
      reduce();
      MP.printFile(DATA_OUT);
    }
    else if (NO_ERROR==MP.parseFile(MP.datafile)) {
      MP.compatibilityOperations();
      reduce();
      MP.printFile(DATA_OUT);
    }
  }

  if (argc == 2) {
    if (NO_ERROR==MP.parseFile(argv[1])) { 
      MP.compatibilityOperations();
      reduce();
      MP.printFile(DATA_OUT);
    }
  }

  if (argc == 3) {
    if (NO_ERROR==MP.parseFile(argv[1]))	{
      MP.compatibilityOperations();
      reduce();
      MP.printFile(argv[2]);
    }
  }

  return 0;
}

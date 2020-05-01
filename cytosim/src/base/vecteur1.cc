//RCS: $Id: vecteur1.cc,v 2.1 2005/02/14 21:50:36 nedelec Exp $
//----------------------------------------------------------------------------
//                            Vecteur1.cc
//----------------------------------------------------------------------------

#include "vecteur1.h"
#include "random.h"


void Vecteur1::setRandomNormed() { 
  XX = RNG.sflip();
}

void Vecteur1::setRandomSphere() { 
  XX = RNG.sreal(); 
}

void Vecteur1::setRandomSphere(const real & s) { 
  XX = s * RNG.sreal(); 
}

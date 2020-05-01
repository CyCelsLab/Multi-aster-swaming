//RCS: $Id: vecteur2.cc,v 2.3 2005/02/14 21:51:32 nedelec Exp $
//----------------------------------------------------------------------------
//                            Vecteur2.cc
//----------------------------------------------------------------------------

#include "vecteur2.h"
#include "random.h"


void Vecteur2::setRandomNormed() {
  real n;
  do {
    XX = RNG.sreal();
    YY = RNG.sreal();
    n = XX*XX + YY*YY;
  } while ( ( n  > 1.0 ) || ( n == 0 ) );
  n = sqrt( n );
  XX /= n;
  YY /= n;
}

void Vecteur2::setRandomSphere() {
  do {
    XX = RNG.sreal();
    YY = RNG.sreal();
  } while ( XX*XX + YY*YY > 1.0 );
}

void Vecteur2::setRandomSphere(const real & s) {
  do {
    XX = RNG.sreal();
    YY = RNG.sreal();
  } while ( XX*XX + YY*YY > 1.0 );
  XX *= s;
  YY *= s;
}

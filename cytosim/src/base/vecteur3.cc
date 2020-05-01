//RCS: $Id: vecteur3.cc,v 2.3 2005/02/14 21:51:53 nedelec Exp $
//----------------------------------------------------------------------------
//                            Vecteur3.cc
//----------------------------------------------------------------------------


#include "vecteur3.h"
#include "random.h"


//hypercube rejection method
void Vecteur3::setRandomNormed() {
  real n;
  do {
    XX = RNG.sreal();
    YY = RNG.sreal();
    ZZ = RNG.sreal();
    n = XX*XX + YY*YY + ZZ*ZZ;
  } while ( ( n  > 1.0 ) || ( n == 0 ) );
  n = sqrt( n );
  XX /= n;
  YY /= n;
  ZZ /= n;
}

void Vecteur3::setRandomSphere() {
  do {
    XX = RNG.sreal();
    YY = RNG.sreal();
    ZZ = RNG.sreal();
  } while (  XX*XX + YY*YY + ZZ*ZZ > 1.0 );
}

void Vecteur3::setRandomSphere(const real & s) {
  do {
    XX = RNG.sreal();
    YY = RNG.sreal();
    ZZ = RNG.sreal();
  } while (  XX*XX + YY*YY + ZZ*ZZ > 1.0 );
  XX *= s;
  YY *= s;
  ZZ *= s;
}


void Vecteur3::computeEulerAngles(real angles[]) const {
  angles[0] = atan2( YY, XX );
  angles[1] = atan2( sqrt( XX * XX + YY * YY ), ZZ );      
}

void Vecteur3::setFromEulerAngles(const real angles[]) {
  XX = cos(angles[0])*sin(angles[1]);
  YY = sin(angles[0])*sin(angles[1]);
  ZZ = cos(angles[1]);
}

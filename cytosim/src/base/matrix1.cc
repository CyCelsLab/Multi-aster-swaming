//RCS: $Id: matrix1.cc,v 2.0 2004/08/16 16:35:58 nedelec Exp $
//--------------------------------matrix1.h---------------------------------
#include "matrix1.h"
#include "random.h"


//-------------------------------------------------------    
void Matrix11::setRotationAroundAxisEuler(const real x, const real, const real) {
  val[0] = ( x > 0 ) ? 1 : -1;
}

//-------------------------------------------------------    
void Matrix11::setRotationFromEulerAngles(const real x, const real, const real) {
  val[0] = ( x > 0 ) ? 1 : -1;
}

//-------------------------------------------------------    
Matrix11 Matrix11::randomRotation() {
  Matrix11 result;
  result.val[0] = RNG.sflip();
  return result;
}


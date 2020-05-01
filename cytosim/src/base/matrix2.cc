//RCS: $Id: matrix2.cc,v 2.3 2004/09/30 14:07:56 nedelec Exp $

#include "matrix2.h"
#include "random.h"
#include "smath.h"
#define SZ 2


//-------------------------------------------------------    
real Matrix22::determinant() const {
  return ( val[0] * val[3] - val[1] * val[2] );
}


//-------------------------------------------------------    
Matrix22 Matrix22::inverse() const {
  Matrix22 result;
  real det = determinant();
  if ( det != 0 ) 
    det = 1.0 / det;
  result.val[0] =  det * val[3];
  result.val[1] = -det * val[1];
  result.val[2] = -det * val[2];
  result.val[3] =  det * val[0];
  return result;
}


//-------------------------------------------------------    
void Matrix22::setRotationFromEulerAngles(const real angle, const real, const real) {
  real ca = cos( angle ), sa = sin( angle );
  val[0+SZ*0] =  ca;
  val[1+SZ*0] =  sa;
  val[0+SZ*1] = -sa;
  val[1+SZ*1] =  ca;
}
  
//-------------------------------------------------------    
void Matrix22::computeEulerAngles(real * angle, real *, real *) {
  * angle = atan2( val[1+SZ*0], val[0+SZ*0] );
}

  
//-------------------------------------------------------    
void Matrix22::rotateAroundAxis( const real angle[2] )
{
  static Matrix22 R;
  R.setRotationFromEulerAngles( angle[0], 0, 0);
  leftMult(R);
}

//-------------------------------------------------------    
Matrix22 Matrix22::randomRotation() {
  Matrix22 result;
  result.setRotationFromEulerAngles( PI*RNG.sreal(), 0, 0);
  return result;
}

#undef SZ


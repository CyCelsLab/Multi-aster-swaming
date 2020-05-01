//RCS: $Id: matrix3.cc,v 2.0 2004/08/16 16:36:11 nedelec Exp $

#include "matrix3.h"
#define SZ 3
#include "random.h"
#include "smath.h"

//-------------------------------------------------------    
real Matrix33::determinant() const {
  return ( + val[0]*val[4]*val[8] + val[2]*val[3]*val[7] + val[1]*val[5]*val[6]
           - val[2]*val[4]*val[6] - val[1]*val[3]*val[8] - val[0]*val[5]*val[7] );
}



//-------------------------------------------------------    
// set this to a rotation of angle a, around one of the main axis X, Y or Z
void Matrix33::setRotationAroundMainAxis(const int axis, const real angle) {
  real ca = cos( angle );
  real sa = sin( angle );
  
  for(int i=0; i<SZ; ++i)
    for(int j=0; j<SZ; ++j) {
      if ( (i!=axis) && (j!=axis) ) {
        if (i == j)
          val[i+SZ*j] = ca;
        else
          val[i+SZ*j] = ((i<j)?-1:1) * ((axis==1)?-1:1) * sa;	      
      } else {
        val[i+SZ*j] = 0;
      }
    }
      val[axis+SZ*axis] = 1.;
}


//-------------------------------------------------------    
void Matrix33::setRotationFromEulerAngles(const real a, const real b, const real c) {
  real ca = cos( a ), sa = sin( a );
  real cb = cos( b ), sb = sin( b );
  real cc = cos( c ), sc = sin( c );
  
  val[0+SZ*0] =  ca*cb;
  val[1+SZ*0] =  sa*cb;
  val[2+SZ*0] = -sb;
  
  val[0+SZ*1] =  ca*sb*sc - sa*cc;
  val[1+SZ*1] =  sa*sb*sc + ca*cc;
  val[2+SZ*1] =  cb*sc;
  
  val[0+SZ*2] =  ca*sb*cc + sa*sc;
  val[1+SZ*2] =  sa*sb*cc - ca*sc;
  val[2+SZ*2] =  cb*cc;
}


//-------------------------------------------------------    
void Matrix33::computeEulerAngles(real * a, real * b, real * c) {
  real cb = sqrt( val[0] * val[0] + val[1] * val[1]);
  
  *b = atan2( -val[2+SZ*0], cb );
  
  if ( cb ) {
    *a = atan2( val[1+SZ*0], val[0+SZ*0] );
    *c = atan2( val[2+SZ*1], val[2+SZ*2] );
  } else {
    *a = 0;
    *c = atan2( -val[0+SZ*1], val[1+SZ*1] );
  }
}


//-------------------------------------------------------    
// set this to a rotation around axis of azimuth b and elevation c
void Matrix33::setRotationAroundAxisEuler(const real a, const real b, const real c) {
  real ca = cos( a ), sa = sin( a ), ca1 = 1 - ca;
  real cb = cos( b ), sb = sin( b );
  real cc = cos( c ), sc = sin( c );
  real sacc        = sa * cc,           sasc        = sa * sc;
  real saccsb      = sacc * sb,         sacccb      = sacc * cb;
  real ccccca1     = cc * cc * ca1,     ccscca1     = cc * sc * ca1;
  real sbccscca1   = sb * ccscca1,      cbccscca1   = cb * ccscca1;
  real cbcbccccca1 = cb * cb * ccccca1, cbsbccccca1 = cb * sb * ccccca1;
  
  val[0+SZ*0] =  cbcbccccca1 + ca;
  val[0+SZ*1] =  cbsbccccca1 - sasc;
  val[0+SZ*2] =  cbccscca1   + saccsb;
  
  val[1+SZ*0] =  cbsbccccca1 + sasc;
  val[1+SZ*1] =  ca - cbcbccccca1 + ccccca1;
  val[1+SZ*2] =  sbccscca1   - sacccb;
  
  val[2+SZ*0] =  cbccscca1 - saccsb;
  val[2+SZ*1] =  sbccscca1 + sacccb;
  val[2+SZ*2] =  1 - ccccca1;
}


//-------------------------------------------------------    
// set this to a rotation around the given axis
void Matrix33::setRotationAroundAxis(const real axis[SZ], const real angle) {
  real aa = atan2( axis[1], axis[0] );
  real ab = atan2( axis[2], sqrt( axis[0]*axis[0] + axis[1]*axis[1] ) );
  setRotationAroundAxisEuler( angle, aa, ab );
}


//-------------------------------------------------------    
// multiply *this on the left by a rotation matrix matrix
void Matrix33::rotateAroundAxis(const real axis[SZ]) {
  static Matrix33 R;    
  real angle = sqrt( axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2] );
  if (angle) {
    real aa = atan2( axis[1], axis[0] );
    real ab = atan2( axis[2], sqrt( axis[0]*axis[0] + axis[1] * axis[1] ) );
    
    R.setRotationAroundAxisEuler( angle, aa, ab );
    leftMult(R);
  }
}

//-------------------------------------------------------    
Matrix33 Matrix33::randomRotation() {
  Matrix33 result;
  //TODO I am not sure that random Euler angles give uniform sampling of rotations... to be verified
  result.setRotationFromEulerAngles( PI*RNG.sreal(), PI*RNG.sreal(), PI*RNG.sreal() );
  return result;
}

#undef SZ


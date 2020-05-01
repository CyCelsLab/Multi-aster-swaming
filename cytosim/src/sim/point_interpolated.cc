//RCS: $Id: point_interpolated.cc,v 2.8 2005/04/10 14:53:36 nedelec Exp $

#include "point_interpolated.h"
#include "point_exact.h"
#include "microtub.h"


void PointInterpolated::clear() 
{
  mPS     = 0;
  mPoint1 = 0;
  mPoint2 = 0;
  mCoef   = 0;
}


Vecteur PointInterpolated::where() const
{ 
  //TODO: optimize PointInterpolated::where() by precalculating dpts = pts[n+DIM]-pts[n]
  Vecteur result;
  const real * pspts = mPS -> getPts();
  result.XX = pspts[ DIM * mPoint1    ] + mCoef * ( pspts[ DIM * mPoint2    ] - pspts[ DIM * mPoint1    ] );
#if ( DIM > 1 )
  result.YY = pspts[ DIM * mPoint1 +1 ] + mCoef * ( pspts[ DIM * mPoint2 +1 ] - pspts[ DIM * mPoint1 +1 ] );
#endif
#if ( DIM > 2 )
  result.ZZ = pspts[ DIM * mPoint1 +2 ] + mCoef * ( pspts[ DIM * mPoint2 +2 ] - pspts[ DIM * mPoint1 +2 ] );
#endif
  return result;
}

Vecteur PointInterpolated::dirInter() const
{ 
  //we could optimize for MT, as we know the distance between points
  //we would need a function dPtsNormalized() in PointSet
  return ( mPS->whereP(mPoint2) - mPS->whereP(mPoint1) ).normalized();
}

real PointInterpolated::length() const
{ 
  return ( mPS->whereP(mPoint2) - mPS->whereP(mPoint1) ).norm(); 
}

bool PointInterpolated::near(const PointExact & ip ) const {
  return (( mPS == ip.getPS() ) &&
          (( mPoint1 == ip.getPoint() ) || ( mPoint2 == ip.getPoint() )));
}

bool PointInterpolated::near(const PointInterpolated & ip ) const {
  return (( mPS == ip.mPS ) &&
          (( mPoint1 == ip.mPoint1 ) || ( mPoint1 == ip.mPoint2 ) ||
           ( mPoint2 == ip.mPoint1 ) || ( mPoint2 == ip.mPoint2 )));
}

int PointInterpolated::looksWrong(const bool forInteraction) const
{
  if ( mPS == 0 )                        return 1;
  if ( mPoint1 < 0 )                     return 2;
  if ( mPoint1 > mPS->lastSegment() )    return 3;
  if ( mPoint2 < 0 )                     return 4;
  if ( mPoint2 > mPS->lastPoint() )      return 5;
  if (( mCoef < 0) || ( mCoef > 1.0 ))   return 6;
  if ( forInteraction ) {
    if ( !mPS ->isLinked() )             return 7;
    if ( mPS -> matIndex() < 0 )         return 8;
  }
  return NO_ERROR;
}

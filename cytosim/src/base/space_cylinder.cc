//RCS: $Id: space_cylinder.cc,v 2.0 2004/08/16 16:38:12 nedelec Exp $

#include "space_cylinder.h"

//-------------------------------------------------------------------------
real  SpaceCylinder::volume() const
{
  return 2 * mLength * PI * mRadius * mRadius;
}

//-------------------------------------------------------------------------
bool  SpaceCylinder::isInside( const real w[] ) const
{
  if ( fabs( w[ 0 ] ) > mLength )  return 0;
  return ( w[1]*w[1]+ w[2]*w[2] <= mRadiusSq );
}


//-------------------------------------------------------------------------
void SpaceCylinder::project( const real w[], real p[] ) const
{
  //calculate the projection on the axis, within boundaries:
  
  int inX = 1;
  
  p[ 0 ] = w[ 0 ];
  p[ 1 ] = w[ 1 ];
  p[ 2 ] = w[ 2 ];
  
  if ( w[ 0 ] >  mLength )
    { p[ 0 ] =  mLength; inX = 0; }
  else
    if ( w[ 0 ] < -mLength )
      { p[ 0 ] = -mLength; inX = 0; }
  
  real n = sqrt( w[1]*w[1]+ w[2]*w[2] );
  int inYZ = ( n <= mRadius );
  
  switch( inX + 2 * inYZ ) {
    case 0:
    case 1:
      n = mRadius / n;
      p[ 1 ] = n * w[ 1 ];
      p[ 2 ] = n * w[ 2 ];
      break;
      
    case 2:
      break;
      
    case 3:
      if ( mLength - fabs( w[ 0 ] ) < mRadius - n ) {
        if ( w[ 0 ] > 0 ) 
          p[ 0 ] =  mLength;
        else
          p[ 0 ] = -mLength;
      } else {
        n = mRadius / n;
        p[ 1 ] = n * w[ 1 ];
        p[ 2 ] = n * w[ 2 ];
      }
      break;
  }
}

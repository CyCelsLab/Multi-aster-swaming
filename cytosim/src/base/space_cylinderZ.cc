//RCS: $Id: space_cylinderZ.cc,v 2.0 2004/08/16 16:38:18 nedelec Exp $

#include "space_cylinderZ.h"

//-------------------------------------------------------------------------
real  SpaceCylinderZ::volume() const
{
  return 2 * mLength * PI * mRadius * mRadius;
}


//-------------------------------------------------------------------------
bool  SpaceCylinderZ::isInside( const real point[] ) const
{
  if ( fabs( point[ 2 ] ) > mLength )  return 0;
  return ( point[0]*point[0] + point[1]*point[1] <= mRadiusSq );
}



//-------------------------------------------------------------------------
void SpaceCylinderZ::project( const real w[], real p[] ) const
{
  //calculate the projection on the axis, within boundaries:
  
  int inZ = 1;
  
  p[ 0 ] = w[ 0 ];
  p[ 1 ] = w[ 1 ];
  p[ 2 ] = w[ 2 ];
  
  if ( w[ 2 ] >  mLength ) { p[ 2 ] =  mLength; inZ = 0; }
  else
    if ( w[ 2 ] < -mLength ) { p[ 2 ] = -mLength; inZ = 0; }
  
  real n = sqrt( w[0]*w[0]+ w[1]*w[1] );
  int inXY = ( n < mRadius );
  
  switch( inZ + 2 * inXY ) {
    case 0:
    case 1:
      n = mRadius / n;
      p[ 0 ] = n * w[ 0 ];
      p[ 1 ] = n * w[ 1 ];
      break;
      
    case 2:
      break;
      
    case 3:
      if ( mLength - fabs( w[ 2 ] ) < mRadius - n )
        {
        if ( w[ 2 ] > 0 ) 
          p[ 2 ] =  mLength;
        else
          p[ 2 ] = -mLength;
        }
      else
        {
        n = mRadius / n;
        p[ 0 ] = n * w[ 0 ];
        p[ 1 ] = n * w[ 1 ];
        }
      break;
  }
}

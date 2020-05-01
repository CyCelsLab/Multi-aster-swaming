//RCS: $Id: space_sphere.cc,v 2.0 2004/08/16 16:39:07 nedelec Exp $


#include "space_sphere.h"

#if (DIM == 1)

real SpaceSphere::volume() const
{
  return 2 * mRadius;
}

bool SpaceSphere::isInside( const real point[] ) const
{
  return point[0] * point[0] <= mRadiusSq;
}

void SpaceSphere::project( const real point[], real proj[] ) const
{
  real n = point[0] * point[0];
  
  if ( n > 0 ) 
    n = mRadius / sqrt( n );
  
  proj[ 0 ] = n * point[ 0 ];
}

#endif


//-------------------------------------------------------------------------


#if (DIM == 2)

real SpaceSphere::volume() const
{
  return PI * mRadius * mRadius;
}

bool SpaceSphere::isInside( const real point[] ) const
{
  return point[0] * point[0] + point[1] * point[1] <= mRadiusSq;
}

void SpaceSphere::project( const real point[], real proj[] ) const
{
  real n = point[0] * point[0] + point[1] * point[1];
  
  if ( n > 0 ) 
    n = mRadius / sqrt( n );
  
  proj[ 0 ] = n * point[ 0 ];
  proj[ 1 ] = n * point[ 1 ];
}

#endif


//-------------------------------------------------------------------------

#if (DIM == 3)

real SpaceSphere::volume() const
{
  return 4/3.0 * PI * mRadius * mRadius * mRadius;
}

bool SpaceSphere::isInside( const real point[] ) const
{
  return point[0] * point[0] + point[1] * point[1] + point[2] * point[2] <= mRadiusSq;
}

void SpaceSphere::project( const real point[], real proj[] ) const
{
  real n = point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
  
  if ( n > 0 ) 
    n = mRadius / sqrt( n );
  
  proj[ 0 ] = n * point[ 0 ];
  proj[ 1 ] = n * point[ 1 ];
  proj[ 2 ] = n * point[ 2 ];
}

#endif

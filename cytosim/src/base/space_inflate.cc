//RCS: $Id: space_inflate.cc,v 1.1 2005/03/30 18:16:56 nedelec Exp $

#include "space_inflate.h"


//--------------------------------------------------------------------
Vecteur     SpaceInflate::getBoundingRect() const
{ 
  Vecteur result = mSpace->getBoundingRect();
  for( int dd = 0; dd < DIM; ++dd ) 
    result[dd] += mRadius;
  return result;
}



//--------------------------------------------------------------------
bool        SpaceInflate::isInside(const real point[]) const
{
  if ( mSpace -> isInside( point ) ) 
    return true;
  
  real proj[DIM];
  mSpace -> project( point, proj );
  
  real n = 0;
  for(int d = 0; d < DIM; ++d)
    n += sqr( point[d] - proj[d] );
  
  return ( n <= mRadiusSq );
}





//--------------------------------------------------------------------
void        SpaceInflate::project(const real point[], real proj[]) const
{
  mSpace -> project(point, proj);
  
  real n = 0, pw[DIM];
  for(int d = 0; d < DIM; ++d) {
    pw[d] = point[d] - proj[d];
    n += pw[d] * pw[d];
  }
  
  //TODO: problem in project-smooth if point is exactly on the box (n==0)
  //if (n==0) we do not know the orthogonal direction to follow. We should
  //take a different position not too far, and project from there.
  
  if ( n > 0 )
    n = ( mSpace->isInside(point) ? -1 : +1 ) * mRadius / sqrt(n);

  for(int d=0; d<DIM; ++d)
    proj[d] += n * pw[d];
}

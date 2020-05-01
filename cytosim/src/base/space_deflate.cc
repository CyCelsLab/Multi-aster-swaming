//RCS: $Id: space_deflate.cc,v 1.1 2005/03/30 18:16:36 nedelec Exp $

#include "space_deflate.h"


//--------------------------------------------------------------------
Vecteur     SpaceDeflate::getBoundingRect() const
{ 
  Vecteur result = mSpace->getBoundingRect();
  for( int dd = 0; dd < DIM; ++dd ) 
    result[dd] -= mRadius;
  return result;
}



//--------------------------------------------------------------------
bool        SpaceDeflate::isInside(const real point[]) const
{
  if ( mSpace -> isInside( point ) ) {
  
    real proj[DIM];
    mSpace -> project( point, proj );
    
    real n = 0;
    for(int d = 0; d < DIM; ++d)
      n += sqr( point[d] - proj[d] );
    
    return ( n >= mRadiusSq );
    
  }
  //case where point was outside:
  return false;
}





//--------------------------------------------------------------------
void        SpaceDeflate::project(const real point[], real proj[]) const
{
  mSpace -> project(point, proj);
  
  real n = 0, pw[DIM];
  for(int d = 0; d < DIM; ++d) {
    pw[d] = point[d] - proj[d];
    n += pw[d] * pw[d];
  }
  
  //TODO: problem in project-smooth if point is exactly on the box (n==0)
  //if (n==0) we do not know the orthogonal direction to follow. We should
  //take another point near by, and project from there.
  if ( n > 0 )
    n = ( mSpace->isInside(point) ? +1 : -1 ) * mRadius / sqrt(n);
  
  for(int d=0; d<DIM; ++d)
    proj[d] += n * pw[d];
  
}

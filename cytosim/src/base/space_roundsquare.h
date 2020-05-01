//RCS: $Id: space_roundsquare.h,v 2.2 2005/04/01 10:17:10 nedelec Exp $

#ifndef SPACE_ROUNDSQUARE_H
#define SPACE_ROUNDSQUARE_H

#include "space.h"

/// square with smooth edges: the cube defined by mSize[0:2] extended by size[3]
class SpaceRoundSquare : public Space
{
private:
  
  ///the three sides of the rectangle
  real mSize[3];
  
  ///the size of the extension
  real mRadius;
  
  ///the square of the radius
  real mRadiusSq;
  
public:
  
  ///creator
  SpaceRoundSquare(const real size_set[6]) {
    //the rounding radius should be positive:
    assert( size_set[3] >= 0 ); 
    mRadius   = fabs( size_set[3] );
    mRadiusSq = mRadius * mRadius;

    for( int dd = 0; dd < DIM; ++dd ) {
      mSize[dd] = size_set[dd];
      //each dimension should be at least the rounding radius:
      if ( mSize[dd] < mRadius )
        mSize[dd] = mRadius;
    }
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_ROUND_SQUARE; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( mSize[0]+mRadius, mSize[1]+mRadius, mSize[2]+mRadius ); }
  
  /// the volume of space inside
  //real        volume()          const;
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
};

#endif

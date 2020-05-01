//RCS: $Id: space_square.h,v 2.0 2004/08/16 15:57:30 nedelec Exp $

#ifndef SPACE_SQUARE_H
#define SPACE_SQUARE_H

#include "space.h"

///a rectangular convex region in space
class SpaceSquare : public Space
{
private:
  
  ///the three sides of the rectangle
  real mSize[3];
  
public:
  
  ///creator
  SpaceSquare(const real size_set[6]) {
    for( int dd = 0; dd < DIM; ++dd )
      mSize[dd] = size_set[dd];
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_SQUARE; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( mSize[0], mSize[1], mSize[2] ); }
  
  /// the volume of space inside
  real        volume()          const;
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
};

#endif


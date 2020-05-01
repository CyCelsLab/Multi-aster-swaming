//RCS: $Id: space_inflate.h,v 1.2 2005/04/22 11:56:52 nedelec Exp $
#ifndef SPACE_INFLATE_H
#define SPACE_INFLATE_H

#include "space.h"

///a Space inflated in all directions, by a uniform length (for convex shapes)
//the method works if the original shape is convex, otherwise no guaranty
class SpaceInflate : public Space
{
private:
  
  ///original space
  const Space * mSpace;
  
  ///length by which the space is extended
  real mRadius;
  
  ///square of the length (optimization, avoids redundant recalculations)
  real mRadiusSq;
  
public:

  ///creator
  SpaceInflate(const Space * space, const real length) {
    assert( space->isConfined() );
    assert( length > 0 );
    mSpace    = space;
    mRadius   = length;
    mRadiusSq = length * length;
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return mSpace->getShape(); }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const;
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
};

#endif

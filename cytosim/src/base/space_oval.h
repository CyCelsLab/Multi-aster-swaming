//RCS: $Id: space_oval.h,v 2.0 2004/08/16 15:57:15 nedelec Exp $

#ifndef SPACE_OVAL_H
#define SPACE_OVAL_H

#include "space.h"

///an oval is a cylinder with two half-sphere caps
class SpaceOval : public Space
{
private:
  
  ///the length of the central cylinder
  real mLength;
  
  ///the radius of the sphere
  real mRadius;
  
  //the square of the radius
  real mRadiusSq;
  
public:
    
  ///creator
  SpaceOval(const real size_set[6]) {
    mLength   = size_set[0];
    mRadius   = size_set[1];
    mRadiusSq = mRadius * mRadius;
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_OVAL; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( mRadius+mLength, mRadius, mRadius ); }
  
  /// the volume of space inside
  real        volume()          const;
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
};

#endif

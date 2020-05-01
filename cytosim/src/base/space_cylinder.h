//RCS: $Id: space_cylinder.h,v 2.0 2004/08/16 16:38:15 nedelec Exp $

#ifndef SPACE_CYLINDER_H
#define SPACE_CYLINDER_H

#include "space.h"

///a cylinder of axis X
class SpaceCylinder : public Space
{
private:
  
  ///the length of the cylinder, in z
  real mLength;
  
  ///the radius
  real mRadius;
  
  ///square of the radius
  real mRadiusSq;
  
public:
  
  ///creator
  SpaceCylinder(const real size_set[6]) {
    mLength   = size_set[0];
    mRadius   = size_set[1];
    mRadiusSq = mRadius * mRadius;
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_CYLINDER; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( mLength, mRadius, mRadius ); }
  
  /// the volume of space inside
  real        volume()          const;
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
};

#endif


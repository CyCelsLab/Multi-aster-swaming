//RCS: $Id: space_sphere.h,v 2.1 2004/08/18 12:49:39 nedelec Exp $

#ifndef SPACE_SPHERE_H
#define SPACE_SPHERE_H

#include "space.h"

///sphere centered at the origin
class SpaceSphere : public Space
{
private:
  
  ///the radius of the sphere
  real mRadius;
  
  //the square of the radius
  real mRadiusSq;
  
public:
  
  ///creator
  SpaceSphere(const real size_set[6]) {
      mRadius   = size_set[0];
      mRadiusSq = mRadius * mRadius;
  }
  
  ///creator with radius
  SpaceSphere(const real radius) {
    mRadius   = radius;
    mRadiusSq = mRadius * mRadius;
  }
  
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_SPHERE; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( mRadius, mRadius, mRadius ); }
  
  /// the volume of space inside
  real        volume()          const;

  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;

};

#endif


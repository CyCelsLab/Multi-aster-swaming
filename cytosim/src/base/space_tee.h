//
// C++ Interface: space_tee
//
// Description: T shape for pombe
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: space_tee.h,v 1.3 2005/03/31 08:56:00 foethke Exp $


#ifndef SPACE_TEE_H
#define SPACE_TEE_H

#include "smath.h"
#include "space.h"
#include "iomessages.h"


//a T is an oval with an extension perpendicular to the symmetry axis
class SpaceTee : public Space
{
private:
  
  ///the length of the central cylinder
  real tLength;
  
  ///the length of the perpendicular part on the cylinder
  real tArmLength;
  
  ///the position of the perpendicular part on the cylinder
  real tJunction;
  
  ///the radius of the caps
  real tCapRadius;
  
  //the square of the radius
  real tCapRadiusSq;
  
  ///the x coordinate of a point relative to tJunction
  mutable real xRel;
  
  ///some angles of a point in space
  mutable real xzAngle, gamma, delta;
  
  ///calculate the mutable vars needed in isInside() and project()
  void        setMutables( const real w[] )     const;
  
public:
    
  ///costructor
  SpaceTee(const real size_set[6]) {
    tLength      = size_set[0];
    tCapRadius   = size_set[1];
    tJunction    = size_set[2];
    tArmLength   = size_set[3];
    tCapRadiusSq = tCapRadius * tCapRadius;
  
    if( (tLength <= 0) || (tArmLength < 0) || (tCapRadius <= 0) ) {
      MSG.error("SpaceTee","The T-shape can't have negative length, arm length or radius.");
      exit(0);
    }
    if( (fabs(tJunction)+tCapRadius) > tLength ) {
      MSG.error("SpaceTee","The position of the branch plus the radius must lie within the length of the T.");
      exit(0);
    }
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_TEE; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( tCapRadius+tLength, tArmLength+2.*tCapRadius, tCapRadius ); }
  
  /// the volume of space inside
  real        volume()          const;
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
};

#endif //SPACE_TEE_H

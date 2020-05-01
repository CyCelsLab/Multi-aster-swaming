//
// C++ Interface: space_boomerang
//
// Description: Can be (re)used to knock out yeast
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: space_boomerang.h,v 1.2 2005/03/17 18:00:05 foethke Exp $


#ifndef SPACE_BOOMERANG_H
#define SPACE_BOOMERANG_H

#include "smath.h"
#include "space.h"
#include "iomessages.h"

///a boomerang with two half-sphere caps
class SpaceBoomerang : public Space
{
private:
  
  ///half the length of the central part of the boomerang
  real bLength;
  
  ///the angle between the arms of the boomerang
  real bAngle;
  
  ///the radius of the boomerang caps
  real bCapRadius;
    
  ///the square of the cap radius
  real bCapRadiusSq;
  
public:
    
  ///constructor
  SpaceBoomerang(const real size_set[6]) {
    bLength      = size_set[0];
    bCapRadius   = size_set[1];
    bAngle       = size_set[2];
    if( (bLength <= 0) || (bCapRadius <= 0) || (bAngle <=0) ) {
      MSG.error("SpaceBoomerang","Boomerangs usually don't have negative length, width or angle.");
      exit(0);
    }
    if( (bLength < tan(bAngle)*bCapRadius) || (bAngle >= PI/2.) ) {
      MSG.error("SpaceBoomerang","An angle that big would make a knot in the boomerang.");
      exit(0);
    }
    bCapRadiusSq = bCapRadius * bCapRadius;
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_BOOMERANG; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( bLength*cos(bAngle)+bCapRadius, bLength*sin(bAngle)+bCapRadius, bCapRadius ); }
  
  /// the volume of space inside
  real        volume()          const;
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
};

#endif // SPACE_BOOMERANG_H

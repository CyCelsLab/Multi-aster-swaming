//
// C++ Interface: space_banana
//
// Description: Yellow and tasty
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: space_banana.h,v 1.3 2005/03/17 14:44:36 foethke Exp $


#ifndef SPACE_BANANA_H
#define SPACE_BANANA_H

#include "smath.h"
#include "space.h"
#include "iomessages.h"

///a banana with two half-sphere caps
class SpaceBanana : public Space
{
private:
  
  ///half the length of the central part of the banana (measured in the middle)
  real bLength;
  
  ///half the normalized arc length of the banana (between 0 and PI)
  real bArcLength;
  
  ///radius of the banana (measured in the middle). Why are bananas curved anyway?
  real bRadius;
  
  ///the radius of the banana caps
  real bCapRadius;
  
  ///the square of the radius
  real bRadiusSq;
  
  ///the square of the cap radius
  real bCapRadiusSq;
  
public:
    
  ///constructor
  SpaceBanana(const real size_set[6]) {
    bLength      = size_set[0];
    bCapRadius   = size_set[1];
    bArcLength   = size_set[2];
    bRadius      = bLength/bArcLength;
    
    //check parameters for sanity
    if( (bLength <= 0) || (bCapRadius <= 0) || (bArcLength <=0) ) {
      MSG.error("SpaceBanana","Bananas usually don't have negative length, width or curvature.");
      exit(0);
    }   
    if( (bRadius - bCapRadius) <= 0 ) {
      MSG.error("SpaceBanana","This banana is too fat, increase the length and/or decrease the arc length.");
      exit(0);
    }
    if( bArcLength > (PI - asin(bCapRadius/bRadius)) ) {
      MSG.error("SpaceBanana","This banana folds on it self, decrease the arc length and/or the cap radius.");
      exit(0);
    }
    
    bRadiusSq    = bRadius    * bRadius;
    bCapRadiusSq = bCapRadius * bCapRadius;
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_BANANA; }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const;
  
  /// the volume of space inside
  real        volume()          const;
    
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
};

#endif // SPACE_BANANA_H

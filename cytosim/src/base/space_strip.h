//RCS: $Id: space_strip.h,v 2.1 2004/12/07 17:41:41 nedelec Exp $

#ifndef SPACE_STRIP_H
#define SPACE_STRIP_H

#include "space.h"

///a rectangular Space with partial periodic boundary conditions
class SpaceStrip : public Space
{
private:
  
  ///the three sides of the rectangle
  real mSize[3];
  
  ///the double of the sizes
  real mSize2[3];
  
public:
    
  ///creator
  SpaceStrip(const real size_set[6]) {
    for( int dd = 0; dd < DIM; ++dd ) {
      mSize[dd]  = size_set[dd];
      mSize2[dd] = 2 * mSize[dd];
    }
    setPeriodic( size_set, 2 );
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_STRIP; }
  
  /// return true if the space has some borders
  bool        isPeriodic()      const { return true; }
    
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


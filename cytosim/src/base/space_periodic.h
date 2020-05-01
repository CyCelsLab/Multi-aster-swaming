
//RCS: $Id: space_periodic.h,v 2.0 2004/08/16 15:57:19 nedelec Exp $

#ifndef SPACE_PERIODIC_H
#define SPACE_PERIODIC_H

#include "space.h"

///a rectangular region in space with periodic boundary conditions
class SpacePeriodic : public Space
{
private:
  
  ///the three sides of the rectangle
  real mSize[3];
  
  ///the double of the sizes
  real mSize2[3];
  
public:
  
  ///creator
  SpacePeriodic(const real size_set[6]) {
    for( int dd = 0; dd < DIM; ++dd ) {
      mSize[dd]  = size_set[dd];
      mSize2[dd] = 2 * mSize[dd];
    }
    setPeriodic( size_set, 1 );
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return SHAPE_PERIODIC; }
  
  /// return true if the space has some borders
  bool        isPeriodic()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return Vecteur( mSize[0], mSize[1], mSize[2] ); }
  
  /// the volume of space inside
  real        volume()          const;
    
};

#endif


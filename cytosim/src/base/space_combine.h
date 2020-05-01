//RCS: $Id: space_combine.h,v 2.1 2004/12/07 17:41:23 nedelec Exp $

#ifndef SPACE_COMBINE_H
#define SPACE_COMBINE_H

#include "space.h"

///SpaceCombine removes an inner Space from another outer Space
class SpaceCombine : public Space
{
private:
  
  ///Space in which objects are contained
  Space * outerSpace;
  
  ///Space in which objects are excluded
  Space * innerSpace;
  
public:
    
  ///creator
  SpaceCombine(Space * space_big, Space * space_small) {
    outerSpace    = space_big;
    innerSpace    = space_small;
    assert( innerSpace -> isConfined() );
  }
  
  /// returns the unique shape Id
  SHAPE_LIST  getShape()        const { return outerSpace->getShape(); }
  
  /// return true if the space has some borders
  bool        isConfined()      const { return true; }
  
  /// returns a rectangle surrounding the space
  Vecteur     getBoundingRect() const { return outerSpace->getBoundingRect(); }
  
  /// true if the point is inside the space
  bool        isInside(const real point[]) const;
  
  /// project point on the closest edge of the space
  void        project(const real point[], real proj[]) const;
  
  ///set point to its periodic representation closest to origin by removing periodic repeats
  void        modulo(real point[]) const { outerSpace->modulo(point); }

};

#endif

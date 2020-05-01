//RCS: $Id: space.h,v 2.12 2005/04/22 13:20:31 nedelec Exp $
//----------------------------------space.h-----------------------------
//  class Space   defines the simulated physical space,
//             its shape, size and an associated grid for divide&conquer
//----------------------------------------------------------------------

#ifndef SPACE_H
#define SPACE_H

//--------------------possible shapes of the box------------------------

enum SHAPE_LIST {
  SHAPE_UNSET     = 0,   SHAPE_SQUARE         = 1,   SHAPE_PERIODIC   = 2,
  SHAPE_SPHERE    = 3,   SHAPE_ELLIPSE        = 4,   SHAPE_BANANA     = 5,
  SHAPE_CYLINDER  = 6,   SHAPE_ROUND_SQUARE   = 7,   SHAPE_OVAL       = 8,
  SHAPE_STRIP     = 9,   SHAPE_OUTSIDE_SPHERE = 10,  SHAPE_CYLINDER_Z = 11,
  SHAPE_BOOMERANG = 12,  SHAPE_TEE            = 13 };

//----------------------------------------------------------------------

#include "vecteur.h"
#include "smath.h"

///Space holds the description of the physical space for cytosim: Derived class include square, sphere, etc.
/** by default, the space is infinite, i.e. has no confinement */
class Space
{
 private:

  ///for the periodic boundary conditions
  static real mPeriod2[3];

  ///for the semi-periodic boundary conditions
  static int  mPeriodic;

 public:

  ///dummy constructor
  Space() {}

  ///destructor
  virtual ~Space() {}

  /// return an int identifying the shape uniquely
  virtual SHAPE_LIST    getShape()                               const { return SHAPE_UNSET; }

  /// true if the space has borders
  virtual    bool       isConfined()                             const { return false; }

  ///returns a rectangle containing the space
  virtual    Vecteur    getBoundingRect()                        const { return Vecteur(10, 10, 10); }

  ///test if w is inside the space, returning false or true, borders are included
  virtual    bool       isInside(const real w[])                 const { return true; }

  ///project(w, p) sets p to the closest point to w on the edge of the box
  virtual    void       project(const real point[], real proj[]) const {
    for( int dd = 0; dd < DIM; ++dd )
      proj[dd] = 0;
  }

  //TODO: a function returning a linearization of the confining force for a box

  //-----------------------------------------------------------------------------
  //the functions below are derived from the ones above:

  ///the volume of the box, calculated via monte-carlo algorithm, but can be overriden
  virtual real       volume()                                    const;

  /// return true if the shape is set
  bool       isSet()                                             const { return getShape() != SHAPE_UNSET; }

  ///true if a sphere (center w, radius) fits in the space, borders included
  bool       isInside(const real center[], const real radius)    const;

  ///project point[] on the space, deflated by radius (should be positive)
  void       project(const real point[], real proj[], const real radius) const;

  ///project, overriding the given vecteur (for convenience)
  void       project(real [])                                    const;

  ///a random place in the volume
  Vecteur    randomPlaceInVolume()                               const;

  ///a random place on the edge
  Vecteur    randomPlaceOnEdge(real radius)                      const;

  ///true if point w is outside the space ( == ! isInside() )
  bool       isOutside(const real w[])                           const { return ! isInside(w); }

  ///true if a sphere (center w, radius) has a portion outside ( == ! isInside() )
  bool       isOutside(const real center[], const real radius)   const { return ! isInside(center, radius); }

  //-----------------PERIODIC BOUNDARY CONDITIONS-------------------
  //all functions are static: only one period for all spaces!

  /// used to set the periodic boundary conditions
  static void       setPeriodic( const real size[], const int semi=0 );

  /// true if the space has periodic boundaries
  static int        isPeriodic()  { return mPeriodic; }

  ///set x to its periodic representation closest to origin by removing periodic repeats
  static void       modulo( real x[] );

  ///calculate the representation of x (1st arg) remove periodic repeats, with o (3rd arg) as origin
  static void       moduloNear( real x[], const real o[] );

  ///set m (1st arg) to be the representation of x (2d arg)  without periodic repeats, with o (3rd arg) as origin
  static void       moduloNear( real m[], const real x[], const real o[] );

  ///calculate the periodic offset to go from y to its representation
  static void       offset( real x[] );

  ///calculate the periodic offset between x and y
  static void       moduloAndOffset( real x[], real off[] );
};


//we provide a gloal "creator" returning Space of any shape:

/// New Space of specific shape and dimensions
Space * newSpace( const int shape, const real dim[6] );

/// New Space with inflation (if >0) or deflation (if < 0) radius
Space * newSpace( const int shape, const real dim[6], const real radius );


#endif // _space_h





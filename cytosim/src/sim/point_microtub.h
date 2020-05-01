//RCS: $Id: point_microtub.h,v 2.11 2005/04/22 11:50:34 nedelec Exp $

//---------------------------------------------------------------------------
//====== The class PointMicrotub defines a position on a microtubule  =======
//---------------------------------------------------------------------------

//it is used in many occasion:
//  - to calculate the interpolated positions on microtubules
//  - in Hand, when the rod / position is needed
//  - to represent rods in the dispatch algorithm for hand attachments

#ifndef POINT_MICROTUB_H
#define POINT_MICROTUB_H

//#define USE_LATTICE

#include "point_interpolated.h"
#include "vecteur.h"

class Microtub;

//Warning: do not change the values for these enums, MT_NOT_END should be zero:
///Enum to specify MT ends, and the origin
enum MTEnd {  
  MT_NOT_END   = 0,   ///< not an end
  MT_MINUS_END = 1,   ///< microtubule minus-end = point 0
  MT_PLUS_END  = 2,   ///< microtubule plus-end = last point
  MT_ORIGIN    = 3    ///< origin of abscissa
};

///\todo Motors should be built on Fiber, not on Microtub
///\todo PointMicrotub should be renamed PositionOnFiber
///represents a position on a Microtub, based on an abscissa
class PointMicrotub : public PointInterpolated
{
  
protected:

  ///flag to see if lattice site on MT has been set
  bool               latticeSiteSet;

  //the hand position at last update of the lattice - needed to clear the lattice before moving
  real               prevHaAbs;
 
  ///the abscissa from the Microtub origin
  real               haAbs;  
  
  ///flag of value MT_PLUS_END or MT_MINUS_END when at end of microtubule
  MTEnd              haEnd;              
   
public:
   
#ifdef TEST_ATTACH
  ///debug field to test attachements
  short              attachCnt;
#endif
  
  ///reset the member variables
  void           clear();
    
  ///the default creator
  PointMicrotub() { clear(); } //FIX the default constructor does not need to reset the variables...
  
  //partial constructor for the attachements
  PointMicrotub( const Microtub * mt, int srod );
  
  ///constructor from abscissa
  PointMicrotub( const Microtub * mt, const real ab, const MTEnd from = MT_ORIGIN ) {
    bindTo( mt, ab, from );
  }
  
  ///empty destructor
  virtual ~PointMicrotub() { }
  
  ///the microtubule to which is attached
  Microtub *     getMT()         const { return (Microtub*)(mPS); }
     
  ///local tangent vector of the Microtub
  Vecteur        dirMT() const;

  ///check haAbs with respect to the boundaries of the Microtub, sets haEnd
  MTEnd          checkBoundaries();
  
  ///update the PointInterpolated, call if haAbs has changed
  void           updateInterpolated();  

  ///sets haAbs from the position described by the PointInterpolated
  void           setFromInterpolated();
  
  ///the distance from the origin of the Microtub 
  real           abscissa() const { return haAbs; }
  
  ///the distance from the specified end/origin of the Microtub
  real           abscissa( const MTEnd from ) const;

  ///this is called after each movement on and attachment to MT
  void           updateLattice(void);

  ///called at detachment from MT
  void           freeFromLattice(void);

  ///Returns flag telling whether it is at the end of a microtubule
  MTEnd          getEnd() const { return haEnd; }

  ///Attaches to a new microtubule
  void           bindTo(const Microtub *, const real ab, const MTEnd from = MT_ORIGIN);
  
  ///attaches to site haLatticeSite on the microtubule lattice
  void           bind(void);

  ///attaches to a different point on the same microtubule
  void           moveTo(const real ab);
  
  ///equivalent to moveTo( haAbs + ab );
  void           moveBy(const real ab)     { moveTo( haAbs + ab ); }
  
  ///Move to a position on a microtubule to which hand is already attached
  void           moveTo(const real ab, const MTEnd from);
  
  //================= function used to prevent double binding to the same rod:

  ///returns the MT segments to which the motor is attached, or 0 if the motor is not attached
  const PointMicrotub * getMTRod() const;
    
  ///calculates the distance of the rod to the point w
  real           distanceToRod( const Vecteur & w );
  
  ///returns true if the given hand binding key prevents does not allow binding
  bool           keyMatch(const BindingKey key) const;
  
  ///debug
  int            looksWrong() const;
};


#endif

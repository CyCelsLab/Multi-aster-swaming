//RCS: $Id: point_interpolated.h,v 2.6 2005/04/10 14:54:14 nedelec Exp $

//---------------------------------------------------------------------------
//======     PointInterpolated defines a position on a microtubule    =======
//---------------------------------------------------------------------------

//it is used in many occasion:
//  - to calculate the interpolated positions on microtubules
//  - in Hand, when the rod / position is needed
//  - to represent rods in the dispatch algorithm for hand attachments

#ifndef POINT_INTERPOLATED_H
#define POINT_INTERPOLATED_H

#include "types.h"
#include "vecteur.h"
#include "simobject.h"

//we use an incomplete declaration:
class PointExact;

///represents an intermediate position between two points on the same SimObject, see also PointExact
class PointInterpolated 
{
  friend class Fiber;
  
protected:

  ///the pointset whose points are being interpolated 
  const SimObject * mPS;

  ///the index of the first point
  int               mPoint1;
 
  ///the index of the second point
  int               mPoint2;
    
  ///the coefficient on the second point, the weight on the first point is 1-coef
  real              mCoef;  ///< mCoef should be between 0 and 1
  
public:

  ///the constructor resets the member variables
  PointInterpolated() { clear(); }
  
  ///constructor
  PointInterpolated( const SimObject * ps, const int srod1, const int srod2, const real c ) {
    mPS     = ps;
    mPoint1 = srod1;
    mPoint2 = srod2;
    mCoef   = c;
  }
  
  ///set to a given position
  void setTo( const SimObject * ps, const int srod1, const int srod2, const real c ) {
    mPS     = ps;
    mPoint1 = srod1;
    mPoint2 = srod2;
    mCoef   = c;
  }
  
  ///resets member variables
  void clear();
  
  ///empty destructor
  virtual ~PointInterpolated()    {}
                       
  ///returns the SimObject as a constant pointer
  const SimObject * getPS() const { return mPS; }

  ///returns the first point index
  int           getPoint1() const { return mPoint1; }

  ///returns the second point index
  int           getPoint2() const { return mPoint2; }

  ///returns the coefficient
  real            getCoef() const { return mCoef; }
  
  ///returns true if the two PointInterpolated share one or more point
  bool near(const PointExact & ip ) const;

  ///returns true if the two PointInterpolated share one or more point
  bool near(const PointInterpolated & ip ) const;

  ///true if getMT() == 0, ie. no SimObject is referenced, or 'attached'
  bool    isFree()        const { return ( mPS == 0 ); }
  
  ///true if getMT() != 0, ie. a SimObject is 'attached'
  bool    isAttached()    const { return ( mPS != 0 ); }
  
  ///the position in space
  Vecteur where() const;
  
  ///the unit vector between the points 1 and 2
  Vecteur dirInter() const;
  
  ///the distance between points 1 and 2
  real   length() const;
  
  ///the index in the big matrix of point 1
  inline int  matIndex1() const    { return mPS->matIndex() + mPoint1; }
  
  ///the index in the big matrix of point 2
  inline int  matIndex2() const    { return mPS->matIndex() + mPoint2; }
  
  ///debug: should return NO_ERROR
  int looksWrong(const bool forInteraction = false) const;
};

#endif

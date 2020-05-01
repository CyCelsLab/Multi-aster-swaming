//RCS: $Id: point_exact.h,v 2.6 2005/04/10 14:53:07 nedelec Exp $

//---------------------------------------------------------------------------
//======  accessory class PointExact, defines a point on a PointSet   =======
//---------------------------------------------------------------------------

//it is used in many occasion:
//  - to calculate the interpolated positions on microtubules
//  - in Hand, when the rod / position is needed
//  - to represent rods in the dispatch algorithm for hand attachments

#ifndef POINT_EXACT_H
#define POINT_EXACT_H

#include "types.h"
#include "vecteur.h"
#include "simobject.h"

///a pointer to a particular point of a PointSet, see also PointInterpolated
class PointExact 
{
protected:

  ///the SimObject containing the point of interest 
  const SimObject * mPS;
  
  ///the index of the point in the SimObject
  int              mPoint;
  
public:

  ///we reset the variable for the default creator
  PointExact() { clear(); }
  
  ///constructor
  PointExact( const SimObject * ps, const int rd ) {
    mPS  = ps;
    mPoint = rd;
  }
  
  ///Set the PointExact to a particular value
  void    setTo( const SimObject * ps, const int rd ) {
    mPS  = ps;
    mPoint = rd;
  }
  
  ///reset the member variables to zero
  void clear();
  
  ///destructor does nothing
  virtual ~PointExact() {}

  ///return the SimObject as a const pointer
  const SimObject * getPS() const { return mPS; }
  
  ///returns the first point index
  int            getPoint() const { return mPoint; }
  
  //true if the this other PointExact is pointing to the same position
  bool         near( const PointExact & pt ) const { 
    return (( mPS == pt.mPS ) && ( mPoint == pt.mPoint ));
  }
  
  ///the position
  Vecteur where() const;
  
  ///the index in the big matrix
  inline int  matIndex() const { return mPS->matIndex() + mPoint; }
  
  ///debug: should return NO_ERROR
  int looksWrong(const bool forInteraction = false) const;
};


#endif

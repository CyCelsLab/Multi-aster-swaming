//RCS: $Id: point_microtub.cc,v 2.15 2005/03/01 10:24:55 clausen Exp $

#include "point_microtub.h"
#include "microtub.h"

//------------------------------------------------------------------------
PointMicrotub::PointMicrotub( const Microtub * mt, const int srod ) 
: PointInterpolated( mt, srod, srod+1, 0)
{
}


//------------------------------------------------------------------------
void PointMicrotub::clear()
{
  PointInterpolated::clear();
  haAbs     = 0;
  prevHaAbs = 0;
  haEnd     = MT_NOT_END;
  latticeSiteSet = false;
}


//------------------------------------------------------------------------
Vecteur PointMicrotub::dirMT() const
{ 
  return getMT() -> dirP( mPoint1 );
}


//------------------------------------------------------------------------
real PointMicrotub::abscissa(const MTEnd from) const
{
  switch( from ) {
    case MT_ORIGIN:    return haAbs;
    case MT_MINUS_END: return haAbs - getMT()->abscissaMinusEnd();
    //in the case of MT_PLUS_END, the argument is measured towards the minus-end
    case MT_PLUS_END:  return getMT()->length() + getMT()->abscissaMinusEnd() - haAbs;
    default: 
      MSG.error("PointMicrotub::abscissa","wrong argument"); return 0;
  }  
}

#ifdef USE_LATTICE
//------------------------------------------------------------------------
/**Takes care of updating getMT() -> lattice. Call after each movement and
attachment to MT. An error is generated downstream from here in in
getMT()->setSiteOccupied( haAbs ) if the lattice site we are trying to find
to is already occupied.
*/
void PointMicrotub::updateLattice(void)
{
  if ( latticeSiteSet ) {
    getMT() -> setSiteFree( prevHaAbs );
  } else {
    latticeSiteSet = true;
  }
  prevHaAbs = haAbs;
  getMT() -> setSiteOccupied( haAbs );
}
#else
void PointMicrotub::updateLattice(void) {}
#endif

//------------------------------------------------------------------------
/**Free lattice point. We use prevHaAbs here as this is the ACTUAL point that
has been set to occupied - haAbs might have been updated by checkBoundaries
without the lattice being updated.
*/
void PointMicrotub::freeFromLattice(void)
{
  getMT() -> setSiteFree( prevHaAbs );
  latticeSiteSet = false;
  prevHaAbs = 0;
}

//------------------------------------------------------------------------
MTEnd PointMicrotub::checkBoundaries()
{
  if ( haAbs <= getMT() -> abscissa(MT_MINUS_END) ) {
    haEnd = MT_MINUS_END;
    haAbs = getMT() -> abscissa(MT_MINUS_END);
    return haEnd;
  }

  if ( haAbs >= getMT() -> abscissa(MT_PLUS_END) ) {
    haEnd = MT_PLUS_END;
    haAbs = getMT() -> abscissa(MT_PLUS_END);
    return haEnd;
  }
  
  haEnd = MT_NOT_END;
  return haEnd;
}


//------------------------------------------------------------------------
//-----------------------------MOVE---------------------------------------
//------------------------------------------------------------------------

void  PointMicrotub::updateInterpolated()
{
  assert ( isAttached() );
  checkBoundaries();
  getMT()->setInterpolation( this, haAbs );
}

//------------------------------------------------------------------------
void  PointMicrotub::setFromInterpolated()
{
  haAbs = getMT() -> abscissa( mPoint1+mCoef );
  checkBoundaries();
}

//------------------------------------------------------------------------
void PointMicrotub::moveTo( const real ab )
{
  assert( isAttached() );

  haAbs = ab;
  updateLattice();
  //updateInterpolated();
}

//------------------------------------------------------------------------
//ab is counted from the given end <from>
void PointMicrotub::moveTo( const real ab, const MTEnd from )
{
  assert( isAttached() );
  haAbs = getMT()->abscissa(ab, from);
  updateLattice();
  updateInterpolated();
}

//------------------------------------------------------------------------
//ab is counted from the given end <from>
void PointMicrotub::bindTo( const Microtub * mt, const real ab, const MTEnd from )
{
  haEnd = mt->setInterpolation( this, ab, from );
  haAbs = mt->abscissa( mPoint1+mCoef );
  updateLattice();
}

//================= function used to prevent double binding to the same rod:

//------------------------------------------------------------------------
const PointMicrotub * PointMicrotub::getMTRod() const
{
  if ( isAttached() ) 
    return getMT() -> getRod( mPoint1 );
  else
    return 0;
}

//------------------------------------------------------------------------
real PointMicrotub::distanceToRod(const Vecteur & w)
{
  assert( isAttached() );
  return getMT() -> distanceToRod(this, w);
}

//------------------------------------------------------------------------
bool PointMicrotub::keyMatch(const BindingKey key) const
{
  //returns TRUE if the motor should not attach!
  return false; //the feature is not used at the moment
}

//------------------------------------------------------------------------
int PointMicrotub::looksWrong() const
{
  if ( mPS == 0 ) return NO_ERROR;
  if ( mPoint2 != mPoint1 + 1 ) return 1;
  if ( haAbs < getMT()->abscissa(MT_MINUS_END) - EPSILON ) return 2;
  if ( haAbs > getMT()->abscissa(MT_PLUS_END) + EPSILON ) return 3;
  return PointInterpolated::looksWrong();
}

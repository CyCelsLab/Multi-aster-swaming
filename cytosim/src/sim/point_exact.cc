//RCS: $Id: point_exact.cc,v 2.6 2005/04/10 14:52:37 nedelec Exp $

#include "point_exact.h"
#include "simobject.h"

void PointExact::clear() 
{
  mPS    = 0;
  mPoint = 0;
}

Vecteur PointExact::where() const
{
  assert( mPS != 0 );
  return mPS -> whereP(mPoint); 
}

int PointExact::looksWrong(const bool forInteraction) const
{
  if ( mPS == 0 )                  return 1;
  if ( mPoint < 0 )                return 2;
  if ( mPoint > mPS->lastPoint() ) return 3;
  if ( forInteraction ) {
    if ( !mPS->isLinked() )        return 4;
    if ( mPS->matIndex() < 0 )     return 5;
  }
  return NO_ERROR;
}

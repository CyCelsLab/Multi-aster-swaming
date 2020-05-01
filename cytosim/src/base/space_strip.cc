//RCS $Id: space_strip.cc,v 2.0 2004/08/16 16:39:16 nedelec Exp $
#include "space_strip.h"





#if (DIM == 1)

real SpaceStrip::volume() const
{
  return 2 * mSize[0];
}

bool  SpaceStrip::isInside( const real point[] ) const
{      
  if ( point[ 0 ] >  mSize[ 0 ] ) return false;
  if ( point[ 0 ] < -mSize[ 0 ] ) return false;
  return true;
}


void SpaceStrip::project( const real point[], real proj[] ) const
{
  if ( point[ 0 ] > 0 )
    proj[ 0 ] =  mSize[ 0 ];
  else
    proj[ 0 ] = -mSize[ 0 ];
}

#endif


//--------------------------------------------------------------------

#if (DIM == 2)

real SpaceStrip::volume() const
{
  return 4 * mSize[0] * mSize[1];
}

bool  SpaceStrip::isInside( const real point[] ) const
{      
  if ( point[ 1 ] >  mSize[ 1 ] ) return false;
  if ( point[ 1 ] < -mSize[ 1 ] ) return false;
  return true;
}


void SpaceStrip::project( const real point[], real proj[] ) const
{
  proj[ 0 ] = point[ 0 ];
  
  if ( point[ 1 ] > 0 )
    proj[ 1 ] =  mSize[ 1 ];
  else
    proj[ 1 ] = -mSize[ 1 ];
}

#endif

//--------------------------------------------------------------------

#if (DIM == 3)

real SpaceStrip::volume() const
{
  return 8 * mSize[0] * mSize[1] * mSize[2];
}

bool  SpaceStrip::isInside( const real point[] ) const
{      
  if ( point[ 2 ] >  mSize[ 2 ] ) return false;
  if ( point[ 2 ] < -mSize[ 2 ] ) return false;
  return true;
}

void SpaceStrip::project( const real point[], real proj[] ) const
{
  proj[ 0 ] = point[ 0 ];
  proj[ 1 ] = point[ 1 ];
  
  if ( point[ 2 ] > 0 )
    proj[ 2 ] =  mSize[ 2 ];
  else
    proj[ 2 ] = -mSize[ 2 ];
}

#endif

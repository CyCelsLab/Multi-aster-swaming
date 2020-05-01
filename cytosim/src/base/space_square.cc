//RCS: $Id: space_square.cc,v 2.0 2004/08/16 16:39:13 nedelec Exp $

#include "space_square.h"

#if (DIM == 1)

real SpaceSquare::volume() const
{
  return 2 * mSize[0];
}

bool  SpaceSquare::isInside( const real w[] ) const
{
  return (  ( w [ 0 ] >= -mSize[ 0 ] )
         && ( w [ 0 ] <=  mSize[ 0 ] ) );
}

void SpaceSquare::project( const real w[], real p[] ) const
{
  if ( w[0] > 0 )
    p[ 0 ] =  mSize[ 0 ];
  else
    p[ 0 ] = -mSize[ 0 ];
}

#endif


//--------------------------------------------------------------------

#if (DIM == 2)

real SpaceSquare::volume() const
{
  return 4 * mSize[0] * mSize[1];
}


bool  SpaceSquare::isInside( const real w[] ) const
{
  return (  ( w [ 0 ] >= -mSize[ 0 ] )
            && ( w [ 0 ] <=  mSize[ 0 ] )
            && ( w [ 1 ] >= -mSize[ 1 ] )
            && ( w [ 1 ] <=  mSize[ 1 ] ) );
}

void SpaceSquare::project( const real w[], real p[] ) const
{
  int out = 0;

  for( int d = 0; d < 2; ++d )
  {
    p[ d ] = w[ d ];
    if ( p[ d ] >  mSize[ d ] ) { p[ d ] =  mSize[ d ]; out = 1; }
    else
    if ( p[ d ] < -mSize[ d ] ) { p[ d ] = -mSize[ d ]; out = 1; }
  }

  if ( out ) return;

  out = (mSize[1] - fabs( w[1] )) < (mSize[0] - fabs( w[0] ));

  if ( w[ out ] > 0 )
    p[ out ] =  mSize[ out ];
  else
    p[ out ] = -mSize[ out ];
}

#endif

//--------------------------------------------------------------------

#if (DIM == 3)

real SpaceSquare::volume() const
{
  return 8 * mSize[0] * mSize[1] * mSize[2];
}



bool  SpaceSquare::isInside( const real w[] ) const
{  
  return (  ( w [ 0 ] >= -mSize[ 0 ] )
            && ( w [ 0 ] <=  mSize[ 0 ] )
            && ( w [ 1 ] >= -mSize[ 1 ] )
            && ( w [ 1 ] <=  mSize[ 1 ] ) 
            && ( w [ 2 ] >= -mSize[ 2 ] )
            && ( w [ 2 ] <=  mSize[ 2 ] ) );
}

void SpaceSquare::project( const real w[], real p[] ) const
{
  int out = 0;

  for(int d = 0; d < 3; ++d )
    {
      p[ d ] = w[ d ];
      if ( p[ d ] >  mSize[ d ] ) { p[ d ] =  mSize[ d ]; out = 1; }
      else
      if ( p[ d ] < -mSize[ d ] ) { p[ d ] = -mSize[ d ]; out = 1; }
    }

  if ( out ) return;

  real l = mSize[ 0 ] - fabs( w[0] );
  
  real u = mSize[ 1 ] - fabs( w[ 1 ] );
  if ( u < l ) { out = 1; l = u; };

  u = mSize[ 2 ] - fabs( w[ 2 ] );
  if ( u < l )  out = 2;

  if ( w[ out ] > 0 )
    p[ out ] =  mSize[ out ];
  else
    p[ out ] = -mSize[ out ];
}

#endif

//RCS: $Id: space_roundsquare.cc,v 2.0 2004/08/16 16:38:36 nedelec Exp $

#include "space_roundsquare.h"



//-------------------------------------------------------------------------
bool  SpaceRoundSquare::isInside( const real w[] ) const
{
  int d;
  real sw[ 3 ];
  
  for( d = 0; d < DIM; ++d ) {
    sw[ d ] = fabs( w[ d ] );      
    if ( sw[ d ] > mSize[ d ] ) 
      return 0;
  }
  
  real n = 0;
  for( d = 0; d < DIM; ++d ) {
    sw[ d ] += mRadius - mSize[ d ];
    if ( sw[ d ] > 0 ) 
      n += sqr( sw[ d ] );
  }
  
  return ( n <= mRadiusSq );
}



//-------------------------------------------------------------------------
void SpaceRoundSquare::project( const real w[], real p[] ) const
{
  int deepin = 1;
  
  //calculate the projection on the inner cube disminished from mRadius
  for( int d = 0; d < DIM; ++d ) {
    p[ d ] = w[ d ];
    
    if ( p[ d ] >  mSize[ d ] - mRadius ) {
      p[ d ] = mSize[ d ] - mRadius;
      deepin = 0;
    }
    
    if ( p[ d ] <  mRadius - mSize[ d ] ) {
      p[ d ] = mRadius - mSize[ d ]; 
      deepin = 0; 
    }
    }
  
  real dist = 0;
  
  if ( deepin ) {
    
    real test;
    int indx = 0;
    dist = mSize[ 0 ] - fabs( w[ 0 ] );
    
    for(int d = 1; d < DIM; ++d ) {
      test = mSize[ d ] - fabs( w[ d ] );
      if ( test < dist ) { 
        indx = d;
        dist = test;
      }
    }
    
    if ( w[ indx ] > 0 )
      p[ indx ] =  mSize[ indx ];
    else
      p[ indx ] = -mSize[ indx ];
    
  } else {
    
    //calculate the distance to the projection:
    for(int d = 0; d < DIM; ++d )
      dist += sqr( w[ d ] - p[ d ] );
    
    //normalize to mRadius, and add to p to get the real projection
    dist = mRadius / sqrt( dist );
    for(int d = 0; d < DIM; ++d )
      p[ d ] += dist * ( w[ d ] - p[ d ] );
    
  }
}

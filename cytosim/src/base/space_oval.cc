//RCS: $Id: space_oval.cc,v 2.4 2005/03/31 13:59:21 nedelec Exp $

#include "space_oval.h"

//-------------------------------------------------------------------------
real SpaceOval::volume() const
{
#if (DIM == 3)
  return 2 * mLength * PI * mRadiusSq + 4/3.0 * PI * mRadius * mRadiusSq;
#else
  return 2 * mLength * 2 * mRadius + PI * mRadiusSq;
#endif
}


//-------------------------------------------------------------------------
bool SpaceOval::isInside( const real w[] ) const
{
  real nrm, x = fabs( w[0] );
  
  if ( x > mLength )
    nrm = ( x - mLength ) * ( x - mLength );
  else
    nrm = 0;
  
  for( int d = 1; d < DIM; ++d )
    nrm += w[d] * w[d];
  
  
  return ( nrm <= mRadiusSq );
}


//-------------------------------------------------------------------------
void SpaceOval::project( const real w[], real p[] ) const
{
  real nrm = 0;
  for(int d = 1; d < DIM; ++d )
    nrm += w[ d ] * w[ d ];
  
  //calculate the projection on the axis, within boundaries:  
  if ( w[0] >  mLength ) {
    
    nrm  += ( w[0] - mLength )*( w[0] - mLength );
    //normalize from this point on the axis
    if ( nrm > 0 ) nrm = mRadius / sqrt( nrm );
    
    p[0] = mLength + nrm * ( w[0] - mLength );
    
    
  } else {
    if ( w[0] < -mLength ) {
      
      nrm  += ( mLength + w[0] )*( mLength + w[0] );
      //normalize from this point on the axis
      if ( nrm > 0 ) nrm = mRadius / sqrt( nrm );
      
      p[0]  = -mLength + nrm * ( w[0] + mLength );
      
    }
    else {
      
      //normalize from this point on the axis
      if ( nrm > 0 ) nrm = mRadius / sqrt( nrm );
      
      p[0] = w[0];
      
    }
  }
  
  if( nrm > 0 )
    for( int d = 1; d < DIM; ++d )
      p[ d ] = nrm * w[ d ];
  else {
    //we project on a arbitrary point on the cylinder
    p[1] = mRadius;
  }
  
  //printf("inside %i\n", insideOval(mLength, mRadius, p) );
}

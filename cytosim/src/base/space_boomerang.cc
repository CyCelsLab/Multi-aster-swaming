//
// C++ Implementation: space_boomerang
//
// Description: Can be (re)used to knock out yeast
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: space_boomerang.cc,v 1.3 2005/03/31 13:40:50 foethke Exp $


#include "space_boomerang.h"


//-------------------------------------------------------------------------
real SpaceBoomerang::volume() const
{
#if (DIM == 3)
  // funny, the volume is the same as for the oval
  return 2.*bLength*PI*bCapRadiusSq + 4./3.*PI*bCapRadiusSq*bCapRadius;
#else
  return 4.*bLength*bCapRadius + PI*bCapRadiusSq;
#endif
}


//-------------------------------------------------------------------------
bool SpaceBoomerang::isInside( const real w[] ) const
{
  real nrm, x;
  real wCurved[DIM];
    
  // we bend the space until the boomerang is straight
  if( w[0] > 0 ) {
    wCurved[0] = w[0]*cos( bAngle ) - w[1]*sin( bAngle );
    #if DIM > 1
    wCurved[1] = w[0]*sin( bAngle ) + w[1]*cos( bAngle );
    #endif
    #if DIM > 2
    wCurved[2] = w[2];
    #endif
  } else if ( w[0] < 0 ) {
    wCurved[0] = w[0]*cos( -bAngle ) - w[1]*sin( -bAngle );
    #if DIM > 1
    wCurved[1] = w[0]*sin( -bAngle ) + w[1]*cos( -bAngle );
    #endif
    #if DIM > 2
    wCurved[2] = w[2];
    #endif
  }
  
  x = fabs( wCurved[0] );
      
  if ( x > bLength )
    nrm = ( x - bLength ) * ( x - bLength );
  else
    nrm = 0;
  
  for( int d = 1; d < DIM; ++d )
    nrm += wCurved[d] * wCurved[d];
  
  
  return ( (nrm <= bCapRadiusSq) );
}


//-------------------------------------------------------------------------
void SpaceBoomerang::project( const real w[], real p[] ) const
{
  real nrm = 0;
  real wCurved[DIM];
  real pCurved[DIM];
  
  // we bend the space until the boomerang is straight
  if( w[0] >= 0 ) {
    wCurved[0] = w[0]*cos( bAngle ) - w[1]*sin( bAngle );
    #if DIM > 1
    wCurved[1] = w[0]*sin( bAngle ) + w[1]*cos( bAngle );
    #endif
    #if DIM > 2
    wCurved[2] = w[2];
    #endif
  } else if( w[0] < 0 ) {
    wCurved[0] = w[0]*cos( -bAngle ) - w[1]*sin( -bAngle );
    #if DIM > 1
    wCurved[1] = w[0]*sin( -bAngle ) + w[1]*cos( -bAngle );
    #endif
    #if DIM > 2
    wCurved[2] = w[2];
    #endif
  }
  
  for(int d = 1; d < DIM; ++d )
    nrm += wCurved[ d ] * wCurved[ d ];
  
  //calculate the projection on the axis, within boundaries:  
  //the right cap:
  if ( (wCurved[0] >  bLength) && (w[0] > 0) ) {
    nrm  += ( wCurved[0] - bLength )*( wCurved[0] - bLength );
    //normalize from this point on the axis
    if ( nrm > 0 ) nrm = bCapRadius / sqrt( nrm );
    
    pCurved[0] = bLength + nrm * ( wCurved[0] - bLength );  
  //the left cap
  } else if ( (wCurved[0] < -bLength) && (w[0] < 0) ) {
      nrm  += ( bLength + wCurved[0] )*( bLength + wCurved[0] );
      //normalize from this point on the axis
      if ( nrm > 0 ) nrm = bCapRadius / sqrt( nrm );
      
      pCurved[0]  = -bLength + nrm * ( wCurved[0] + bLength );
  } else {
    
    #if DIM == 2
    //points in the cone above the boomerang are projected to x=0
    if( (wCurved[0] < -tan(bAngle)*bCapRadius) && (w[0] >= 0) ) {
      wCurved[0] = -tan(bAngle)*bCapRadius;
    } else if( (wCurved[0] > tan(bAngle)*bCapRadius) && (w[0] < 0) ) {
      wCurved[0] =  tan(bAngle)*bCapRadius;
    }

    //points inside, in the cone above the angle are projected to x=0
    if( (wCurved[0] < tan(bAngle)*bCapRadius) && (wCurved[1] < 0) && (w[0] >= 0) ) {
      wCurved[0] =  tan(bAngle)*bCapRadius;
    } else if( (wCurved[0] > -tan(bAngle)*bCapRadius) && (wCurved[1] < 0) && (w[0] < 0) ) {
      wCurved[0] = -tan(bAngle)*bCapRadius;
    }
    #endif
    
    #if DIM == 3
    //in 3D we need to find the polar coordinates of the point in the yz-plane
    real yzRad, yzAngle;
    yzRad   = sqrt( wCurved[1]*wCurved[1] + wCurved[2]*wCurved[2] );
    yzAngle = atan2( wCurved[2], wCurved[1] );
    //points in the cone above the boomerang are projected to x=0
    if( (wCurved[0] < -tan(bAngle)*bCapRadius*cos(yzAngle)) && (w[0] >= 0) ) {
      wCurved[0] = -tan(bAngle)*bCapRadius*cos(yzAngle);
    } else if( (wCurved[0] > tan(bAngle)*bCapRadius*cos(yzAngle)) && (w[0] < 0) ) {
      wCurved[0] =  tan(bAngle)*bCapRadius*cos(yzAngle);
    }

    //points inside, in the cone above the angle are projected to x=0
    if( (wCurved[0] < tan(bAngle)*bCapRadius*cos(yzAngle)) && (wCurved[1] < 0) && (yzRad < bCapRadius) && (w[0] >= 0) ) {
      wCurved[0] =  tan(bAngle)*bCapRadius*cos(yzAngle);
    } else if( (wCurved[0] > -tan(bAngle)*bCapRadius*cos(yzAngle)) && (wCurved[1] < 0) && (yzRad < bCapRadius) && (w[0] < 0) ) {
      wCurved[0] = -tan(bAngle)*bCapRadius*cos(yzAngle);
    }
    #endif
    
    //normalize from this point on the axis
    if ( nrm > 0 ) nrm = bCapRadius / sqrt( nrm );
    
    pCurved[0] = wCurved[0];    
  }
  
  for( int d = 1; d < DIM; ++d )
    pCurved[ d ] = nrm * wCurved[ d ];
  
  // now we have to bend back the boomerang
  if( w[0] >= 0 ) {
    p[0] = pCurved[0]*cos( -bAngle ) - (pCurved[1])*sin( -bAngle );
    #if DIM > 1
    p[1] = pCurved[0]*sin( -bAngle ) + (pCurved[1])*cos( -bAngle );
    #endif
    #if DIM > 2
    p[2] = pCurved[2];
    #endif
  } else if( w[0] < 0 ) {
    p[0] = pCurved[0]*cos( bAngle ) - (pCurved[1])*sin( bAngle );
    #if DIM > 1
    p[1] = pCurved[0]*sin( bAngle ) + (pCurved[1])*cos( bAngle );
    #endif
    #if DIM > 2
    p[2] = pCurved[2];
    #endif
  }
}

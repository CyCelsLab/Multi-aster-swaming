//
// C++ Implementation: space_banana
//
// Description: Yellow and tasty
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: space_banana.cc,v 1.3 2005/03/31 13:40:50 foethke Exp $


#include "space_banana.h"


//-------------------------------------------------------------------------
real SpaceBanana::volume() const
{
#if (DIM == 3)
  // the volume of a torus is 2*PI^2*bRadius*bCapRadius^2
  return 2.*PI*PI*bRadius*bCapRadiusSq / ( bArcLength/2.*PI ) \
         + 4./3.*PI*bCapRadiusSq*bCapRadius;
#else
  return 2.*PI*((bRadiusSq + bCapRadiusSq)) / ( bArcLength/2.*PI ) \
         + PI*bCapRadiusSq;
#endif
}

Vecteur SpaceBanana::getBoundingRect() const
{
  if( bArcLength < PI/2. )
    return Vecteur( bRadius*sin(bArcLength)+bCapRadius, bRadius*(1.-cos(bArcLength))+bCapRadius, bCapRadius );
  else
    return Vecteur( bRadius+bCapRadius, bRadius*(1.-cos(bArcLength))+bCapRadius, bCapRadius );
}

//-------------------------------------------------------------------------
bool SpaceBanana::isInside( const real w[] ) const
{
  real center[DIM];      // a point on the backbone of the banana
  real nrm, wrad, angle;
    
  //calculate polar coordinates of point w relative to the origin
  //where the banana is bent around (0, -bRadius, 0)
  wrad  = w[0]*w[0];
  #if DIM > 1
  wrad += (-bRadius - w[1])*(-bRadius - w[1]);
  #endif
  wrad  = sqrt(wrad);
  
  //calculate the angle
  angle = acos( ( w[1] + bRadius )/wrad );
  if( w[0] < 0 ) angle = -angle;
  
  nrm = 0;
  if( angle > bArcLength ) {
    center[0] = bRadius*sin(bArcLength);
    #if DIM > 1
    center[1] = bRadius*cos(bArcLength) -bRadius;
    #endif
    #if DIM > 2
    center[2] = 0.;
    #endif
  }
  else if( angle < -bArcLength ) {
    center[0] = -bRadius*sin(bArcLength);
    #if DIM > 1
    center[1] = bRadius*cos(bArcLength) -bRadius;
    #endif
    #if DIM > 2
    center[2] = 0.;
    #endif
  }
  else {
    center[0] = bRadius*sin(angle);
    #if DIM > 1
    center[1] = bRadius*cos(angle) -bRadius;
    #endif
    #if DIM > 2
    center[2] = 0.;
    #endif
  }
  
  for( int ii = 0; ii < DIM; ii++ )
    nrm += (w[ii] - center[ii]) * (w[ii] - center[ii]);
  
  return ( (nrm <= bCapRadiusSq) );
}
  
//-------------------------------------------------------------------------   
void SpaceBanana::project( const real w[], real p[] ) const
{
  real center[DIM];      // a point on the backbone of the banana
  real dist[DIM];        // a vector pointing from center to w
  real nrm, wrad, angle;
    
  //calculate polar coordinates of point w relative to the origin
  //where the banana is bent around (0, -bRadius, 0)
  wrad  = w[0]*w[0];
  #if DIM > 1
  wrad += (-bRadius - w[1])*(-bRadius - w[1]);
  #endif
  wrad  = sqrt(wrad);
  
  //calculate the angle
  angle = acos( ( w[1] + bRadius )/wrad );
  if( w[0] < 0 ) angle = -angle;  
  
  //determine the point on the backbone of the banana
  //where w will be projected on
  if( angle > bArcLength ) {
    center[0] = bRadius*sin(bArcLength);
    #if DIM > 1
    center[1] = bRadius*cos(bArcLength) - bRadius;
    #endif
    #if DIM > 2
    center[2] = 0.;
    #endif
  }
  else if( angle < -bArcLength ) {
    center[0] = -bRadius*sin(bArcLength);
    #if DIM > 1
    center[1] = bRadius*cos(bArcLength) -bRadius;
    #endif
    #if DIM > 2
    center[2] = 0.;
    #endif
  }
  else {
    center[0] = bRadius*sin(angle);
    #if DIM > 1
    center[1] = bRadius*cos(angle) -bRadius;
    #endif
    #if DIM > 2
    center[2] = 0.;
    #endif
  }

  //calculate the distance between center and w
  dist[0] = w[0] - center[0];
  #if DIM > 1
  dist[1] = w[1] - center[1];
  #endif
  #if DIM > 2
  dist[2] = w[2];
  #endif
  
  nrm = 0;
  for( int ii = 0; ii < DIM; ii++ )
    nrm += dist[ii] * dist[ii];
      
  if ( nrm > 0 ) nrm = bCapRadius / sqrt( nrm );
  
  for( int ii = 0; ii < DIM; ii++ )
    p[ii] = center[ii] + dist[ii] * nrm;
}

//
// C++ Implementation: space_tee
//
// Description: 
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: space_tee.cc,v 1.12 2005/04/08 10:42:02 foethke Exp $


#include "space_tee.h"
#include "poly_roots.h"


//-------------------------------------------------------------------------
real SpaceTee::volume() const
{
///todo: calculate correct volume for T-shape
#if (DIM == 3)
  return 2 * tLength * PI * tCapRadiusSq + 4/3.0 * PI * tCapRadius * tCapRadiusSq;
#else
  return 2 * tLength * 2 * tCapRadius + PI * tCapRadiusSq;
#endif
}


//-------------------------------------------------------------------------
void SpaceTee::setMutables( const real w[] ) const
{
  xRel = (w[0] - tJunction);
  
  #if DIM == 2
  xzAngle = 0;
  delta   = atan2( fabs(xRel), w[1] );
  gamma   = PI/4.;
  #endif
  #if DIM == 3
  xzAngle = atan2( w[2], xRel );
  delta   = atan2( sqrt( xRel*xRel + w[2]*w[2] ), w[1] );
//  gamma   = PI/4. + atan2( sin(xzAngle)*sin(xzAngle), fabs(cos(xzAngle)) ) / 2.;
  gamma = atan2( tCapRadius, tCapRadius*fabs(cos(xzAngle)) );
  #endif  
}


//-------------------------------------------------------------------------
bool SpaceTee::isInside( const real w[] ) const
{
  real nrm;
  real x = fabs( w[0] );
  
  setMutables( w );
  
  if( (delta <= gamma) && (w[1] > tCapRadius*fabs(cos(xzAngle))) ) {
    if( w[1] > tArmLength+tCapRadius )
      nrm = (w[1] - (tArmLength+tCapRadius)) * (w[1] - (tArmLength+tCapRadius));
    else
      nrm = 0;
    
    nrm += xRel * xRel;
    #if (DIM > 2)
    nrm += w[2] * w[2];
    #endif
  } else {
    if ( x > tLength )
      nrm = (x - tLength) * (x - tLength);
    else
      nrm = 0;
  
    for( int d = 1; d < DIM; ++d )
      nrm += w[d] * w[d];  
  }
  
  return ( nrm <= tCapRadiusSq );
}


//-------------------------------------------------------------------------
void SpaceTee::project( const real w[], real p[] ) const
{
  real nrm = 0;
  real y;
    
  setMutables( w );
  
  if( (delta <= gamma) && (w[1] > tCapRadius*fabs(cos(xzAngle))) ) {
    //the point will be projected somewhere on the arm
    
    nrm += xRel * xRel;
    #if (DIM > 2)
    nrm += w[2] * w[2];
    #endif        
    
    if( w[1] > tArmLength+tCapRadius ) {
      y    = w[1] - (tArmLength+tCapRadius);
      nrm += y*y;
      
      if( nrm > 0 ) nrm = tCapRadius / sqrt( nrm );
    
      p[1] = tArmLength+tCapRadius + nrm*y;
    } else {
      if( nrm > 0 ) nrm = tCapRadius / sqrt( nrm );
      p[1] = w[1];
    }
    
    if( nrm > 0 ) p[0] = nrm*xRel   + tJunction;
    else          p[0] = tCapRadius + tJunction;
    #if (DIM > 2)
    p[2] = nrm*w[2];
    #endif    
  
  } else {
    //outside points are projected on the base cylinder
    //inside points are handled below
    if ( w[0] >  tLength ) {
      for(int d = 1; d < DIM; ++d )
        nrm += w[ d ] * w[ d ];
      nrm  += ( w[0] - tLength )*( w[0] - tLength );
      
      if ( nrm > 0 ) nrm = tCapRadius / sqrt( nrm );
    
      p[0] = tLength + nrm * ( w[0] - tLength );
        
      if ( nrm > 0 )
        for( int d = 1; d < DIM; ++d )
          p[ d ] = nrm * w[ d ];
      else p[1] = tCapRadius;
      
    } else if ( w[0] < -tLength ) {
      
      for(int d = 1; d < DIM; ++d )
      nrm += w[ d ] * w[ d ];
      
      nrm  += ( tLength + w[0] )*( tLength + w[0] );
      //normalize from this point on the axis
      if ( nrm > 0 ) nrm = tCapRadius / sqrt( nrm );
      
      p[0]  = -tLength + nrm * ( w[0] + tLength );

      if ( nrm > 0 )
        for( int d = 1; d < DIM; ++d )
          p[ d ] = nrm * w[ d ];
      else p[1] = tCapRadius;
      
    } else {
      
      #if( DIM == 2 )
      nrm += w[ 1 ] * w[ 1 ];
      //points inside, under the arm of the T are projected to the corners, well that's not the whole story...
      if( (w[0] > (tJunction-tCapRadius)) && (w[0] <= tJunction) && (w[1] > 0) && (w[1] <= tCapRadius) )
        if( w[1] >= ((-tCapRadius - xRel)*(-tCapRadius - xRel)/(4.*tCapRadius)) ) {
          //we project to the corner
          p[0] = tJunction-tCapRadius;
          p[1] = w[1];
        } else {
          //actually we are closer to the bottom of the base cylinder and project there
          p[0] = w[0];
          p[1] = -w[1]; // this is arbitrary, but conserves the norm
        }
      else if( (w[0] > tJunction) && (w[0] < tJunction+tCapRadius) && (w[1] > 0) && (w[1] <= tCapRadius) )
        if( w[1] >= ((tCapRadius - xRel)*(tCapRadius - xRel)/(4.*tCapRadius)) ) {
          //we project to the corner
          p[0] = tJunction+tCapRadius;
          p[1] = w[1];        
        } else {
          //actually we are closer to the bottom of the base cylinder and project there
          p[0] = w[0];
          p[1] = -w[1];
        }
      else {
        // we are not in the intersection area but project on the base cylinder
        p[0] = w[0];
        p[1] = w[1];
      }
  
      //normalize from this point on the axis
      if( nrm > 0 ) nrm  = tCapRadius / sqrt( nrm );
      if( nrm > 0 ) p[1] = nrm * p[1];
      else          p[1] = tCapRadius;
      #endif
      
      #if( DIM == 3 )
      if( (w[0] > (tJunction-tCapRadius)) && (w[0] <= tJunction+tCapRadius) && (w[1] > 0) && (w[1] <= tCapRadius) && isInside( w ) ) {
        //we are inside the intersection area and have to decide what to do
        real yzAngle    = atan2( w[1], fabs(w[2]) );
        real yzAngleMax = atan2( fabs(xRel), sqrt(tCapRadiusSq - xRel*xRel) );
        real yMin       = tCapRadius*fabs(cos(xzAngle));
                
        if( yzAngle <= yzAngleMax ) {
          //we project on the base cylinder
          for(int d = 1; d < DIM; ++d )
            nrm += w[ d ] * w[ d ];
          p[0] = w[0];
          if( nrm > 0 ) nrm = tCapRadius / sqrt( nrm );      
          for( int d = 1; d < DIM; ++d )
            p[ d ] = nrm * w[ d ];
        } else if( w[1] > yMin ) {
          //we project on the arm
          //we should almost never get here because these points are intercepted already above
          printf("Do we ever get here?\n");
          printf("delta: %e, gamma: %e, gamma-delta: %e\n", delta, gamma, gamma-delta);
          printf("w[1]: %e, yMin: %e, w[1]-yMin: %e\n", w[1], yMin, w[1]-yMin);
          nrm += xRel*xRel;
          nrm += w[2]*w[2];
          p[1] = w[1];
          if( nrm > 0 ) nrm = tCapRadius / sqrt( nrm );      
          p[0] = nrm*xRel + tJunction;
          p[2] = nrm*w[2];
        } else {
          //we are in the black hole and project on the intersection line
          if( w[2] != 0 ) {
            
            real xTurned, xTurnedSq;     // x position of turned point (see below)
            real zsq = w[2]*w[2];        // square of the points z coordinate
            real asq = 2.*tCapRadiusSq;  // square of the long half axis of the intersection ellipse
            real a, b, c, d, e ;         // coefficients of the quartic
            double sol[4];               // solutions of the quartic
            int numSol;                  // number of real solutions
            real tSol, xSol, xSolTurned; // the correct solutions of the quartic and of x
            real tCapRadiusQuad = tCapRadiusSq*tCapRadiusSq;  // tCapRadius^4
                                    
            // we turn the point, so that the intersection ellipse is in the xz-plane
            if( xRel >= 0 ) {
              xTurned   = (xRel + w[1]) / sqrt(2.);
              xTurnedSq = (xRel + w[1])*(xRel + w[1]) / 2.;
            } else {
              xTurned   = (xRel - w[1]) / sqrt(2.);
              xTurnedSq = (xRel - w[1])*(xRel - w[1]) / 2.;
            }
            
            // set the coefficients of the quartic
            a =   1.;
            b =   6.*tCapRadiusSq;
            c = (13.*tCapRadiusQuad                -    tCapRadiusSq*(2.*xTurnedSq + zsq));
            d = (12.*tCapRadiusQuad*tCapRadiusSq   - 4.*tCapRadiusQuad*(xTurnedSq + zsq));
            e =   4.*tCapRadiusQuad*tCapRadiusQuad - 2.*tCapRadiusQuad*tCapRadiusSq*(xTurnedSq + 2.*zsq);
            
            // solve the quartic
            numSol = PolyRoots::solveQuartic( a, b, c, d, e, sol );
            if( numSol < 0 ) {
              MSG.error("space_tee.cc","Failed to solve quartic for the intersection area.");
              exit(0);
            }

            // find the biggest t, (this is the shortest distance)
            tSol = sol[0];
            for( int ii = 1; ii < numSol; ii++ )
              if( sol[ii] > tSol ) tSol = sol[ii];
            
            // calculate x from t
            xSolTurned = asq*xTurned/(tSol + asq);
            
            // turn the point back to it's original position
            xSol = xSolTurned / sqrt(2.);
            p[0] = xSol + tJunction;
            p[1] = fabs(xSol);
            if( w[2] > 0 ) p[2] =  sqrt( tCapRadiusSq - xSol*xSol );
            else           p[2] = -sqrt( tCapRadiusSq - xSol*xSol );
          } else {
            if( xRel >= 0 ) p[0] = tJunction + tCapRadius;
            else            p[0] = tJunction - tCapRadius;
            p[1] = tCapRadius;
            p[2] = w[2];
          }
        } 
      } else {
        // we are not in the intersection area but project on the base cylinder
        for(int d = 1; d < DIM; ++d )
          nrm += w[ d ] * w[ d ];
        
        p[0] = w[0];
        
        //normalize from this point on the axis
        if( nrm > 0 ) nrm = tCapRadius / sqrt( nrm );
        
        if ( nrm > 0 )
          for( int d = 1; d < DIM; ++d )
            p[ d ] = nrm * w[ d ];
        else
          p[1] = tCapRadius;
      }
      #endif
    }
  }
}

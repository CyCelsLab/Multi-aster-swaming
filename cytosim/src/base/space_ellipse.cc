//RCS: $Id:"

#include "space_ellipse.h"

// void projectEllipse3( real mSize[3], real w[3], real p[3] )
//
// find the closest point (p) on an ellipse of mSize (mSize) to a given point (w)
//
// Analytically, it's a order 4 polynom to solve
// here we just iteratively find a point close on the ellipse in 2D...

#if (DIM == 1)

real SpaceEllipse::volume() const
{
  return 2 * mSize[0];
}

bool  SpaceEllipse::isInside( const real w[] ) const
{
  return (( w [ 0 ] >= -mSize[ 0 ] ) && ( w [ 0 ] <=  mSize[ 0 ] ));
}

void SpaceEllipse::project( const real w[], real p[] ) const
{
  if ( w[0] > 0 )
    p[ 0 ] =  mSize[ 0 ];
  else
    p[ 0 ] = -mSize[ 0 ];
}

#endif


#if (DIM == 2)


real SpaceEllipse::volume() const
{
  return PI * mSize[0] * mSize[1];
}


bool  SpaceEllipse::isInside( const real w[] ) const
{
  return ( sqr( w[ 0 ] / mSize[ 0 ] ) + sqr( w[ 1 ] / mSize[ 1 ] ) <= 1 );
}



void SpaceEllipse::project( const real w[2], real p[2] ) const
{
  real a = atan2( w[1]*mSize[0], w[0]*mSize[1] );
  real h = 0;
  real da;
  real ca, sa, t[2], tt, xw[2];
  
  const real coef = mSize[0] * mSize[1];
  //--------------iterate to refine the projection:
  for( int ii = 0; ii < 20; ++ii ) {
    
    ca = cos( a );
    sa = sin( a );
    // q is the current point on the ellipse, defined by angle a
    p[0] = mSize[0] * ca;
    p[1] = mSize[1] * sa;
    // nn is the tangent to the ellipse in p
    t[0] = -mSize[0] * sa;
    t[1] =  mSize[1] * ca;
    
    xw[0] = w[0] - p[0];
    xw[1] = w[1] - p[1];
    
    tt = t[0] * t[0] + t[1] * t[1];
    
    h = ( xw[0]*t[1] - xw[1]*t[0] ) / tt;
    if ( h > 0 )
      da = ( xw[0] * t[0] + xw[1] * t[1] ) / ( tt + h * coef );
    else
      da = ( xw[0] * t[0] + xw[1] * t[1] ) / tt;
    
    a += da * 0.5;
    if ( fabs(da) < EPSILON ) break;
  }
}

#endif


#if (DIM == 3)

real SpaceEllipse::volume() const
{
  return 4/3.0 * PI * mSize[0] * mSize[1] * mSize[2];
}

bool SpaceEllipse::isInside( const real w[] ) const
{  
  return  sqr( w[0]/mSize[0] ) + sqr( w[1]/mSize[1] ) + sqr( w[2]/mSize[2]) <= 1;
}

void SpaceEllipse::project( const real w[3], real p[3] ) const
{
  real ca, cb, sa, sb;

  ca = w[0] / mSize[0];
  sa = w[1] / mSize[1];
  sb = w[2] / mSize[2];
  
  real a = atan2( sa, ca );
  real b = atan2( sb, sqrt( ca*ca + sa*sa ) );
  real h = 0;
  
  Vecteur3 Ta, Tb, N, XW;
  real Ta2, Tb2, TaTb;
  real num, da, db;
  
  //const real minSize = minT( mSize[0], mSize[1], mSize[2] );
  const real coef = mSize[0] * mSize[1] * mSize[2];

  for( int ii=1; ii < 50; ++ii ) {
    
    ca = cos( a ); sa = sin( a );
    cb = cos( b ); sb = sin( b );
    
    p[0] = mSize[0] * ca * cb;
    p[1] = mSize[1] * sa * cb;
    p[2] = mSize[2] * sb;
    //P.set(   mSize[0] * ca * cb, mSize[1] * sa * cb, mSize[2] * sb );
    
    Ta.set( -mSize[0] * sa * cb, mSize[1] * ca * cb, 0 );
    Ta2 = Ta.normSquare();

    Tb.set( -mSize[0] * ca * sb, -mSize[1] * sa * sb, mSize[2] * cb );
    Tb2 = Tb.normSquare();
    
    XW = Vecteur( w[0]-p[0], w[1]-p[1], w[2]-p[2] );
    h = XW.norm();
    
    //the two tangent are not perpendicular:
    TaTb = Ta * Tb;
    num = Ta2 * Tb2 - TaTb * TaTb + h * coef;
    da = (( XW * Ta ) * Tb2 - ( XW * Tb ) * TaTb ) / num;
    db = (( XW * Tb ) * Ta2 - ( XW * Ta ) * TaTb ) / num;
    
    a += da;// * 0.5;
    b += db;// * 0.5;
    
    if ( fabs(da) + fabs(db) < EPSILON ) break;
    
    //if ( ii > 40 )
    //  printf("|XW| %9.3f   da %9.2e   db %9.2e\n", h, da, db );
    
  }
  //printf("\n");
}
#endif

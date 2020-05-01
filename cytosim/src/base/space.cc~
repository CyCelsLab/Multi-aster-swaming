//RCS: $Id: space.cc,v 2.13 2005/04/22 18:21:48 nedelec Exp $

#include "space.h"

#include "exceptions.h"
#include "space_square.h"
#include "space_periodic.h"
#include "space_sphere.h"
#include "space_oval.h"
#include "space_banana.h"
#include "space_roundsquare.h"
#include "space_strip.h"
#include "space_ellipse.h"
#include "space_cylinder.h"
#include "space_cylinderZ.h"
#include "space_combine.h"
#include "space_boomerang.h"
#include "space_tee.h"

#include "space_inflate.h"
#include "space_deflate.h"

#include "iomessages.h"

int Space::mPeriodic = 0;
real Space::mPeriod2[3] = { 0, 0, 0 };

//-------------------------------------------------------------------------
Space * newSpace( const int shape, const real size_set[6] )
{
  switch( shape ) {

    case SHAPE_UNSET:
      return new Space();

    case SHAPE_SQUARE:
      return new SpaceSquare(size_set);

    case SHAPE_PERIODIC:
      return new SpacePeriodic(size_set);

    case SHAPE_SPHERE:
      return new SpaceSphere(size_set);

    case SHAPE_OVAL:
      return new SpaceOval(size_set);

    case SHAPE_BANANA:
      return new SpaceBanana(size_set);

    case SHAPE_ROUND_SQUARE:
      return new SpaceRoundSquare(size_set);

    case SHAPE_STRIP:
      return new SpaceStrip(size_set);

    case SHAPE_ELLIPSE:
      return new SpaceEllipse(size_set);

    case SHAPE_CYLINDER:
      return new SpaceCylinder(size_set);

    case SHAPE_CYLINDER_Z:
      return new SpaceCylinderZ(size_set);

    case SHAPE_OUTSIDE_SPHERE: {
      Space * inS  = new SpaceSquare( size_set );
      Space * outS = new SpaceSphere( size_set[4] );
      return new SpaceCombine( inS, outS );
    }

    case SHAPE_BOOMERANG:
      return new SpaceBoomerang(size_set);

    case SHAPE_TEE:
      return new SpaceTee(size_set);

    default:
      return new Space();

  }
}

//-------------------------------------------------------------------------
Space * newSpace( const int shape, const real sizes[6], const real inflation )
{
  if ( inflation > 0 )
    return ( new SpaceInflate( newSpace( shape, sizes ), +inflation) );
  if ( inflation < 0 )
    return ( new SpaceDeflate( newSpace( shape, sizes ), -inflation) );

  return newSpace( shape, sizes );
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
// The volume is estimated in a simple monte-carlo. This function should be
// overwritten by an analytic formula when possible...
real Space::volume() const
{
  static real result = 0;

  static bool done = false;
  if ( done ) return result;
  done = true;

  const long CNT = 1<<18;
  long in = 0;
  real fac = 1.0, w[3];

  Vecteur rect = getBoundingRect();
  for( int d = 0; d < DIM; ++d )
    fac *= 2 * rect[ d ];

  for( long cnt = 0; cnt < CNT; ++cnt )
    {
      w[0] = rect.XX * RNG.sreal();
      w[1] = rect.YY * RNG.sreal();
      w[2] = rect.ZZ * RNG.sreal();
      in += isInside( w );
    }
  result = fac * real( in ) / real( CNT );
  MSG(4, "Monte-carlo estimated volume = %.2f +/- %.2f\n",
	result, sqrt( fac * result / CNT ));

  return result;
}


//-------------------------------------------------------------------------
//We provide an uniform random distribution in the volume by Monte-Carlo
//we draw a random point in a bounding box, until one pass the isInside() test.
Vecteur Space::randomPlaceInVolume() const
{
  assert( this );
  unsigned long ouf = 0;
  Vecteur result;

  do {

    result = getBoundingRect();

    result[ 0 ] *= RNG.sreal();
    if ( DIM > 1 ) result[ 1 ] *= RNG.sreal();
    if ( DIM > 2 ) result[ 2 ] *= RNG.sreal();
    if ( ++ouf > MAX_TRIALS_BEFORE_STUCK )
      throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::randomPlaceInVolume()");

    } while( ! isInside( result ) );

  return result;
}

//-------------------------------------------------------------------------
void Space::project( real w[] ) const
{
  real c[ 3 ];
  for( int d = 0; d < DIM; ++d)
    c[ d ] = w[ d ];
  project( c, w );
}


//-------------------------------------------------------------------------
//test if the sphere centered at w, of radius 'radius' entirely fits in the box
bool Space::isInside(const real center[], const real radius ) const
{
  assert( radius > 0 );
  //if the center is outside, we are outside for sure
  if ( ! isInside( center ) ) return false;

  //we project, to get the closest point W, on the surface of the space
  real p[ 3 ];
  project( center, p );

  //we calculate the distance to P (we could make DIM-specific code)
  real n = 0;
  for( int d = 0; d < DIM; ++d)
    n += ( center[d] - p[d] ) * ( center[d] - p[d] );

  //we are inside if the distance is greater than the given radius
  return ( n >= radius * radius );
}

//--------------------------------------------------------------------
//this code is equivalent to project() in SpaceDeflate
void Space::project(const real point[], real proj[], const real radius) const
{
  assert( radius > 0 );

  project(point, proj);

  real n = 0, pw[DIM];
  for(int d = 0; d < DIM; ++d) {
    pw[d] = point[d] - proj[d];
    n += pw[d] * pw[d];
  }

  ///\todo: problem in project() with radius if point is exactly on the box (n==0)
  //if (n==0) we do not know the orthogonal direction to follow. We should
  //take another point near by, and project from there.
  if ( n > 0 )
    n = ( isInside(point) ? +1 : -1 ) * radius / sqrt(n);

  for(int d=0; d<DIM; ++d)
    proj[d] += n * pw[d];
}

//-------------------------------------------------------------------------
Vecteur Space::randomPlaceOnEdge( real radius ) const
{
  unsigned long ouf = 0;
  if ( radius == 0 ) radius = 0.01;
  Vecteur result;
  do {
    result = randomPlaceInVolume();
    assert( isInside( result ) );
    if (++ouf>MAX_TRIALS_BEFORE_STUCK)
      throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in Space::randomPlaceOnEdge()");
  } while ( isInside( result, radius ) );
  return result;
}

//--------------------------------------------------------------------
void Space::setPeriodic( const real size[], const int semi )
{
  mPeriodic = semi;

  for(int ii = 0; ii < DIM; ++ii )
    mPeriod2[ii] = 2 * size[ii];
}

//--------------------------------------------------------------------


void  Space::modulo( real point[] )
{
  if ( mPeriodic == 0 )
    return;

  if ( mPeriodic == 2 ) {
    //semi-periodic, ie not in the last dimension
#if ( DIM > 1 )
    point[0] = drem( point[0], mPeriod2[0] );
#endif
#if ( DIM > 2 )
    point[1] = drem( point[1], mPeriod2[1] );
#endif

  } else {
    //periodic in all dimensions
    point[0] = drem( point[0], mPeriod2[0] );
#if ( DIM > 1 )
    point[1] = drem( point[1], mPeriod2[1] );
#endif
#if ( DIM > 2 )
    point[2] = drem( point[2], mPeriod2[2] );
#endif

  }
}


//-------------------------------------------------------------------------
//this makes modulo around the center o
void Space::moduloNear( real x[], const real o[] )
{
  for( int dd = 0; dd < DIM; ++dd )
    x[dd] -= o[dd];

  modulo( x );

  for( int dd = 0; dd < DIM; ++dd )
    x[dd] += o[dd];
}

//-------------------------------------------------------------------------
//this makes modulo around the center o
void Space::moduloNear( real m[], const real x[], const real o[] )
{
  for( int dd = 0; dd < DIM; ++dd )
    m[dd] = x[dd] - o[dd];

  modulo( m );

  for( int dd = 0; dd < DIM; ++dd )
    m[dd] += o[dd];
}

//-------------------------------------------------------------------------
//calculate both the integral part (in div) and the reminder (in x)
void Space::moduloAndOffset( real x[], real div[] )
{
  for( int dd = 0; dd < DIM; ++dd )
    div[dd] = x[dd];

  modulo(x);

  for( int dd = 0; dd < DIM; ++dd )
    div[dd] -= x[dd];
}

//-------------------------------------------------------------------------
//calculate both the integral part (in div) and the reminder (in x)
void Space::offset( real x[] )
{
  static real y[DIM];
  for( int dd = 0; dd < DIM; ++dd )
    y[dd] = x[dd];

  modulo(y);

  for( int dd = 0; dd < DIM; ++dd )
    x[dd] -= y[dd];
}



//RCS: $Id: fiber.cc,v 1.2 2005/04/22 12:12:40 nedelec Exp $

#include "fiber.h"
#include "point_interpolated.h"
#include "point_microtub.h"
#include "assert_macro.h"
#include "iowrapper.h"
#include "iomessages.h"
#include "exceptions.h"
#include "sim_param.h"

#include "fiber_solve.cc"

//===================================================================
//fiberConstructorBase() should be called only by the constructors
void Fiber::fiberConstructorBase()
{
  lagrangeMult = 0;
  constructorProjection();
  reset();
}

//-------------------------------------------------------------------
void Fiber::reset()
{
  mtType         = 0;
  mtcut          = 0;
  mtcutbest      = MP.mtrodlength;
  mtabminus      = 0;
  mtrigidkm      = 0;  //set in setMobility()
  lagrangeValid  = false;
#ifdef CUT_WITH_CURVATURE
  mtcuterror     = 0;
  ///giving a different 'seed' desynchronize the optimalCut() operations
  mtcuterrorindx = RNG.pint_exc(MTCUT_PERIOD);
#endif
}

//-------------------------------------------------------------------
int Fiber::allocate(const int size)
{
  int psize = PointSet::allocate(size);
  if ( psize ) {
    if ( lagrangeMult ) delete[] lagrangeMult;
    
    lagrangeMult = new real[ psize ];
    
    if ( lagrangeMult == 0 ) {
      fprintf(stderr, "Fiber::allocate() memory allocation failed\n");
      exit(1);
    }
  }
  return psize;
}

//-------------------------------------------------------------------
real Fiber::getLagrange(const int ii) const
{
  if ( lagrangeMult && lagrangeValid )
    return lagrangeMult[ii];
  return 0;
}


//-------------------------------------------------------------------
void Fiber::resetLagrange()
{
  lagrangeValid = false;
}

//-------------------------------------------------------------------
void Fiber::setNbPoints(int size)
{
  PointSet::setNbPoints(size);
  //we reset the values of all the Rods
  mtrods.setSize(size);
  for( int ii=0; ii < size; ++ii )
    mtrods[ ii ].setTo( this, ii, ii+1, 0 );
}

//-------------------------------------------------------------------
Fiber::~Fiber()
{
  mtrods.erase();
  deallocateProjection();
  if ( lagrangeMult ) delete[] lagrangeMult;
  lagrangeMult = 0;
}


//-------------------------------------------------------------------
PointMicrotub * Fiber::getRod( int pos )
{
  assert(( pos >= 0 ) && ( pos < nbSegments()) );
  assert( mtrods.size() >= nbSegments() );
  return &( mtrods[ pos ] );
}


//===================================================================

//-------------------------------------------------------------------
void Fiber::setStraight(const Vecteur where, const Vecteur dir)
{
  //we normalize the direction for more safety:
  assert( dir.norm() > EPSILON );
  Vecteur dpts = dir * ( mtcut / dir.norm() );
  for(int p = 0 ; p < nbPoints(); ++p )
    setPoint( p, where + p * dpts );
}

//-------------------------------------------------------------------
void Fiber::setStraight(const Vecteur where, const Vecteur dir, const MTEnd start)
{
  switch( start ) {
    
    case MT_MINUS_END: 
      setStraight( where, dir );
      break;
      
    case MT_PLUS_END:  
      setStraight( where + dir*length(), -dir );
      break;
      
    case MT_ORIGIN:                  //in fact the middle of the Microtub
      setStraight( where - 0.5*dir*length(), dir );
      break;
      
    default:
      MSG.error("Fiber::setStraight()","wrong MTEnd");
  }
}

//-------------------------------------------------------------------
void Fiber::setStraight(const Vecteur where, const Vecteur dir, const real newLength)
{
  assert( mtcutbest > EPSILON );
  int maxrod = maxT(1, int( newLength / mtcutbest ));
  setNbPoints( maxrod + 1 );
  mtcut = newLength / real( maxrod );
  setStraight( where, dir );
}

//-------------------------------------------------------------------
void Fiber::setStraight(const Vecteur where, const Vecteur dir, const real newLength, const MTEnd start)
{
  switch( start ) {

  case MT_MINUS_END: 
    setStraight( where, dir, newLength);
    break;

  case MT_PLUS_END:  
    setStraight( where + dir*newLength, -dir, newLength);
    break;

  case MT_ORIGIN:                  //in fact the middle of the Microtub
    setStraight( where - 0.5*dir*newLength, dir, newLength);
    break;

  default:
    MSG.error("Fiber::setStraight()","wrong MTEnd");
  }
}

//===================================================================
//  - - - - - - - - - - - CONSTRUCTORS - - - - - - - - - - - - - - -
//===================================================================

//set microtubules as a line, starting at w, of length alen
Fiber::Fiber(Vecteur w, Vecteur dir, real len, MTEnd end)
{
  fiberConstructorBase();
  setStraight(w, dir, len, end);
}

//===================================================================
//===================================================================

//move the points to satisfy the length constraints
//while conserving the center of gravity of the tubule
void Fiber::reshape()
{
  int d, p, q, ii, jj = 0;
  real rightleft, dl, dd, dy[3];

  for( p=0; p < nbSegments(); ++p )  {
    
    ii = jj;
    jj = DIM * ( p + 1 );
    
    for(dd=0, d=0; d < DIM; ++d ) {
      dy[ d ] = pspts[ jj + d ] - pspts[ ii + d ];
      dd += dy[ d ] * dy[ d ];
    }
    
    dd = sqrt( dd );
    
    if ( dd > 0 ) {
      dl = ( nbSegments() - p ) * ( dd - mtcut ) / ( dd * nbPoints() ) ;
      rightleft = real( p + 1 ) / real( nbSegments() - p );
      
      for( d = 0; d < DIM; ++d ) {
        dy[ d ] *= dl;
        for( q = 0; q <= p; ++q ) 
          pspts[ DIM * q + d ] += dy[ d ];
                
        dy[ d ] *= rightleft;
        for( ; q < nbPoints(); ++q ) 
          pspts[ DIM * q + d ] -= dy[ d ];
      }
    } else {
      /** Two consecutive points overlap: this is definitely a bad thing !
      we do nothing, hoping that the Brownian motion will push them appart
      */
    }
  }
}

/*
// ------------  standard ( unoptimized ) version:

void Fiber::reshape()
  //move the points to satisfy the length constraints
  //while conserving the center of gravity of the tubule
{
  int d, p;
  Vecteur a, b;
  real distance;
  Vecteur correction;
  real left, right;

  for( p = 0; p < nbSegments(); ++p )
    {
      a = whereP( p );
      b = whereP( p+1 );

      distance   = ( b - a ).norm();
      correction = ( mtcut / distance - 1 ) * ( a - b );

      left  = real( nbSegments()-p ) / real( nbPoints() );

      for( d = 0; d <= p; ++d )
	       movePoint( d, left * correction );

      right = real( p+1 ) / real( nbPoints() );

      for( d = p+1; d < nbPoints(); ++d )
	       movePoint( d, -right * correction );
    }
}
*/



//------------------------------------------------------------------------
//------------- DISTANCE FROM A POINT TO A SECTION OF TUBE ---------------
//------------------------------------------------------------------------
/** the only consideration for calculating the distance of a point to a rod is
the case where the projection falls outside the rods, i.e. the point is
closest to one of the ending point of that rod.

the function distance() only considers binding to the ends of the
microtubules (+/- ends). Binding to the ends of inners rods is forbiden.
*/

//------------------------------------------------------------------------
//caculates the distance from this rod to position w, returning the result in *si
real Fiber::distanceToRod( PointInterpolated * si, const Vecteur & w ) const
{
  //The PointInterpolated si should have two values set: mPS and mPoint1
  //to define on which MT segment the distance will be calculated
  //the value si->mCoef will be set here as the projection of w on the segment
  assert( si->mPS == this );
  assert( si->mPoint1 >= 0 );
  assert( si->mPoint1 < nbSegments() );
  
  Vecteur dx = dpts( si->mPoint1 );
  Vecteur aw = w - whereP( si->mPoint1 );
  Space::modulo( aw );

  si->mCoef = ( aw * dx );

  if ( si->mCoef < 0 ) {
    if ( si->mPoint1 > 0 ) 
      return 10000;      //this rod is not the minus end
    si->mCoef  = 0;
    return aw.norm();
  }
  
  if ( si->mCoef > mtcut*mtcut ) {
    if ( si->mPoint1 < lastSegment() ) 
      return 10000;      //this rod is not the plus-end
    si->mCoef = 1.0;
    aw -= dx;
    Space::modulo( aw );
    return aw.norm();
  }

  //we could optimize here, as mtcut is constant:
  si->mCoef /= mtcut * mtcut;

  return ( aw - si->mCoef * dx ).norm();
}

//------------------------------------------------------------------------
//  distanceToRodEnd calculates the real distance, considering the end points
//  for any rod 
//------------------------------------------------------------------------

real Fiber::distanceToRodWithEnds( PointInterpolated * si, const Vecteur & w ) const
{
  assert( si->mPS == this );
  assert( si->mPoint1 >= 0 );
  assert( si->mPoint1 < nbSegments() );
  
  Vecteur dx = dpts( si->mPoint1 );
  Vecteur aw = w - whereP( si->mPoint1 );
  Space::modulo( aw );

  si->mCoef = ( aw * dx ) / ( mtcut * mtcut );

  if ( si->mCoef < 0 ) {
    si->mCoef = 0;
    return aw.norm();
  }

  if ( si->mCoef > 1.0 ) {
    si->mCoef = 1.0;
    aw = w - whereP( si->mPoint1 + 1 );
    Space::modulo( aw );
    return aw.norm();
  }

  return ( aw - si->mCoef * dx ).norm();
}


//------------------------------------------------------------------------
// distanceToMT considers all the rod in the given microtubule, finding
// the closest rod... not optimized... use tryToAttach() instead
//------------------------------------------------------------------------

real Fiber::distanceToMT( PointInterpolated * si, const Vecteur & w ) const
{
  assert( si->mPS == this );
  
  PointInterpolated IP(this, 0, 1, 0), result;
    
  real test, closest = distanceToRodWithEnds( &IP, w );
  *si = IP;
  
  //try all the rods one by one:
  for( int ii = 1; ii < nbSegments(); ++ii ) {
    IP.setTo( this, ii, ii+1, 0 );
    test = distanceToRodWithEnds( &IP, w );
    if ( test < closest ) {
      closest = test;
      *si = IP;
    }
  }
  
  return closest;
}




//------------------------------------------------------------------------
//  distanceToRodFast is an optimized version of distanceToRod
//------------------------------------------------------------------------

#if ( DIM == 1 )

real Fiber::distanceToRodFast( PointInterpolated * si, const Vecteur & w ) const
{
  assert( si->mPS == this );
  assert(( si->mPoint1 >= 0 ) && ( si->mPoint1 < nbSegments() ));
  
  real dx = pspts[ si->mPoint1 + 1 ] - pspts[ si->mPoint1 ] ;
  real d  = w.XX - pspts[ si->mPoint1 ];

  //Space::modulo( dx );
  assert( Space::isPeriodic() == false );

  si->mCoef = d / dx;
  
  if ( si->mCoef < 0 ) {
    if ( si->mPoint1 > 0 )      
      return 10000;  //rod is not at the end of mt
    si->mCoef = 0;
    return fabs( d );
  }
  
  if ( si->mCoef > 1.0 ) {
    if ( si->mPoint1 < lastSegment() )
      return 10000;  //rod is not at the end of mt
    si->mCoef = 1.0;
    return fabs( d - dx );
  }
  
  return 0;
}


#endif


#if ( DIM == 2 )

real Fiber::distanceToRodFast( PointInterpolated * si, const Vecteur & w ) const
{
  assert( si->mPS == this );
  assert(( si->mPoint1 >= 0 ) && ( si->mPoint1 < nbSegments() ));

  real dx = pspts[ DIM * si->mPoint1 + DIM ]     - pspts[ DIM * si->mPoint1 ] ;
  real dy = pspts[ DIM * si->mPoint1 + DIM + 1 ] - pspts[ DIM * si->mPoint1 + 1 ] ;

  real ax = w.XX - pspts[ DIM * si->mPoint1 ];
  real ay = w.YY - pspts[ DIM * si->mPoint1 + 1 ];

  //Space::modulo( ax, ay );
  assert( Space::isPeriodic() == false );

  si->mCoef = ( ax * dx + ay * dy ) / ( mtcut * mtcut );

  if ( si->mCoef < 0 ) {
    if ( si->mPoint1 > 0 ) 
      return 10000;      //rod is not an end
    si->mCoef = 0;
    return sqrt( ax * ax + ay * ay );
  }
  
  if ( si->mCoef > 1.0 ) {
    if ( si->mPoint1 < lastSegment() ) return 10000; //rod is not an end
    si->mCoef = 1.0;
    return ( w - whereEnd(MT_PLUS_END) ).norm();
  }
  
  return fabs( dx * ay - dy * ax ) / mtcut;
}
  
#endif        //   DIM == 2
  
  
#if ( DIM == 3 )

real Fiber::distanceToRodFast( PointInterpolated * si, const Vecteur & w ) const
{
  assert( si->mPS == this );
  assert(( si->mPoint1 >= 0 ) && ( si->mPoint1 < nbSegments() ));
    
  real dx = pspts[ DIM * si->mPoint1 + DIM ]     - pspts[ DIM * si->mPoint1 ] ;
  real dy = pspts[ DIM * si->mPoint1 + DIM + 1 ] - pspts[ DIM * si->mPoint1 + 1 ] ;
  real dz = pspts[ DIM * si->mPoint1 + DIM + 2 ] - pspts[ DIM * si->mPoint1 + 2 ] ;
    
  real ax = w.XX - pspts[ DIM * si->mPoint1 ];
  real ay = w.YY - pspts[ DIM * si->mPoint1 + 1 ];
  real az = w.ZZ - pspts[ DIM * si->mPoint1 + 2 ];

  //Space::modulo( ax, ay, az );
  assert( Space::isPeriodic() == false );
  
  si->mCoef = ( ax * dx + ay * dy + az * dz ) / ( mtcut * mtcut );
    
  if ( si->mCoef < 0 ) {
    if ( si->mPoint1 > 0 )
      return 10000;      //rod is not an end
    si->mCoef = 0;
    return sqrt( ax * ax + ay * ay + az * az );
  }
  
  if ( si->mCoef > 1.0 ) {
    if ( si->mPoint1 < lastSegment() ) 
      return 10000;  //rod is not an end
    si->mCoef = 1.0;
    return ( w - whereEnd(MT_PLUS_END) ).norm();
  }
  
  ax -= si->mCoef * dx;
  ay -= si->mCoef * dy;
  az -= si->mCoef * dz;
  return sqrt( ax * ax + ay * ay + az * az );
}
#endif


//========================================================================
//=====================GROWING/SHRINKING==================================
//========================================================================

/** 
When a microtubule grows, the distance between points increases. 
When the distance (mtcut) is more than twice as long than the ideal cut 
(mtcutbest), a new point is added. 
*/

void Fiber::recut(int nbpts)
{
  if ( nbpts == nbPoints() ) return;

  assert( nbpts > 1 );
  real mtcut_new = length() / real( nbpts-1 );

  MSG(12,"Fiber::recut M%lx : points %i -> %i,  section %.3f -> %.3f\n",
      name, nbPoints(), nbpts, mtcut, mtcut_new );

  //we allocate a tmp array to hold the new points:
  static real * tmp = 0;
  static int allocated = 0;
  if ( nbpts > allocated ) {
    if ( tmp ) delete[] tmp;
    allocated = nbpts + 2;
    tmp = new real[ DIM * allocated ];
    if ( tmp == 0 ) {
      fprintf(stderr, "Fiber::recut(int) memory allocation failed\n");
      exit(1);
    }
  }

  Vecteur w;

  // TODO: 2d order interpolation in Microtub recut
  // calculate intermediate points into tmp[]:
  for(int p = 1; p < nbpts-1; ++p ) {
      w = where( p * mtcut_new, MT_MINUS_END );
      tmp[ DIM*p   ] = w.XX;
#if ( DIM >= 2 )
      tmp[ DIM*p+1 ] = w.YY;
#endif
#if ( DIM >= 3 )
      tmp[ DIM*p+2 ] = w.ZZ;
#endif
    }

  // copy the position of plus-end into tmp[]:
  int pp = DIM*(nbpts-1);
  int lp = DIM*lastPoint();
  for(int d = 0 ; d < DIM; ++d )
    tmp[ pp + d ] = pspts[ lp + d ];

  setNbPoints( nbpts );

  // copy calculated points back into pspts
  for(int p = DIM; p < DIM*nbpts; ++p )
    pspts[ p ] = tmp[ p ];
  
  mtcut = mtcut_new;
}

//----------------------------------------------------------------------------
/** cutting to minimize abs( mtcut - mtcutbest ):
    this ensures      2/3  < mtcut/mtcutbest <  4/3  */


real Fiber::minCosRodAngles()
{
  real s, result = mtcut * mtcut;
  Vecteur dir1 = dpts(0), dir2;
  
  for( int p = 1; p < nbSegments() ; ++p ) {
    if ( p % 2 )
      dir2 = dpts(p);
    else
      dir1 = dpts(p);
    s = dir1 * dir2;
    if ( s < result ) result = s;
  }
  return result / ( mtcut * mtcut );
}


#ifdef CUT_WITH_CURVATURE


// the MT section length dependent on the curvature, not the total fiber length!
// we only recut once every time step, to allow the solver to equilibrate the position
void Fiber::optimalCut()
{
  //we override constantly mtcutbest, to allow refining live from the player
  //mtcutbest = MP.mtrodlength;  //live refining
  
  if ( nbPoints() < 3 ) {
    if ( length() > 1.333 * mtcutbest )
      recut( 3 );
    return;
  }
  
  //very short Fibers only get two points
  if (( nbPoints() == 3 ) && ( length() < 1.333 * mtcutbest )) {
    recut( 2 );
    return;
  }
  
  //accumulate the error in variable mtcuterror, to smoothen the Brownian influence
  mtcuterror += mtcut * sqrt( 1.0 - minCosRodAngles() );
  //increment the counter
  ++mtcuterrorindx;

  if ( mtcuterrorindx < MTCUT_PERIOD ) return;
  
  //we scale the error accumulated over the last PERIOD calls
  mtcuterror /= ( MTCUT_PERIOD * MP.mtcuterror );
  
  //cut more finely if the error is large
  if ( mtcuterror > 1.0 ) {
    if ( mtcut > MP.mtcutmin )
      recut( nbPoints() + 1 );
  } else {
    //cut less finely if the error is small
    if (( mtcuterror < 0.4 ) && ( nbPoints() > 3 )) {
      recut( nbPoints() - 1 );
    }
  }
  
  //reset the counter and error accumulator
  mtcuterror     = 0;
  mtcuterrorindx = 0;
}

#else

//the sectionning of fibers is done only as a function of their length
void Fiber::optimalCut()
{
  mtcutbest = MP.mtrodlength;      //to allow recutting live from the player

  while( mtcutbest * nbPoints() < ( mtcut - EPSILON ) * ( nbPoints()-0.5 ))
    recut( nbPoints() + 1 );

  while(( nbPoints() > 2 )
	 && (( mtcut + EPSILON ) * ( nbPoints()-1.5 ) < mtcutbest * ( nbPoints()-2 )) )
    recut( nbPoints() - 1 );
}

#endif



//____________________________________________________________________________
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
/**
 grow1(lon) and grow2(lon) increase length of microtubule by 'lon', at the minus- 
 or the plus-end. We do not interpolate, which works for small increments and
 fairly straight Microtub.
*/

//----------------------------------------------------------------------------
void Fiber::grow1(const real lon)
{
  real a = -lon / length();

  if ( a < 0 ) {

    int p = 0, q = nbSegments();
    Vecteur dp0 = dpts( 0 ), dp1 = dp0;
    movePoint( p, ( a * q ) * dp0 );
    while( ++p < nbSegments() ) {
      --q;
      if ( p % 2 == 1 ) {
        dp0 = dpts( p );
        movePoint( p, ( a * q ) * dp1 );
      } else {
        dp1 = dpts( p );
        movePoint( p, ( a * q ) * dp0 );
      }
    }
    
  } else {
    for( int p = 0, q = nbSegments() ; p < nbSegments() ; ++p, --q )
      movePoint( p, ( a * q ) * dpts( p ));
  }
  
  mtcut += lon / real( nbSegments() );
  mtabminus -= lon;
}


//----------------------------------------------------------------------------
void Fiber::grow2(const real lon)
{
  real a = lon / length();

  if ( a > 0 ) {
  
    int p = nbSegments();
    Vecteur dp0 = dpts( p - 1 ), dp1 = dp0;
    movePoint( p, ( a * p ) * dp0 );
    while( --p > 0 ) {
      if ( p % 2 == 1 ) {
        dp0 = dpts( p - 1 );
        movePoint( p, ( a * p ) * dp1 );
      } else {
        dp1 = dpts( p - 1 );
        movePoint( p, ( a * p ) * dp0 );
      }
    }

  } else
    for( int p = nbSegments() ; p > 0 ; --p )
      movePoint( p, ( a * p ) * dpts( p - 1 ));

  mtcut += lon / real( nbSegments() );
}


//----------------------------------------------------------------------------
void Fiber::growAtEnd( const MTEnd which, const real lon)
{
  assert(( which == MT_PLUS_END ) || ( which == MT_MINUS_END ));
  
  if ( which == MT_MINUS_END ) {
    if ( lon != 0 ) grow1( lon ); 
  } else {
    if ( lon != 0 ) grow2( lon );
  }
}


//----------------------------------------------------------------------
//---------------------     convenience functions     ------------------
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//difference of two consecutive points, p and p+1
Vecteur Fiber::dpts( const int p ) const
{
  assert(( p >= 0 ) && ( p < nbSegments() ));
  
#if ( DIM == 1 )
  return Vecteur( pspts[ p+1 ] - pspts[ p ], 0, 0 );
#endif
#if ( DIM == 2 )
  return Vecteur( pspts[ 2*p+2 ] - pspts[ 2*p ], pspts[ 2*p+3 ] - pspts[ 2*p+1 ], 0 );
#endif
#if ( DIM == 3 )
  return Vecteur( pspts[ 3*p+3 ] - pspts[ 3*p ], pspts[ 3*p+4 ] - pspts[ 3*p+1 ], pspts[ 3*p+5 ] - pspts[ 3*p+2 ] );
#endif
  /*
  Vecteur result;
  for(int d = 0; d < DIM; ++d )
    result[ d ] = pspts[ DIM*p + d + DIM ] - pspts[ DIM*p + d ];
  return result;
   */
}


//----------------------------------------------------------------------------
// tangent vector to the tubule, ie. the normalized difference
Vecteur Fiber::dirP( const int p ) const
{
  assert(( p >= 0 ) && ( p < nbSegments() ));
  //here we divide by mtcut, which should be the distance between points.
  //we could otherwise normalize, but that is more costly
  Vecteur result;     
  for(int d=0; d < DIM; ++d )
    result[d] = ( pspts[ DIM*(p+1) + d ] - pspts[ DIM*p + d ] ) / mtcut;
  return result;
}

//----------------------------------------------------------------------------
real Fiber::abscissa( const MTEnd end ) const
//return the abs with respect to the MT_ORIGIN.
{
  switch( end ) {
    case MT_ORIGIN:    return 0;
    case MT_MINUS_END: return mtabminus;
    case MT_PLUS_END:  return mtabminus + mtcut * nbSegments();
    default:
      MSG.error("abscissa(MTEnd)","wrong argument"); return 0;
  }
}


//----------------------------------------------------------------------------
real Fiber::abscissa( const MTEnd end, const MTEnd from ) const
{
  switch( from ) {
    
    case MT_ORIGIN:
      switch( end ) {
        case MT_ORIGIN:    return 0;
        case MT_MINUS_END: return mtabminus;
        case MT_PLUS_END:  return mtabminus + mtcut * nbSegments();
        default:           MSG.error("abscissa", "wrong argument");
      }
      
    case MT_MINUS_END:
      switch( end ) {
        case MT_ORIGIN:    return -mtabminus;
        case MT_MINUS_END: return 0;
        case MT_PLUS_END:  return mtcut * nbSegments();
        default:           MSG.error("abscissa", "wrong argument");
      }
      
    case MT_PLUS_END:                  //taken towards the minus-end
      switch( end ) {
        case MT_ORIGIN:    return mtabminus + mtcut * nbSegments();
        case MT_MINUS_END: return mtcut * nbSegments();
        case MT_PLUS_END:  return 0;
        default:           MSG.error("abscissa", "wrong argument");
      }
      
    default:
      MSG.error("abscissa", "wrong argument");
  }
  return 0;
}

//----------------------------------------------------------------------------
real Fiber::abscissa( const real rod ) const
{  
  return  mtcut * rod + mtabminus;
}

//----------------------------------------------------------------------------
real Fiber::abscissa( const real rod, const MTEnd from ) const
{  
  switch( from ) {
    
    case MT_ORIGIN:
      return  mtcut * rod + mtabminus;
      
    case MT_MINUS_END:
      return  mtcut * rod;
      
    case MT_PLUS_END:                  //taken towards the minus-end
      return  mtcut * ( nbSegments() - rod );
      
    default:
      MSG.error("abscissa", "wrong argument");
  }
  return 0;
}

//----------------------------------------------------------------------------
void Fiber::setInterpolation( PointInterpolated * si, const MTEnd which ) const
{
  switch( which ) {
    case MT_MINUS_END:
      si->setTo( this, 0, 1, 0.0 );
      break;
    case MT_PLUS_END:
      si->setTo( this, lastSegment(), lastPoint(), 1.0 );
      break;
    default:
      MSG.error("Fiber::setInterpolation(MTEnd)", "wrong argument");
    }
}

//----------------------------------------------------------------------------
/** converts the abscisse ab into a rod position r, and a relative abscisse a
    the point X is defined by the distance <ab> counted from the specified end
    (r,a) is calculated to fullfill X = P(r) * (1-a) + P(r+1) * a
   r is an integer 0 < r < maxrod ,   and    0 < a < 1
*/

MTEnd Fiber::setInterpolation( PointInterpolated * si, const real ab, const MTEnd from ) const
{
  MTEnd result = MT_NOT_END;
  double n, co;
  int rd;

  switch( from ) {

    case MT_ORIGIN:
      assert( ab >= abscissa( MT_MINUS_END ) - EPSILON );
      assert( ab <= abscissa( MT_PLUS_END  ) + EPSILON );
      co = modf(( ab - mtabminus ) / mtcut, &n );
      rd = int(n);
      break;
    
    case MT_MINUS_END:
      assert( ab >= - EPSILON );
      assert( ab <= length() + EPSILON );
      co = modf( ab / mtcut, &n );
      rd = int(n);
      break;
      
    case MT_PLUS_END:         //this is counted from the plus towards the minus end
      assert( ab >= - EPSILON );
      assert( ab <= length() + EPSILON );
      co = 1 - modf( ab / mtcut, &n );
      rd = lastSegment() - int( n );
     break;
      
    default:
      rd = 0; co = 0.0;  //not necessary, but avoid a warning
      MSG.error("setInterpolation","wrong argument");
  }

  //bound the values, to produce a valid output:
  if ( rd < 0 ) {
    rd = 0;
    co = 0.0;
    result = MT_MINUS_END;
  }
  else
  if ( rd > lastSegment() ) { 
    rd = lastSegment();
    co = 1.0;
    result = MT_PLUS_END;
  }
  
  assert(( co >= 0.0 ) && ( co <= 1.0 ));
  si->setTo( this, rd, rd+1, co );
  return result;
}

//----------------------------------------------------------------------------
Vecteur Fiber::whereEnd( MTEnd which ) const
{
  switch( which ) {
    case MT_MINUS_END: return whereP( 0 );
    case MT_PLUS_END:  return whereP( lastPoint() );
    default: MSG.error("Fiber::whereEnd(MTEnd)", "wrong argument");
      return VZERO;
  }
}


//----------------------------------------------------------------------------
real Fiber::forceOnEnd( MTEnd which ) const
{
  switch( which ) {
    case MT_MINUS_END: return ( -detForcesP(0)           * dirEnd(MT_MINUS_END) );
    case MT_PLUS_END:  return (  detForcesP(lastPoint()) * dirEnd(MT_PLUS_END) );
    default: MSG.error("Fiber::forceOnEnd(MTEnd)", "wrong argument");
  }
  return 0;
}


//----------------------------------------------------------------------------
Vecteur Fiber::dirEnd( MTEnd which ) const
{
  switch( which ) {
    case MT_MINUS_END: return dirP( 0 );
    case MT_PLUS_END:  return dirP( lastSegment() );
    default: MSG.error("Fiber::dirEnd(MTEnd)", "wrong argument");
      return VZERO;
  }
}
//----------------------------------------------------------------------------
Vecteur Fiber::where( const PointInterpolated * si ) const
{
  assert( si->mPS == this );
  assert( si->looksWrong() == NO_ERROR );

  if ( si->mPoint1 > lastSegment() ) {
    //assert( si->mCoef == 0 );
    return whereEnd(MT_PLUS_END);
  }

  //TODO: optimize by precalculating dpts = pts[n+DIM]-pts[n]
  Vecteur result;
  for(int d = 0, p = DIM * si->mPoint1; d < DIM; ++d, ++p )
    result[ d ] = pspts[ p ] + si->mCoef * ( pspts[ p + DIM ] - pspts[ p ] );
  return result;
}


//----------------------------------------------------------------------------
Vecteur Fiber::dirMT( const PointInterpolated * si ) const
{
  assert( si != 0 );
  assert( si->mPS == this );
  assert( si->looksWrong() == NO_ERROR );

  if ( si->mPoint1 > lastSegment() ) {
    //assert( si->mCoef == 0 );
    return dirP(MT_PLUS_END);
  }

  return dirP( si->mPoint1 );
}

//----------------------------------------------------------------------------
Vecteur Fiber::where( const real ab, const MTEnd from ) const
{
  PointInterpolated si;
  setInterpolation( &si, ab, from );
  return where( & si );
}


//----------------------------------------------------------------------------
Vecteur Fiber::dirMT( const real ab, const MTEnd from ) const
{
  PointInterpolated si;
  setInterpolation( &si, ab, from );
  return dirMT( & si );
}

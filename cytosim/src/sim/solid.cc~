//RCS: $Id: solid.cc,v 2.19 2005/04/22 12:13:22 nedelec Exp $
//---------------------------------solid.cc----------------------------------

#include "assert_macro.h"
#include "solid.h"
#include "aster.h"
#include "clapack.h"
#include "exceptions.h"
#include "grafted.h"
#include "iowrapper.h"
#include "sim.h"

#if ( DIM == 3 )
  #include "quaternion.h"
  #include "matrix3.h"
#endif


// the Space in which solids live:
Space * Solid::space = 0;

//==============================METHODS======================================
void Solid::diffuse()
{
  //beware: this is only a translation diffusion
  assert( nbPoints() == 1 );
  static Vecteur w, v;
  w = whereP(0);
  w.addRandom(MP.cxdiff_dt[0]);
  if ( space->isOutside(w) ) {
    space->project(w, v);
    w = 2 * v - w;
  }
  setPoint(0, w);
}

//---------------------------------------------------------------------------
// Base Constructor,  to be called by all contructors:
// reset essential variables, 
void Solid::solidConstructorBase()
{
  soType          = 0;
  soShape         = 0;
  soShapeSize     = 0;
  soShapeMomentum = 0;
}


//---------------------------------------------------------------------------
// 1. allocate pspts[] in PointSet,
// the returns value can be zero if nothing was done, 
// or the new size of pspts[] which was allocated.
int Solid::allocate(const int size)
{
  int psize = PointSet::allocate(size);
  
  //if PointSet allocated for pspts[], we allocate the same for soShape[]
  if ( psize ) {
    // allocate a new array of the right size:
    real  *  soShape_new = new real [ DIM*psize ];
    
    if ( soShape_new == 0 ) {
      fprintf(stderr, "Solid::allocate(int) memory allocation failed\n");
      exit(1);
    }
    
    // copy the current values in the new array:
    if ( soShape ) {
      for ( int p = 0; p < DIM * nbPoints(); ++p )
        soShape_new[ p ] = soShape[ p ];
      // delete the 'current' array:
      delete[] soShape;
    }
    // the 'new' array becomes the 'current' one:
    soShape = soShape_new;
        
    MSG(12, "Solid::allocate %lx ask %d done %d\n", this, size, psize);
  }
  
  // return the size allocated (following PointSet::allocate()
  return psize;
}

//---------------------------------------------------------------------------
void Solid::deallocate()
{
  if ( soShape )
    delete[] soShape;
  soShape = 0;
}
  

//---------------------------------------------------------------------------
//the constructor called to create an initial configuration,
//creates a Solid, according to parameters MP.so????
Solid::Solid(int)
{
  //MSG(12, "Solid::Solid(int) %lx\n", this);

  //reset the member variables:
  solidConstructorBase();

  Vecteur w, x;

  //find a position in the Solid space, according to initial code:  
  w = initPosition( space, MP.soinit );
  
  //allocate the requested number of points (MP.soptsmax)
  //a Grafted motor is attached at every point, of a type chosen
  //randomly according to the ratios defined in MP.soghmax:

  switch( MP.soptsmax ) {
    
    case 0:
      MSG.error("Solid::Solid(int)", "MP.soptsmax == 0");

    case 1:
      //special case: the point is set at w:
      addPointWithGrafted( w, RNG.pint_ratio(MP.MAX, MP.soghmax) );
      break;

    case 2:
      //special case: we put the two point symmetrically around w:
      x = MP.sosize[0] * Vecteur::randNormed();
      addPointWithGrafted( w + x, RNG.pint_ratio(MP.MAX, MP.soghmax));
      addPointWithGrafted( w - x, RNG.pint_ratio(MP.MAX, MP.soghmax));
      break;

    default:  
      //general case: the points are put randomly around w:
      for( int n = 0; n < MP.soptsmax; ++n ) {
        switch( MP.soshape ) {
          case 0: x = MP.sosize[0] * Vecteur::randSphere(); break;
          case 1: x = MP.sosize[0] * Vecteur::randNormed(); break;
          case 2: x = MP.sosize[0] * Vecteur::random(); break;
          default: MSG.error("Solid::Solid()", "illegal value for MP.soshape");
        }
        addPointWithGrafted( w + x, RNG.pint_ratio(MP.MAX, MP.soghmax));
      }
      break;
    }

  //store the current shape, (see fixShape() below):
  fixShape();
}


//---------------------------------------------------------------------------
Solid::~Solid()
{
  MSG(13, "Solid::destructor k%lx\n", name);
  deallocate();
}


//mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
//---------------------------------------------------------------------------
//wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

void Solid::fixShape()
// copy the current shape in the reference array soShape,
// also calculates the moment of inertia.
// this reference will be used by 'reshape()', which will restore this shape,
// while preserving current position and orientation of the solid (see below)
{
  calculatePosition();  // calculate the center of gravity in pscenter[]


  //set reference to current shape translated by -G (center) :
  //calculate the momentum soShapeMomentum needed in rescale():
  soShapeMomentum = 0;
  for(int p = 0, d = 0; p < DIM * nbPoints(); ++p, d = (d+1) % DIM ) {
    soShape[ p ] = pspts[ p ] - pscenter[ d ];
    soShapeMomentum += soShape[ p ] * soShape[ p ];
  }
  //we store the current number of points, for checking in reshape()
  soShapeSize = nbPoints();
}


//---------------------------------------------------------------------------
void Solid::rescale()
  // rescale the current positions to have the same 'size' as the
  // reference shape. size is measured as sum( ( x - g )^2 )
  // this is a 'sloppy' but faster version of 'reshape',
  // correcting only for a specific kind of global deformation.
{
  // for zero or one point we do nothing:
  if ( nbPoints() < 2 ) 
    return;

  int p, d;

  // calculate the center of gravity:
  calculatePosition();

  // we translate the object to subtract the center of gravity.
  // at the same time, we calculate the momentum of the current position:
  real current_size = 0;
  for( p = 0, d = 0; p < DIM * nbPoints(); ++p, d = (d+1) % DIM )
  {
    pspts[ p ] -= pscenter[ d ]; 
    current_size += pspts[ p ] * pspts[ p ];
  }

  if ( current_size <= 0. ) return;

  // calculate the scaling factor to restore the size to 'soShapeMomentum':
  real scale = sqrt( soShapeMomentum / current_size );

  // restore the shape, scaled from the center of gravity:
  for( p = 0, d = 0; p < DIM * nbPoints(); ++p, d = (d+1) % DIM )
    pspts[ p ] = scale * pspts[ p ] + pscenter[ d ];
}


//---------------------------------------------------------------------------
// reshape_really() finds the best ( rotation + translation ) to bring
// the reference shape stored in soShape[] onto the current shape in pspts[],
// it then set pspts[] to the transformed soShape[]. The result is to restore
// the shape in soShape[] while preserving place and orientation of the object
// 
// the best translation is the ones that conserves the center of gravity,
// the best rotation is obtained differently in 2D and 3D, see below.

#if ( DIM == 1 )

void Solid::reshape_really()
{
  // for zero or one point we do nothing:
  if ( nbPoints() < 2 ) 
    return; 

  //we check that the number of points is the same as when fixShape() was called.
  if ( soShapeSize != nbPoints() )
    MSG.error("Solid::reshape_really()", "Number of points missmatch between soShape and current");
  
  // we just have to translate soShape[] to the current center of gravity
  // no rotation.

  calculatePosition();
  
  for( int pp = 0; pp < nbPoints(); ++pp )
    pspts[ pp ] = soShape[ pp ] + pscenter[ 0 ];
}

#elif ( DIM == 2 )

void Solid::reshape_really()
{
  //  for zero or one point we do nothing:
  if ( nbPoints() < 2 ) 
    return;

  //we check that the number of points is the same as when fixShape() was called.
  if ( soShapeSize != nbPoints() )
    MSG.error("Solid::reshape_really()", "Number of points missmatch between soShape and current");
  
  //  calculate the center of gravity:
  calculatePosition();

  //  the best rotation is obtained by simple math on the cross products
  // and vector products of soShape[] and pspts[]: (see it on paper)

  real a = 0, b = 0;

  for( int pp = 0; pp < nbPoints(); ++pp ) {
      a += pspts[ DIM*pp     ] * soShape[ DIM*pp     ]
        +  pspts[ DIM*pp + 1 ] * soShape[ DIM*pp + 1 ];
      
      b += soShape[ DIM*pp     ] * pspts[ DIM*pp + 1 ]
        -  soShape[ DIM*pp + 1 ] * pspts[ DIM*pp     ];
    }
  
  real n = sqrt( a * a + b * b );

  // cosine and sinus of the rotation:
  real c = 1, s = 0;
  if ( n > EPSILON ) {
    c = a / n;
    s = b / n;
  }

  //printf(" n %8.3f, c %8.3f, s %8.3f norm = %8.3f\n", n, c, s, c*c + s*s);

  // apply transformation = rotation + translation:

  for( int pp = 0; pp < nbPoints(); ++pp ) {
    pspts[DIM*pp ]   = c * soShape[ DIM*pp ] - s * soShape[ DIM*pp+1 ]
      + pscenter[0];
    pspts[DIM*pp+1 ] = s * soShape[ DIM*pp ] + c * soShape[ DIM*pp+1 ]
      + pscenter[1];
  }
}

#elif ( DIM == 3 )

// We follow the procedure described by Berthold K.P. Horn in
// "Closed-form solution of absolute orientation using unit quaternions"
// Journal of the optical society of America A, Vol 4, Page 629, April 1987
// this is quite a heavy computation.
void Solid::reshape_really()
{
  // for zero or one point we do nothing:
  if ( nbPoints() < 2 )
    return;

  //we check that the number of points is the same as when fixShape() was called.
  if ( soShapeSize != nbPoints() )
    MSG.error("Solid::reshape_really()", "Number of points missmatch between soShape and current");
  
  // Calculate the center of gravity:
  calculatePosition();

  int pp, dd, ee;

  static real S[ 3 * 3 ], N[ 4 * 4 ];

  for( ee = 0; ee < 9; ++ee )
    S[ ee ] = 0;

  for( pp = 0; pp < nbPoints(); ++pp )
    for( dd = 0; dd < DIM; ++dd )
      for( ee = 0; ee < DIM; ++ee )
        S[ dd + DIM*ee ] += soShape[ DIM*pp + dd ] * pspts[ DIM*pp + ee ];

  N[ 0 + 4 * 0 ] = S[ 0 + DIM * 0 ] + S[ 1 + DIM * 1 ] + S[ 2 + DIM * 2 ];
  N[ 0 + 4 * 1 ] = S[ 1 + DIM * 2 ] - S[ 2 + DIM * 1 ];
  N[ 0 + 4 * 2 ] = S[ 2 + DIM * 0 ] - S[ 0 + DIM * 2 ];
  N[ 0 + 4 * 3 ] = S[ 0 + DIM * 1 ] - S[ 1 + DIM * 0 ];
  N[ 1 + 4 * 1 ] = S[ 0 + DIM * 0 ] - S[ 1 + DIM * 1 ] - S[ 2 + DIM * 2 ];
  N[ 1 + 4 * 2 ] = S[ 0 + DIM * 1 ] + S[ 1 + DIM * 0 ];
  N[ 1 + 4 * 3 ] = S[ 2 + DIM * 0 ] + S[ 0 + DIM * 2 ];
  N[ 2 + 4 * 2 ] = S[ 1 + DIM * 1 ] - S[ 0 + DIM * 0 ] - S[ 2 + DIM * 2 ];
  N[ 2 + 4 * 3 ] = S[ 1 + DIM * 2 ] + S[ 2 + DIM * 1 ];
  N[ 3 + 4 * 3 ] = S[ 2 + DIM * 2 ] - S[ 1 + DIM * 1 ] - S[ 0 + DIM * 0 ];


  //we call lapack to find the largest engeinvalue, and corresponding vector,
  //which is the quaternion corresponding to the rotation we want:

  int nbvalues;
  real eigenvalue;
  static Quaternion quat;
  static real work[ 8*4 ];
  static int iwork[ 5*4 ];
  static int ifail[4];
  int info = 0;

  lapack_xsyevx('V','I','U', 4, N, 4, 0, 0, 4, 4, EPSILON,
                &nbvalues, &eigenvalue, quat, 4, work, 8*4, iwork, ifail, &info );

  if ( info ) {
    MSG.warning("Solid::reshape_really()", "lapack_xsyevx() failed with code %i", info); 
    return;
  }
  //MSG("optimal LWORK = %i\n", work[0] );
  //MSG("eigen value %6.2f,", eigenvalue);
  //quat.println();

  //get the rotation matrix corresponding to the quaternion:
  quat.setThisMatrix33( S );
  //Matrix33( S ).print();

  //apply the transformation = rotation + translation:
  for( pp = 0; pp < nbPoints(); ++pp ) {
    blas_xcopy(DIM, pscenter, 1, pspts+DIM*pp, 1);
    blas_xgemv('N', DIM, DIM, 1.0, S, DIM, soShape+DIM*pp, 1, 1.0,
               pspts+DIM*pp, 1);
  }
}
#endif    // ( DIM == 3 )


//---------------------------------------------------------------------------
// Solid::reshape() is a combination of reshape_really() and rescale()
// reshape_really() is CPU intensive, and is called every 10 steps.
// we desinchronize the call by comparing with the name of the solid,
// to get a smooth simulation
void Solid::reshape()
{
  if ( sim.iterationCnt() % 10 == name % 10 )
    reshape_really();
  else
    rescale();
}


//===========================================================================
//zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
//===========================================================================


void Solid::setMobility()
{
  //TODO: make solid with constant mobility
  psMobility = nbPoints() / ( 6 * PI * MP.visc * MP.sosize[0] );
  //psMobility *= 5;   //TODO: increase mobility of solids ?
}


//---------------------------------------------------------------------------
#if ( DIM == 1 )

//The projection in 1D is just summing all the forces, and distributing
//equally to all the points:

void Solid::prepareProjection() {}

void Solid::setProjectedForces(const real * X, real * Y) const
{
  real fo = 0;
  for( int p = 0; p < nbPoints(); ++p )
    fo += X[ p ];
  fo /= real( nbPoints() );
  for( int p = 0; p < nbPoints(); ++p )
    Y[ p ] = fo;
}

#elif ( DIM == 2 )

//---------------------------------------------------------------------------
// To project in 2D or 3D, we calculate the resulting tensor by summing all
// the forces on all points, reducing it at the center of gravity.
// from this, we can deduce the forces compatible with solid motion, which
// are only translation + rotation.

void Solid::prepareProjection()
{
  if ( nbPoints() == 1 ) 
    return;

  // we have a special case with two points
  if ( nbPoints() == 2 ) 
    {
      soAxis = Vecteur( pspts+DIM ) - Vecteur( pspts );
      soMom[ 0 ] = 1.0 / ( 2 * soAxis.normSquare() );
      return;
    }
  
  // calculate center of gravity:
  calculatePosition();

  // calculate the rotational momentum of inertia:
  soMom[ 0 ] = 0;
  for(int p = 0; p < DIM * nbPoints(); ++p )
    soMom[ 0 ] += pspts[ p ] * pspts[ p ];
  soMom[ 0 ] -= nbPoints() * pscenter.normSquare();
}


void Solid::setProjectedForces(const real * X, real * Y) const
{
  int p;

  if ( nbPoints() == 1 ) {
    for( p = 0; p < DIM; ++p )
      Y[ p ] = X[ p ];
    return;
  }

  Vecteur f;           //resulting force
  
  //special case with two points:
  if ( nbPoints() == 2 ) {
      real l; 
      l = soMom[ 0 ] * ( ( X[ DIM   ] - X[ 0 ] ) * soAxis.XX +
			 ( X[ DIM+1 ] - X[ 1 ] ) * soAxis.YY );

      f.XX = l * soAxis.XX;
      f.YY = l * soAxis.YY;

      Y[ 0 ]       = ( X[ 0 ]       + f.XX );
      Y[ 1 ]       = ( X[ 1 ]       + f.YY );
      Y[ 0 + DIM ] = ( X[ 0 + DIM ] - f.XX );
      Y[ 1 + DIM ] = ( X[ 1 + DIM ] - f.YY );
      return;
    }

  f = VZERO;                 //resulting force
  CrossProduct  m = MZERO;   //resulting torque

  for( p = 0; p < nbPoints(); ++p ) {
    f.XX += X[ p*DIM     ];
    f.YY += X[ p*DIM + 1 ];
    m += pspts[ p*DIM ] * X[ p*DIM+1 ] - pspts[ p*DIM+1 ] * X[ p*DIM ];
  }
  
  m += f ^ pscenter;                  //reduce to the center of gravity
  m /= soMom[ 0 ];                    //final rotation term
  f /= real( nbPoints() );
  f += pscenter ^ m;
  
  for( p = 0; p < nbPoints(); ++p ) {
    Y[ p*DIM     ] = f.XX - m * pspts[ p*DIM + 1 ];
    Y[ p*DIM + 1 ] = f.YY + m * pspts[ p*DIM     ];
  }
}  

#elif ( DIM == 3 )

//---------------------------------------------------------------------------
void Solid::prepareProjection()
{
  if ( nbPoints() == 1 ) 
    return;
  if ( nbPoints() == 2 ) 
    {
      soAxis = Vecteur( pspts+DIM ) - Vecteur( pspts );
      soMom[ 0 ] = 1.0 / ( 2 * soAxis.normSquare() );
      return;
    }

  //Calculate center of gravity:
  calculatePosition();

  //Calculate the rotational momentum of inertia:
  for(int ii = 0; ii < DIM * DIM; ++ii )
    soMom[ ii ] = 0;

  for(int p = 0; p < nbPoints(); ++p ) {
    real * pos = pspts + DIM * p;
    soMom[ 0 + DIM * 0 ] += pos[ 0 ] * pos[ 0 ];
    soMom[ 0 + DIM * 1 ] += pos[ 0 ] * pos[ 1 ];
    soMom[ 0 + DIM * 2 ] += pos[ 0 ] * pos[ 2 ];
    soMom[ 1 + DIM * 1 ] += pos[ 1 ] * pos[ 1 ];
    soMom[ 1 + DIM * 2 ] += pos[ 1 ] * pos[ 2 ];
    soMom[ 2 + DIM * 2 ] += pos[ 2 ] * pos[ 2 ];
  }

  //transform it to the center of gravity (cf. equations of motion):

  real diag = soMom[0] + soMom[4] + soMom[8] - nbPoints() * pscenter.normSquare() ;

  soMom[ 0 + DIM * 0 ] = diag - soMom[ 0 + DIM * 0 ] + nbPoints() * pscenter[ 0 ] * pscenter[ 0 ];
  soMom[ 0 + DIM * 1 ] =      - soMom[ 0 + DIM * 1 ] + nbPoints() * pscenter[ 0 ] * pscenter[ 1 ];
  soMom[ 0 + DIM * 2 ] =      - soMom[ 0 + DIM * 2 ] + nbPoints() * pscenter[ 0 ] * pscenter[ 2 ];
  soMom[ 1 + DIM * 1 ] = diag - soMom[ 1 + DIM * 1 ] + nbPoints() * pscenter[ 1 ] * pscenter[ 1 ];
  soMom[ 1 + DIM * 2 ] =      - soMom[ 1 + DIM * 2 ] + nbPoints() * pscenter[ 1 ] * pscenter[ 2 ];
  soMom[ 2 + DIM * 2 ] = diag - soMom[ 2 + DIM * 2 ] + nbPoints() * pscenter[ 2 ] * pscenter[ 2 ];


  //Matrix33(soMom).print();
  // for a 3D object, the matrix should be symmetrix positive definite
  // try to compute the cholesky factorization:
  int info;
  lapack_xpotf2('U', DIM, soMom, DIM, &info);
  if ( info ) {
    MSG.error("solid.cc","solid has aligned points"); //TODO deal with this special case
  }

  /*
  //calculate matrix inverse
  lapack_xpotri('U', DIM, soMom, DIM, &info);
  assert( info == 0 );
  */
}

//---------------------------------------------------------------------------
void Solid::setProjectedForces(const real * X, real * Y) const
  // Y <- P * X
  // P is the projection associated with the constraints of motion without
  // deformation (solid object), 
  // we calculate the total force and momentum in zero, and distribute
  // it according to solid motion mechanics.
{
  int p;

  if ( nbPoints() == 1 ) {
    for( p = 0; p < DIM; ++p ) 
      Y[ p ] = X[ p ];
    return;
  }

  // in 3D, with two points, the matrix of momentum soMom is singular:
  // we have a special case:
  if ( nbPoints() == 2 ) {
    real l = 0;
    for( p = 0; p < DIM; ++p )
      l += ( X[ p+DIM ] - X[ p ] ) * soAxis[ p ];
    
    Vecteur internal = ( soMom[ 0 ] * l ) * soAxis;
    
    for( p = 0; p < DIM; ++p ) {
      Y[ p ]       = ( X[ p ]       + internal[ p ] );
      Y[ p + DIM ] = ( X[ p + DIM ] - internal[ p ] );
    }
    return;
  }
 
  Vecteur       f = VZERO;           //resulting force
  CrossProduct  m = MZERO;           //resulting momentum

  for( p = 0; p < nbPoints(); ++p ) {
    f.XX += X[ p*DIM     ];
    f.YY += X[ p*DIM + 1 ];
    f.ZZ += X[ p*DIM + 2 ];
    m.XX += pspts[ p*DIM+1 ] * X[ p*DIM+2 ] - pspts[ p*DIM+2 ] * X[ p*DIM+1 ];
    m.YY += pspts[ p*DIM+2 ] * X[ p*DIM   ] - pspts[ p*DIM   ] * X[ p*DIM+2 ];
    m.ZZ += pspts[ p*DIM   ] * X[ p*DIM+1 ] - pspts[ p*DIM+1 ] * X[ p*DIM   ];
  }

  m += f ^ pscenter; //reduce to the center of gravity
		     
  //solve the 3x3 linear system, using the factorization stored in somon:
  int info;
  lapack_xpotrs('U', DIM, 1, soMom, DIM, m, DIM, &info );
  assert( info == 0 );
  f /= real( nbPoints() );
  f += pscenter ^ m;

  for( p = 0; p < nbPoints(); ++p ) {
    Y[ p*DIM   ] = f.XX + m.YY * pspts[ p*DIM+2 ] - m.ZZ * pspts[ p*DIM+1 ];
    Y[ p*DIM+1 ] = f.YY + m.ZZ * pspts[ p*DIM   ] - m.XX * pspts[ p*DIM+2 ];
    Y[ p*DIM+2 ] = f.ZZ + m.XX * pspts[ p*DIM+1 ] - m.YY * pspts[ p*DIM   ];
  }
}

#endif


void  Solid::prepareProjectionDiff(const real * forces)        {}
void  Solid::addProjectionDiff(const real * X, real * Y) const {}

//***************************************************************************
//***************************************************************************
//***************************************************************************

void Solid::addPointWithGrafted(Vecteur w, int ty)
{
  linkGrafted( new Grafted( ty, this, addPoint( w ) ) );
}


//---------------------------------------------------------------------------
void Solid::write()
{
  assert( name > 0 );
  IO.writeRecordTag("so");
  IO.writeUInt32( name );
  IO.writeUInt8( soType );
  IO.writeUInt8( 0 );                 // reserved for later use
  IO.writeReal32( 0 );               // reserved for later use
  PointSet::write();
}


//---------------------------------------------------------------------------
void Solid::read()
{
  setFlag( sim.frameInBuffer() );
  try {

    soType     = IO.readUInt8();
    IO.readUInt8();
    IO.readReal32();
    PointSet::read();

  } catch( IOException e ) {

    //MSG("S%lx %i ", name, r);
    e.addBeforeMessage("Solid::read : ");
    clearPoints();
    throw e;

  }

  fixShape();
}

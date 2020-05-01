//RCS: $Id: fiber_solve.cc,v 1.5 2005/04/25 13:30:48 nedelec Exp $
//-----------------------------fiber_solve.cc------------------------
//--------------------   functions used in solve   ------------------
//-------------------------------------------------------------------

#include "cblas.h"
#include "clapack.h"
#include "matrix.h"

#define CRAZY_OPTIMIZATION

void Fiber::setMobility()
{
  psMobility = nbPoints() * HYDRO / ( 4 * PI * MP.visc * length() );
  //All the microtubules have the same rigidity, defined by MP.mtrigid
  //the scaling of the bending elasticity depends on the length of the segments:
  mtrigidkm  = MP.mtrigid / ( cub( mtcut ) * MP.km );
  
  //FIX EXPERIMENTAL: increased rigidity of KINETOCHORE microtubules (Vale & Gohta, Jan 2005)
  if ( getType() == MAGIC_TYPE )
    mtrigidkm *= 10;
}


//-----------------------------------------------------------------------

#if ( DIM > 1 )

#ifdef PROJECT_WITH_MATRIX
  #include "fiber_projectmat.cc"
#else
  #include "fiber_project.cc"
#endif

#else   // DIM == 1 is a trivial case:

void Fiber::constructorProjection() {}
void Fiber::prepareProjection()     {}
void Fiber::deallocateProjection()  {}

void Fiber::setProjectedForces( const real * X, real * Y ) const
{
  real sum = 0;
  for(int ii = 0; ii < nbPoints(); ++ii )
    sum += X[ ii ];

  sum *= 1.0 / real( nbPoints() );
  for(int ii = 0; ii < nbPoints(); ++ii )
    Y[ ii ] = sum;
}

void Fiber::prepareProjectionDiff(const real*) {}
void Fiber::addProjectionDiff(const real*, real*) const {}

#endif


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//                        RIGIDITY CONTRIBUTION
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

/*
void Fiber::setRigidityUp(Matrix & mB, const int offset) const
  //only the upper diagonal is set
{
  for( int ii = 0; ii < nbPoints() - 2 ; ++ii )
    {
      mB( offset + ii   ,  offset + ii   ) += -1 * mtrigidkm; 
      mB( offset + ii   ,  offset + ii+1 ) += +2 * mtrigidkm; 
      mB( offset + ii   ,  offset + ii+2 ) += -1 * mtrigidkm; 
      mB( offset + ii+1 ,  offset + ii+1 ) += -4 * mtrigidkm; 
      mB( offset + ii+1 ,  offset + ii+2 ) += +2 * mtrigidkm; 
      mB( offset + ii+2 ,  offset + ii+2 ) += -1 * mtrigidkm; 
    }
}
*/

//set rigidity matrix elements, only the upper diagonal is set
void Fiber::setRigidityUp(Matrix & mB, const int offset) const
{
  int sz = nbPoints();
  if ( sz < 3 ) return;

  int ii;
  int s = offset;
  int e = s + sz;

  for (ii = s+2;  ii < e-2 ; ++ii )  mB( ii, ii   ) = -6*mtrigidkm;
  for (ii = s+1;  ii < e-2 ; ++ii )  mB( ii, ii+1 ) = +4*mtrigidkm;
  for (ii = s;    ii < e-2 ; ++ii )  mB( ii, ii+2 ) =   -mtrigidkm;

  mB( s  , s   ) = -mtrigidkm;
  mB( e-1, e-1 ) = -mtrigidkm;
  if ( sz == 3 )
    mB( s+1, s+1 ) = -4*mtrigidkm;
  else
    {
      mB( s+1, s+1 ) = -5*mtrigidkm;
      mB( e-2, e-2 ) = -5*mtrigidkm;
    }
  mB( s  , s+1 ) = 2*mtrigidkm;
  mB( e-2, e-1 ) = 2*mtrigidkm;
}


//==========================================================================
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//==========================================================================


#ifndef CRAZY_OPTIMIZATION

void Fiber::addRigidity(const real * X, real * Y) const
// Y <- Y + ( Rigidity Matrix ) * X
{
  if ( DIM < 2 ) return;
  if ( nbPoints() < 3 ) return;
  int kk, jj, d;

  //----- we allocate an array tmp of size DIM * nbPoints() at least:
  static real * tmp = 0;
  static int allocated = 0;
  if ( nbPoints() > allocated ) {
    if ( tmp ) delete[] tmp;
    allocated = nbPoints() + 2;
    tmp = new real[ DIM * allocated ];
    if ( tmp == 0 ) {
      fprintf(stderr, "Fiber::addRigidity(const real *, real*) memory allocation failed\n");
      exit(1);
    }    
  }
  

  //calculate the second-differential of points in tmp:
  //scale by the rigidity term:
  for( kk = 0, jj = 0; jj < lastSegment(); ++jj )
    for( d = 0; d < DIM; ++d ) {
      assert( kk == DIM*jj + d );
      tmp[ kk ] = mtrigidkm *( X[ kk ] - 2 * X[ kk+DIM ] + X[ kk+2*DIM ]);
      ++kk;
    }

  //printf("\n");
  //PrintVect( tmp, DIM * nbPoints() );

  //add back to the vector Y:
  for( jj = 0, kk=0; jj < lastSegment(); ++jj )
    for( d = 0; d < DIM; ++ d ) {
      assert( kk == DIM*jj + d );
      Y[ kk         ] -=   tmp[ kk ];
      Y[ kk + DIM   ] += 2*tmp[ kk ];
      Y[ kk + 2*DIM ] -=   tmp[ kk ];
      ++kk;
    }
}

#else //CRAZY_OPTIMIZATION

#if ( DIM == 1 )
void Fiber::addRigidity(const real * X, real * Y) const {}
#else
void Fiber::addRigidity(const real * X, real * Y) const
{
  if ( nbPoints() < 3 ) return;
  
  //----- we allocate an array tmp of size DIM * nbPoints() at least:
  static real * tmp = 0;
  static int allocated = 0;
  if ( nbPoints() > allocated ) {
    if ( tmp ) delete[] tmp;
    allocated = nbPoints() + 2;
    tmp = new real[ DIM * allocated ];
    if ( tmp == 0 ) {
      fprintf(stderr, "Fiber::addRigidity(const real *, real*) memory allocation failed\n");
      exit(1);
    }    
  };
  
  //calculate the second-differential of points in tmp:
  real * pW   = tmp;
  const real * pX   = X;
  const real * pend = pX + DIM * nbSegments();
  
  register real dx0, bufdx0 = pX[DIM+0] - pX[0];
  register real dx1, bufdx1 = pX[DIM+1] - pX[1];
#if ( DIM == 3 )
  register real dx2, bufdx2 = pX[DIM+2] - pX[2];
#endif
  pX += DIM;
  
  while( pX < pend ) {
    dx0 = pX[DIM+0] - pX[0];
    dx1 = pX[DIM+1] - pX[1];
#if ( DIM == 3 )
    dx2 = pX[DIM+2] - pX[2];
#endif

    pX += DIM;
    
    pW[0] = dx0 - bufdx0;
    pW[1] = dx1 - bufdx1;
#if ( DIM == 3 )
    pW[2] = dx2 - bufdx2;
#endif
   pW += DIM;
   bufdx0 = dx0;
   bufdx1 = dx1;
#if ( DIM == 3 )
   bufdx2 = dx2;
#endif

  }
  pW[0] = 0;
  pW[1] = 0;
#if ( DIM == 3 )
  pW[2] = 0;
#endif

  //printf("\n");
  //PrintVect( tmp, DIM * nbSegments() );
  
  //scale by the rigidity term:
  blas_xscal( DIM * lastSegment(), mtrigidkm, tmp, 1 );
  
  //add the rigidity terms to Y:
  real * pY = Y + DIM;
  pW   = tmp + DIM;
  pend = tmp + DIM * lastPoint();
  
  register real w0 = tmp[0], w3;
  register real w1 = tmp[1], w4;
#if ( DIM == 3 )
  register real w2 = tmp[2], w5;
#endif
  
  Y[0] -= w0;
  Y[1] -= w1;  
#if ( DIM == 3 )
  Y[2] -= w2;  
#endif
  
  bufdx0 = pY[0] + w0;
  bufdx1 = pY[1] + w1;
#if ( DIM == 3 )
  bufdx2 = pY[2] + w2;
#endif
  

  while( pW < pend ) {
    w3  = pW[0];
    w4  = pW[1];
#if ( DIM == 3 )
    w5  = pW[2];
#endif
    pW += DIM;
    
    dx0 = w3 - w0;
    w0  = w3;
    dx1 = w4 - w1;
    w1  = w4;
#if ( DIM == 3 )
    dx2 = w5 - w2;
    w2  = w5;
#endif

    pY[0] = bufdx0 - dx0;
    pY[1] = bufdx1 - dx1;
#if ( DIM == 3 )
    pY[2] = bufdx2 - dx2;
#endif

    pY += DIM;
    
    bufdx0 = pY[0] + dx0;
    bufdx1 = pY[1] + dx1;
#if ( DIM == 3 )
    bufdx2 = pY[2] + dx2;
#endif

  }
  pY[0] = bufdx0;
  pY[1] = bufdx1;
#if ( DIM == 3 )
  pY[2] = bufdx2;
#endif

  //PrintVect( Y, DIM * nbPoints() );
}
#endif
#endif 


#ifdef PROJECT_WITH_BLAS

//this is old unused code, but it might be useful some day. 
//Check before using!

void Fiber::addRigidity(const real * X, real * Y) const
// Y <- Y + ( Rigidity Matrix ) * X
{
  if ( DIM < 2 ) return;
  if ( nbPoints() < 3 ) return;

  //----- we allocate an array tmp of size DIM * nbPoints() at least:
  static real * tmp = 0;
  static int allocated = 0;
  if ( nbPoints() > allocated ) {
    if ( tmp ) delete[] tmp;
    allocated = nbPoints() + 2;
    tmp = new real[ DIM * allocated ];
    if ( tmp == 0 ) {
      fprintf(stderr, "Fiber::addRigidity(const real *, real*) memory allocation failed\n");
      exit(1);
    }    
  }

  //calculate the second-differential of points in tmp:
  blas_xcopy( DIM*nbPoints(),    X, 1, tmp, 1 );
  blas_xaxpy( DIM*nbSegments(),  -1.0, tmp+DIM, 1, tmp, 1 );
  blas_xaxpy( DIM*lastSegment(), -1.0, tmp+DIM, 1, tmp, 1 );
  
  //scale by the rigidity term:
  blas_xscal( DIM*lastSegment(), mtrigidkm, tmp, 1 );

  //add back to the vector Y:
  blas_xaxpy( DIM*lastSegments(), -1.0, tmp, 1, Y, 1 );
  blas_xaxpy( DIM*lastSegments(), +2.0, tmp, 1, Y+DIM, 1 );
  blas_xaxpy( DIM*lastSegments(), -1.0, tmp, 1, Y+2*DIM, 1 );
}

#endif  //PROJECT_WITH_BLAS

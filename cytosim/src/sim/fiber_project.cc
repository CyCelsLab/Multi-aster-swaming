//RCS: $Id: fiber_project.cc,v 1.4 2005/04/25 13:30:40 nedelec Exp $
//----------------------------------------------------------------------------
//------------------- projections for constrained dynamics  ------------------
//----------------------------------------------------------------------------

void Fiber::constructorProjection()
{
  //reset all variables for the projections:
  allocatedP        = 0;
  mtdpts            = 0;  
  mtJJt             = 0;
  mtJJtiJforce      = 0;
  mtJJtiJforce_null = true;
}

//----------------------------------------------------------------------------
void Fiber::allocateProjection()
{
  int nbv = DIM * nbPoints();
  if ( allocatedP < nbv ) {
      MSG(12, "Fiber::allocateP M%lx %d\n", name, nbv );
      if ( mtdpts )       delete[] mtdpts;
      if ( mtJJt )        delete[] mtJJt;
      if ( mtJJtiJforce ) delete[] mtJJtiJforce;
      
      allocatedP   = nbv + 2*DIM;  //we use a chunck-size of 2 points
      
      mtdpts       = new real[ allocatedP ];
      mtJJt        = new real[ 2*( allocatedP / DIM ) ];
      mtJJtiJforce = new real[ allocatedP / DIM ];
      
      if (( mtdpts == 0 ) || ( mtJJt == 0 ) || ( mtJJtiJforce == 0 )) {
        fprintf(stderr, "Fiber::allocateProjection() memory allocation failed\n");
        exit(1);
      }      
    }
}

//----------------------------------------------------------------------------
void Fiber::deallocateProjection()
{
  MSG(12, "Fiber::deallocateProjection M%lx\n", name );
  if ( mtdpts )       delete[] mtdpts;
  if ( mtJJt )        delete[] mtJJt;
  if ( mtJJtiJforce ) delete[] mtJJtiJforce;
  mtdpts       = 0;
  mtJJt        = 0;
  mtJJtiJforce = 0;
}

//----------------------------------------------------------------------------
// prepare variables for the projection
void Fiber::prepareProjection()
{
  //------------first allocate the matrix if needed:
  allocateProjection();
 
  //----------- calculate the vectors dpts, vPd and vPe used in the inversion
  register real x, y[DIM], dxx, dxy;

  for(int d = 0; d < DIM ; ++d )
    y[ d ] = 0.0;

  for(int kk = 0, jj = 0; jj < nbSegments() ; ++jj ) {
    dxx = 0.0;
    dxy = 0.0;

    kk = DIM * jj;
    for(int d = 0; d < DIM ; ++d )	{
      x = pspts[ kk + d + DIM ] - pspts[ kk + d ];
      mtdpts[ kk + d ] = x;
      dxx += x * x;
      dxy += x * y[ d ];
      y[ d ] = x;
    }

    //set the diagonal and off-diagonal of J*J'
    mtJJt[ jj ]               = 2 * dxx;
    mtJJt[ jj + nbPoints() ]  =   - dxy;  //the first value for jj=0 is discarded
  }

  int info = 0;
  lapack_xpttrf( nbSegments(), mtJJt, mtJJt + nbPoints()+1, &info );
  if ( info ) 
    MSG.error("Fiber::prepareProjection()", "lapack_pttrf() failed");
}


//----------------------------------------------------------------------------
#ifndef CRAZY_OPTIMIZATION
void Fiber::setProjectedForces(const real * X, real * Y) const
{
   int nbc = nbSegments();  //number of constraints

   //----- we allocate a temporary array tmp of size nbc at least:
   static real * tmp = 0;
   static int allocated = 0;
   if ( nbc > allocated ) {
     if ( tmp ) delete[] tmp;
     allocated = nbc + 2;
     tmp = new real[ allocated ];
     if ( tmp == 0 ) {
       fprintf(stderr, "Fiber::setProjectedForces(const real *, real*) memory allocation failed\n");
       exit(1);
     }     
   }

   //-------   tmp <- J * X
   for(int jj = 0, kk = 0; jj < nbc; ++jj, kk+=DIM ) {
#if ( DIM == 2 )
     tmp[ jj ] = mtdpts[ kk   ] * ( X[ kk+DIM   ] - X[ kk   ] )
               + mtdpts[ kk+1 ] * ( X[ kk+DIM+1 ] - X[ kk+1 ] );
#elif ( DIM == 3 )
     tmp[ jj ] = mtdpts[ kk   ] * ( X[ kk+DIM   ] - X[ kk   ] )
               + mtdpts[ kk+1 ] * ( X[ kk+DIM+1 ] - X[ kk+1 ] )
               + mtdpts[ kk+2 ] * ( X[ kk+DIM+2 ] - X[ kk+2 ] );
#endif
   }
   
   // tmp <- inv( J * Jt ) * tmp to find the lagrangeMult multipliers
   int info = 0;
   lapack_xpttrs(nbc, 1, mtJJt, mtJJt+nbPoints()+1, tmp, nbc, &info);
   if ( info )
     MSG.error("Fiber::setProjectedForces()", "lapack_pttrs() failed");

   //printf("lagrangeMult "); vectPrint( nbc, tmp );

   // Y <- X + Jt * tmp :
   // Y = X + mtdpts * tmp[ /DIM] - mtdpts * tmp[ /DIM]
   
   register real x, buf[DIM];
   for(int d = 0; d < DIM; ++d ) 
     buf[d] = 0;
   
   int kk = 0;
   for(int jj = 0; jj < nbc; ++jj ) {
     for(int d = 0; d < DIM; ++d ) {
       x      = mtdpts[kk] * tmp[jj];
       Y[kk]  = X[kk] - buf[d] + x;
       buf[d] = x;
       ++kk;
     }
   }
   for(int d = 0; d < DIM; ++d, ++kk ) 
     Y[kk] = X[kk] - buf[d];
   assert( kk == DIM*nbPoints() );
   //printf("Y  "); vectPrint( DIM * nbPoints(), Y );
 }



#else  //CRAZY_OPTIMIZATION



void Fiber::setProjectedForces(const real * X, real * Y) const
{
  int nbc = nbSegments();  //number of constraints
  
  //----- we allocate a temporary array tmp of size nbc at least:
  static real * tmp = 0;
  static real * dpx = 0;
  static int allocated = 0;
  if ( nbc > allocated ) {
    if ( tmp ) delete[] tmp;
    allocated = nbc + 2;
    tmp = new real[ allocated ];
    dpx = new real[ DIM*allocated ];
    if ( tmp == 0 ) {
      fprintf(stderr, "Fiber::setProjectedForces(const real *, real*) memory allocation failed\n");
      exit(1);
    }     
  }
  

  real * pT = tmp;
  const real * pX = X + DIM;
  const real * pM = mtdpts;
  register real x3, x0 = X[0];
  register real x4, x1 = X[1];
#if ( DIM == 3 )
  register real x5, x2 = X[2];
#endif
  
  for(int jj = 0; jj < nbc; ++jj) {
    x3 = pX[0];
    x4 = pX[1];
#if ( DIM == 2 )
    pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1);
#elif ( DIM == 3 )
    x5 = pX[2];
    pT[0] = pM[0] * (x3 - x0) + pM[1] * (x4 - x1) + pM[2] * (x5 - x2);
    x2 = x5;
#endif
    ++pT;
    pX += DIM;
    pM += DIM;
    x0 = x3;
    x1 = x4;
   }
  
  // tmp <- inv( J * Jt ) * tmp to find the lagrangeMult multipliers
  int info = 0;
  lapack_xpttrs(nbc, 1, mtJJt, mtJJt+nbPoints()+1, tmp, nbc, &info);
  if ( info )
    MSG.error("Fiber::setProjectedForces()", "lapack_pttrs() failed");
  
  //printf("lagrangeMult "); vectPrint( nbc, tmp );
  
  
  real * pY = Y;
  pX = X + DIM;
  pM = mtdpts;
  
  x3 = X[0];
  x4 = X[1];
#if ( DIM == 3 )
  x5 = X[2];
#endif
  
  for(int kk = 0; kk < nbc; ++kk) {

    x0  = pM[0] * tmp[kk];
    x1  = pM[1] * tmp[kk];
#if ( DIM == 3 )
    x2  = pM[2] * tmp[kk];
#endif
    pM += DIM;
    
    pY[0]  = x3    + x0;
    x3     = pX[0] - x0;
    
    pY[1]  = x4    + x1;
    x4     = pX[1] - x1;
    
#if ( DIM == 3 )
    pY[2]  = x5    + x2;
    x5     = pX[2] - x2;
#endif
    
    pX += DIM;
    pY += DIM;
  }

  pY[0] = x3;
  pY[1] = x4;
#if ( DIM == 3 )
  pY[2] = x5;
#endif  
  //printf("Y  "); vectPrint( DIM * nbPoints(), Y );
}
#endif  //CRAZY_OPTIMIZATION


//----------------------------------------------------------------------------
//========================= PROJECTION DERIVATIVE ============================
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv



//----------------------------------------------------------------------------
void Fiber::prepareProjectionDiff( const real * force )
{
  int jj, kk;
  int nbc = nbSegments();                 //number of constraints

  //-------   mtJJtiJforce <- J * X
  for(jj = 0; jj < nbc; ++jj ) {
    kk = DIM * jj;
    mtJJtiJforce[ jj ] = mtdpts[ kk ] * ( force[ kk + DIM ] - force[ kk ] );
    for( int dd = 1; dd < DIM; ++dd ) //we start from 1 here!
      mtJJtiJforce[ jj ] += mtdpts[ kk+dd ] * ( force[ kk+dd + DIM ] - force[ kk+dd ] );
  }
  
  // mtJJtiJforce <- inv( J * Jt ) * mtJJtiJforce (lagrangeMult multipliers)
  int info = 0;
  lapack_xpttrs(nbc, 1, mtJJt, mtJJt+nbPoints()+1, mtJJtiJforce, nbc, &info);
  if ( info ) 
    MSG.error("Fiber::prepareProjectionDiff()", "lapack_pttrs() failed");
  
  
  //printf("all  lagrangeMult "); vectPrint( nbc, mtJJtiJforce );
  blas_xcopy( nbc, mtJJtiJforce, 1, lagrangeMult, 1 );
  lagrangeValid = true;

  //----- we remove compressive forces ( lagrangian < 0 )

  //experimental: instead of 0, where have a small positive threshold
  //which was numerically hand-determinated: we found that if one of
  //the lagrangeMult multiplier was above 0.05 / km * dt, it was unstable
  //so we keep here all above 0.01... gain -20% CPU in some situations

  mtJJtiJforce_null = true;
  //const real stab = 0.01 / ( MP.km * MP.dt );  //TODO: verify that optimization
  const real stab = 0.0;                         //this is the safer choice (a bit slower)
  for( jj = 0; jj < nbc; ++jj )
    if ( mtJJtiJforce[ jj ] < stab )
      mtJJtiJforce[ jj ] = 0.0;
    else
      mtJJtiJforce_null = false;

  /*
  if ( !mtJJtiJforce_null ) {
    printf("\ngood Lagrange multipliers ");
    for( jj = 0; jj < nbc; ++jj )
      printf( "%6.3f ", mtJJtiJforce[jj]*MP.km*MP.dt);
    mtJJtiJforce_null = true;
  }
   */
}


//----------------------------------------------------------------------------
#ifndef CRAZY_OPTIMIZATION

void Fiber::addProjectionDiff( const real * X, real * Y ) const
{
  if ( mtJJtiJforce_null )
    return;
  
  //straightforward implementation:
  for( int jj = 0; jj < nbSegments(); ++jj ) {
    if ( mtJJtiJforce[ jj ] ) {
      int ll = DIM*jj;
      for( int d = 0; d < DIM; ++d ) {
        real w = mtJJtiJforce[ jj ] * ( X[ ll+d+DIM ] - X[ ll+d ] );
        Y[ ll+d     ] += w;
        Y[ ll+d+DIM ] -= w;
      }
    }
  }
}

#else //CRAZY_OPTIMIZATION

void Fiber::addProjectionDiff( const real * X, real * Y ) const
{
  if ( mtJJtiJforce_null )
    return;
  
  const real * pX = X;
  real * pY = Y;
  register real w0, w1;
#if ( DIM == 3 )
  register real w2;
#endif
  
  int last = nbSegments(); 
  for( int jj = 0; jj < last; ++jj ) {
    if ( mtJJtiJforce[ jj ] ) {
      w0 = mtJJtiJforce[ jj ] * ( pX[ 0+DIM ] - pX[0] );
      w1 = mtJJtiJforce[ jj ] * ( pX[ 1+DIM ] - pX[1] );
#if ( DIM == 3 )
      w2 = mtJJtiJforce[ jj ] * ( pX[ 2+DIM ] - pX[2] );
#endif
      pX += DIM;

      pY[ 0     ] += w0;
      pY[ 0+DIM ] -= w0;
      pY[ 1     ] += w1;
      pY[ 1+DIM ] -= w1;
#if ( DIM == 3 )
      pY[ 2     ] += w2;
      pY[ 2+DIM ] -= w2;
#endif
      pY += DIM;
    } else {
      pX += DIM;
      pY += DIM;
    }
  }
}
#endif  //CRAZY_OPTIMIZATION


//this is old unused code, but it might be useful some day. 
//Check before using!

#ifdef PROJECT_WITH_BLAS
//----------------------------------------------------------------------------
void Fiber::setProjectedForces( const real * X, real * Y ) const
  // Y <- P * X = X + internal forces
{
  assert( DIM > 1 );
  //----- we allocate an array tmp of size DIM * nbPoints() at least:
  static real * tmp1 = 0;
  static real * tmp2 = 0;
  static int allocated = 0;
  if ( nbPoints() > allocated ) {
    if ( tmp1 ) delete[] tmp1;
    if ( tmp2 ) delete[] tmp2;
    allocated = nbPoints() + 2;
    tmp1 = new real[ DIM * allocated ];
    tmp2 = new real[ DIM * allocated ];
    if (( tmp1 == 0 ) || ( tmp2 == 0 )) {
      fprintf(stderr, "Fiber::setProjectedForces(const real *, real*) memory allocation failed\n");
      exit(1);
    }
  }

  blas_xcopy( DIM*nbSegments(), X, 1, tmp1, 1 );
  blas_xaxpy( DIM*nbSegments(), -1.0, X+DIM, 1, tmp1, 1 );

  assert( allocatedP > 0 );
  blas_xsbmv('U', nbSegments(), 0, -1.0, mtdpts, DIM, tmp1, DIM, 0.0, tmp2, 1 );
  for(int d = 1; d < DIM; ++d )
    blas_xsbmv('U', nbSegments(), 0, -1.0, mtdpts+d, DIM, tmp1+d, DIM, 1.0, tmp2, 1);
  
  int info = 0;
  lapack_xpttrs(nbSegments(), 1, mtJJt, mtJJt+nbPoints()+1, tmp2, nbSegments(), &info);
  assert( info == 0 );

  for(int d = 0; d < DIM; ++d )
    blas_xsbmv('U', nbSegments(), 0, -1.0, mtdpts+d, DIM, tmp2, 1, 0.0, tmp1+d, DIM);

  blas_xcopy( DIM*nbPoints(), X, 1, Y, 1 );
  blas_xaxpy( DIM*nbSegments(), -1.0, tmp1, 1, Y, 1 );
  blas_xaxpy( DIM*nbSegments(),  1.0, tmp1, 1, Y+DIM, 1 );
}

#endif  //PROJECT_WITH_BLAS

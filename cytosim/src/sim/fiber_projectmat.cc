//RCS: $Id: fiber_projectmat.cc,v 1.1 2005/04/10 14:42:16 nedelec Exp $
//------------------------------------------------------------------
//---------  projections performed with explicit matrices  ---------
//------------------------------------------------------------------

void Fiber::constructorProjection()
{
  //reset all variables for the projections:
  allocatedM     = 0;
  mtP            = 0;
  mtDiffP        = 0;
  mtJJtiJ        = 0;
}

//------------------------------------------------------------------
void Fiber::allocateProjection()
{
  int nbv = DIM * nbPoints();         //number of variables needed
  
  if ( allocatedM < nbv ) {
    MSG(12, "Fiber::allocateProjection MAT M%lx %d\n", name, nbv );
    if ( mtP )          delete[] mtP;
    if ( mtDiffP )      delete[] mtDiffP;
    if ( mtJJtiJ )      delete[] mtJJtiJ;
    allocatedM   = nbv;
    mtP          = new real[ allocatedM * allocatedM ];
    mtDiffP      = new real[ allocatedM * allocatedM ];
    mtJJtiJ      = new real[ allocatedM * nbSegments() ];
    
    if (( mtP == 0 ) || ( mtDiffP == 0 ) || ( mtJJtiJ == 0 )) {
      fprintf(stderr, "Fiber::allocateProjection() memory allocation failed\n");
      exit(1);
    }    
  }
}

//------------------------------------------------------------------
void Fiber::deallocateProjection()
{
  MSG(12, "Fiber::deallocateProjection MAT M%lx\n", name );
  if ( mtP )          delete[] mtP;
  if ( mtDiffP )      delete[] mtDiffP;
  if ( mtJJtiJ )      delete[] mtJJtiJ;
  mtP          = 0;
  mtDiffP      = 0;
  mtJJtiJ      = 0;
}


//------------------------------------------------------------------
// computes the derivation of the projection matrix P = I - J' ( J J' )^-1 J
// the constraints are  norm( point(p+1) - point(p) ) - mtcut^2 = 0
// we drop the factor 2 when we derive, it does not change anything
void Fiber::computeProjectionMat()
{
  int nbc = nbSegments();             //number of constraints
  int nbv = DIM * nbPoints();         //number of variables

  allocateProjection();

  //----- we allocate the static arrays needed:
  static real * J = 0;
  static real * JJt = 0;
  static int allocated = 0;

  if ( nbc > allocated ) {
    if ( J )   delete[] J;
    if ( JJt ) delete[] JJt;
    allocated = nbc;
    J   = new real[ nbv * nbc ];
    JJt = new real[  2  * nbc ];
    if (( J == 0 ) || ( JJt == 0 )) {
      fprintf(stderr, "Fiber::computeProjectionMat() memory allocation failed\n");
      exit(1);
    }
  }

 //------------compute the projection matrix
  real x;
  Vecteur v, w, dv, dw;
  int ofs, jj, kk;
  int info=0;

  blas_xzero( nbv*nbc, J, 1 );

  //set up the Jacobian matrix J and the diagonals of J * Jt	
  w  = whereP( 0 );
  for( jj = 0; jj < nbc ; ++jj )
  {
    //set J:
    for( ofs = 0; ofs < DIM ; ++ofs )
    {
      kk = DIM * jj + ofs;
      x = pspts[ kk + DIM ] - pspts[ kk ];
      J[ jj + nbc * kk         ] = -x;
      J[ jj + nbc * ( kk+DIM ) ] =  x;
    }

    //set the diagonal and off-diagonal term of JJt:
    v  = w;
    w  = whereP( jj + 1 );
    dv = dw;
    dw = w - v;
    JJt[ jj ] = 2 * dw.normsq();       //diagonal term
    JJt[ jj + nbc ] = - dv * dw;       //off-diagonal term
  }

  // JJtiJ <- J
  blas_xcopy( nbc * nbv, J, 1, mtJJtiJ, 1 );

  // JJtiJ <- inv( JJt ) * J
  lapack_xptsv(nbc, nbv, JJt, JJt+nbc+1, mtJJtiJ, nbc, &info );
  if ( info ) MSG.error("Fiber::computeProjectionMat()", "lapack_ptsv() failed");

  // matP <-  - Jt * JJtiJ
  blas_xgemm('T', 'N', nbv, nbv, nbc, -1., J, nbc,
	     mtJJtiJ, nbc, 0., mtP, nbv );

  // matP <- matP + I
  for( jj = 0; jj < nbv*nbv; jj += nbv+1 )
    mtP[ jj ] += 1.0;
  //printf("M %lx\n",name ); matPrint( nbv, nbv, mtP );
  //TODO: non-isotropic mobility corrections : MT move easier along axis
}




//------------------------------------------------------------------
void Fiber::prepareProjection()
{
  computeProjectionMat();
}


//------------------------------------------------------------------
void Fiber::setProjectedForces( const real * X, real * Y ) const
  // Y <- P * X = X + internal forces
{     
  int nbv = DIM * nbPoints();
  blas_xsymv('U', nbv, 1.0, mtP, nbv, X, 1, 0.0, Y, 1);
}



//------------------------------------------------------------------
//========================= PROJECTION DERIVATIVE ==================
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

void Fiber::prepareProjectionDiff(const real * force)
{
  int nbc = nbSegments();             //number of constraints
  int nbv = DIM * nbPoints();         //number of variables

  //----- we allocate the static arrays needed:
  static real * JJtiJforce = 0;
  static int allocated = 0;

  if ( nbc > allocated ) {
    if ( JJtiJforce ) delete[] JJtiJforce;
    allocated = nbc;
    JJtiJforce = new real[ nbc ];
  }

  //caculate the lagrangian coefficients:
  blas_xgemv('N',nbc,nbv, 1., mtJJtiJ, nbc, force, 1, 0., JJtiJforce, 1);
  
  blas_xcopy( nbc, JJtiJforce, 1, lagrangeMult, 1 );
  lagrangeValid = true;
  //printf("all  lagrangeMult "); vectPrint( nbc, lagrangeMult );
 
  //we remove compressive forces ( lagrangian < 0 )
  for( int ii = 0; ii < nbc; ++ii )
    if ( JJtiJforce[ ii ] < 0 )
      JJtiJforce[ ii ] = 0;

  //printf("L ");vectPrint( nbc, JJtiJforce );

  //set up the first term in the derivative of J with respect to variable x[ii]
  //set up term  P * (DJ)t (JJti) J force:
  for( int jj = 0; jj < nbv; ++jj ) {
    real * coljj = mtDiffP + nbv * jj;
    blas_xzero( nbv, coljj, 1 );
    int lin = jj / DIM;
    if ( lin > 0 ) {
      coljj[ jj-DIM ] = +JJtiJforce[ lin-1 ];
      coljj[ jj     ] = -JJtiJforce[ lin-1 ];
    }
    if ( lin < nbc ) {
      coljj[ jj     ] += -JJtiJforce[ lin ];
      coljj[ jj+DIM ]  = +JJtiJforce[ lin ];
    }
    //printf("projectionDiff M%lx\n", name); matPrint( nbv, nbv, mtDiffP );
    //the final matrix is symmetric, for any force, 
    //as can be seen from the above relations to set its columns
  }
}

//------------------------------------------------------------------
void Fiber::addProjectionDiff( const real * X, real * Y ) const
{
  int nbv = DIM * nbPoints();
  blas_xsymv('U', nbv, 1., mtDiffP, nbv, X, 1, 1., Y, 1);
}


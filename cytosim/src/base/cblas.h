//RCS: $Id: cblas.h,v 2.3 2005/04/24 22:46:50 nedelec Exp $
//---------------------------------cblas.h-----------------------------------
// This contains C front-ends to the standard fortran linear algebra package BLAS
// see http://www.netlib.org/blas
//
// functions are renamed : eg. blas_xcopy calls scopy if SINGLE is defined,
//                                           or dcopy if real == double.
// we added a function blas_xzero for convenient initialization of vectors,
// and blas_xnrm8(), to infinite norm

#ifndef BLAS_H 
#define BLAS_H

#include "types.h"
#include "smath.h"

extern "C" {

#undef F_NAME
#ifdef SINGLE
  #define F_NAME(x) s##x##_
  #define iF_NAME(x) is##x##_
#else
  #define F_NAME(x) d##x##_
  #define iF_NAME(x) id##x##_
#endif

  
  
/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

inline real blas_xdot( int N, const real * X, int incX, const real * Y, int incY) 
{
    real F_NAME(dot)( int *, const real*, int *, const real*, int *);
    return F_NAME(dot)(&N, X, &incX, Y, &incY);
}

inline real blas_xnrm2( int N, const real *X, int incX)
{
  real F_NAME(nrm2)( int *, const real *, int *);
  return F_NAME(nrm2)(&N, X, &incX);
}

inline real blas_xasum( int N, const real *X, int incX)
{
  real F_NAME(asum)( int *, const real *, int *);
  return F_NAME(asum)(&N, X, &incX);
}

inline real blas_xsum( int N, const real *X, int incX)
{
  real F_NAME(sum)( int *, const real *, int *);
  return F_NAME(sum)(&N, X, &incX);
}

inline int blas_ixamax( int N, const real *X, int incX)
{
  int iF_NAME(amax)( int *, const real *, int *);
  return ( iF_NAME(amax)(&N, X, &incX) - 1 );
}


//Our little addition: the infinite norm:

inline real blas_xnrm8( int N, const real * X, int incX ) {
  int indx = blas_ixamax( N, X, incX );
  return fabs( X[ indx ] ); 
}


inline int blas_ixmax( int N, const real *X, int incX) 
{
  int iF_NAME(max)( int *, const real *, int *);
  return (iF_NAME(max)(&N, X, &incX) - 1);
}

inline int blas_ixamin( int N, const real *X, int incX) 
{
  int iF_NAME(amin)( int *, const real *, int *);
  return (iF_NAME(amin)(&N, X, &incX) - 1);
}

inline int blas_ixmin( int N, const real *X, int incX) 
{
  int iF_NAME(min)( int *, const real *, int *);
  return (iF_NAME(min)(&N, X, &incX) - 1);
}

inline void blas_xswap( int N, real *X, int incX, real *Y, int incY) 
{
  void F_NAME(swap)( int *, real *, int *, real *, int *);
  F_NAME(swap)(&N, X, &incX, Y, &incY);
}

inline void blas_xcopy( int N, const real *X, int incX, real *Y, int incY) 
{
  void F_NAME(copy)( int *, const real *, int *, real *, int *);
  F_NAME(copy)(&N, X, &incX, Y, &incY);
}

///another addition: set the vector to zero
inline void blas_xzero( int N, real * X, const int incX )
{
  
  const real zero = 0.0;
  blas_xcopy( N, &zero, 0, X, incX );
  /*
   if ( incX == 1 )
     for( int ii = 0; ii < N; ++ii ) X[ ii ] = 0;
   else	
     for( int ii = 0; ii < N; ++ii ) X[ ii*incX ] = 0;
  */
}

inline void blas_xaxpy( int N, real alpha, const real *X, int incX, real *Y, int incY) 
{
  void F_NAME(axpy)( int *, real *, const real *, int *, real *, int *);
  F_NAME(axpy)(&N, &alpha, X, &incX, Y, &incY );
}

inline void blas_xaxpby( int N, real alpha, const real *X, int incX, real beta, real *Y, int incY)
{
  void F_NAME(axpby)( int *, real *, const real *, int *, real *, real *, int *);
  F_NAME(axpby)( &N, &alpha, X, &incX, &beta, Y, &incY );
}

inline void blas_xrotg(real *a, real *b, real *c, real *s) 
{
  void F_NAME(rotg)(real *, real *, real *, real *);
  F_NAME(rotg)(a, b, c, s);
}

inline void blas_xrotmg(const real *d1, const real *d2, const real *b1,	real b2, real *P) 
{
  void F_NAME(rotmg)(const real *, const real *, const real *, real *, real *);
  F_NAME(rotmg)(d1, d2, b1, &b2, P);
}

inline void blas_xrot( int N, real *X, int incX, real *Y, int incY, real c, real s) 
{
  void F_NAME(rot)( int *, real *, int *, real *, int *, real *, real *);
  F_NAME(rot)(&N, X, &incX, Y, &incY, &c, &s);
}

inline void blas_xrotm( int N, real *X, int incX,	real *Y, int incY, real *P) 
{
  void F_NAME(rotm)( int *, real *, int *, real *, int *, real *);
  F_NAME(rotm)(&N, X, &incX, Y, &incY, P);
}

inline void blas_xscal( int N, real alpha, real *X, int incX) 
{
  void F_NAME(scal)( int *, real *, real *, int *);
  F_NAME(scal)( &N, &alpha, X, &incX);
}

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

inline void blas_xgemv(char TransA, int M, int N,
		       real alpha, const real *A, int lda,
		       const real *X, int incX, real beta,
		       real *Y, int incY)
{
  void F_NAME(gemv)(char *, int *, int *, real *, 
		    const real *, int *, const real *, int *, 
		    real *, real *, int *);
  F_NAME(gemv)(&TransA, &M, &N, &alpha, A, &lda, X, &incX, &beta, Y, &incY);  
}

inline void blas_xtrmv( char Uplo, char TransA, char Diag, int N, 
			const real *A, int lda, real *X, int incX)
{
  void F_NAME(trmv)( char *, char *, char *, int *, 
		     const real *, int *, real *, int *);
  F_NAME(trmv)(&Uplo, &TransA, &Diag, &N, A, &lda, X, &incX);
  
}

inline void blas_xtrsv( char Uplo, char TransA, char Diag,
			int N, const real *A, int lda, real *X,
			int incX);

inline void blas_xgbmv(char TransA, int M, int N, int Kl, 
		       int Ku, real alpha, const real *A, 
		       int lda, const real *X, int incX, 
		       real beta, real *Y, int incY);

inline void blas_xtbmv( char Uplo, char TransA, char Diag,
			int N, int K, const real *A, int lda,
			real *X, int incX);

inline void blas_xtbsv( char Uplo, char TransA, char Diag,
			int N, int K, const real *A, int lda,
			real *X, int incX);

inline void blas_xtpsv( char Uplo, char TransA, char Diag,
			int N, const real *Ap, real *X, int incX);

inline void blas_xtpmv( char Uplo, char TransA, char Diag,
			int N, const real *Ap, real *X, int incX);

inline void blas_xger( int M, int N, real alpha, 
		       const real *X, int incX, const real *Y,
		       int incY, real *A, int lda)
{
  void F_NAME(ger)( int *, int *, real * alpha, 
		    const real *, int *, const real *,
		    int *, real *, int *);
  F_NAME(ger)(&M, &N, &alpha, X, &incX, Y, &incY, A, &lda);
}


inline void blas_xsymv(char Uplo, int N, real alpha, 
		       const real *A, int lda, const real *X, 
		       int incX, real beta, real *Y, 
		       int incY)
{
  void F_NAME(symv)( char *, int *, real *, const real *,
		     int *, const real *, int *,
		     real *, real *, int *);
  F_NAME(symv)(&Uplo,&N,&alpha,A,&lda,X,&incX,&beta,Y,&incY);
}

inline void blas_xsbmv( char Uplo, int N, int K, real alpha,
			const real *A, int lda, const real *X, 
			int incX, real beta, real *Y, int incY)
{
  void F_NAME(sbmv)( char *, int *, int *, real *,
		     const real *, int *, const real *, 
		     int *, real *, real *, int *);
  F_NAME(sbmv)(&Uplo,&N,&K,&alpha,A,&lda,X,&incX,&beta,Y,&incY);
}

inline void blas_xspmv( char Uplo, int N, real alpha, 
			const real *A, const real *X, int incX, 
			real beta, real *Y, int incY)
{
  void F_NAME(spmv)( char *, int *, real *, 
		     const real *, const real *, int *, 
		     real *, real *, int *);
  F_NAME(spmv)(&Uplo,&N,&alpha,A,X,&incX,&beta,Y,&incY);
}

inline void blas_xsyr( char Uplo, int N, real alpha, const real *X, 
		       int incX, real *A, int lda)
{
  void F_NAME(syr)( char *, int *, real *, 
		    const real *, int *,
		    real *, int *);
  F_NAME(syr)( &Uplo, &N, &alpha, X, &incX, A, &lda );
}

inline void blas_xsyr2( char Uplo, int N, real alpha, const real *X, 
			int incX, const real *Y, int incY, real *A, 
			int lda);

inline void blas_xspr( char Uplo, int N, real alpha, const real *X, 
		       int incX, real *Ap)
{
  void F_NAME(spr)( char *, int *, real *,
		    const real *, int *, real *);
  F_NAME(spr)( &Uplo, &N, &alpha, X, &incX, Ap);
}

inline void blas_xspr2( char Uplo, int N, real alpha, const real *X, 
			int incX, const real *Y, int incY, real *A);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

inline void blas_xgemm( char TransA, char TransB, int M, int N,
			int K, real alpha, const real *A,
			int lda, const real *B, int ldb,
			real beta, real *C, int ldc)
{
  void F_NAME(gemm)( char *, char *, int *, int *,
		     int *, real *, const real * ,
		     int *, const real *, int *,
		     real *, real *, int *);
  F_NAME(gemm)(&TransA, &TransB,&M,&N,&K,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

inline void blas_xsymm( char Side, char Uplo, int M, int N,
			real alpha, const real *A, int lda,
			const real *B, int ldb, real beta,
			real *C, int ldc)
{
  void F_NAME(symm)( char *, char *, int *, int *,
		     real *, const real *, int *,
		     const real *, int *, real *,
		     real *, int *);
  F_NAME(symm)(&Side, &Uplo,&M,&N,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}

inline void blas_xsyrk( char Uplo, char Trans, int N, int K,
			real alpha, const real *A, int lda,
			real beta, real *C, int ldc);

inline void blas_xsyr2k( char Uplo, char Trans, int N, int K,
			 real alpha, const real *A, int lda,
			 const real *B, int ldb, real beta,
			 real *C, int ldc);

inline void blas_xtrmm( char Uplo, char TransA, char Diag, 
			int M, int N,
			real alpha, const real *A, int lda,
			real *B, int ldb);

inline void blas_xtrsm( char side, char uplo, char transA, char diag, 
			int M, int N,
			real alpha, const real *A, int lda,
			real *B, int ldb)
{
  void F_NAME(trsm)( char *, char *, char *, char *,
		     int *, int *, real *, const real *, int *,
		     real *, int *);
  F_NAME(trsm)(&side, &uplo, &transA, &diag, &M, &N, &alpha, A, &lda,
	       B, &ldb);
}
}


#endif /* BLAS_H */

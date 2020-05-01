//RCS: $Id: clapack.h,v 2.0 2004/08/16 16:35:01 nedelec Exp $
// This contains C-front ends to some of the Fortran routines of LAPACK.
// see http://www.netlib.org/lapack
//
// LAPACK contains more than 1000 linear algebra functions in FORTRAN, 
// here we just implement the functions needed in Cytosim.

#ifndef LAPACK_H
#define LAPACK_H


#include "cblas.h"

extern "C" {
  
  inline void lapack_xptsv( int N, int NHRS, real* D, real* E, 
                            real* B, int LDB, int* INFO)
  {
    void F_NAME(ptsv)( int*, int*, real*, real*, real*, int*, int* );
    F_NAME(ptsv)(&N, &NHRS, D, E, B, &LDB, INFO );
  }
  
  
  inline void lapack_xptsvx( char fact, int N, int NHRS, 
                             const real* D, const real* E, real * DF, real * EF,
                             const real* B, int LDB, real * X, int LDX,
                             real * RCOND, real * FERR, real * BERR,
                             real * work, int * INFO)
  {
    void F_NAME(ptsvx)( char *, int *, int *, 
                        const real *, const real *, real *, real *,
                        const real *, int *, real *, int *,
                        real *, real *, real *,
                        real *, int * );
    F_NAME(ptsvx)(&fact, &N, &NHRS, D, E, DF, EF, B, &LDB, X, &LDX, 
                  RCOND, FERR, BERR, work, INFO );
  }
  
  inline void lapack_xpttrf( int N, real * D, real * E, int * INFO)
  {
    void F_NAME(pttrf)( int *, real *, real *, int *);
    F_NAME(pttrf)( &N, D, E, INFO );
  }
  
  inline void lapack_xpttrs( int N, int NHRS, 
                             const real * D, const real * E, real * B, 
                             int LDB, int * INFO )
  {
    void F_NAME(pttrs)( int *, int *,
                        const real *, const real *, real *, int *, int *);
    F_NAME(pttrs)( &N, &NHRS, D, E, B, &LDB, INFO );
  }
  
  
  
  inline void lapack_xposv( char uplo, int N, int NHRS, 
                            real* A, int LDA,
                            real* B, int LDB, int * INFO)
  {
    void F_NAME(posv)( char *, int*, int*, real*, int*, real*, int*, int*);
    F_NAME(posv)(&uplo, &N, &NHRS, A, &LDA, B, &LDB, INFO );
  }
  
  inline void lapack_xpotrf( char uplo, int N,
                             real* A, int LDA, int * INFO)
  {
    void F_NAME(potrf)( char *, int*, real*, int*, int*);
    F_NAME(potrf)(&uplo, &N, A, &LDA, INFO );
  }
  
  inline void lapack_xpotrs( char uplo, int N, int NHRS, 
                             const real* A, int LDA,
                             real* B, int LDB, int * INFO)
  {
    void F_NAME(potrs)( char *, int*, int*,
                        const real*, int*, real*, int*, int*);
    F_NAME(potrs)(&uplo, &N, &NHRS, A, &LDA, B, &LDB, INFO );
  }
  
  inline void lapack_xpotf2( char uplo, int N,
                             real* A, int LDA, int * INFO)
  {
    void F_NAME(potf2)( char *, int*, real*, int*, int*);
    F_NAME(potf2)(&uplo, &N, A, &LDA, INFO );
  }
  
  inline void lapack_xpotri( char uplo, int N,
                             real* A, int LDA, int * INFO)
  {
    void F_NAME(potri)( char *, int*, real*, int*, int*);
    F_NAME(potri)(&uplo, &N, A, &LDA, INFO );
  }
  
  inline void lapack_xpptrf( char uplo, int N,
                             real* A, int * INFO)
  {
    void F_NAME(pptrf)( char *, int*, real*, int*);
    F_NAME(pptrf)(&uplo, &N, A, INFO );
  }
  
  inline void lapack_xpptrs( char uplo, int N, int NHRS, 
                             const real* A,
                             real* B, int LDB, int * INFO)
  {
    void F_NAME(pptrs)( char *, int*, int*,
                        const real*, real*, int*, int*);
    F_NAME(pptrs)(&uplo, &N, &NHRS, A, B, &LDB, INFO );
  }
  
  inline void lapack_xpptri( char uplo, int N,
                             real* A, int * INFO)
  {
    void F_NAME(pptri)( char *, int*, real*, int*);
    F_NAME(pptri)(&uplo, &N, A, INFO );
  }
  
  inline void lapack_xtrtrs( char uplo, char trans, char diag, int N, int NHRS,
                             const real* A, int LDA,
                             real* B, int LDB, int * INFO)
  {
    void F_NAME(trtrs)( char *, char *, char *, int*, int*,
                        const real*, int *, real*, int*, int*);
    F_NAME(trtrs)(&uplo, &trans, &diag, &N, &NHRS, A, &LDA, B, &LDB, INFO );
  }
  
  inline void lapack_xgesv( int N, int NRHS, 
                            real * A, int LDA, int * IPIV, 
                            real * B, int LDB, int * INFO)
  {
    void F_NAME(gesv)( int *, int *, real *, int *, int *, 
                       real *, int *, int *);
    F_NAME(gesv)(&N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
  }
  
  
  inline void lapack_xgetrf( int M, int N, real * A, int LDA, 
                             int * IPIV, int * INFO)
  {
    void F_NAME(getrf)( int *, int *, real *, int *, int *, int *);
    F_NAME(getrf)(&M, &N, A, &LDA, IPIV, INFO);
  }
  
  inline void lapack_xgetri( int N, real * A, int LDA, 
                             const int * IPIV, real * WORK, int LWORK,
                             int * INFO)
  {
    void F_NAME(getri)( int *, real *, int *,
                        const int *, real *, const int *, int *);
    F_NAME(getri)(&N, A, &LDA, IPIV, WORK, &LWORK, INFO);
  }
  
  inline void lapack_xgetrs( char trans, int N, int NRHS, 
                             const real * A, int LDA, 
                             const int * IPIV, real * B, int LDB, 
                             int * INFO)
  {
    void F_NAME(getrs)( char *, int *, int *,
                        const real *, int *, 
                        const int *, real *, int *, int *);
    F_NAME(getrs)(&trans, &N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
  }
  
  
  inline void lapack_xsyevx( char jobz, char range, char uplo, int N,
                             real * A, int LDA, real VL, real VU,
                             int IL, int IU, real ABSTOL, int * M, 
                             real * W, real * Z, int LDZ, real * WORK, 
                             int LWORK, const int * IWORK,
                             int * IFAIL, int * INFO)
  {
    
    void F_NAME(syevx)( char *, char *, char *, int *,
                        real *, int *, real *, real *,
                        int *, int *, real *, int *, 
                        real *, real *, int *, real * ,
                        int *, const int *, int *, int *);
    
    F_NAME(syevx)(&jobz, &range, &uplo, &N,
                  A, &LDA, &VL, &VU, &IL, &IU,
                  &ABSTOL, M, W, Z,
                  &LDZ, WORK, &LWORK,
                  IWORK, IFAIL, INFO);
  }
  
  
  
  
  
  inline void F_NAME(gtsv)( int*, int*, real*, real*, real*,
                            real*, int*, int*);
  
  inline void F_NAME(ysv)( char*, int*, int*, real*, int*, int*, 
                           real*, int*, real*, int*, int*);



  
  inline void lapack_xgeqrf( int M, int N, real * A, int LDA, 
                             real * TAU, real * WORK, int LWORK,
                             int * INFO)
    {
    void F_NAME(geqrf)( int *, int *, real *, int *, 
                        real *, real *, const int *, int *);
    F_NAME(geqrf)(&M, &N, A, &LDA, TAU, WORK, &LWORK, INFO);
    }
  
  inline void lapack_xormqr( char side, char trans, int M, int N, int K,
                             const real * A, int LDA, const real * TAU,
                             real * C, int LDC,
                             real * WORK, int LWORK, int * INFO)
    {
    void F_NAME(ormqr)( char *, char *, int *, int *, int *,
                        const real *, int *, const real *,
                        real *, int *, real *, const int *, int *);
    F_NAME(ormqr)(&side, &trans, &M, &N, &K, A, &LDA, TAU, C, &LDC,
                  WORK, &LWORK, INFO);
    }
  
  inline void lapack_xgels( char trans, int M, int N, int NRHS,
                            const real * A, int LDA, const real * B, int LDB,
                            real * WORK, int LWORK, int * INFO)
    {
    void F_NAME(gels)( char *, int *, int *, int *,
                       const real *, int *, const real *, int *,
                       real *, const int *, int *);
    F_NAME(gels)(&trans, &M, &N, &NRHS, A, &LDA, B, &LDB,
                 WORK, &LWORK, INFO);
    }
  
}  //extern("C")



#endif //LAPACK_H

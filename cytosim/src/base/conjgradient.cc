//RCS: $Id: conjgradient.cc,v 2.7 2005/03/14 18:12:17 nedelec Exp $
//-----------------------------------------------------------------------------
//               conjugate gradient and related iterative methods
//           to solve linear systems: http://www.netlib.org/templates
//-----------------------------------------------------------------------------

#include "conjgradient.h"
#include "cblas.h"
#include "clapack.h"
#include "smath.h"

//The linear system is only defined by the functions, which calculate
//the product of the matrix by a given vector. These functions are defined
//in simsolve.cc and passed to the solvers as an argument.
//The preconditionner P is defined in the same way.

#ifndef NORM
   #define NORM blas_xnrm8
#endif

//------------------ the vector-vector dot product used in the algorithm

#ifndef SINGLE
  #define DOT blas_xdot
#else
//we sum in double precision to improve accuracy...
double DOT( int size, real * X, int inc_x, real * Y, int inc_y )
{
  double result = 0;
  if ( ( inc_x == 1 ) && ( inc_y == 1 ) ) {
    
    for( int ii = 0; ii < size; ++ii )
      result += double( X[ ii ] ) * double( Y[ ii ] );
    
  } else {
    
    for( int ii = 0; ii < size; ++ii )
      result += double( X[ inc_x * ii ] ) * double( Y[ inc_y * ii ] );
  }
  return result;
}
#endif



//-------all purpose allocation function:

void allocateCG(int size, int& allocated,
                real** vec1=0, real** vec2=0, real** vec3=0, real** vec4=0,
                real** vec5=0, real** vec6=0, real** vec7=0, real** vec8=0)
{
    // allocate the working-arrays for the solvers 
    
    if ( size > allocated ) {
      if ( vec1 && *vec1 ) delete[]( *vec1 ); 
      if ( vec2 && *vec2 ) delete[]( *vec2 );
      if ( vec3 && *vec3 ) delete[]( *vec3 );
      if ( vec4 && *vec4 ) delete[]( *vec4 );
      if ( vec5 && *vec5 ) delete[]( *vec5 );
      if ( vec6 && *vec6 ) delete[]( *vec6 );
      if ( vec7 && *vec7 ) delete[]( *vec7 );
      if ( vec8 && *vec8 ) delete[]( *vec8 );
      
      if ( vec1 ) { 
        *vec1 = new real[ size ];
        if (*vec1 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }
      
      if ( vec2 ) { 
        *vec2 = new real[ size ];
        if (*vec2 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }

      if ( vec3 ) { 
        *vec3 = new real[ size ];
        if (*vec3 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }

      if ( vec4 ) { 
        *vec4 = new real[ size ];
        if (*vec4 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }

      if ( vec5 ) { 
        *vec5 = new real[ size ];
        if (*vec5 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }

      if ( vec6 ) { 
        *vec6 = new real[ size ];
        if (*vec6 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }

      if ( vec7 ) { 
        *vec7 = new real[ size ];
        if (*vec7 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }
      
      if ( vec8 ) { 
        *vec8 = new real[ size ];
        if (*vec8 == 0) {
          fprintf(stderr, "ConjugateGradient::allocate() memory allocation failed\n");
          exit(1);
        }
      }
      
      allocated = size;
    }
}



//=============================================================================
//              Conjugate Gradient, no Preconditionning
//=============================================================================


void Solver::conjGradient(int size, const real* b, real* x,
                          void (*matVect)( const real*, real* ),
                          int & max_iter, real & tol)
{
  static real *d = 0, *s=0, *r=0, *p=0, *q=0;
  static int alc = 0;
  
  allocateCG( size, alc, &d, &s, &r, &p, &q );
  
  real alpha, beta, dold, dnew;
  
  blas_xcopy( size, b, 1, r, 1 );
  matVect( x, s );
  blas_xaxpy( size, -1, s, 1, r, 1);            //   r <- b - A * x
  blas_xcopy( size, r, 1, d, 1 );               //   d <- r 
  dnew = DOT(size, r, 1, r, 1);
  
  real normb = NORM( size, b, 1 );
  real limit = tol * normb;
  real test = NORM( size, r, 1 );
  
  int ii;
  for (ii=0; (ii <= max_iter) && ( test > limit ); ++ii ) 
    {
      matVect( d, q );                          //   q = A * d
      
      alpha = dnew / DOT(size, d, 1, q, 1);
      blas_xaxpy( size,  alpha, d, 1, x, 1 );   //   x += alpha * d
      blas_xaxpy( size, -alpha, q, 1, r, 1 );   //   r -= alpha * q
      
      dold = dnew;
      dnew = DOT(size, r, 1, r, 1);
      beta = dnew / dold;
      blas_xscal( size, beta, d, 1 );
      blas_xaxpy( size, 1, r, 1, d, 1 );        //   d = beta * d + r
      
      test = NORM( size, r, 1 );
    }
  
  tol = test / normb;
  max_iter = ii;
}



//=============================================================================
//              Conjugate Gradient, with Preconditioning
//=============================================================================



void Solver::conjGradientPrecond(int size, const real * b, real * x, 
                                 void (*matVect)( const real*, real* ),
                                 void (*precond)( const real*, real* ),
                                 int & max_iter, real & tol)
{
  static real *d = 0, *s=0, *r=0, *p=0, *q=0; 
  static int alc = 0;

  allocateCG( size, alc, &d, &s, &r, &p, &q );
  real alpha, beta, dold, dnew;
  
  blas_xcopy( size, b, 1, r, 1 );
  matVect( x, s );
  blas_xaxpy( size, -1, s, 1, r, 1);             //   r = b - M * x
  
  precond( r, d );                               //   d <- inv(M) * r
  
  dnew = DOT(size, r, 1, d, 1);
  
  real normb = NORM( size, b, 1 );
  real limit = tol * normb;                      //  limit <-  | b |
  real test = NORM(size, r, 1 );
  
  int ii;
  for (ii=0; (ii <= max_iter) && ( test > limit ); ++ii ) 
    {
      matVect( d, q );                           //   q = M * d
      
      alpha = dnew / DOT(size, d, 1, q, 1);
      blas_xaxpy( size,  alpha, d, 1, x, 1 );    //   x += alpha * d
      blas_xaxpy( size, -alpha, q, 1, r, 1 );    //   r -= alpha * q
            
      precond( r, s );                           //   s = inv(M) * r;
      
      dold = dnew;
      dnew = DOT(size, r, 1, s, 1);
      beta = dnew / dold;
      blas_xscal( size, beta, d, 1 );
      blas_xaxpy( size, 1, s, 1, d, 1 );         //   d = beta * d + s

      test = NORM( size, r, 1 );
    }
  tol = test / normb;
  max_iter = ii;
}



//=============================================================================
//                      Bi-Conjugate Gradient
//=============================================================================



void Solver::biConjGradient(int size, const real * b, real * x, 
                            void (*matVect)( const real*, real* ),
                            void (*matVectTrans)( const real*, real* ),
                            int & max_iter, real & tol)
{
  static real *r=0, *rb=0, *p=0, *pb=0, *q=0, *qb=0;
  static int alc = 0;
  
  allocateCG( size, alc, &r, &rb, &p, &pb, &q, &qb );
  
  real alpha, beta, dold, dnew;
  
  blas_xcopy( size, b, 1, r, 1 );
  matVect( x, rb );
  blas_xaxpy( size, -1, rb, 1, r, 1);            //   r = b - A * x
  
  blas_xcopy( size, r, 1, p, 1 );
  blas_xcopy( size, r, 1, rb, 1 );
  blas_xcopy( size, r, 1, pb, 1 );
  
  dnew = DOT(size, rb, 1, r, 1);
  
  real normb = NORM( size, b, 1 );
  real limit = tol * normb;                      //  limit <-  | b |
  real test = NORM( size, r, 1 );
  
  int ii;
  for (ii=0; (ii <= max_iter) && ( test > limit ); ++ii ) 
    {
      matVect( p, q );                           //   q = A * p
      matVectTrans( pb, qb );                    //   qb = A' * pb
      
      alpha = dnew / DOT(size, pb, 1, q, 1);
      blas_xaxpy( size,  alpha, p, 1, x, 1 );    //   x  += alpha * p
      blas_xaxpy( size, -alpha, q, 1, r, 1 );    //   r  -= alpha * q
      blas_xaxpy( size, -alpha, qb, 1, rb, 1 );  //   rb -= alpha * qb
      
      dold = dnew;
      dnew = DOT(size, r, 1, rb, 1);
      beta = dnew / dold;
      blas_xscal( size, beta, p, 1 );
      blas_xaxpy( size, 1, r, 1, p, 1 );         //   p  = beta * p  + r
      blas_xscal( size, beta, pb, 1 );
      blas_xaxpy( size, 1, rb, 1, pb, 1 );       //   pb = beta * pb + rb
      
      test = NORM( size, r, 1 );
    }
  tol = test / normb;
  max_iter = ii;
}



//=============================================================================
//                 Bi-Conjugate Gradient Stabilized
//=============================================================================



int Solver::biCGstab(int size, const real * b, real * x,
                     void (*matVect)( const real*, real* ),
                     int & iter, real & tol, real normb)
{
  double rho_1, rho_2 = 1, alpha = 0, beta, omega = 1.0;
  static real *r=0, *rtilde=0, *p=0, *t=0,  *v=0;
  static int alc = 0;
  
  allocateCG( size, alc, &r, &rtilde, &p, &t, &v );
  
  blas_xcopy( size, b, 1, r, 1 );
  matVect( x, rtilde );
  blas_xaxpy( size, -1.0, rtilde, 1, r, 1);       // r = b - A * x
  blas_xcopy( size, r, 1, rtilde, 1 );
  
  if ( normb == 0.0 ) normb = NORM( size, b, 1 );
  if ( normb == 0.0 ) normb = 1.0;
  
  real   max_tol  = tol;
  int    max_iter = iter;
  
  for(iter = 1; iter <= max_iter; ++iter) {
    rho_1 = DOT(size, rtilde, 1, r, 1);
    
    if (rho_1 == 0.0) {
      tol = NORM(size, r, 1) / normb;
      return 2;
    }
    
    if (iter == 1)
      blas_xcopy(size, r, 1, p, 1 );          // p = r;
    else {
      beta = (rho_1/rho_2) * (alpha/omega);
      blas_xaxpy(size, -omega, v, 1, p, 1);
      blas_xscal(size, beta, p, 1);
      blas_xaxpy(size, 1.0, r, 1, p, 1);      // p = r + beta*(p-omega*v);
    }
    
    matVect( p, v );                          // v = A * p;
    alpha = rho_1 / DOT(size, rtilde, 1,  v, 1);
    blas_xaxpy(size, -alpha, v, 1, r, 1);     // r = r - alpha * v;
    blas_xaxpy(size,  alpha, p, 1, x, 1);     // x = x + alpha * p;
    
    tol = NORM(size, r, 1) / normb;
    if ( tol < max_tol ) return 0;
    
    matVect( r, t );                          // t = A * s;
    
    omega = DOT(size, t, 1, r, 1) / DOT(size, t, 1, t, 1);
    blas_xaxpy(size,  omega, r, 1, x, 1);     // x = x + omega * r;
    blas_xaxpy(size, -omega, t, 1, r, 1);     // r = r - omega * t;
    
    tol = NORM(size, r, 1) / normb;
    if ( tol < max_tol ) return 0;
    if ( omega == 0.0 )  return 3;
    rho_2 = rho_1;
    }
  return 1;
}



//=============================================================================
//        Bi-Conjugate Gradient Stabilized with Preconditionning
//=============================================================================


int Solver::biCGstabPrecond(int size, const real * b, real * x, 
                            void (*matVect)( const real*, real* ),
                            void (*precond)( const real*, real* ),
                            int & iter, real & tol, real normb)
{
  double rho_1, rho_2 = 1, alpha = 0, beta, omega = 1.0, delta;
  static real *r=0, *rtilde=0, *p=0, *phat=0, *shat=0, *t=0, *v=0;
  static int alc = 0;

  allocateCG( size, alc, &r, &rtilde, &p, &t, &v, &phat, &shat );

  blas_xcopy( size, b, 1, r, 1 );
  matVect( x, rtilde );
  blas_xaxpy( size, -1.0, rtilde, 1, r, 1);         //   r = b - A * x
  
  blas_xcopy( size, r, 1, rtilde, 1 );              //   r_tilde = r

  if ( normb == 0.0 ) normb = NORM( size, b, 1 );
  if ( normb == 0.0 ) normb = 1.0;

  real   max_tol  = tol;
  int    max_iter = iter;

  for (iter = 1; iter <= max_iter; ++iter) {
    rho_1 = DOT(size, rtilde, 1, r, 1);
    
    //we should rather test if fabs(rho_1) is greater than an EPSILON
    //since we will divide by rho_1 at the next time step (rho_2 <- rho_1)
    if (rho_1 == 0.0) {
      tol = NORM(size, r, 1) / normb;
      return 2;
    }
    
    if (iter == 1)
      blas_xcopy(size, r, 1, p, 1 );            // p = r;
    else {
      beta = (rho_1/rho_2) * (alpha/omega);
      //we should test here the value of beta, which is scalar
      blas_xaxpy(size, -omega, v, 1, p, 1);
      blas_xscal(size, beta, p, 1);
      blas_xaxpy(size, 1.0, r, 1, p, 1);        // p = r + beta*(p-omega*v);
    }
    
    precond( p, phat );                         // phat = inv(M) * p;
    matVect( phat, v );                         // v = M * phat;
    
    //added test for failure detected by D. Foethke, Jan 2005
    delta = DOT(size, rtilde, 1,  v, 1);
    if (delta == 0.0) {
      tol = NORM(size, r, 1) / normb; 
      return 4;
    }
    alpha = rho_1 / delta;
    blas_xaxpy(size, -alpha, v, 1, r, 1);       // r = r - alpha * v;
    blas_xaxpy(size,  alpha, phat, 1, x, 1);    // x = x + alpha * phat;
    
    tol = NORM(size, r, 1) / normb;
    if ( tol < max_tol ) return 0;
    
    precond( r, shat );                         // shat = inv(M) * r
    matVect( shat, t );                         // t = M * shat
    
    omega = DOT(size, t, 1, r, 1) / DOT(size, t, 1, t, 1);
    blas_xaxpy(size,  omega, shat, 1, x, 1);    // x = x + omega * shat
    blas_xaxpy(size, -omega, t, 1, r, 1);       // r = r - omega * t
    
    tol = NORM(size, r, 1) / normb;
    if ( tol < max_tol ) return 0;
    if ( omega == 0.0 )  return 3;
    rho_2 = rho_1;
  }
  return 1;
}


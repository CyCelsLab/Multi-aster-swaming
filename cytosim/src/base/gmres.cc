//RCS: $Id: gmres.cc,v 2.3 2005/03/14 18:12:25 nedelec Exp $
//-----------------------------------------------------------------------------
//               GMRES and related iterative methods
//                     to solve linear systems
//-----------------------------------------------------------------------------

#include <cstdio>
#include "gmres.h"
#include "cblas.h"
#include "clapack.h"
#include "smath.h"

//The linear system is only defined by the functions, which calculate
//the product of the matrix by a given vector. These functions are defined
//in simsolve.cc and passed to the solvers as an argument.
//The preconditionner P is defined in the same way.


//the NORM used to measure convergence can be modified
#ifndef NORM
#define NORM blas_xnrm8
#endif

//-------all purpose allocation function:

void allocateSolver(int size, int& allocated,
                    real** vec1=0, real** vec2=0, real** vec3=0, real** vec4=0,
                    real** vec5=0, real** vec6=0, real** vec7=0, real** vec8=0,
                    real** mat1=0, real** mat2=0)
{
    // allocate the working-arrays for the solvers 
    
    if ( size > allocated ) {
        if ( vec1 && *vec1 ) delete[] *vec1;
        if ( vec2 && *vec2 ) delete[] *vec2;
        if ( vec3 && *vec3 ) delete[] *vec3;
        if ( vec4 && *vec4 ) delete[] *vec4;
        if ( vec5 && *vec5 ) delete[] *vec5;
        if ( vec6 && *vec6 ) delete[] *vec6;
        if ( vec7 && *vec7 ) delete[] *vec7;
        if ( vec8 && *vec8 ) delete[] *vec8;
        if ( mat1 && *mat1 ) delete[] *mat1;
        if ( mat2 && *mat2 ) delete[] *mat2;
        
        if ( vec1 ) *vec1 = new real[ size + 1];
        if ( vec2 ) *vec2 = new real[ size + 1];
        if ( vec3 ) *vec3 = new real[ size + 1];
        if ( vec4 ) *vec4 = new real[ size ];
        if ( vec5 ) *vec5 = new real[ size ];
        if ( vec6 ) *vec6 = new real[ size ];
        if ( vec7 ) *vec7 = new real[ size ];
        if ( vec8 ) *vec8 = new real[ size ];
        if ( mat1 ) *mat1 = new real[ (size+1)*(size+1) ];
        if ( mat2 ) *mat2 = new real[ (size+1)*(size+1) ];
        
        allocated = size;
    }
}


void Solver::FOM(int size, const real* b, real* dx,
		 void (*matVec)( const real*, real* ),
		 int & max_iter)
{
  // this solver uses pure FOM (Full Orthogonalization Method)

  int mBroken = max_iter;
  real beta;
  
  int info;

  static real *V=0, *H=0, *Ym=0;
  static int alc = 0;

  allocateSolver(size, alc, &Ym, 0, 0, 0, 0, 0, 0, 0, &V, &H);
  
  // calculate r0
  // (we are starting with x0 = 0, so r0 = b and there's nothing to do)
  
  // calculate beta and v1
  beta = blas_xnrm2(size,b,1);
  blas_xcopy(size,b,1,V,1); 
  blas_xscal(size,1./beta,V,1);

  //
  for(int i = 0; i < max_iter*max_iter; i++) {
    H[i] = 0.;
  }
  
  // do the loop

  //bool secondOrtho = 0;
  //real normAv;

  for(int j = 0; j < max_iter; j++) {
    //if(!secondOrtho) {
      matVec(V+j*size, V+(j+1)*size);
      //  normAv = blas_xnrm2(size, V+(j+1)*size, 1);
      //}
    //vectPrint(size, V+(j+1)*size, 4);
    //printf("normAv %d: %e\n", j+1, blas_xnrm2(size,V+(j+1)*size,1));
    for(int i = 0; i <= j; i ++) {
      H[i + j*max_iter] = blas_xdot(size,V+(j+1)*size,1,V+i*size,1);
      blas_xaxpy(size,-1.*H[i + j*max_iter],V+(i*size),1,V+(j+1)*size,1);
      //vectPrint(size, V+(j+1)*size, 4);
      //printf("norm %d: %e\n", j+1, blas_xnrm2(size,V+(j+1)*size,1));      
    }
    H[(j+1) + j*max_iter] = blas_xnrm2(size,V+(j+1)*size,1);
    //printf("Final Norm: H[%d + %d*(max_iter+1)]: %e\n", j+1, j, H[(j+1) + j*max_iter]);

    // check for breakdown
    if(H[(j+1) + j*max_iter] < 1.e-10) {
      printf("breakdown after %d of %d steps!\n", j, max_iter);
      mBroken = j+1;
      break;
    }
//     else {
      
//       printf("Cancellation grade: %e\n", H[(j+1) + j*max_iter]/normAv);
//       if(H[(j+1) + j*max_iter]/normAv < 1.e-8) {
// 	printf("Severe cancellation occured during orthogonalization!\n");
// 	if(secondOrtho) {
// 	  printf("Second orthogonalization failed as well - giving up!\n");
// 	  exit(0);
// 	  secondOrtho = 0;
// 	  blas_xscal(size,1./H[(j+1) + j*max_iter],V+(j+1)*size,1);      
// 	}
// 	else {
// 	  secondOrtho = 1;
// 	  blas_xscal(size,1./H[(j+1) + j*max_iter],V+(j+1)*size,1);      
// 	  j -= 1;
// 	}
//       }
//       else {
// 	secondOrtho = 0;
// 	blas_xscal(size,1./H[(j+1) + j*max_iter],V+(j+1)*size,1);      
//       }
//     }
    
//     printf("norm %d: %f\n", j+1, blas_xnrm2(size,V+(j+1)*size,1));
//     printf("\n");
    blas_xscal(size,1./H[(j+1) + j*max_iter],V+(j+1)*size,1);      
  }


//   printf("FOM with %d steps (m=%d):\n", mBroken, max_iter);
//   //  matPrint(size, mBroken, V, 5);
//   printf("\n");
//   for(int j = 0; j < mBroken; j++) {
//     for(int i = 0; i < j; i++) {
//       printf("%d:%d %e ", i+1, j+1, blas_xdot(size, V+i*size, 1, V+j*size, 1));
//     }
//     printf("\n");
//   }

  // calculate ym
  for(int i = 0; i < mBroken; i++) {
    Ym[i] = 0.;
  }
  Ym[0] = beta;
  int* ipiv = new int[mBroken];
  lapack_xgetrf(mBroken,mBroken,H,max_iter,ipiv,&info);
  lapack_xgetrs('n',mBroken,1,H,max_iter,ipiv,Ym,max_iter,&info);
  if(info) {
    printf("Lapack-Error: Could not solve system! Info Code: %d\n", info);
    matPrint(max_iter, mBroken, H);
    exit(0);
  }
  delete[] ipiv;

  // calculate xm = x0 + Vm*ym (x0 = 0)
  blas_xgemv('n',size,mBroken,1.,V,size,Ym,1,0.,dx,1);

//   real *test = new real[size];
//   matVec(dx, test);
//   vectPrint(size, b, 5);
//   vectPrint(size, test, 5);
//   blas_xaxpy(size, -1., b, 1, test, 1);
//   vectDump(size, test);
//   printf("\n");
//   delete[] test;

}

void Solver::GMRES(int size, const real *b, real *dx,
		   void (*matVec)( const real*, real* ),
		   int &max_iter)
{
  // this solver uses pure GMRES (Generalized Minimum Residual)

  int mBroken = max_iter;
  real beta;

  int info;
  real *tau  = new real[size];
  int lwork  = size*16;
  real *work = new real[lwork];
  
  static real *V=0, *H=0, *Ym=0;
  static int alc = 0;

  allocateSolver(size, alc, &Ym, 0, 0, 0, 0, 0, 0, 0, &V, &H);

  // calculate r0
  // (we are starting with x0 = 0, so r0 = b and there's nothing to do)
  
  // calculate beta and v1
  beta = blas_xnrm2(size,b,1);
  blas_xcopy(size,b,1,V,1); 
  blas_xscal(size,1./beta,V,1);

  //
  for(int i = 0; i < (max_iter+1)*max_iter; i++) {
    H[i] = 0.;
  }
  
  // do the loop
  for(int j = 0; j < max_iter; j++) {
    matVec(V+j*size, V+(j+1)*size);
    for(int i = 0; i <= j; i ++) {
      H[i + j*(max_iter+1)] = blas_xdot(size,V+(j+1)*size,1,V+i*size,1);
      blas_xaxpy(size,-H[i + j*(max_iter+1)],V+i*size,1,V+(j+1)*size,1);
    }
    H[(j+1) + j*(max_iter+1)] = blas_xnrm2(size,V+(j+1)*size,1);
    
    // check for breakdown
    if(H[(j+1) + j*(max_iter+1)] < 1.e-15) {
      printf("breakdown after %d of %d steps!\n", j, max_iter);
      mBroken = j+1;
      break;
    }

    blas_xscal(size,1./H[(j+1) + j*(max_iter+1)],V+(j+1)*size,1);
  }
  
  // calculate ym
  for(int i = 0; i < mBroken+1; i++) {
    Ym[i] = 0.;
  }
  Ym[0] = beta;
  lapack_xgeqrf(mBroken+1,mBroken,H,max_iter+1,tau,work,lwork,&info);
  if(info) {
    printf("Lapack-Error in geqrf: Info Code: %d\n", info);
    exit(0);
  }
  lapack_xormqr('l','t',mBroken+1,1,mBroken,H,max_iter+1,tau,Ym,max_iter+1,work,lwork,&info);
  if(info) {
    printf("Lapack-Error in ormqr: Info Code: %d\n", info);
    exit(0);
  }
  blas_xtrsm('l','u','n','n',mBroken,1,1.,H,max_iter+1,Ym,max_iter+1);
  
  // calculate xm = x0 + Vm*ym (x0 = 0)
  blas_xgemv('n',size,mBroken,1.,V,size,Ym,1,0.,dx,1);
  
  delete[] tau;
  delete[] work;
}

void Solver::GMRESHH(int size, const real *b, real *dx,
		     void (*matVec)( const real*, real* ),
		     int& max_iter, real& tol, real normb)
{
  real dummy;
  real sum;
  real sigma;
  int info;
  
  static real *z=0, *normsq=0, *g=0, *c=0, *s=0, *W=0, *H=0;
  static int alc = 0;
  
  // in the worst case (max_iter = size) this amount of memory is needed:
  // z, normsq, g: n+1, c, s: n, W, H: n*(n+1)+1
  allocateSolver(size, alc, &z, &normsq, &g, &c, &s, 0, 0, 0, &W, &H);

  // calculate r0
  // (we are starting with x0 = 0, so r0 = b and there's nothing to do)
  
  // load the starting vector r0 = b into z
  blas_xcopy(size,b,1,z,1);
  z[size] = 0; // this is only needed for technical reasons if max_iter=size;
  
  if ( normb == 0 ) normb = blas_xnrm2(size, b, 1);
  if ( normb == 0 ) normb = 1.0;
  
  // compute the Householder unit vectors
  for(int j = 0; j < max_iter+1; j++) {
    
    // calculate the Householder vectors w_j, their norm and P_jz
    for(int i = 0; i < j; i++) {
      //W[i + j*size] = 0.;
      H[i + j*size] = z[i];
    }
    sum = 0.;
    for(int i = j; i < size; i++) {
      sum += z[i]*z[i];
      W[i + j*size] = z[i];
      //H[i + j*size] = 0;
    }
    W[j + j*size] += sign(z[j])*sqrt(sum);
    H[j + j*size]  = z[j] - W[j + j*size];
    normsq[j]      = 2*(sum + fabs(z[j])*sqrt(sum));

    if( j > 0 ) {
      // apply old rotations to new column of H
      for(int i = 1; i < j; i++) {
	dummy = H[(i-1) + j*size];
	H[(i-1) + j*size] =  c[i-1]*dummy + s[i-1]*H[i + j*size];
	H[  i   + j*size] = -s[i-1]*dummy + c[i-1]*H[i + j*size];
      }

      // calculate c_i and s_i for the new plane-rotation
      sigma = sqrt(H[(j-1) + j*size]*H[(j-1) + j*size]
		   + H[j + j*size]*H[j + j*size]);
      s[j-1] = H[  j   + j*size]/sigma;
      c[j-1] = H[(j-1) + j*size]/sigma;
      
      // apply new rotation to new column of H and to g
      H[(j-1) + j*size] = sigma;
      H[  j   + j*size] = 0;
      dummy             = g[j-1];
      g[j-1]            =  c[j-1]*dummy;
      g[ j ]            = -s[j-1]*dummy;
      
      // check the residual
//       printf("Residual:                  %e\n", fabs(g[j]));
//       printf("Residual devided by b:     %e\n", fabs(g[j])/blas_xnrm2(size,b,1));
//       if(normb > 0) {
// 	printf("Residual devided by normb: %e normb: %e\n", fabs(g[j])/normb, normb);
//       }
      if (fabs(g[j])/normb < tol) {
        tol      = fabs(g[j])/normb;
        max_iter = j;
        break;
      }
    }
    else g[0] = H[0];

    if(j < max_iter) {
    // calculate v and store it in dx
    for(int i = 0; i < j; i++) {
      dx[i] = 0.;
    }
    dx[j] = 1. - 2.*W[j + j*size]*W[j + j*size]/normsq[j];
    for(int i = j+1; i < size; i++) {
      dx[i] = -2.*W[i + j*size]*W[j + j*size]/normsq[j];
    }
    
    for(int i = j-1; i >= 0; i--) {
      sigma = 0.;
      for(int k = i; k < size; k++) {
	// calculate dot product of w_i and P_i...P_je_j
	sigma += W[k + i*size]*dx[k];
      }
      sigma *= 2./normsq[i];
      for(int k = i; k < size; k++) {
	dx[k] -= W[k + i*size]*sigma;
      }
    }
    
    // calculate P_j...P_1Av
      matVec(dx, z);
      for(int i = 0; i <= j; i++) {
	sigma = 0.;
	for(int k = i; k < size; k++) {
	  // calculate dot product of w_i and P_i-1...P_1v
	  sigma += W[k + i*size]*z[k];
	}
	sigma *= 2./normsq[i];
	for(int k = i; k < size; k++) {
	  z[k] -= W[k + i*size]*sigma;
	}
      }
    }
  }
  
  // calculate ym
  lapack_xtrtrs('u', 'n', 'n', max_iter, 1, H+size, size, g, size, &info);
  if(info) {
    printf("Lapack-Error in trtrs: Info Code: %d\n", info);
    exit(0);
  }
  
  // Calculate xm = x0 + Vm*ym (x0 = 0)
  for(int i = 0; i < size; i++) {
    dx[i] = 0.;
  }
  for(int j = max_iter - 1; j >= 0; j--) {
    dx[j] += g[j];
    
    // now apply P_j
    sigma = 0.;
    for(int k = j; k < size; k++) {
      // calculate dot product of w_j and (eta_j e_j + z)
      sigma += W[k + j*size]*dx[k];
    }
    sigma *= 2./normsq[j];
    for(int k = j; k < size; k++) {
      dx[k] -= W[k + j*size]*sigma;
    }
  }
  
//   real *test = new real[size];
//   matVec(dx, test);
//   //vectPrint(size, b, 10);
//   //vectPrint(size, test, 10);
//   blas_xaxpy(size, -1., b, 1, test, 1);
//   //vectDump(size, test);
//   printf("Test: |bm - b| = %e\n", blas_xnrm2(size, test, 1));
//   delete[] test;
}

void Solver::GMRESHHlapack(int size, const real *b, real *dx,
			   void (*matVec)( const real*, real* ),
			   int &max_iter, real& tol, real normb)
{
  real dummy;
  real sum;
  real sigma;
  
  int info;
  real *tau  = new real[size];
  int lwork  = size*64;
  real *work = new real[lwork];
  
  static real *z=0, *normsq=0, *g=0, *c=0, *s=0, *W=0, *H=0;
  static int alc = 0;
  
  // in the worst case (max_iter = size) this amount of memory is needed:
  // z, normsq, g: n+1, c, s: n, W, H: n*(n+1)+1
  allocateSolver(size, alc, &z, &normsq, &g, &c, &s, 0, 0, 0, &W, &H);
  real *HBackup = new real[(size+1)*(size+1)];
  //  real *VVectors = new real[size*(size+1)];
  
  // calculate r0
  // (we are starting with x0 = 0, so r0 = b and there's nothing to do)
  
  // load the starting vector r0 = b into z
  blas_xcopy(size,b,1,z,1);
  z[size] = 0; // this is only needed for technical reasons if max_iter=size;
  
  if ( normb == 0 ) normb = blas_xnrm2(size, b, 1);
  if ( normb == 0 ) normb = 1.0;

  // compute the Householder unit vectors
  for(int j = 0; j < max_iter+1; j++) {
    
    // calculate the Householder vectors w_j, their norm and P_jz
    for(int i = 0; i < j; i++) {
      W[i + j*size]       = 0.;
      H[i + j*size]       = z[i];
      HBackup[i + j*size] = z[i];
    }
    sum = 0.;
    for(int i = j; i < size; i++) {
      sum += z[i]*z[i];
      W[i + j*size] = z[i];
      H[i + j*size] = 0;
      HBackup[i + j*size] = 0;
    }
    W[j + j*size]       += sign(z[j])*sqrt(sum);
    H[j + j*size]        = z[j] - W[j + j*size];
    HBackup[j + j*size]  = z[j] - W[j + j*size];
    normsq[j]            = 2*(sum + fabs(z[j])*sqrt(sum));

    if( j > 0 ) {
      // apply old rotations to new column of H
      for(int i = 1; i < j; i++) {
	dummy = H[(i-1) + j*size];
	H[(i-1) + j*size] =  c[i-1]*dummy + s[i-1]*H[i + j*size];
	H[  i   + j*size] = -s[i-1]*dummy + c[i-1]*H[i + j*size];
      }

      // calculate c_i and s_i for the new plane-rotation
      sigma = sqrt(H[(j-1) + j*size]*H[(j-1) + j*size]
		   + H[j + j*size]*H[j + j*size]);
      s[j-1] = H[  j   + j*size]/sigma;
      c[j-1] = H[(j-1) + j*size]/sigma;
      
      // apply new rotation to new column of H and to g
      H[(j-1) + j*size] = sigma;
      H[  j   + j*size] = 0;
      dummy             = g[j-1];
      g[j-1]            =  c[j-1]*dummy;
      g[ j ]            = -s[j-1]*dummy;
      
      // check the residual
//       printf("Residual:                  %e\n", fabs(g[j]));
//       printf("Residual devided by b:     %e\n", fabs(g[j])/blas_xnrm2(size,b,1));
//       if(normb > 0) {
// 	printf("Residual devided by normb: %e normb: %e\n", fabs(g[j])/normb, normb);
//       }
      if(fabs(g[j])/normb < tol) {
        tol      = fabs(g[j])/normb;
        max_iter = j;
        break;
      }
      //matPrint(size, max_iter+1, H);
    }
    else g[0] = H[0];

    if(j < max_iter) {
      // calculate v and store it in dx
      for(int i = 0; i < j; i++) {
	dx[i] = 0.;
      }
      dx[j] = 1. - 2.*W[j + j*size]*W[j + j*size]/normsq[j];
      for(int i = j+1; i < size; i++) {
	dx[i] = -2.*W[i + j*size]*W[j + j*size]/normsq[j];
      }
      
      for(int i = j-1; i >= 0; i--) {
	sigma = 0.;
	for(int k = i; k < size; k++) {
	  // calculate dot product of w_i and P_i...P_je_j
	  sigma += W[k + i*size]*dx[k];
	}
	sigma *= 2./normsq[i];
	for(int k = i; k < size; k++) {
	  dx[k] -= W[k + i*size]*sigma;
	}
      }
      
//       blas_xcopy(size, dx, 1, VVectors+j*size, 1);
//       printf("GMRESHH:\n");
//       //vectPrint(size, VVectors+j*size, 5);
//       //printf("Norm: %f\n", blas_xnrm2(size,VVectors+j*size,1));
//       for(int jj = 1; jj <= j; jj++) {
// 	for(int ii = 0; ii < jj; ii++) {
// 	  printf("%d:%d %e ", ii+1, jj+1, blas_xdot(size, VVectors+ii*size, 1, VVectors+jj*size, 1));
// 	}
// 	printf("\n");
//       }
      
      // calculate P_j...P_1Av
      matVec(dx, z);
      for(int i = 0; i <= j; i++) {
	sigma = 0.;
	for(int k = i; k < size; k++) {
	  // calculate dot product of w_i and P_i-1...P_1v
	  sigma += W[k + i*size]*z[k];
	}
	sigma *= 2./normsq[i];
	for(int k = i; k < size; k++) {
	  z[k] -= W[k + i*size]*sigma;
	}
      }
    }
  }
  
  //---------------------------------------------------------------------------

  // now solve the least square problem for ym
  for(int i = 0; i < max_iter+1; i++) {
    z[i] = 0.;
  }
  z[0] = HBackup[0];
  lapack_xgeqrf(max_iter+1,max_iter,HBackup+size,size,tau,work,lwork,&info);
  //printf("Needed size for work array: %f\n", work[0]);
  if(info) {
    printf("Lapack-Error in geqrf: Info Code: %d\n", info);
    exit(0);
  }
  lapack_xormqr('l','t',max_iter+1,1,max_iter,HBackup+size,size,tau,z,max_iter+1,work,lwork,&info);
  if(info) {
    printf("Lapack-Error in ormqr: Info Code: %d\n", info);
    exit(0);
  }
  blas_xtrsm('l','u','n','n',max_iter,1,1.,HBackup+size,size,z,max_iter+1);

  // this lapack driver can be used instead of the individual functions above
//   lapack_xgels('n', max_iter+1, max_iter, 1,
// 	       HBackup+size, size, z, max_iter+1, work, lwork, &info);

  // Calculate xm = x0 + Vm*ym (x0 = 0)
  for(int ii = 0; ii < size; ii++) {
    dx[ii] = 0.;
  }
  for(int j = max_iter - 1; j >= 0; j--) {
    dx[j] += z[j];
    
    // now apply P_j
    sigma = 0.;
    for(int k = j; k < size; k++) {
      // calculate dot product of w_j and (eta_j e_j + z)
      sigma += W[k + j*size]*dx[k];
    }
    sigma *= 2./normsq[j];
    for(int k = j; k < size; k++) {
      dx[k] -= W[k + j*size]*sigma;
    }
  }
  
  //---------------------------------------------------------------------------
  // This can be used instead of the above block
  // Instead of solving the least square problem, FOM is used to find the
  // solution. For this to work you have to comment in the "VVectors" lines!

//   for(int i = 0; i < max_iter; i++) {
//     z[i] = 0.;
//   }
//   z[0] = HBackup[0];
//   int* ipiv = new int[max_iter];
//   lapack_xgetrf(max_iter,max_iter,HBackup+size,size,ipiv,&info);
//   lapack_xgetrs('n',max_iter,1,HBackup+size,size,ipiv,z,max_iter,&info);
//   if(info) {
//     printf("Lapack-Error: Could not solve system! Info Code: %d\n", info);
//     matPrint(size, max_iter, HBackup);
//     exit(0);
//   }
//   delete[] ipiv;

//   // calculate xm = x0 + Vm*ym (x0 = 0)
//   blas_xgemv('n',size,max_iter,1.,VVectors,size,z,1,0.,dx,1);

  //---------------------------------------------------------------------------

//   real *test = new real[size];
//   matVec(dx, test);
//   //vectPrint(size, b, 10);
//   //vectPrint(size, test, 10);
//   blas_xaxpy(size, -1., b, 1, test, 1);
//   //vectDump(size, test);
//   printf("Test: |bm - b| = %e\n", blas_xnrm2(size, test, 1));
//   delete[] test;
  
  delete[] tau;
  delete[] work;
  delete[] HBackup;
  //  delete[] VVectors;
}

//RCS: $Id: conjgradient.h,v 2.1 2005/04/22 11:52:47 nedelec Exp $
//-----------------------------------------------------------------------------
//               conjugate gradient and related iterative methods
//           to solve linear systems: http://www.netlib.org/templates
//-----------------------------------------------------------------------------

#include "cblas.h"

/// iterative methods to solve large sparse linear systems: Conjugate-Gradient
namespace Solver
{
        
//=============================================================================
///              Conjugate Gradient, no Preconditionning
//=============================================================================

void conjGradient(int size, const real * b, real * x,
                  void (*matVect)( const real *, real * ),
                  int & max_iter, real & tol);


//=============================================================================
///              Conjugate Gradient, with Preconditionning
//=============================================================================

void conjGradientPrecond(int size, const real * b, real * x, 
                         void (*matVect)( const real*, real* ),
                         void (*precond)( const real*, real* ),
                         int & max_iter, real & tol);


//=============================================================================
///                      Bi-Conjugate Gradient
//=============================================================================

void biConjGradient(int size, const real * b, real * x, 
                    void (*matVect)( const real*, real* ),
                    void (*matVectTrans)( const real*, real* ),
                    int & max_iter, real & tol);


//=============================================================================
///                 Bi-Conjugate Gradient Stabilized
//=============================================================================

int biCGstab(int size, const real * b, real * x,
             void (*matVect)( const real*, real* ),
             int & iter, real & tol, real normb = 0);


//=============================================================================
///        Bi-Conjugate Gradient Stabilized with Preconditionning
//=============================================================================

int biCGstabPrecond(int size, const real * b, real * x, 
                    void (*matVect)( const real*, real* ),
                    void (*precond)( const real*, real* ),
                    int & iter, real & tol, real normb = 0);

};


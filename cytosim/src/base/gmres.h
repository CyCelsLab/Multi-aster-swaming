//RCS: $Id: gmres.h,v 2.1 2005/04/22 11:52:52 nedelec Exp $
//-----------------------------------------------------------------------------
//               GMRES and related iterative methods
//                     to solve linear systems
//-----------------------------------------------------------------------------

#include "types.h"

/// iterative methods to solve large sparse linear systems: GMRES
namespace Solver
{
    
//=============================================================================
///        this solver uses pure FOM (Full Orthogonalization Method)
//=============================================================================

void FOM(int size, const real* b, real* dx,
         void (*matVec)( const real*, real* ),
         int & max_iter);


//=============================================================================
///       this solver uses pure GMRES (Generalized Minimum Residual)
//=============================================================================

void GMRES(int size, const real *b, real *dx,
           void (*matVec)( const real*, real* ),
           int & max_iter);


//=============================================================================
///           this solver uses GMRES (Generalized Minimum Residual)
///                   with Householder orthogonalisation
//=============================================================================

void GMRESHH(int size, const real *b, real *dx,
             void (*matVec)( const real*, real* ),
             int& max_iter, real& tol, real normb = 0);

void GMRESHHlapack(int size, const real *b, real *dx,
                   void (*matVec)( const real*, real* ),
                   int& max_iter, real& tol, real normb = 0);

};

//RCS: $Id: simobject.cc,v 1.1 2005/04/10 14:39:33 nedelec Exp $

#include "simobject.h"
#include "exceptions.h"
#include "cblas.h"
#include "iowrapper.h"
#include "iomessages.h"
#include "sim_param.h"

//----------------------------------------------------------------------------
SimObject::SimObject()
{
  psForces       = 0;
  psDetForces    = 0;
  //the initializations below are not necessary, but nicer when debugging:
  psMobility     = 0;
  psMatrixIndex  = -7;  //this is just a recognizable value for debugging
}


//----------------------------------------------------------------------------
// V <- s*psMobility*V
void SimObject::multMobilityVect( real * V, const real s ) const
{
  blas_xscal(DIM*nbPoints(), s*psMobility, V, 1);      
}


//----------------------------------------------------------------------------
// M <- psMobility*M
void SimObject::multMobilityMat( real * M, const real s ) const
{
  int vecDim = DIM*nbPoints();
  
  for(int ii = 0; ii < vecDim; ++ii) {
    multMobilityVect( M + ii*vecDim, s );
  }
}


//----------------------------------------------------------------------------
/** force <- force + Brownian forces
and return an estimate of Brownian motion to the rhs. */
real SimObject::addBrownianForces( real * force ) const
{
  //we set the forces, it will scaled by MP.km * MP.dt later
  real brownian = sqrt( 2 * MP.kT / ( MP.dt * psMobility )) / MP.km;
  
  for(int jj = 0; jj < DIM*nbPoints(); ++jj )
    force[ jj ] += brownian * RNG.gauss();
  
  //we estimate the contribution of this to vRHS in SIM::solve().
  return ( psMobility * brownian );
}


//----------------------------------------------------------------------------
/** rhs <- rhs + Brownian motion
and return an estimate of Brownian motion to the rhs. */
real SimObject::addBrownianMotion( real * rhs ) const
{
  /** By default we do nothing here!
     This function is used as an alternative addBrownianForces.
     It adds the Brownian motion AFTER projecting the Forces and multiplying with mu.
     If this is used, addBrownianForces should do nothing but return -1! */
  return -1;
}

//****************************************************************************
//****************************************************************************


//----------------------------------------------------------------------------
void SimObject::setProjectedForces(const real * X, real * Y) const
{
  blas_xcopy( DIM * nbPoints(), X, 1, Y, 1 );
}


//----------------------------------------------------------------------------
void SimObject::setSpeedsFromForces(const real * X, real * Y, const real sc) const
{
  assert( X != Y );
  setProjectedForces( X, Y );
  multMobilityVect( Y, sc );
}


void SimObject::setSpeedsWithCorrection(const real * X, real * Z, real * Y, const real sc) const
{
  assert( X != Y );
  assert( X != Z );
  assert( Y != Z );
  addProjectionDiff( X, Z );
  setSpeedsFromForces( Z, Y, sc );
}



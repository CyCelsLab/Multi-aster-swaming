//RCS: $Id: simobject.h,v 1.1 2005/04/10 14:39:29 nedelec Exp $


#ifndef SIMOBJECT_H
#define SIMOBJECT_H

#include "types.h"
#include "pointset.h"

class Matrix;

///Set of points that can be simulated. Parent of Microtub, Nucleus and Solid. 
class SimObject : public PointSet
{
  
protected:

  ///const pointer to an array allocated in SIM::solve(), containing the forces on the mt
  const real   *  psForces; 

  ///const pointer to an array allocated in SIM::solve(), containing the deterministic forces (without Brownian)
  const real   *  psDetForces;
  
  ///mobility (speed = mobilily * force). 
  ///We generally assume that all points have the same mobility.
  real            psMobility;

  ///position in the vectors used to solve the linear system, see SIM::solve().
  ///index for the X-coordinate of the first point. we have then X1, Y1, Z1, X2, Y2...
  int             psMatrixIndex;

  //--------------------------------------------------------------

public:

  //The constructor resets the pointers
  SimObject();
      
  ///Destructor does nothing, as memory is allocated in SIM::solve()
  virtual ~SimObject() {}

  //--------------------------------------------------------------

  ///stores the index where the coordinates are located in SIM::solve()
  void          setMatIndex(const int index)  { psMatrixIndex = index; }
  
  ///index for the X-coordinate of the first point. we have then X1, Y1, Z1, X2, Y2...
  inline int    matIndex()              const { return psMatrixIndex; }

  ///index for the X-coordinate of the first point. we have then X1, Y1, Z1, X2, Y2...
  inline int    matIndexDIM()           const { return DIM * psMatrixIndex; }

  //--------------------------------------------------------------

  ///set the pointer to array of forces, allocated in SIM::solve()
  void          setForcesAddr(const real * f) { psForces = f; }
  
  ///array of forces
  const real *  getForcesAddr()         const { return psForces; }
  
  ///thee force at the given point, as calculated in the last SIM::solve()
  Vecteur       forcesP(const int p)    const { if ( psForces ) return Vecteur( psForces + DIM*p ); return VZERO; }
  
  ///set a pointer to real allocated in solve(), containing the deterministic forces (without Brownian)
  void          setDetForcesAddr(const real * f) { psDetForces = f; }
  
  ///the deterministic part of the forces at the given point, calculate in the last SIM::solve()
  Vecteur       detForcesP(const int p) const { if ( psDetForces) return Vecteur( psDetForces + DIM*p ); return VZERO; }

  //--------------------------------------------------------------

  ///returns the diagonal element of the mobility-matrix (speed = mobility * force)
  virtual real  getMobility(const int d) const { return psMobility; }
  
  ///should update the mobility parameter, as a function of size (called in solve() )
  virtual void  setMobility() = 0;
  
  ///adds the Brownian forces to argument (before projection), returning a magnitude estimate
  virtual real  addBrownianForces(real * force) const;
  
  ///adds the displacements due to the Brownian motion to vRHS (after projection and multiplication with mu*dt)
  virtual real  addBrownianMotion(real * rhs) const;

  //--------------------------------------------------------------

  ///tells SIM::solve() to use preconditionning on this block or not
  virtual bool  precondition()          const { return ( nbPoints() < 10 ); }
  
  //--------------------------------------------------------------
  //rigidity: this is force acting internally to the points of the objects
  //for example, the flexural rigidity of Fibers. There are two independent
  //mechanism to include them: either in the matrix, or directly in the vectors
  
  ///add a rigidity term Y <- Y + Rigidity * X, used e.g. for Microtub flexibility
  virtual void  addRigidity(const real * X, real * Y) const {}
  
  ///add a rigidity term to upper part of matrix mB, used e.g. for Microtub flexibility
  virtual void  setRigidityUp(Matrix & mB, const int offset) const {}
  
  //--------------------------------------------------------------
  //---functions associated with the object's internal constraints
  
  ///re-establish the internal constraints of the object by moving the points
  ///this is needed to correct second-order errors due to linearization, called by SIM::solve()
  virtual void  reshape() = 0;
  
  ///should set the variables needed to call setProjectedForces()
  virtual void  prepareProjection() {}
  
  ///setProjectedForces( X, Y ) calculates Y <- ProjectionMatrix * X
  virtual void  setProjectedForces(const real * X, real * Y) const;
  
  ///multiplies a vector V with the diag. mobility-matrix mu, and with a scalar s if given
  virtual void  multMobilityVect(real* V, const real s=1.0) const;
    
  ///calculates the speed of points in Y, for the forces given in X
  virtual void  setSpeedsFromForces(const real * X, real * Y, const real sc=1.0) const;

  ///multiplies a matrix M with the diag. mobility-matrix mu, and with a scalar s if given
  void          multMobilityMat(real* M, const real s=1.0) const;

  //---Correction term for the projection
  ///set the Projection correction term, P', from the given forces
  virtual void  prepareProjectionDiff(const real * forces) {}
  
  ///addProjectionDiff(real *X, real *Y) should perform Y <- Y + P' * X;
  virtual void  addProjectionDiff(const real * X, real * Y) const {}

  ///calculates the speed of points in Y, for the forces given in X, adding the correction term P'
  virtual void  setSpeedsWithCorrection(const real * X, real * Z, real * Y, const real sc=1.0) const;
      
};

#endif

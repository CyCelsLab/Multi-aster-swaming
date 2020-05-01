//RCS: $Id: sim_solve.cc,v 2.34 2005/04/24 22:49:51 nedelec Exp $
//--------------------------------------------------------------------------
//                            simsolve.cc
//--------------------------------------------------------------------------
//                    solve the equations of motion
//                  implicit integration, sparse matrix
//--------------------------------------------------------------------------

///\todo We would like to know the Lagrange multipliers along the constraints,
///because they are physically relevant. We could for example display them to
///see the compression/traction forces along the Microtub axes.
///We should try Lagrangian Dynamics.

#include <time.h>
#include "clapack.h"
#include "matblock.h"
#include "matsym.h"
#include "matsparse.h"
#include "matsparsesym.h"
#include "matsparsesym2.h"
//#include "matsparsebandsym.h"
#include "conjgradient.h"
//#include "gmres.h"
#include "point_exact.h"

///a namespace to contain the 'global' variables of simsolve.cc
namespace SimSolveGlobals 
{
  
  /**
  Solves the motion of point-objects defined by the dynamic equation:
   d vPTS/dt = mobility * km * mP * ( vBAS + ( mB + mC + mdiffP ) * vPTS )
   
   vPTS is a vector containing all the coordinates (x, y, z) of all the points
   defining the filaments. size(V) = DIM * nbpts;
   
   mB is the isotropic part of the linearization of the forces,
   i.e. it has the same values on the different dimension x, y and z, 
   it is square of size nbpts, symmetric and sparse.
   
   mC is the non-isotropic part of the linearization of the forces
   it can include the confinement forces, or and other interactions
   is is square of size DIM*nbpts, symmetric and sparse
   
   mP is a matrix of projection on the direction allowed by the contraints
   of constant length on the MT segments. it is diagonal by blocks, each
   block corresponding to a PointSet. P is symmetric of size DIM*nbpts
   it is not calculated here: each block is done by the corresponding PointSet
   
   mdiffP is the term comming from the derivative of the projection P,
   we put it under the projection, for reasons of numerical stability
   
   vBAS is a vector of same size as vPTS, 
   it includes the positions of grafted hands, calibrated random forces
   to simulate brownian motion and offsets for periodic boundary conditions;
   */

/// CPU_CHECK_INTERVAL defines the interval (in time steps) at which 
/// Preconditionned and un-preconditionned solvers are compared for cpu-time, 
/** ( the fastest is used until the next cpu-speed comparison is performed )
The iterative solver converges faster with preconditionning, but each
iteration is more expensive, so actually comparing the CPU time is probably
the only safe way to find which one is fastest */
const int CPU_CHECK_INTERVAL = 250;      //(range 100 to 1000)

/// INCLUDE_RIGIDITY Toggles between including the rigidity terms in mB, or
/// calculating them on the fly, which is faster.
/// if enabled, use an appropriate matrix type for mB, for example sparseBand
//#define INCLUDE_RIGIDITY              //(normally OFF)

//include the correction terms for stability into the preconditionner
//BCG should converge a little better, but with a small additional CPU cost
//#define INCLUDE_CORRECTION_IN_PRECONDITIONNER
  
//reset the initial guess for the iterative method at every step
//this option seems to converge faster, and makes more sense:
//the objects in the simulation are mixed are every time step,
//and therefore swap indices in the matrix. Hence the solution at a 
//particular time is mostly uncorrelated with the next time step.
//Conclusion: you should absolutely enable this option!!!
#define INITIAL_GUESS_IS_ZERO
  
// numerical integration scheme
// SOLVER_IMPLICIT is the most stable, the choice is clear!

//#define SOLVER_EULER
#define SOLVER_IMPLICIT

///\todo implement MONITOR_STABILITY,
///automatically detect numerical instabilities (not implemented)
//#define  MONITOR_STABILITY
  
//--------------------------------------------------------------------------

Array<SimObject *> objs; ///<the list of objects containing points to simulate

int  nbpts = 0;         ///<the total number of points in the system
int  maxinblock = 0;    ///<the maximum nb of points encountered in pointsets

//---------------------------------------------------------vectors:
// all vectors below will be of size DIM * nbpts

real * vPTS = 0;         ///< position of the points
real * vSOL = 0;         ///< the points at the previous steps
real * vBAS = 0;         ///< base points of forces
real * vRHS = 0;         ///< right hand side of the final system
real * vTMP = 0;         ///< intermediate of calculus, TEST & DEBUG
real * vFOR = 0;         ///< the calculated forces, available to the simulation
real * vDET = 0;         ///< the calculated deterministic forces, without the Brownian component

//---------------------------------------------------------matrices:
MatrixBlock             mD;      ///< block-diagonal preconditionner.
MatrixSparseSymmetric2  mB;      ///< isotropic symmetric part of the dynamic 
MatrixSparseSymmetric2  mC;      ///< non-isotropic symmetric part of the dynamic 
                                 //MatrixSparseSymmetric or MatrixSparseSymmetric2 are compatible

///we use flags to use only the necessary matrices:
bool  use_mB;   ///< the matrix mB is non-zero and should be used
bool  use_mC;   ///< the matrix mC is non-zero and should be used

bool  dump_now = false;  ///< if true, dump all the matrix for matlab debugging
};

//we make this namespace available in this file only:
using namespace SimSolveGlobals;

#include "sim_interactions.cc"

//----------------------------------------------------------------------------------
void SIM::addTestForces()
{
    Solid * soi, * soj;

  //========    Coulomb interactions between solid points:
  int ii, jj;

  for(soi=firstSolid(); soi ; soi=soi->next())
  for(soj=soi, soj=soj->next(); soj ; soj=soj->next() )
    for( ii = 0; ii < soi->nbPoints(); ++ii )
      for( jj = 0; jj < soj->nbPoints(); ++jj )
      interactionCoulomb( soi->whereP(jj) - soj->whereP(ii), 
                          soi->matIndex() + ii, soj->matIndex() + jj );

}


//==========================================================================
//                                ALLOCATE
//==========================================================================


void SIM::solveAllocate()
{
  static int nbpts_alc = 0;

  //------------  allocate the block matrix for the preconditionner:
  mD.allocate( objs.size() ); 

  //------------  set up the variables for each object:
  nbpts = 0;
  maxinblock = 0;

  for(int ii = 0; ii < objs.size(); ++ii) {
    //The third argument sets the corresponding block in mD 
    //to behave as the Identity if precondition() return 0:
    mD.setBlockSize( ii, DIM * objs[ii]->nbPoints(), !objs[ii]->precondition() );
    objs[ii]->setMatIndex( nbpts );    //starting point in vect. V, W, A
    nbpts        += objs[ii]->nbPoints();
    maxinblock    = maxT( maxinblock, objs[ii]->nbPoints() );
  }
  assert( nbpts > 0 );

  ///\todo catch bad_alloc errors

  //the second argument of the allocation is only for sparseBand matrix
  mB.allocate( nbpts, 2 );
  mC.allocate( DIM * nbpts, 2*DIM-2 );

  if ( nbpts > nbpts_alc ) {
    
      nbpts_alc = nbpts;

      if ( vBAS )   delete[] vBAS;
      if ( vPTS )   delete[] vPTS;
      if ( vSOL )   delete[] vSOL;
      if ( vTMP )   delete[] vTMP;
      if ( vRHS )   delete[] vRHS;
      if ( vFOR )   delete[] vFOR;
      if ( vDET )   delete[] vDET;

      vBAS   = new real[ DIM * nbpts ];
      vPTS   = new real[ DIM * nbpts ];
      vSOL   = new real[ DIM * nbpts ];
      vTMP   = new real[ DIM * nbpts ];
      vRHS   = new real[ DIM * nbpts ];
      vFOR   = new real[ DIM * nbpts ];
      vDET   = new real[ DIM * nbpts ];

      if ((vBAS == 0) || (vPTS == 0) || (vSOL == 0) || (vTMP == 0) || (vRHS == 0) || (vFOR == 0) || (vDET == 0)) {
        fprintf(stderr, "SIM::solveAllocate() memory allocation failed\n");
        exit(1);
      }
      
#ifndef INITIAL_GUESS_IS_ZERO
      //we reset vSOL if the size has changed, because it is used as an initial
      //guess for the biCGstab, and therefore not reset at every time step
      blas_xzero( DIM*nbpts, vSOL, 1 );
#endif
  }
}



//==========================================================================
//                          MAT_VECT + PROJECT...
//==========================================================================


void addForces( const real * X, real * Y )
  // Just the linear part of forces :  Y <- Y + ( mB + mC ) * X;
{
  assert( X != Y );

#if ( DIM > 1 )
#ifndef INCLUDE_RIGIDITY
  for(int ii = 0; ii < objs.size(); ++ii) {
    objs[ii]->addRigidity( X + objs[ii]->matIndexDIM(), Y + objs[ii]->matIndexDIM() );
  }
#endif
#endif
      
  if ( use_mB )        // Y <- Y + mB * X
    mB.multiplyVectorAddIsotropic( X, Y );
  
  if ( use_mC )        // Y <- Y + mC * X
    mC.multiplyVectorAdd( X, Y );
}


//----------------------------------------------------------------------------------
// calculate the forces / MP.km for all points in Y
void calculateForces( const real * X, real * Y )
{
  assert(( X != Y ) && ( vBAS != X ) && ( vBAS != Y ));
  blas_xcopy( DIM*nbpts, vBAS, 1, Y, 1 );      // Y <- BAS
  addForces( X, Y );
}



//----------------------------------------------------------------------------------
void setSpeedsFromForces( const real * X, real * Y, const real sc )
{
  assert( X != Y );
  for(int ii = 0; ii < objs.size(); ++ii) {
    objs[ii]->setSpeedsFromForces( X + objs[ii]->matIndexDIM(), Y + objs[ii]->matIndexDIM(), sc );
  }
}

//----------------------------------------------------------------------------------
// this calculates the matrix product for the implicit integration
// it is called in the conjugate gradient algorithm
// Y <- ( I + dt * P ( mB + mC )) * X;
// uses vTMP as a temporary vector for the calculation !
void matVect( const real * X, real * Y )
{
  assert(( X != Y ) && ( X != vTMP ) && ( Y != vTMP ));

  blas_xzero( DIM*nbpts, vTMP, 1 ); 
  addForces( X, vTMP );

  for(int ii = 0; ii < objs.size(); ++ii) {
    objs[ii]->setSpeedsWithCorrection( X + objs[ii]->matIndexDIM(), vTMP + objs[ii]->matIndexDIM(), Y + objs[ii]->matIndexDIM(), -MP.dt*MP.km );
    //objs[ii]->addProjectionDiff( X + objs[ii]->matIndexDIM(), vTMP + objs[ii]->matIndexDIM() );
    //objs[ii]->setSpeedsFromForces( vTMP + objs[ii]->matIndexDIM(), Y + objs[ii]->matIndexDIM(), -MP.dt*MP.km );
  }
  
  blas_xaxpy( DIM*nbpts, 1.0, X, 1, Y, 1 );
}


//----------------------------------------------------------------------------------
void matVectTrans( const real *, real *)
{
  MSG.error("matVectTrans","disabled: the system is non symmetric");
}


//==========================================================================
//======================   PRECONDITIONNING   ==============================
//==========================================================================

void duplicateMat( int ps, const real * X, real * Y )
  //duplicates X into Y, DIM times, and symmetrize it
{
  real xx = 0;
  int d, kk, ii, jj, ll;
  int bs = DIM * ps;

  blas_xzero( bs*bs, Y, 1 );

  for( ii = 0; ii < ps; ++ii ) {
    xx = X[ ii + ps * ii ];
    
    kk = ( bs+1 ) * DIM * ii;
    for( d = 0; d < DIM; ++d, kk += bs+1 ) 
      Y[ kk ] = xx;
    
    for( jj = ii+1; jj < ps; ++jj ) {
	  xx = X[ ii + ps * jj ];
	  kk = DIM * ( ii + bs * jj );
	  ll = DIM * ( jj + bs * ii );
	  for( d = 0; d < DIM; ++d, kk += bs+1, ll += bs+1 ) {
        Y[ kk ] = xx;
        Y[ ll ] = xx;
      }
	}
  }
  /*
  if ( MP.debug == 98 ) {
      printf("\nmB:\n");
      matPrint( ps, ps, X );
      matPrint( bs, bs, Y );
    }
  */
}


//----------------------------------------------------------------------------------
// return the block corresponding to mt on the diagonal of the total matrix
// i.e. I + P ( B + C + P' ). P, C and D are symmetric, but not their product !

///\todo P is not longer symmetric for the nucleus!!
///check if this is a problem anywhere in simsolve!

void SIM::getBlock( const SimObject * mt, real * DB )
{
  int ps = mt->nbPoints();
  int bs = DIM * ps;
  
  //----- we allocate working arrays needed:
  static real * tmp1 = 0;
  static real * tmp2 = 0;
  static int allocated = 0;
  if ( ps > allocated ) {
    if ( tmp1 ) delete[] tmp1;
    if ( tmp2 ) delete[] tmp2;
    allocated = ps;
    tmp1 = new real[ maxT( ps * ps, bs ) ];
    tmp2 = new real[ bs * bs ];
    if ((tmp1 == 0) || (tmp2 == 0)) {
      fprintf(stderr, "SIM::getBlock() memory allocation failed\n");
      exit(1);
    }
  }

  blas_xzero( ps*ps, tmp1, 1 );

#if ( DIM > 1 )
#ifndef INCLUDE_RIGIDITY
  MatrixSymmetric adapterMatrix( ps, tmp1 );
  mt->setRigidityUp( adapterMatrix, 0 );
#endif
#endif
  
  if ( use_mB )
    mB.addTriangularBlock( tmp1, mt->matIndex(), ps );

  duplicateMat( ps, tmp1, tmp2 );
  
  if ( use_mC )
    mC.addDiagonalBlock( tmp2, DIM * mt->matIndex(), bs );

  //printf("dynamic block:\n"); matPrint( bs, bs, tmp2 );

#ifdef INCLUDE_CORRECTION_IN_PRECONDITIONNER
  //include the stretch correction P' in preconditioner, vector by vector:
  blas_xzero( bs, tmp1, 1 );
  for( int ii = 0; ii < bs; ++ii ) {
    tmp1[ ii ] = 1;
    mt->addProjectionDiff( tmp1, tmp2 + bs * ii );
    tmp1[ ii ] = 0;
  }
#endif
  
  //printf("dynamic with P'\n"); matPrint( bs, bs, tmp2 );

  //compute the projection, by applying it to each column vector:
  for( int ii = 0; ii < bs; ++ii )
    mt->setSpeedsFromForces(tmp2 + bs * ii, DB + bs * ii, -MP.dt );
  
  //add the identity matrix: 
  for( int ii=0; ii < bs*bs; ii += bs+1)
    DB[ ii ] += 1.0;
  //real one = 1.0; blas_xaxpy( bs, 1.0, &one, 0, DB, bs+1 ); //mac bug here
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//the preconditionner is obtained by inverting each diagonal block
int SIM::computePreconditionner()
{
  int  bs, info;
  static int allocated = 0;
  static int * iPivot = 0;
  static int workSize = 0;
  static real * work  = 0;

  for(int ii = 0; ii < objs.size(); ++ii)
    if ( ! mD.treatBlockLikeIdentity(ii) ) {

    real      * DB = mD.getBlock( ii );
    assert( mD.getBlockSize(ii) == DIM*objs[ii]->nbPoints() );
    assert( DB != 0 );
 
    //we calculate the block corresponding to this microtubule:
    getBlock( objs[ii], DB );

    //allocate the tmp arrays needed for the inversion
    bs = DIM * objs[ii]->nbPoints();

    if ( bs > allocated ) {
      allocated = bs;
      if ( iPivot ) {
        delete[] iPivot;
        delete[] work;
      }
      //ask lapack for the optimal size of the array:
      real optimalSize;
      lapack_xgetri( bs, 0, bs, 0, &optimalSize, -1, &info );
      if ( info ) return 0;      //failed to query the optimal size !!!
      //printf("SIM::computePreconditionner() lapack::getri query: size %i optimalSize %.0f\n", bs, optimalSize);
      workSize = (int)optimalSize;
      iPivot = new int[ allocated ];
      work   = new real[ workSize ];
      if ((work == 0) || (iPivot == 0)) {
        fprintf(stderr, "SIM::computePreconditionner() memory allocation failed\n");
        exit(1);
      }      
    }

    //invert the matrix DB by LU factorization:
    info = 0;
    
    lapack_xgetrf( bs, bs, DB, bs, iPivot, &info );
    if ( info ) return 0;      //failed to factorize matrix !!!
    
    lapack_xgetri( bs, DB, bs, iPivot, work, workSize, &info );
    if ( info ) return 0;      //failed to invert matrix !!!
  }
  return 1;                    //successfull
}



//----------------------------------------------------------------------------
// apply the preconditionner   Y <- inv(M) * X
void preconditionner( const real * X, real * Y)
{
  mD.multiplyVector( X, Y );
}


//==========================================================================
//=============================   COPY    ==================================
//==========================================================================


void SIM::copyPoints( int dir )
{
  SimObject * ps = 0;
  
  switch( dir ) {
    
    case 0:         //----------------------------- objects -> vPTS
      
      for(int ii = 0; ii < objs.size(); ++ii) {
        ps = objs[ ii ];
        blas_xcopy( DIM*ps->nbPoints(), ps->getPts(), 1, vPTS + ps->matIndexDIM(), 1);
      }
      break;
      
    case 1:          //----------------------------- objects <- vPTS
      
      for(int ii = 0; ii < objs.size(); ++ii) {
        ps = objs[ ii ];
        blas_xcopy( DIM*ps->nbPoints(), vPTS + ps->matIndexDIM(), 1, ps->pts(), 1);
        ps->setForcesAddr( vFOR + ps->matIndexDIM() );
        ps->setDetForcesAddr( vDET + ps->matIndexDIM() );
        ps->reshape();
      }
      break;
      
    case 2:          //----------------------------- objects <- objects + vSOL
      
      for(int ii = 0; ii < objs.size(); ++ii) {
        ps = objs[ ii ];
        blas_xaxpy(DIM*ps->nbPoints(), 1., vSOL + ps->matIndexDIM(), 1, ps->pts(), 1);
        ps->setForcesAddr( vFOR + ps->matIndexDIM() );
        ps->setDetForcesAddr( vDET + ps->matIndexDIM() );
        ps->reshape();
      }
      break;
  }
}

//==========================================================================
//++++++++++++++++++++++++++  PREPARE  +++++++++++++++++++++++++++++++++++++
//==========================================================================

void SIM::solveSetMatrices()
{  
  //---- reset matrix B and C

  mB.setToZero();
  mC.setToZero();
  
  //---- initialize vBAS to zero
  blas_xzero( DIM*nbpts, vBAS, 1 );
  
  //---- prepare for the projections
  int ii;
  SimObject * ps = 0;
  Vecteur w, b;

  for(ii = 0; ii < objs.size(); ++ii) {
    
    ps = objs[ ii ];
    ps->prepareProjection();
    ps->setMobility();
            
#if ( DIM > 1 )
#ifdef INCLUDE_RIGIDITY
      ps->setRigidityUp( mB, ps->matIndex() );  //include the rigidity term in matrix mB
#endif
#endif

    //--------add the contributions due to confinements
    ///\todo fix: confinements can lead to numerical instabilities, for example
    /// in a Dogterom setup with fast MT_SHRINKING tubes ( dt 0.01, vshrink -2 )
    /// or if MP.km is very high ( 1000 ). This is due to the fact that in one
    /// step, and outside point is brought back inside, turning off the force
    /// at the next step...
    
    if ( ps->isConfined() ) {
      PointExact pte;
      for(int pp = 0; pp < ps->nbPoints(); ++pp ) {
        if ( ! ps->insidePoint( pp )) {
          ps->projectPoint( pp, b );
          pte.setTo( ps, pp );
          ///\todo use interactionLongClamp() instead of interactionPlane more systematically
          /// if shape gives a center and radius to approximate the local boundaries
          /// we should do that in OVAL at least.
          if ( ps->getSpace()->getShape() == SHAPE_SPHERE )
            interactionLongClamp( pte, VZERO, MP.boxsize[0], MP.boxkmratio );
          else
            interactionPlane( pte, ps->whereP( pp ), b, MP.boxkmratio );
        }
      }
    }
    
  }

  //addTestForces();             //DEBUG option
        
  //--------add the contributions of the motor complexes:
  
  for( Complex * co = firstBridgeComplex(); co ; co=co->next() ) {
    //normal complex of length zero
    if ( MP.cxlength[ co->getType() ] == 0 )
      interactionLink( co->getHand1(), co->getHand2(), co->getStiffness() );
    else
      interactionAsymmetricSideLink( co->getHand1(), co->getHand2(), MP.cxlength[co->getType()]);
  }  
  
  //--------add the contributions of the grafted hands:
      
  for( Grafted * gh = firstBoundGrafted(); gh ; gh=gh->next() ) {
    if ( gh->isLink() ) {
      //interactionSideLink( gh->getInterpolation(), gh->getBase(), MP.cxlength[0] );
      interactionLink( gh->getHand(), gh->getBase() );
    } else {
      //interactionStiffClamp( gh->getInterpolation(), gh->whereGrafted(), MP.cxlength[0]);
      interactionClamp( gh->getHand(), gh->whereGrafted(), gh->getStiffness());
    }
  }
  
  //--------add the contributions of the aster microtubule-solid links:
  PointInterpolated pti1, pti2, pti3, pti4;
  PointExact        pte, ptc;
 
  for( Aster * as = firstAster(); as ; as=as->next() ) {
    if ( as->getSolid() ) {
      for( int ii = 0 ; ii < as->maxNbMicrotub(); ++ii ) {
        
        if ( as -> setClamp1( &pte, &ptc, ii ))
          interactionLink( pte, ptc, MP.askmratio );
        
        if ( as -> setClamp2( &pte, &pti1, ii ))
          interactionLink( pti1, pte, MP.askmratio );
      }
    }
  }
  
  //--------add the contributions of the nucleus microtubule-sphere links
  //        and of the links between bundled microtubules
    
  for( Nucleus * nu = firstNucleus(); nu ; nu=nu->next() ) {
    for( int nuindx = 0 ; nuindx < nu->maxNbMicrotub(); nuindx++ ) {
      
      //add links between points on the sphere and mitcrotubules
      if ( nu -> setPointClamp( &pte, &ptc, nuindx ))
        interactionLink( pte, ptc, MP.nukmratio[0] );
            
      //add links between bundled microtubules
      if ( nu -> setMTClamp( &pti1, &pti2, &pti3, &pti4, nuindx )) {
        interactionLink( pti1, pti2, MP.nukmratio[1] );
        interactionLink( pti3, pti4, MP.nukmratio[1] );
      }
    
      //add long links to push out microtubules that went inside the nucleus
      if ( nu->hasMicrotub( nuindx ) ) {
        Microtub* nuMT = nu->getMicrotub( nuindx );
        pte.setTo( nu, 0 );
        for ( int kk=0; kk < nuMT->nbPoints(); kk++ ) {
          if ( (nuMT->whereP( kk ) - nu->getPosition()).normSquare() < nu->radius()*nu->radius() ) {
            ptc.setTo( nuMT, kk );
            interactionLongLink( pte, ptc, nu->radius(), MP.nukmratio[2] );
          }
        }
      }
    }
  }
  //------------------------------------------------------------------------
  //Microtub * mt = findMicrotub(1);
  //if ( mt ) {
  //  PointExact pta( mt, 0 );
  //  interactionClamp( pta, VZERO );
  // }
  //------------------------------------------------------------------------
  
    
  use_mB = mB.nonZero();
  use_mC = mC.nonZero();
  
  //the calls to optimizeForMultiply() are optional:
  mB.optimizeForMultiply();
  mC.optimizeForMultiply();
}



//==========================================================================

void monitorCPUToSetPreconditionning( bool & precond_flag )
{
  //to decide between preconditionning and non-preconditionning:
  //we try both, compare the CPU time spent, and retain the fastest.
  //this is done every CPU_CHECK_INTERVAL iterations,
  //and we pool 3 sucessive steps together

  static clock_t cpu, time[2];

  switch( sim.iterationCnt() % CPU_CHECK_INTERVAL ) 
    {
    case 1:
      cpu = clock();                         // take time after step 1
      break;
    case 4:
      time[ precond_flag ] = clock() - cpu;  // record time for steps 2,3,4
      cpu = clock();                         // take time after step 4
      precond_flag = 1 - precond_flag;       // switch to other method
      break;
    case 7:
      time[ precond_flag ] = clock() - cpu;  // record time for steps 5,6,7
      precond_flag = ( time[0] > time[1] );  // choose fastest method
      
      MSG(5, "CPU at step %8i: no preconditionner %9.2f ms; with precond. %9.2f ms: using P=%i\n",
          sim.iterationCnt(), float(1000*time[0]) / CLOCKS_PER_SEC, float(1000*time[1]) / CLOCKS_PER_SEC, precond_flag);
      break;
    default:
      break;
    }
}

//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&       SOLVE        &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



void SIM::solve()
{
  static bool precond_flag = true;
  //the initial value applies to the first three iterations
  //after which precond is set automatically by monitorCPU().

  //------------  register all objects with points in the list -----------

  objs.clear();

  for(Microtub * mt=firstMicrotub();   mt ; mt=mt->next() )
    objs.pushFront( mt );

  for(Solid *    so=firstSolid();      so ; so=so->next() )
    objs.pushFront( so );

  for(Nucleus *  nu=firstNucleus();    nu ; nu=nu->next() )
    objs.pushFront( nu );

  if ( objs.size() == 0 ) return;

  //-----------------------------------------------------------------------

  solveAllocate();                           // initialization & allocation
  copyPoints( 0 );                           // vPTS <- coordinates of mts/pts
  solveSetMatrices();                        // set all matrices and vectors

  calculateForces( vPTS, vDET );             // vDET = deterministic forces before constraints
  blas_xcopy( DIM*nbpts, vDET, 1, vFOR, 1 ); // vFOR = vDET
  
  //we make an upper estimate of the Brownian contribution to vRHS
  //the system will be solved with a residual MP.tolerance * normRHS_estimate
  //when MP.tolerance is smaller than 1, that should work well.
  real normRHS_estimate = 0;
  
  // add the thermal random forces to vFOR
  // and get the estimate of Brownian motions in vRHS:
  for(int ii = 0; ii < objs.size(); ++ii) {
    real th = objs[ii]->addBrownianForces( vFOR + objs[ii]->matIndexDIM() );
    if ((( normRHS_estimate == 0 ) || ( normRHS_estimate > th )) && ( th > 0 ))
      normRHS_estimate = th;
  }


#ifdef SOLVER_EULER
  //warning: that is broken for Nucleus! Missing the brownian noise.
  setSpeedsFromForces( vFOR, vTMP, MP.dt * MP.km );
  blas_xaxpy( DIM*nbpts, 1.0, vTMP, 1, vPTS, 1 ); 
  copyPoints( 1 );
  return;
#endif
  
  const real tolerance_ask = MP.tolerance;
  real tolerance = tolerance_ask;
    
#ifdef SOLVER_IMPLICIT

  //----------------- OPTIONS FOR PRECONDITIONNING -----------------
  //automatically choose the fastest preconditionning method:
  monitorCPUToSetPreconditionning( precond_flag );
  
  //this is an option to manually control preconditionning from the player:
  //precond_flag = user_control_keys[1];
  
  //to test the preconditionning in sim, you may enable the line below:
  //precond_flag = iterationCnt() % 2;
  
  //With monitorCPUToSetPreconditionning(), the CPU time influences the program. 
  //This input is not controllable by the user. Running two simulations with the
  //same parameters can produce difference output. To make cytosim perfectly 
  //deterministic, the preconditionning should always be the same:
  //precond_flag = true;
  //----------------------------------------------------------------
  
  //we can solve in two equivalent ways:
  // OPT1 :  (I - hPA) ( Xnew - Xold ) = hP * force( Xold )
  // OPT2 :  (I - hPA) Xnew = Xold + hP vBAS
  // with force( X ) = PA * X + vBAS
  // the matrix is the same, just the vectors change.
  // We use here OPT1

  //set the differential of the projection with the current forces:
  for(int ii = 0; ii < objs.size(); ++ii) {
    objs[ii] -> prepareProjectionDiff( vFOR + objs[ii]->matIndexDIM() );
  }
  
  //calculate the final right-hand side of the system in vRHS:
  setSpeedsFromForces( vFOR, vRHS, MP.dt * MP.km );
  
  //we also have to rescale the estimate we have made from vFOR 
  normRHS_estimate *= MP.km * MP.dt;

  //add projected Brownian motion for objects that didn't do
  //that in addBronianForces() above (so far only the nucleus)
  for(int ii = 0; ii < objs.size(); ++ii) {
    real th = objs[ii]->addBrownianMotion( vRHS + objs[ii]->matIndexDIM() );
    if ((( normRHS_estimate == 0 ) || ( normRHS_estimate > th )) && ( th > 0 ))
      normRHS_estimate = th;
  }
  
  // biConjugate Gradient is said (in exact arithmetic) to converge at most
  // in a number of steps equal to the size of the system. We use this limit:
  const int nbiter_ask = DIM*nbpts;
  int nbiter = nbiter_ask;

#ifdef INITIAL_GUESS_IS_ZERO
  // set the initial guess for the solution of the system (Xnew - Xold):
  // we can use the solution at the previous step, or a vector of zeros.
  // Using the previous solution would work only if the objects are not 
  // considered in a random order to build the linear system (our case).
  blas_xzero( DIM*nbpts, vSOL, 1 );
#endif
  
  int code = 0;
  // We now solve the system MAT * vSOL = vRHS  by an iterative method:
  ///\todo implement another iterative method, e.g. gmres
  
  // the tolerance is in scaled to the contribution of Brownian
  // motions contained in vRHS, assuming that the error is equally spread 
  // along all degrees of freedom, this should work for tolerance << 1
  // here a printf() can be used to check that the estimate is correct:
  //printf("normRHS_estimated = %8.2e   variance( vRHS ) / estimate = %8.4f\n", 
  //       normRHS_estimate, blas_xnrm2(DIM*nbpts, vRHS, 1) / (normRHS_estimate * sqrt( DIM*nbpts)) );

  // TEST: the tolerance to solve the system should be such that the solution
  // found does not depend on the initial guess.

  /*
  ///\todo check GMRES, which should perform better than biCGstab
  Solver::GMRESHH( DIM*nbpts, vRHS, vSOL, matVect,
                   nbiter, tolerance, normRHS_estimate );
  */
 
  //------- call the iterative solver:
  
  if ( precond_flag && computePreconditionner() )
      code = Solver::biCGstabPrecond( DIM*nbpts, vRHS, vSOL, matVect, preconditionner,
                                      nbiter, tolerance, normRHS_estimate );
  else
      code = Solver::biCGstab( DIM*nbpts, vRHS, vSOL, matVect,
                               nbiter, tolerance, normRHS_estimate );

  
  //------- in case the solver did not converge, we try other methods:
  if ( tolerance > tolerance_ask ) {
    
    MSG("biCGstab%i failure %i after %i iter. tol %.2e ", 
        precond_flag, code, nbiter, tolerance);

    //---we first try Zero as a different initial guess:
    blas_xcopy( DIM*nbpts, vRHS, 1, vSOL, 1 );
    
    //---reset the desired tolerance and max nb of iterations:
    tolerance = tolerance_ask;
    nbiter    = nbiter_ask;
    
    //---try the same method again:
    if ( precond_flag ) {
      Solver::biCGstabPrecond( DIM*nbpts, vRHS, vSOL, matVect, preconditionner,
                               nbiter, tolerance, normRHS_estimate );
    } else {
      Solver::biCGstab( DIM*nbpts, vRHS, vSOL, matVect,
                        nbiter, tolerance, normRHS_estimate );
    } 
    
    //check again for convergence:
    if ( tolerance <= tolerance_ask ) 
      MSG("...rescued by changing seed: %i iter. tol %.2e\n", nbiter, tolerance);
    else {
        
      //---use zero as an initial guess:
      blas_xzero( DIM*nbpts, vSOL, 1 );

      //---reset the desired tolerance and max nb of iterations:
      tolerance = tolerance_ask;
      nbiter    = nbiter_ask;

      //try the other method:
      if ( precond_flag ) {
        Solver::biCGstab( DIM*nbpts, vRHS, vSOL, matVect,
                          nbiter, tolerance, normRHS_estimate );
      } else {
        if ( computePreconditionner() )
          Solver::biCGstabPrecond( DIM*nbpts, vRHS, vSOL, matVect, preconditionner,
                                   nbiter, tolerance, normRHS_estimate );
        else
          MSG("...failed to compute precondionner");
      }
    
      //check again for convergence:
      if ( tolerance <= tolerance_ask ) 
        MSG("...rescued with other method: %i iter. tol %.2e\n", nbiter, tolerance);
      else {
        //no solver could converge... this is really bad!
        MSG("...dead %i iter. tol %.2e\n", nbiter, tolerance);
        finish_now = true;
        dump_now = true;
        return;
      }
    }
  }
  
  //this dumping enables verification of the matrix operations in matlab
  if ( dump_now ) {
    MSG("simsolve.cc: dumping the matrices and finishing simulation");
    dumpAll();
    finish_now = true;
  }
    
  copyPoints( 2 );   ///< copy back vPTS to the PointSets
  
  
  ///\todo: in sim_solve.cc: warning if the number of iterations is higher than 200
  
  //report on the matrix type and size, sparsity, and the number of iterations
  if (( iterationCnt() % 50 == 0 ) || ( MP.debug > 10 )) {
      MSG(4, "%6i  size %5i * %i", iterationCnt(), nbpts, DIM );
      if ( use_mB ) MSG(4, " mB: %s ", mB.shortDescription() );
      if ( use_mC ) MSG(4, " mC: %s ", mC.shortDescription() );
      MSG(4, " biCGstab%i iter. %4i res. %.2e", precond_flag, nbiter, tolerance );
#ifdef INCLUDE_RIGIDITY
      MSG(4, " i");
#endif
      MSG(4, "\n");
    }

#endif
}






//==========================================================================
//                            DEBUG  -  DUMP
//==========================================================================



//produce a output of matrices which can be imported in matlab for verification
void SIM::dumpAll()
{
#ifndef INCLUDE_RIGIDITY
  MSG.error("dumpAll","incorrect dump if rigidity not included in matB");
#endif
  printf("entering dumpAll() in simsolve.cc\n");
  
  real * tmp1 = new real[ DIM * maxinblock ];
  real * tmp2 = new real[ sqr( DIM * maxinblock ) ];

  if ((tmp1 == 0) || (tmp2 == 0)) {
    fprintf(stderr, "SIM::dumpAll() memory allocation failed\n");
    exit(1);
  }
  
  FILE * f = fopen("dump_DIM", "w"); 
  if ((f==0) || ferror(f)) return;
  fprintf(f,"%% DIM, nbpts, dt, tolerance\n");
  fprintf(f,"%i\n%i\n%.8e\n%.8e\n", DIM, nbpts, MP.dt, MP.tolerance );
  fclose( f );

  f = fopen("dump_mobility", "w"); if ((f==0) || ferror(f)) return;
  fprintf(f,"%% microtubules&solids : name, coord, mobility\n");
  for(int ii = 0; ii < objs.size(); ++ii ) {
    for(int d=0; d < DIM * objs[ii]->nbPoints(); ++d ) 
      fprintf(f,"%5i %16.8e %16.8e\n", ii+1, objs[ii]->coord(d), objs[ii]->getMobility(d) );
  }
  fclose(f);

  f = fopen("dump_sol", "w"); 
  if ((f==0) || ferror(f)) return;
  vectDump(DIM*nbpts, vPTS, f);
  fclose(f);

  f = fopen("dump_rhs", "w");
  if ((f==0) || ferror(f)) return;
  vectDump(DIM*nbpts, vRHS, f);
  fclose(f);

  f = fopen("dump_matB", "w"); 
  if ((f==0) || ferror(f)) return;
  mB.printSparse(f);
  fclose(f);

  f = fopen("dump_matC", "w"); 
  if ((f==0) || ferror(f)) return;
  mC.printSparse(f);
  fclose(f);


  f = fopen("dump_project", "w");
  if ((f==0) || ferror(f)) return;
  for(int oo = 0; oo < objs.size(); ++oo ) {
    int bs = DIM * objs[oo]->nbPoints();

    blas_xzero( bs*bs, tmp2, 1 );
    blas_xzero( bs, tmp1, 1 );
    //compute the matrix, by applying it to each column vector:
    for(int ii = 0; ii < bs; ++ii ) {
      tmp1[ ii ] = 1;
      objs[oo]->setSpeedsFromForces( tmp1, tmp2 + bs * ii );
      tmp1[ ii ] = 0;
    }

    matDumpOffset( bs, bs, tmp2, objs[oo]->matIndexDIM(), f );
  }
  fclose(f);


  f = fopen("dump_projectdiff", "w");
  if ((f==0) || ferror(f)) return;
  for(int oo = 0; oo < objs.size(); ++oo ) {
    int bs = DIM * objs[oo]->nbPoints();

    blas_xzero( bs*bs, tmp2, 1 );
    blas_xzero( bs, tmp1, 1 );
    //compute the matrix, by applying it to each column vector:
    for( int ii = 0; ii < bs; ++ii ) {
      tmp1[ ii ] = 1;
      objs[oo]->addProjectionDiff( tmp1, tmp2 + bs * ii );
      tmp1[ ii ] = 0;
    }
    matDumpOffset( bs, bs, tmp2, objs[oo]->matIndexDIM(), f );
  }
  fclose(f);

  computePreconditionner();
  f = fopen("dump_precond", "w"); 
  if ((f==0) || ferror(f)) return;
  mD.printSparse(f);
  fclose(f);

  f = fopen("dump_diagblock", "w"); 
  if ((f==0) || ferror(f)) return;
  for(int oo = 0; oo < objs.size(); ++oo ) {
    int bs = DIM * objs[oo]->nbPoints();
    getBlock( objs[oo], tmp2 );
    matDumpOffset( bs, bs, tmp2, objs[oo]->matIndexDIM(), f );
  }
  fclose(f);

  delete[] tmp1;
  delete[] tmp2;
}

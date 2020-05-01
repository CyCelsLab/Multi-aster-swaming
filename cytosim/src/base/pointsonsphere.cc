//RCS: $Id: pointsonsphere.cc,v 2.7 2005/04/19 22:09:30 nedelec Exp $
 
#include "pointsonsphere.h"

//-----------------------------------------------------------
PointsOnSphere::PointsOnSphere() 
{ 
  mCoord = 0;
}  
//-----------------------------------------------------------
PointsOnSphere::PointsOnSphere( int nbp ) 
{ 
  mCoord = 0;
  distributePoints( nbp ) ; 
}   
//-----------------------------------------------------------
PointsOnSphere::~PointsOnSphere( ) 
{ 
  if (mCoord) {
      delete[] mCoord;
      mCoord = 0; 
    }
}   
//-----------------------------------------------------------
void PointsOnSphere::copyCoordinatesOfPoint( real x[3], int ii )
{
  x[ 0 ] = mCoord[ 3*ii + 0 ];
  x[ 1 ] = mCoord[ 3*ii + 1 ];
  x[ 2 ] = mCoord[ 3*ii + 2 ];
}

//-----------------------------------------------------------
void PointsOnSphere::copyCoordinatesOfPoint( real * x, real * y, real * z, int ii )
{
  *x = mCoord[ 3*ii + 0 ];
  *y = mCoord[ 3*ii + 1 ];
  *z = mCoord[ 3*ii + 2 ];
}

//-----------------------------------------------------------
void PointsOnSphere::copyPositionsForAllPoints( real x[] )
  //method to copy the array points coordinates onto the given array x
{
  for( int ii = 0; ii < 3*mNumberPts; ++ii )
    x[ ii ] = mCoord[ ii ];
}

//-----------------------------------------------------------
void PointsOnSphere::reportConvergence(FILE * file) 
{
  fprintf(file, "PointOnSphere converged after %i steps to an energy = %f\n",
          mNumberSteps, mEnergy);
}

//-----------------------------------------------------------
void PointsOnSphere::printAllPositions( FILE * file ) 
{
  for( int ii = 0; ii < mNumberPts; ++ii )
    fprintf( file, "%f %f %f\n", mCoord[3*ii], mCoord[3*ii+1], mCoord[3*ii+2]);
}

//-----------------------------------------------------------
void PointsOnSphere::setRandom( real * P )
  /** hypercube rejection method, calling the Random Number Generator RNG */
{
  real n;
  do {
    P[0] = RNG.sreal();
    P[1] = RNG.sreal();
    P[2] = RNG.sreal();
    n = P[0]*P[0] + P[1]*P[1] + P[2]*P[2];
    if ( n == 0 ) {
      //this solves a strange condition only occuring on macs...
      fprintf(stderr, "PointOnSphere: RNG returned 3 zeros: re-seeding\n");
      RNG.seedTimer(); //exit( EXIT_FAILURE);
    }
  } while (( n  > 1.0 ) || ( n == 0.0 ));
  n = sqrt( n );
  P[0] /= n;
  P[1] /= n;
  P[2] /= n;
}


//-----------------------------------------------------------
real PointsOnSphere::distance( real P[], real Q[] )
{
  return sqrt( (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]) );
}

//-----------------------------------------------------------
real PointsOnSphere::distanceSquare( real P[], real Q[] )
{
  return (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]);
}

//-----------------------------------------------------------
real PointsOnSphere::minimumDistance() 
{
  real dist, result = 2.;
  for(int ii = 1; ii < mNumberPts; ++ii )
    for(int jj = 0; jj < ii; ++jj ) {
      dist = distanceSquare( mCoord + 3 * ii, mCoord + 3 * jj );
      if ( dist < result ) 
        result = dist;
    }
  return sqrt(result);
}

//-----------------------------------------------------------
void PointsOnSphere::scaleUniformly(const real factor) 
{
  for(int ii = 0; ii < 3*mNumberPts; ++ii )
    mCoord[ ii ] *= factor;
}


//-----------------------------------------------------------
real PointsOnSphere::coulombEnergy( real P[] ) 
{
  real dist, result = 0;
  for(int ii = 1; ii < mNumberPts; ++ii )
    for(int jj = 0; jj < ii; ++jj ) {
      dist = distance( P + 3 * ii, P + 3 * jj );
      if ( dist > 0 ) result += 1.0 / dist;
    }
  return result;
}

//-----------------------------------------------------------
void PointsOnSphere::calculateCoulombForces( real forces[] ) 
{
  int ii, jj;
  short dd;
  real dx[ 3 ];
  real dist;
  mMinDistance = 4;
  
  //--------- reset forces:
  for( ii = 0; ii < 3 * mNumberPts; ++ii )
    forces[ ii ] = 0;
  
  //--------- calculate coulomb pair interactions:
  // first particle is ii, second one is jj:
  for( ii = 1; ii < mNumberPts; ++ii ) {
    for( jj = 0; jj < ii; ++jj ) {
      
      //calculate vector and distance^2 between from jj to ii
      for( dist = 0, dd = 0; dd < 3 ; ++dd ) {
        dx[ dd ] = mCoord[ 3*ii + dd ] - mCoord[ 3*jj + dd ];
        dist += dx[ dd ] * dx[ dd ];
      }
      
      if ( dist > 0 ) {      //if they do not overlap:
      
        if ( dist < mMinDistance ) mMinDistance = dist;
        
        //force = vector / r^3, but here dist = r^2
        dist = 1.0 / ( dist * sqrt( dist ));
        //update forces for jj and ii:
        for( dd = 0 ; dd < 3; ++dd ) {
          dx[ dd ] *= dist;
          forces[ 3*ii + dd ] += dx[ dd ];
          forces[ 3*jj + dd ] -= dx[ dd ];
	      }
        
      } else {   //if ii and jj overlap, we use a random force
        for( dd = 0 ; dd < 3; ++dd ) {
          dx[ dd ] = 0.1 * RNG.sreal();
          forces[ 3*ii + dd ] += dx[ dd ];
          forces[ 3*jj + dd ] -= dx[ dd ];
        }
      }
    }
  }
  
  mMinDistance = sqrt( mMinDistance );
  //------------- remove centripede contribution of forces:
  // assuming here that points are already on the sphere (norm=1)
  // ( the algorithm converge even without this, but slower )
  
  for( ii = 0; ii < mNumberPts; ++ii ) {
    dist = 0;
    for( dd = 0; dd < 3; ++dd )
      dist += mCoord[ 3*ii + dd ] * forces[ 3*ii + dd ];
    
    for( dd = 0; dd < 3; ++dd )
      forces[ 3*ii + dd ] -= dist * mCoord[ 3*ii + dd ];
  }
}
//-----------------------------------------------------------
// move the points in the directions set by the forces,
// but scaling by the factor S:
void PointsOnSphere::movePointsAccordingToForces( real Pnew[], real Pold[], real forces[], real S )
{
  short dd;
  real norm;
  for(int ii = 0; ii < mNumberPts; ++ii ) {
    //------------- move the point, and calulate norm:
    norm = 0;
    for( dd = 0; dd < 3; ++dd ) {
	  Pnew[ 3*ii + dd ] = Pold[ 3*ii + dd ] + S * forces[ 3*ii + dd ];
	  norm += Pnew[ 3*ii + dd ] * Pnew[ 3*ii + dd ];
    }
      
    //------------- normalize (project on the sphere):
      norm = 1.0/ sqrt( norm );
    for( dd = 0; dd < 3; ++dd )
      Pnew[ 3*ii + dd ] *= norm;
  }
}


//-----------------------------------------------------------
// creates a relatively even distribution of nbp points on the sphere
// the coordinates are stored in real array mCoord
void PointsOnSphere::distributePoints( int nbp, real precision ) 
{
  // with no points, we should get upset, but we just return.
  if ( nbp <= 0 ) return;
  
  // with N points on the sphere, they should each occupy
  // an area of 4*PI / N, the distance should be sqrt(4PI/N)
  // the precision is rescaled with that expected distance:
  real precision_scaled = precision * sqrt( 4 * 3.14159 / nbp );
  
  //reallocate the array if needed:
  if ( mCoord && ( nbp != mNumberPts )) {
    delete[] mCoord;
    mCoord = 0; 
  }
  //fprintf(stderr, "PointsOnSphere::distributePoints(%i)...", nbp);
  
  int nb_good_steps = 0;
  mNumberPts = nbp;
  
  //make an initial guess for the step size:
  real step_size = 10 * precision / mNumberPts, energy_new;
  
  //allocate the array of coordinates if needed:
  if ( mCoord == 0 )
    mCoord         = new real[ mNumberPts * 3 ];
  
  //allocate also the forces, these will be deleted below:
  real * forces    = new real[ mNumberPts * 3 ];
  real * coord_new = new real[ mNumberPts * 3 ];
  
  if (( mCoord == 0 ) || ( forces == 0 ) || (coord_new == 0 )) {
    fprintf(stderr, "PointsOnSphere::distributePoints():: memory allocation failed\n");
    exit(1);
  }
  
  //------------ distribute the points randomly on the sphere:
  for(int ii = 0; ii < mNumberPts; ii++ )
    setRandom( mCoord + ii * 3 );
  
  //--------- for one point only, we return:
  if ( mNumberPts < 2 ) {
    mEnergy = 0;
    mNumberSteps = 0;
    return; 
  }
  
  //------------ calculate the initial energy:
  mEnergy = coulombEnergy( mCoord );
  
  //------------ calculate forces:
  calculateCoulombForces( forces );
  
  for(mNumberSteps = 0; mNumberSteps < max_nb_iterations; ++mNumberSteps ) {
    
    //------------ move the mCoord:
    movePointsAccordingToForces( coord_new, mCoord, forces, step_size );
    
    //------------ calculate the new potential energy
    energy_new = coulombEnergy( coord_new );
    
    //printf("%3i : step %5i : mEnergy = %18.8f   step_size = %8.5f %s\n",
    //     nbp, mNumberSteps, mEnergy, step_size, (energy_new<mEnergy?"yes":"no"));      
    
    //------------ accept of reject new configuration:
    // accept if the new configuration has a lower energy: better repartition
    if ( energy_new < mEnergy ) {
      //---- values for 'magic' below were tested in term of convergence 
      //     few trials seemed to agree for magic = 7... really magic!
      
      // basically, if we have done 'magic' successful moves at a given step size,
      // then we try to increase the step size by a factor 2, 
      if ( ++nb_good_steps >= magic ) {
        step_size *= 3.1415;   //this value is just for fun, 3 would work
        nb_good_steps = 0;
      }
	
      //---- accept the new positions, by just swapping the arrays
      real * t = mCoord;
      mCoord = coord_new;
      coord_new = t;
      
      mEnergy = energy_new;
      //---- re-calculate the forces:
      calculateCoulombForces( forces );
    
    } else {
      //here the new configuration has higher energy:
      //it is rejected, and we try a smaller step size:
      nb_good_steps = 0;
      step_size /= 2;
      //if we are below our precision, we stop trying
      if ( step_size < precision_scaled ) break;
    }
  }
  delete[] coord_new;
  delete[] forces;
  //fprintf(stderr, "done\n", nbp);
}
//-----------------------------------------------------------

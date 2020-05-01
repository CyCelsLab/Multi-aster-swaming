//RCS: $Id: pointsonsphere.h,v 2.4 2005/04/19 22:09:42 nedelec Exp $
//=========================== pointsonsphere.h ===============================
//



#ifndef POINTSONSPHERE_H
#define POINTSONSPHERE_H

#include "smath.h"
#include "random.h"
#include "types.h"


//a limitation of icc obliges us to declare this outside the class
///default precision
static const real  default_precision = 1e-3;        

///  A class to distribute points on the unit sphere, as uniformly as possible
/**  A class to distribute points as uniformly as possible on the unit sphere,
for ANY number of points.
see http://mathworld.wolfram.com/SphericalCode.html

Algorithm:

1. The points are first distributed randomly on the sphere
2. a 1/r^3 repulsive force is assumed for all points, 
to which corresponds a certain potential energy in 1/r^2
3. new positions are calculated form the current one, the forces
and an adjustable scaling factor:  dx = scale * ( force at x )
4. the potential energy of the new configuration is calculated
5. the move is accepted if its total energy is lower,
the move is rejected if the energy is higher,
if the move is rejected, the scaling factor is reduced,
if the move is accepted, the scaling factor is increased
6. steps 2 to 5 are repeated until convergence.

The main method is the class constructor, or distributePoints()
which take the number of points as argument and calls the calculation:
points coordinates can be retreived by the functions
copyPositionsForAllPoints()
or
copyCoordinatesOfPoint()

F. Nedelec, created August 2002, last modified October 10, 2002

\todo generalize PointsOnSphere to any surface defined by isInside() & project()

*/
class PointsOnSphere
{
  ///This really magic number that affect convergence speed, but not the answer
  static const int  magic = 7;

  ///max. number of iterations
  static const int  max_nb_iterations = 50000;        
  
  /// number of point on the sphere
  int  mNumberPts;

  /// coordinates of the points in a array
  /** in the array all the coordinates are together (x,y,z) point 1, (x,y,z) point 2, etc.
    so the coordinates for the first point are:
       x = mCoord[0], y = mCoord[1], z = mCoord[3]
    the coordinates of point ii are:
       x = mCoord[3*i+0], y = mCoord[3*i+1], z = mCoord[3*i+2] 
  */
  real * mCoord;

  /// number of steps in the minimation calculation
  int  mNumberSteps;

  /// Coulomb energy that must be minimize 
  real mEnergy;
  
  /// Smallest distance between 2 points in the configuration
  real mMinDistance;
  
private:
    
  /// set coordinates point P randomly on the sphere
  void setRandom( real * P );
  
  /// Calculate distance between point given their coordinates P and Q (3-dim)
  real distance( real P[], real Q[] );
  
  /// Calculate distance between point given their coordinates P and Q (3-dim)
  real distanceSquare( real P[], real Q[] );

  /// coulomb energy
  real coulombEnergy( real P[] );
  
  /// coumomb force
  void calculateCoulombForces( real forces[] );
  
  /// move point from old to new coordinates
  void movePointsAccordingToForces( real Pnew[], real Pold[], real forces[], real S );
  
  
public:

  /// default constructor, does nothing
  PointsOnSphere();

  /// constructor that also calls distributePoints(), 
  PointsOnSphere( int nbp );

  /// default destructor
  ~PointsOnSphere();

  /// number of points in the configuration
  int    nbPoints()                     { return mNumberPts;  }
  
  /// number of interation performed in the last distributePoints()
  int    nbIterationsToConverge()       { return mNumberSteps; }
  
  /// the 'virtual' total energy of the configuration
  real   finalEnergy()                  { return mEnergy; }
  
  /// minimum distance in the configuration, in 3D space
  real   minimumDistance();

  /// scale all the points by a given factor (call after the calculation)
  void scaleUniformly(const real factor);
  
  /// address where the coordinates for point ii are  
  real * addressForThisPoint(int ii)    { return mCoord + 3 * ii; }

  /// copy the coordinates from point ii onto the given 3-dim array x
  void copyCoordinatesOfPoint( real x[3], int ii );
      
  /// copy the coordinates from point ii onto x,y,z
  void copyCoordinatesOfPoint( real * x, real * y, real * z, int ii );

  /// copy the array points coordinates onto the given array x
  void copyPositionsForAllPoints( real x[] );
  
  /// report results on the file
  void reportConvergence( FILE * file = stdout );

  /// write points coordinates
  void printAllPositions( FILE * file = stdout );
    
  /// distribute the nbp points on the sphere and store their coordinates
  void distributePoints( int nbp, real precision = default_precision );

};

#endif

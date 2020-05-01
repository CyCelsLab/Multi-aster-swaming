//RCS: $Id: quaternion.h,v 2.3 2005/03/31 20:57:48 foethke Exp $
// Class Quaternion:
// Copyright F. Nedelec, EMBL, Oct 2002   Email: nedelec@embl.de
//

#ifndef QUATERNION_H
#define QUATERNION_H

#include "types.h"
#include "smath.h"

/** Quaternions are extensions of complex number in dimension 4.
    A quaternion is =  q[0] + i * q[1] + j * q[2] + k * q[3],
    where q[?] are four real numbers. q[0] is the real part, and
    the other are imaginary.
    ATTENTION: other implementations may have different conventions!

    i, j and k are the three unit quaternions extending the real axis.
    They multiply as: i*i = -1, i*j = k and permutations of these rules.
    Quaternion-multiplication is ANTI- commutative: i*j = k and j*i = -k
    
    Unit quaternions are mostly useful to represent rotations in space:
    The group of rotation in 3D space has three degrees of freedom.
    Symmetric matrices can represent them using 6 scalar numbers, while
    unit quaternion represent them using only 4 scalar numbers.
    besides the economy, using quaternion reduces numerical drifts due
    to round-off errors. It is also easier to normalise a quaternion 
    than a matrix, to restore a rotation which has been abused by
    numerical errors, or to interpolate two rotations.
    
    The class Quaternion below implements the standard mathematical
    operations, plus conversions to and from 3x3 real matrices, 
    and to 4x4 transformation matrices used in OpenGL.
*/

/// Quaternion is a mathematical concept similar to complex number, in dimension four
class Quaternion 
{
 private:
  
  //quaternions have four coordinates
  real q[4]; //!<  q[] represents the quaternion  =  q[0] + i * q[1] + j * q[2] + k * q[3]
  
 public:
  
  /// The default constructor does not reset any value
  Quaternion() {}
  
  /// Constructor which can be used to convert from a real
  Quaternion(const real a, const real b=0, const real c=0, const real d=0) ;
  
  /// default destructor 
  virtual ~Quaternion() {}

  /// setting the values from Cartesian coordinates
  void set(const real a, const real b=0, const real c=0, const real d=0) ;  

  /// set all coordinates to zero
  void reset();

  /// access to a modifiable coordinate
  real & operator [] ( const int n );
  
  /// access to a non-modifiable coordinate 
  real operator [] ( const int n ) const;

  /// conversion operator to a "real array"
  operator real * ();

  /// Opposition: change sign in all coordinates 
  Quaternion operator - () const;
    
  /// multiply by a real value
  Quaternion operator * ( const real f ) const;
    
  /// divide by a real value
  Quaternion operator / ( const real f ) const;

  /// add a real value to your quaternion
  void  operator += ( const real f );

  /// subtract a real value to your quaternion
  void  operator -= ( const real f );

  /// multiply for a real value your cute little quaternion
  void  operator *= ( const real f );

  /// divide by a real value your hugly fat quaternion
  void  operator /= ( const real f );

  /// sum two quaternions
  Quaternion  operator + ( const Quaternion &a ) const;

  /// subtract two quaternions
  Quaternion  operator - ( const Quaternion &a ) const;
    
  /// add another quaternion to your quaternion
  void  operator += ( const Quaternion & a );

  /// subtract a quaternion to your quaternion
  void  operator -= ( const Quaternion & a );
    
  ///multiplication from the right side
  void  operator *= ( const Quaternion & a );

  /// divide your quaternion by another quaternion
  void  operator /= ( const Quaternion & a );

  /// multiplication between quaternions
  Quaternion  operator * ( const Quaternion & a ) const;

  /// division between quaternions
  Quaternion  operator / ( const Quaternion & a ) const;

  /// extract the square of the norm, i.e. norm*norm
  real normSquare() const; 

  /// extract the norm 
  real norm() const;

  /// return the normalized quaternion
  Quaternion normalized() const;

  /// normalize your quaternion
  void normalizeit();

  /// conjugated quaternion +, -, -, -
  Quaternion conjugated() const;

  /// conjugate your quaternion +, -, -, -
  void conjugateit();

  /// inversed quaternion:  1/*this
  Quaternion inverted() const;

  /// inverse your quaternion
  void inverseit();

  /// the opposed quaternion:  -*this
  Quaternion opposed() const;
  
  /// oppose your quaternion
  void opposeit();
  
  /// squared quaternion: (*this)*(*this)
  Quaternion squared() const;

  ///  square your quaternion
  void squareit();

  /// multiplication from the right side
  void rightMult( const Quaternion & a );

  ///multiplication from the left side
  void leftMult( const Quaternion & a );

  /// multiplication from the right side, different implementation
  void rightMult_fast( const Quaternion & a );
  
  ///multiplication from the left side, different implementation
  void leftMult_fast( const Quaternion & a );

  /// generate the associated 3x3 rotation matrix
  void setThisMatrix33( real m[3*3] ) const;

  /// set *this from the given rotation matrix
  void setFromMatrix33( const real m[3*3] );

  /// generate an OpenGL compatible transformation matrix
  // we assume here that type GLdouble will be the normal double
  void setThisMatrix16( double m[16], const real trans[]=0 ) const;
  //gcc 3.2 (RedHat 8) doesn't like the array dimension to be specified
  //void setThisMatrix16( double m[16], const real trans[3]=0 ) const;
  
  /// calculate the polar coordinates
  void setThisPolar( real * r, real * phi, real * theta, real * psi );
  
  /// set quaternion from polar coordinates
  void setFromPolar( const real r, const real phi, const real theta, const real psi );
  
  /// create a new quaternion with polar coordinates
  static Quaternion newFromPolar( const real r, const real phi, const real theta, const real psi);

  /// set *this to represent the rotation of axis v, and angle 'angle'
  void setFromAxisAngle( const real v[], const real angle );
  //gcc 3.2 (RedHat 8) doesn't like the array dimension to be specified
  //void setFromAxisAngle( const real v[3], const real angle );
  
  /// set *this to represent the rotation of axis X, Y or Z (0,1,2) and angle 'angle'
  void setFromMainAxisAngle( const int axis, const real angle );
  
  /// set *this to represent the rotation of axis v, of angle = v.norm();
  void setFromAxisAngle( const real v[] );
  //gcc 3.2 (RedHat 8) doesn't like the array dimension to be specified
  //void setFromAxisAngle( const real v[3] );
  
  /// compute the axis and angle of the rotation described by 'this'
  void computeAxisAngle( real v[3], real * angle ) const;

  /// get angle of the rotation described by 'this'
  real getAngle() const;
   
  /// Linear interpolation between two rotations 'a' and 'b'.
  Quaternion slerp( const Quaternion &a, const Quaternion &b, const real t );
  
  /// printf
  void print( FILE * out = stdout );
  
  /// printf with a new-line
  void println( FILE * out = stdout );
};

#endif

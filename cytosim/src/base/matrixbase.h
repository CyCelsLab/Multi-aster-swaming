//RCS: $Id: matrixbase.h,v 2.1 2005/01/10 17:27:47 foethke Exp $

#include <cstdio>
#include "smath.h"
#include "types.h"


/// Fortran-style matrices of small dimensions: 1, 2 and 3
/** 
The matrice's elements are stored in a one dimensional array to be
compatible with FORTRAN routines. 
Matrices for small dimensions: 1 2 or 3 can be used to represent
transformations of the Vecteur class */
class MATRIX
{
  
private:

  /// the array of value, column-major
  real val[ SZ*SZ ];

public:

  /// The default creator does not initialize
  MATRIX() {}
  
  /// The default destructor does nothing
  virtual ~MATRIX() {}
  
  
  /// Copy-creator from a C-style array
  MATRIX(real v[SZ][SZ]) {
    for(int i=0; i<SZ; ++i)
      for(int j=0; j<SZ; ++j)
        val[i + SZ * j] = v[i][j]; 
  }
  
  /// Copy-creator from a Fortran-style array
  MATRIX(real v[SZ*SZ]) {	
    for(int ii = 0; ii < SZ*SZ; ++ii)
      val[ ii ] = v[ ii ];
  }
  
  /// set all to zero
  void setToZero() {
    for(int ii = 0; ii < SZ*SZ; ++ii)
      val[ ii ] = 0;
  }  
  
  // set to the Identity
  void setToOne() {
    setToZero();
    for(int ii = 0; ii < SZ*SZ; ii += SZ+1)
      val[ ii ] = 1.;
  }
  
  /// conversion to a real array
  operator real * () {
    return val; 
  }

  /// access to modifiable elements of the matrix by (line, column)
  real & operator()(const int ii, const int jj) { 
    return val[ ii + SZ * jj ];
  }

  /// access to modifiable elements by index in the array
  const real & operator[](const int ii) const {
    return val[ ii ]; 
  }
  
  /// return the transposed matrix
  MATRIX transposed() const {
    MATRIX result; 
    for(int ii = 0; ii < SZ; ++ii)
      for(int jj = 0; jj < SZ; ++jj)
        result.val[ ii + SZ*jj ] = val[ jj + SZ*ii ];
    return result;
  }
  
  /// transpose this matrix
  void transposeit() {
    real x;
    for(int ii = 0; ii < SZ; ++ii)
      for(int jj = 0; jj < ii; ++jj) {
        x = val[ ii + SZ*jj ];
        val[ ii + SZ*jj ] = val[ jj + SZ*ii ];
        val[ jj + SZ*ii ] = x;
      }
  }
    
  /// determinant of the matrix
  real determinant() const;
    
  /// the inverse of the matrix
  MATRIX inverse() const;
  
  /// replace by the inverse
  void inverseit() {
    *this = inverse();
  }
  
  
  /// maximum of all fabs( elements )
  real maxNorm() const {
    real result = 0;
    for( int ii = 0; ii < SZ*SZ; ++ii ) {
      if ( fabs( val[ ii ] ) > result ) 
        result = fabs( val[ ii ] );
    }
    return result;
  }
    

    
  /// simple printf
  void print(FILE * out = stdout) const {
    for( int ii = 0; ii < SZ; ++ii ) {
      for( int jj = 0; jj < SZ; ++jj )
        fprintf(out, " %9.4f", val[ii+SZ*jj] );
      fprintf(out, "\n");
    }
  }
  
  
  /// set this to a rotation around axis of azimuth b and elevation c
  void setRotationAroundAxisEuler(const real angle, const real b, const real c);
 
  /// compute the Euler angles of the matrix
  void computeEulerAngles(real * a, real * b, real * c);

  /// multiply *this on the left by a rotation matrix
  void rotateAroundAxis(const real axis[SZ]);  

  // set *this from the azimuth a, elevation b and angle c
  void setRotationFromEulerAngles(const real a, const real b, const real c);

  /// set from a vector containing the Euler angles
  void setRotationFromEulerAngles(const real a[3]) {
    setRotationFromEulerAngles( a[0], a[1], a[2] );
  }
  
  static MATRIX randomRotation();
  
#if ( SZ == 3 )
    
  /// set this to a rotation of angle a, around one of the main axis X, Y or Z
  void setRotationAroundMainAxis(const int axis, const real angle);
  
  /// set this to a rotation around the given axis
  void setRotationAroundAxis(const real axis[SZ], const real angle);
  
#endif
  
  
  /// calculate a distance to the subspace of rotation
  real maxDeviationFromRotation() const {
    MATRIX mt = transposed();
    mt.rightMult( *this );
    for( int ii = 0; ii < SZ; ++ii )
      mt( ii, ii ) -= 1.0;
    return mt.maxNorm();
  }
  
  /// returns the matrix multiplied by a real scalar
  friend MATRIX operator * (const MATRIX & a, const real b) {
    MATRIX result;
    for(int ii = 0; ii < SZ*SZ; ++ii )
      result.val[ ii ] = a.val[ ii ] * b;
    return result;
  }

  /// returns the matrix multiplied by a real scalar
  friend MATRIX operator * (const real a, const MATRIX & b) {
    MATRIX result;
    for(int ii = 0; ii < SZ*SZ; ++ii )
      result.val[ ii ] = a * b.val[ ii ];
    return result;
  }
  
  /// multiplication by a real scalar
  void operator *=(const real a) {
    for(int ii = 0; ii < SZ*SZ; ++ii)
      val[ ii ] *= a;
  }
  
  /// division by a real scalar
  void operator /=(const real a) {
    for(int ii = 0; ii < SZ*SZ; ++ii)
      val[ ii ] /= a;
  }
  
	/// addition of another matrix
  void operator += (const MATRIX m) {	
    for(int ii = 0; ii < SZ*SZ; ++ii )
      val[ ii ] += m.val[ ii ];
  }

  /// subtraction of another matrix
  void operator -= (const MATRIX m) {	
    for(int ii = 0; ii < SZ*SZ; ++ii )
      val[ ii ] -= m.val[ ii ];
  }

  /// addition of two matrices
  friend MATRIX operator + (const MATRIX & a, const MATRIX & b) {
    MATRIX result;
    for(int ii = 0; ii < SZ*SZ; ++ii )
      result.val[ ii ] = a.val[ ii ] + b.val[ ii ];
    return result;
  }
  
  /// substraction of two matrices
  friend MATRIX operator - (const MATRIX & a, const MATRIX & b) {
    MATRIX result;
    for(int ii = 0; ii < SZ*SZ; ++ii )
      result.val[ ii ] = a.val[ ii ] - b.val[ ii ];
    return result;
  }
  
  /// multiplication of two matrices
  friend MATRIX operator * (const MATRIX & a, const MATRIX & b) { 
    MATRIX result; 
    result.setToZero();
    for(int ii = 0; ii < SZ; ++ii)
      for(int jj = 0; jj < SZ; ++jj)
        for(int kk = 0; kk < SZ; ++kk)
          result.val[ ii + SZ*jj ] += a.val[ ii + SZ*kk ] * b.val[ kk + SZ*jj ];
    
    return result;
  }  

	/// multiplication on the right by another matrix
  void rightMult(const MATRIX & b) { 
    *this = (*this) * b;
  }
  
  /// multiplication on the left by another matrix
  void leftMult(const MATRIX & b) { 
    *this = b * ( *this );
  }

  /// multiplication of a vecteur on the right
  void rightMult(real output[SZ], real input[SZ]) const {
    for( int ii = 0; ii < SZ; ++ii ) {
      output[ ii ] = 0;
      for( int jj = 0; jj < SZ; ++jj )
        output[ ii ] += val[ ii + SZ*jj ] * input[ jj ];
    }
  }
  
  /// multiplication of a vecteur on the right, overridding the input
  void rightMultOverride(real input[SZ]) const {
    static real copy[SZ];
    for( int ii = 0; ii < SZ; ++ii ) 
      copy[ ii ] = input[ ii ];
    rightMult( input, copy );
  }

};

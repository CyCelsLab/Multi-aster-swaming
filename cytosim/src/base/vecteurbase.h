//RCS: $Id: vecteurbase.h,v 2.8 2005/03/30 12:29:29 nedelec Exp $

//vecteurbase.h defined the common interface followed by vecteur1,2,3
//this is achieved by including vecteur_base.h from the header files
//For example in vecteur2.h:
//#define VECTEUR vecteur2
//#include "vecteur_base.h"
//This way, we garantee that the three classes share the same interface.

#include <cstdio>
#include "smath.h"
#include "types.h"
#include "random.h"
#include "assert_macro.h"

// we define short-cuts to access each coordinate:

#ifndef XX
  #define XX mXYZ[0]
  #define YY mXYZ[1]
  #define ZZ mXYZ[2]
#endif

/// Vecteur1 Vecteur2 and Vecteur3 are vectors with 1, 2 and 3 components.
/** VECTEUR is defined as Vecteur1, Vecteur2, Vecteur3 before compilation. 
The file vecteur_base.h is thus an interface that all these classes obey. 
*/
class VECTEUR
{

public:

  /// the coordinates of the vectors are stored in a (public) array  
  real mXYZ[3];
  
public:
  
  /// default constructor: coordinates are not initialized
  VECTEUR() { }

  /// dummy destructor does nothing
  virtual ~VECTEUR() { }

  /// constructor from individual coordinates
  inline VECTEUR(const real & xx, const real & yy, const real & zz);
  
  /// constructor from an array of coordinates
  inline VECTEUR(const real b[]);
    
  /// copy of another vecteur that's passed as an array of real
  inline void operator =(const real b[]);

  /// return the address of the coordinate vector
  const real * getXYZ()      const     { return mXYZ; }

  /// return X value
  real getX()                const     { return mXYZ[0]; }

  /// return Y value
  real getY()                const     { return mXYZ[1]; }

  /// return Z value
  real getZ()                const     { return mXYZ[2]; }

  /// implicit conversion to a modifiable real[] array
  operator real*()                     { return mXYZ; }

  /// implicit conversion to a const real[] array
  operator const real*()     const     { return mXYZ; }
  
  /// modifiable access to individual coordinates
  real & operator[](const int ii)       { assert((ii>=0)&&(ii<3)); return mXYZ[ii]; }
  
  /// constant access to individual coordinates
  real   operator[](const int ii) const { assert((ii>=0)&&(ii<3)); return mXYZ[ii]; }  

  /// set values by giving coordinates
  inline void set(const real & xx, const real & yy, const real & zz);
  
  /// set values by giving an array of coordinates
  inline void set(const real b[]);

  /// set to zero
  inline void clear();
  
  /// change signs of all coordinates
  inline void oppose();

  //------------------------------------------------------------------
  
  /// returns the square of the standard norm
  inline real normSquare()  const;
  
  /// returns the standard norm
  inline real norm()    const;
  
  /// normalize to norm=1
  inline void normalize();
  
  /// normalize to norm=n
  inline void normalizeAndScale(const real & n);

  /// normalize to norm=1/n
  inline void normalizeAndDivide(const real & n);

  /// normalize and scale vertically
  inline void normalizeAndVmult(const real & x, const real a[]);
  
  /// normalized vecteur
  inline VECTEUR normalized() const;
  //------------------------------------------------------------------

  /// set all values independantly randomly in [-1,+1]
  inline void setRandom();
  
  /// set all values independently randomly in [-s,+s]
  inline void setRandom(const real & s);
  
  /// equivalent to this += random(s);
  inline void addRandom(const real & s);
  
  /// set randomly in the sphere (norm < 1)
  void setRandomSphere();
  
  /// set randomly in the sphere (norm < s)
  void setRandomSphere(const real & s);

  /// set randomly on the sphere (norm=1)
  void setRandomNormed();
    
  /// returns a random vecteur in [-1,+1]
  inline static VECTEUR random();
  
  /// returns a random vecteur in [-s,+s]
  inline static VECTEUR random(const real & s);
  
  /// returns a random vecteur in the sphere (norm < 1 )
  static VECTEUR randSphere() {
    VECTEUR result; result.setRandomSphere(); return result; 
  }

  /// returns a random vecteur in the sphere (norm < s)
  static VECTEUR randSphere(const real & s) {
    VECTEUR result; result.setRandomSphere(s); return result; 
  }

  /// returns a random vecteur on the sphere (norm == 1)
  static VECTEUR randNormed() {
    VECTEUR result; result.setRandomNormed(); return result; 
  }
  //------------------------------------------------------------------

  /// get the Euler angles
  void computeEulerAngles(real angles[]) const;
  
  /// set from Euler angles
  void setFromEulerAngles(const real angles[]);
  
  /// multiplication coordinate by coordinates
  inline void vmult(const real b[]);
  
  /// division coordinates by coordinates
  inline void vdiv(const real b[]);

  /// linear interpolation: returns a + x * b
  inline friend VECTEUR interpolate(const VECTEUR & a, const real & x, const VECTEUR & b);
  
  /// square of the distance between two points == (a-b).normSquare(), but faster
  inline friend real distanceSquare(const VECTEUR & a, const VECTEUR & b);
  
  /// distance between two points == (a-b).norm(), but faster
  inline friend real distance(const VECTEUR & a, const VECTEUR & b);  
  
  //------------------------------------------------------------------
  
  /// addition of two vecteurs
  inline friend VECTEUR operator +(const VECTEUR & a, const VECTEUR & b);
  
  /// subtraction of two vecteurs
  inline friend VECTEUR operator -(const VECTEUR & a, const VECTEUR & b);
  
  /// unary + operator does nothing
  inline friend VECTEUR operator +(const VECTEUR & b)   { return b; }

  /// opposition of a vecteur
  inline friend VECTEUR operator -(const VECTEUR & b);
  
  
#ifdef CROSS_PRODUCT_IS_REAL
  
  // in dimension 1 and 2, the Z axis corresponds to simple reals
  // we define the cross-product with a real, i.e. to a vector along Z
  
  /// cross product of two vecteurs
  inline friend real operator ^(const VECTEUR & a, const VECTEUR & b);
  /// cross product of two vecteurs
  inline friend VECTEUR operator ^(const VECTEUR & a, const real & b);
  /// cross product of two vecteurs
  inline friend VECTEUR operator ^(const real & a, const VECTEUR & b);
  
#else
  
  // Dimension = 3 : the cross product of two vecteurs is a vecteur
  /// cross product of two vecteurs
  inline friend VECTEUR operator ^(const VECTEUR & a, const VECTEUR & b);

  inline friend VECTEUR operator ^(const VECTEUR & a, const real * b);
  
#endif
  
  /// scalar product of two vecteurs
  inline friend real operator *(const VECTEUR & a, const VECTEUR & b);
  
  /// multiplication by a real scalar
  inline friend VECTEUR operator *(const VECTEUR & a, const real & s);

  /// mutiplication by a real scalar
  inline friend VECTEUR operator *(const real & s, const VECTEUR & a);
  
  /// division by a real scalar
  inline friend VECTEUR operator /(const VECTEUR & a, const real & s);
  
  /// addition of another vecteur
  inline void operator +=(const VECTEUR & b);
  
  /// addition of another vecteur that's passed as an array of real
  inline void operator +=(const real * b);
  
  /// subtraction of another vecteur
  inline void operator -=(const VECTEUR & b);
  
  /// multiplication by a scalar
  inline void operator *=(const real & b);
  
  /// division by a scalar
  inline void operator /=(const real & b);  
  
  //------------------------------------------------------------------

  /// equality test
  inline friend int operator ==(const VECTEUR & a, const VECTEUR & b);
  
  /// different test
  inline friend int operator !=(const VECTEUR & a, const VECTEUR & b);
  //------------------------------------------------------------------
  
  
  /// print to a file
  inline void print( FILE * out = stdout );
  
  /// print to a file with a new line
  inline void println( FILE * out = stdout );

};

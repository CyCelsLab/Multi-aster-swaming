//RCS: $Id: smath.h,v 2.14 2005/04/07 20:21:58 nedelec Exp $
//These are the basic mathematical functions for cytosim.
//If the intel compiler is installed, we use the optimized intel math library
//instead of <cmath>. Never include <cmath> directly, but include this file.


#ifndef SMATH_H
#define SMATH_H

//use the optimized intel math library for the intel compiler
#ifdef __INTEL_COMPILER
    #include <mathimf.h>
//    #include "mathimf.h"
#else
    #include <cmath>
    //"cmath" is a bit weird. It undefines some c99 macros that were already
    //defined in "math.h" and redefines them in a seperate namespace,
    //called "__gcc_cxx". If you want to use functions like isfinite, isnan,
    //etc, you have to write __gcc_cxx::isfinite.
#endif

#include <cstdio>
#include <cstdlib>

#ifndef REAL
   #include "types.h"
#endif

///a famous mathematical number, troncating all the precious digits:
const real PI=3.14159265358979323846264338327950288;


//----------------------- some templated macros:

///maximum value
template<class T> T maxT(const T & a, const T & b)
{
  return (a > b) ? a : b;
}

///minimum value
template<class T> T minT(const T & a, const T & b)
{
  return (a > b) ? b : a;
}

///maximum of three arguments
template<class T> T maxT(const T & a, const T & b, const T & c)
{
  if ( a > b )
    return (a > c) ? a : c;
  else
    return (b > c) ? b : c;
}

///minimum of three arguments
template<class T> T minT(const T & a, const T & b, const T & c)
{
  if ( a > b )
    return (b < c) ? b : c;
  else
    return (a < c) ? a : c;
}

///the sign of a number
template<class T> int sign(const T a)
{
  return (a>0)?1:((a<0)?-1:0);
}

///square of a number
template<class T> T  sqr(const T a)
{
  return a * a;
}

///cube of a number
template<class T> T  cub(const T a)
{
  return a * a * a;
}

///extract a digit form a number:
//1st digit is really the first one, we do not start at zero!
template<class T> int digit(T x, int p)
{
  for(int q=1; q < p; ++q ) x /= 10;
  return x % 10;
}


//returns 1 if w[] is inside the
//ellipse of axis   (1,0,0), (0,1,0), (0,0,1)
//and of half width  d[0],    d[1],    d[2]
int insideEllipse(const real * w, const real * d);

//cross(a,b,c,d) returns 1 if the segments [a,b] and [c,d] are crossing
//setting the x and y to the abscisse of the intersection on [a,b] and [c,d]
int cross( const real a[], const real b[],
           const real c[], const real d[],
           real * x, real * y);

//this is needed under windows:
/*
#include "assert_macro.h"
inline real drem( real a, real b )
{
  assert( b > 0 );
  int p = (int)floor( 0.5 + a / b );
  if ( p )
     return a - p * b;
  else
     return a;
}
*/

//copy bytes
void copyBytes( void * dest, const void * src, const int cnt);

///print the bit values of the given adress...
void printBits(const void * a, const int size);

///check array ap of size size for nans
//int isnan_array( const real* ap, const int size );

///check array ap of size size for infs
//int isinf_array( const real* ap, const int size );

///check array ap of size size for nans and infs
//int isnotfinite_array( const real* ap, const int size );


/// Some functions usefull for output/debugging, BLAS style:

void vectPrint( int m, const real * X, FILE * f = stdout, int digits = 3 );
void vectPrint( int m, const real * X, int digits, FILE * f = stdout );
void vectDump( int m, const real * X, FILE * f = stdout, int digits = 8 );
void vectDump( int m, const real * X, int digits, FILE * f = stdout );
void matPrint( int m, int n, const real * X, FILE * f = stdout, int digits = 3 );
void matPrint( int m, int n, const real * X, int digits, FILE * f = stdout );
void matDump( int m, int n, const real * X, FILE * f = stdout, int digits = 8 );
void matDump( int m, int n, const real * X, int digits, FILE * f = stdout );
void matDumpOffset( int m, int n, const real * X, int off, FILE * f = stdout, int digits = 8 );
void matDumpOffset( int m, int n, const real * X, int off, int digits, FILE * f = stdout );

#endif //#ifdef SMATH_H

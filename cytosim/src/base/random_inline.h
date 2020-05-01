//RCS: $Id: random_inline.h,v 2.1 2005/01/10 17:27:47 foethke Exp $

// Modified by F. Nedelec, October 2002 and later for the needs of Cytosim
// 29 March 2003: added gaussian random numbers
// August 2003: added signed integer arguments to pint_()
// Sept 2003: added random integer according to given distribution
// END modications by F. nedelec

// MersenneTwister.h
// Mersenne Twister random number generator -- a C++ class MTRand
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v0.8  24 March 2002  rjwagner@writeme.com

// The Mersenne Twister is an algorithm for generating random numbers.  It
// was designed with consideration of the flaws in various other generators.
// The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
// are far greater.  The generator is also fast; it avoids multiplication and
// division, and it benefits from caches and pipelines.  For more information
// see the inventors' web page at http://www.math.keio.ac.jp/~matumoto/emt.html

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

// Copyright (C) 2002  Richard J. Wagner
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// 
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// The original code included the following notice:
//
//     Copyright (C) 1997, 1999 Makoto Matsumoto and Takuji Nishimura.
//     When you use this, send an email to: matumoto@math.keio.ac.jp
//     with an appropriate reference to your work.
//
// It would be nice to CC: rjwagner@writeme.com and Cokus@math.washington.edu
// when you write.

#ifndef RANDOM_H
#define RANDOM_H

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own MTRand object)

#include <climits>
#include <ctime>
#include <cstdio>
#include "assert_macro.h"
#include "smath.h"

class Random {
//=================================  Data  =================================
public:
	typedef unsigned int uint32;  // unsigned integer type, min 32 bits
	typedef   signed int  int32;

	
	enum { N = 624 };              // length of state vector
	enum { SAVE = N + 1 };         // length of array for save()

protected:
	enum { M = 397 };              // period parameter
	enum { MAGIC = 0x9908b0dfU };  // magic constant
	
	uint32 state[N];  // internal state
	uint32 *pNext;    // next value to get from state
	int left;         // number of values left before reload needed


//================================ Methods =================================
public:


//---------------------------------reload-----------------------------------

void reload()
{
	// Generate N new values in state
	// Made clearer and faster by Matthew Bellew (matthew.bellew@home.com)
	register uint32 *p = state;
	register int i;
	for( i = N - M; i--; ++p )
		*p = twist( p[M], p[0], p[1] );
	for( i = M; --i; ++p )
		*p = twist( p[M-N], p[0], p[1] );
	*p = twist( p[M-N], p[0], state[0] );

	left = N, pNext = state;
}


uint32 hash( time_t t, clock_t c )
{
	// Get a uint32 from t and c
	// Better than uint32(x) in case x is floating point in [0,1]
	// Based on code by Lawrence Kirby (fred@genesis.demon.co.uk)

	static uint32 differ = 0;  // guarantee time-based seeds will change

	uint32 h1 = 0;
	unsigned char *p = (unsigned char *) &t;
	for( size_t i = 0; i < sizeof(t); ++i )
	{
		h1 *= UCHAR_MAX + 2U;
		h1 += p[i];
	}
	uint32 h2 = 0;
	p = (unsigned char *) &c;
	for( size_t j = 0; j < sizeof(c); ++j )
	{
		h2 *= UCHAR_MAX + 2U;
		h2 += p[j];
	}
	return ( h1 + differ++ ) ^ h2;
}


//--------------------------creators-------------------------


Random( const uint32  oneSeed ) { seed(oneSeed); }
Random( uint32 *const bigSeed ) { seed(bigSeed); }
Random() { seed(); }



//---------------------------reals----------------------------

double preal()                        // real number in [0,1)
{ 
  return double( pint() ) * 2.3283064365386962890625e-10; // 1/2^32
}

double preal( const double& n )       // real number in [0,n)
{ 
  return preal() * n; 
}

double preal_exc()                   // real number in (0,1)
{
  return (1.0+double(pint())) * 2.328306435996595202945965527e-10;  
  // 1/( 2^32+1 )
}

double preal_exc( const double& n )  // real number in (0,n)
{ 
  return preal_exc() * n; 
}

double sreal()                       // signed real number in (-1,1)
{ 
  return double( sint() ) * 4.656612873077392578125e-10;    // 1/2^31
}

double real_range(const double low, const double high)
{
  return low + preal() * ( high - low );
}

double gauss()        // signed real number, following a normal law N(0,1)
{
  //a trick from Numerical Recipe

  static double buffer_value = 0;
  static int    buffer_flag = 0;

  if ( buffer_flag ) { 
    buffer_flag = 0;
    return buffer_value;
  } else {
    double x, y, normsq, fac;
    do {
       x = sreal();
       y = sreal();
       normsq = x * x + y * y;
    } while (( normsq >= 1.0 ) || ( normsq == 0.0 ));
    fac = sqrt( -2.0 * log( normsq ) / normsq );
    buffer_value = fac * x;
    buffer_flag  = 1;
    return fac * y;
  }
}

double gauss2()        // signed real number, following a normal law N(0,1)
{
  //const double PI = 3.14159265358979323846264338327950288;

  static double buffer_value = 0;
  static int    buffer_flag = 0;

  if ( buffer_flag ) {
    buffer_flag = 0;
    return buffer_value;
  } else {
    double angle = double( pint() ) * 1.46291807926715968105133780430979e-9;
    //the constant is 2*pi/2^32
    double norm  = sqrt( -2.0 * log( preal_exc() ) );
    buffer_value = norm * cos( angle );
    buffer_flag = 1;
    return norm * sin( angle );
  }
}

//---------------------------integers-----------------------------

uint32 pint()                        // integer in [0,2^32-1]
{
  if( left == 0 ) reload();
  --left;
  
  register uint32 s1;
  s1 = *pNext++;
  s1 ^= (s1 >> 11);
  s1 ^= (s1 <<  7) & 0x9d2c5680U;
  s1 ^= (s1 << 15) & 0xefc60000U;
  return ( s1 ^ (s1 >> 18) );
}

uint32 pint_exc( const uint32 n )   // unsigned integer in [0,n-1] for n < 2^32
{
  return uint32( preal() * n );
}

int32 pint_exc( const int32 n )     // integer in [0,n-1] for n < 2^31
{
  assert( n >= 0 );
  return int32( preal() * n );
}

uint32 pint_inc( const uint32 n )   // unsigned integer in [0,n] for n < 2^32
{
  return uint32( preal() * (n+1) );
}

int32 pint_inc( const int32 n )     // integer in [0,n] for n < 2^32
{
  assert( n >= 0 );
  return int32( preal() * (n+1) );
}

uint32 pint_inc_true( const uint32 n )    // integer in [0,n] for n < 2^32
{
  // Find which bits are used in n
  uint32 used = ~0;
  for( uint32 m = n; m; used <<= 1, m >>= 1 ) {}
  used = ~used;
  
  // Draw numbers until one is found in [0,n]
  uint32 i;
  do
    i = pint() & used;  // toss unused bits to shorten search
  while( i > n );
  return i;
}


uint32 pint_exc_true( const uint32 n )       // integer in [0,n) for n < 2^32
{
  if ( n == 0 ) return 0;
  // Find which bits are used in n
  uint32 used = ~0;
  for( uint32 m = n-1; m; used <<= 1, m >>= 1 ) {}
  used = ~used;
  
  // Draw numbers until one is found in [0,n]
  uint32 i;
  do
    i = pint() & used;  // toss unused bits to shorten search
  while( i >= n );
  return i;
}


int32 sint()                            // signed integer;
  {
    if( left == 0 ) reload();
    --left;
    
    register uint32 s1;
    s1 = *pNext++;
    s1 ^= (s1 >> 11);
    s1 ^= (s1 <<  7) & 0x9d2c5680U;
    s1 ^= (s1 << 15) & 0xefc60000U;
    return ( s1 ^ (s1 >> 18) );
  }
 

int32 sint_inc(const int32 N)          // integer in [-N, N];
{
  assert( N >= 0 );
  return pint_inc( 2 * N ) - N;
}

int32 sint_exc(const int32 N)          // integer in (-N, N);
{
  assert( N >= 0 );
  return pint_inc( 2 * ( N - 1 ) ) - N + 1;
}

bool test( const double & p ) 
{
  return ( preal() < p ); 
}

bool flip()                             // 0  or  1
{
  return ( pint() % 2 ); 
}

int32 int_range(const int32 low, const int32 high)  // in [low, high]
{ 
  return low + pint_inc( high - low );
}


uint32 pint_ratio(const uint32 n, const int ratio[]) 
     //returns an integer in [0 n[, with the ratios given in the array of ints
{
  int sum = 0;
  uint32 ii;
  for( ii = 0; ii < n; ++ii ) sum += ratio[ ii ];
  // sum==0 denotes a careless use of the function, with wrong arguments.
  // we return here 0, but in harder times, we should call an MSG.error()
  if ( sum == 0 ) return 0;
  sum = (int) floor( sum * preal() );
  ii = 0;
  while ( sum >= ratio[ ii ] )
    sum -= ratio[ ii++ ];
  return ii;
}

//-------------------------- uniform choice among the values given:

double choose(const double x, const double y)
{
  if ( flip() ) return x; else return y; 
}

double choose(const double x, const double y, const double z)
{
  switch( pint_exc(3) ) {
  case 0:  return x;
  case 1:  return y;
  default: return z;
  }
}

double choose(const double x, const double y, const double z, 
	      const double t)
{
  switch( pint_exc(4) ) {
  case 0:  return x;
  case 1:  return y;
  case 2:  return z;
  default: return t;
  }
}

double choose(const double x, const double y, const double z, 
	      const double t, const double u)
{
  switch( pint_exc(5) ) {
  case 0:  return x;
  case 1:  return y;
  case 2:  return z;
  case 3:  return t;
  default: return u;
  }
}

double choose(const double x, const double y, const double z, 
	      const double t, const double u, const double v)
{
  switch( pint_exc(6) ) {
  case 0:  return x;
  case 1:  return y;
  case 2:  return z;
  case 3:  return t;
  case 4:  return u;
  default: return v;
  }
}

//--------------------------seeding------------------------------------

uint32 seed( uint32 oneSeed )
     //included a call to the timer, if argument is zero 2002/Oct/3, FNedelec 
{
  // last attempt to get a correct seeding:
  if ( oneSeed == 0 ) oneSeed = 1;
  
  // Seed the generator with a simple uint32
  register uint32 *s;
  register int i;
  for( i = N, s = state;
       i--;
       *s = oneSeed & 0xffff0000,
	 *s++ |= ( (oneSeed *= 69069U)++ & 0xffff0000 ) >> 16,
	 (oneSeed *= 69069U)++ ) {}  // hard to read, but fast
  reload();

  return oneSeed;
}


void seed( uint32 *const bigSeed )
{
  // Seed the generator with an array of 624 uint32's
  // There are 2^19937-1 possible initial states.  This function allows
  // any one of those to be chosen by providing 19937 bits.  The lower
  // 31 bits of the first element, bigSeed[0], are discarded.  Any bits
  // above the lower 32 in each element are also discarded. Theoretically,
  // the rest of the array can contain any values except all zeroes.
  // Just call seed() if you want to get array from /dev/urandom
  register uint32 *s = state, *b = bigSeed;
  register int i = N;
  for( ; i--; *s++ = *b++ & 0xffffffff ) {}
  reload();
}


void seed()
{
  // Seed the generator with an array from /dev/urandom if available
  // Otherwise use a hash of time() and clock() values
  
  // First try getting an array from /dev/urandom
  FILE* urandom = fopen( "/dev/urandom", "rb" );
  if( urandom )
    {
      register uint32 *s = state;
      register int i = N;
      register bool success = true;
      while( success && i-- )
	{
	  success = fread( s, sizeof(uint32), 1, urandom );
	  *s++ &= 0xffffffff;  // filter in case uint32 > 32 bits
	}
      fclose(urandom);
      if( success )
	{
	  // There is a 1 in 2^19937 chance that a working urandom gave
	  // 19937 consecutive zeroes and will make the generator fail
	  // Ignore that case and continue with initialization
	  reload();
	  return;
	}
    }
  
  // Was not successful, so use time() and clock() instead
  seed( hash( time(NULL), clock() ) );
}

uint32 seedTimer()
{
  uint32 s = hash( time(NULL), clock() );
  seed( s );
  return s;
}

	
//-------------- Saving and loading generator state-------------------


void save( uint32* saveArray ) const  // to array of size SAVE
{
  register uint32 *sa = saveArray;
  register const uint32 *s = state;
  register int i = N;
  for( ; i--; *sa++ = *s++ ) {}
  *sa = left;
}


void load( uint32 *const loadArray )  // from such array
{
  register uint32 *s = state;
  register uint32 *la = loadArray;
  register int i = N;
  for( ; i--; *s++ = *la++ ) {}
  left = *la;
  pNext = &state[N-left];
}

protected:

uint32 hiBit( const uint32& u ) const { return u & 0x80000000U; }
uint32 loBit( const uint32& u ) const { return u & 0x00000001U; }
uint32 loBits( const uint32& u ) const { return u & 0x7fffffffU; }
uint32 mixBits( const uint32& u, const uint32& v ) const
{ return hiBit(u) | loBits(v); }
uint32 twist(const uint32& m, const uint32& s0, const uint32& s1) const
{ return m ^ (mixBits(s0,s1)>>1) ^ (loBit(s1) ? MAGIC : 0U); }

};


extern Random RNG;


#endif  //MERSENNETWISTER_H

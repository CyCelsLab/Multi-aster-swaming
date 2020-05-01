//RCS: $Id: random.h,v 2.16 2005/03/15 09:38:06 nedelec Exp $

// Formated Maria Mora Corral, Jan 2004
/** Modified by F. Nedelec, October 2002 and later for the needs of Cytosim
// 29 March 2003: added gaussian random numbers
// August 2003: added signed integer arguments to pint_()
// Sept 2003: added random integer according to given distribution
// END modications by F. nedelec

// MersenneTwister.h
// Mersenne Twister random number generator -- a C++ class Random
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
*/
#ifndef RANDOM_H
#define RANDOM_H

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own Random object)

#include <climits>
#include <ctime>
#include <cstdio>
#include "smath.h"
#include "types.h"
#include "assert_macro.h"

/// a random number generator class based on MersenneTwister
class Random {
//=================================  Data  =================================
public:
  
  /// type for unsigned integers
  typedef unsigned int uint32;   //!< unsigned integer type, min 32 bits
  
  ///	type for signed integers
  typedef   signed int  int32;   //!< signed integer type, min 32 bits

  /// type for big unsigned integers
  typedef unsigned long long uint64;

protected:

  /// mersenne twister internal variable
  enum { N = 624 };              //!<  length of state vector
  
  /// mersenne twister internal variable
  enum { SAVE = N + 1 };         //!<  length of array for save()
  
  /// mersenne twister internal variable
  enum { M = 397 };              //!<  period parameter
  
  /// mersenne twister internal variable
  enum { MAGIC = 0x9908b0dfU };  //!<  magic constant
  
  /// the last bit in a 32-bits integer
  enum { BIT31 = 1U << 31, BITS = 0x7FFFFFFFU, EXPON32 = 127U << 23 };

  /// exponent for a double precision float
  enum { BIT63 = 1ULL << 63, EXPON64 = 1023ULL << 52 };
  
  /// mersenne twister internal variable
  uint32 state[N];               //!<  internal state
  
  /// mersenne twister internal variable
  uint32 *pNext;                 //!<  next value to get from state
  
  /// mersenne twister internal variable
  int left;                      //!<  number of values left before reload needed


//================================ Methods =================================
  public:
    
    /// Default constructor
    Random();
    
    /// Constructor which also seeds the generator
    Random( const uint32  oneSeed );
    
    /// Constructor with seed = array of integers
    Random( uint32 *const bigSeed );


    /// Refill the random number buffer
    void reload();

    /// produce a correct seed from the time data
    uint32 hash( time_t t, clock_t c);
    
    /// unsigned integer in [0,2^32-1]
    uint32 pint() {
      if( left == 0 ) reload();
      --left;
      
      register uint32 s1;
      s1 = *pNext++;
      s1 ^= (s1 >> 11);
      s1 ^= (s1 <<  7) & 0x9d2c5680U;
      s1 ^= (s1 << 15) & 0xefc60000U;
      return ( s1 ^ (s1 >> 18) );
    }
    
    /// unsigned integer in [0,n-1] for n < 2^32
    uint32 pint_exc( const uint32 n ) {
      return uint32( preal() * n );
    }
    
    /// integer in [0,n-1] for n < 2^31
    int32 pint_exc( const int32 n ) {
      assert( n >= 0 ); //the cast would not work for negative numbers
      return int32( preal() * n );
    }
    
    /// unsigned integer in [0,n] for n < 2^32
    uint32 pint_inc( const uint32 n ) {
      return uint32( preal() * (n+1) );
    }
    
    /// integer in [0,n] for n < 2^32
    int32 pint_inc( const int32 n ) {
      assert( n >= 0 );
      return int32( preal() * (n+1) );
    }

    /// integer in [0,n] for n < 2^32, integer based algorithm using modulo
    uint32 pint_inc_true( const uint32 n );

    /// integer in [0,n[ for n < 2^32, integer based algorithm using modulo
    uint32 pint_exc_true( const uint32 n );

    /// signed integer in [-2^31+1,2^31-1]; inlined for speed
    int32 sint() {
      if ( left == 0 ) reload();
      --left;
      
      register uint32 s1;
      s1 = *pNext++;
      s1 ^= (s1 >> 11);
      s1 ^= (s1 <<  7) & 0x9d2c5680U;
      s1 ^= (s1 << 15) & 0xefc60000U;
      return ( s1 ^ (s1 >> 18) );
    }
    
    /// integer in [-N, N], boundaries included
    int32 sint_inc(const int32 n) {
      return pint_inc( 2 * n ) - n;
    }
    
    /// integer in [1-N, N-1], i.e. in ]-N,N[ with boundaries excluded
    int32 sint_exc(const int32 n) {
      return pint_inc( 2 * ( n - 1 ) ) - n + 1;
    }
    
    /// random integer in [low, high], boundaries included
    int32 int_range(const int32 low, const int32 high) {
      assert( high >= low );
      return low + int32( preal() * ( high - low + 1 ));
    }    
    
    /// returns 1 with a probability p, 0 with probability (1-p)
    bool test( const real & p ) {
      ///\todo: we could optimize by precalculating the probabilities in sim_param.cc
      static const real F = 4294967296.0; //constant factor = 2^32, biggest uint32
      //negative probabilities are not good, we should issue a warning here
      assert( p >= 0 );
      if (  p <= 0. ) return false; 
      register real tmp = p * F;          //F is a power of 2, that should be fast
      if ( tmp >= F ) return true;
      return ( pint() < uint32(tmp) );    //comparison is faster for integers than for floats
    }
    
    /// returns 1 with a probability p, 0 with probability (1-p)
    bool test_old( const real & p ) {
      return ( preal() < p ); 
    }
        
    /// returns  0  or  1 with equal chance
    bool flip() {
      return ( pint() & 1 ); 
    }
  
    /// returns -1  or  1 with equal chance
    int sflip() {
      return ((pint() & 1) ? -1 : 1); 
    }
    
    ///returns an integer in [0 n], with the ratios given in the array of ints
    uint32 pint_ratio(const uint32 n, const int ratio[]) ;
    
    /// fast (dirty) random float in [0,1[, requires IEEE Standard 754
    float pfloat_fast() {
      //This assumes IEEE Standard 754 Floating point numbers
      //32 bits: 1 for sign, 8 for exponents, 23 for fraction
      uint32 result = EXPON32 | ( pint() >> 9 );
      //we get this way a random number between 1 and 2. We substract 1.0,
      //but that drops the lower bits, reducing the precision in small numbers
      return *(float*)&result - 1.0;
    }
    
    /// random float in [0,1[, requires IEEE Standard 754 
    float pfloat();    
    
    /// random float in ]-1,1[, requires IEEE Standard 754
    float sfloat();    
    
    /// slow random double in [0,1[, using two uint32 to set all the fraction bits, requires IEEE Standard 754
    double pdouble();    
    
    /// slow random double in ]-1,1[, using two uint32 to set all the fraction bits, requires IEEE Standard 754
    double sdouble();
    
    /// inline positive real number in (0,1), boundaries included
    real preal() { 
      return real( pint() ) * 2.3283064365386962890625e-10; 
      // the constant scaling factor is 1/2^32
      //TODO: the factor could be moved in the (precalculated) probabilities being tested
    }
    
    /// inline signed real number in (-1,1), boundaries included
    real sreal() { 
      return real( sint() ) * 4.656612873077392578125e-10;
      // the constant scaling factor is 1/2^31
    }
        
    /// real number in [0,n], boundaries included
    real preal( const real& n );
    
    /// non-zero real number in ]0,1]: cannot be zero, be can be very small...
    real preal_exc(); 
    
    /// non-zero real number in ]0,n]
    real preal_exc( const real& n );
    
    /// real number uniformly between low and high
    real real_uniform_range(const real low, const real high);

    /// real number between low and high, linearly biased towards high
    real real_increasing_range(const real low, const real high);

    /// signed real number, following a normal law N(0,1)
    real gauss();
    
    /// signed real number, following a normal law N(0,1), slightly different implementation
    real gauss2(); 
    
    /// positive real with exponential law prob(u) = p * exp( -p * u )
    real exponential(const real p) {
      return -log( preal() ) / p;
    }
    
    ///uniform choice among the 2 values given:  x,y 
    real choose(const real x, const real y);

    ///uniform choice among the 3 values given:  x,y,z
    real choose(const real x, const real y, const real z);

    ///uniform choice among 4 values given
    real choose(const real x, const real y, const real z, 
                const real t);

    ///uniform choice among 5 values given
    real choose(const real x, const real y, const real z, 
                const real t, const real u);

    ///uniform choice among 6 values given
    real choose(const real x, const real y, const real z, 
                const real t, const real u, const real v);

    ///uniform choice among 7 values given
    real choose(const real x, const real y, const real z, 
                const real t, const real u, const real v,
                const real w);

    ///uniform choice among 8 values given
    real choose(const real x, const real y, const real z, 
                const real t, const real u, const real v,
                const real w, const real a);
    
    ///uniform choice among 9 values given
    real choose(const real x, const real y, const real z, 
                const real t, const real u, const real v,
                const real w, const real a, const real b);
    
    //uniform shuffling of the array T[], from knuth's The Art of Programming, Vol2 chp. 3.4.2
    template<class T> void mixWell(T val[], int size) {
      int  jj = size, kk;
      while( jj > 1 ) {
        kk = pint_exc( jj ); 
        --jj;
        T tmp = val[ jj ];
        val[ jj ] = val[ kk ];
        val[ kk ] = tmp;
      }
    }
    
    /// seed from uint32
    void seed( uint32 oneSeed );

    /// seed with a bigSeed
    void seed( uint32 *const bigSeed );

    /// seed 
    void seed();

    /// seed using time as seed-number
    uint32 seedTimer();

    /// Saving generator state to array of size SAVE
    void save( uint32* saveArray ) const;

    /// Loading generator state from such array
    void load( uint32 *const loadArray );


  protected:
    ///Mersenne Twister internal function
    uint32 hiBit( const uint32& u ) const                                    { return u & 0x80000000U; }
    ///Mersenne Twister internal function
    uint32 loBit( const uint32& u ) const                                    { return u & 0x00000001U; }
    ///Mersenne Twister internal function
    uint32 loBits( const uint32& u ) const                                   { return u & 0x7fffffffU; }
    ///Mersenne Twister internal function
    uint32 mixBits( const uint32& u, const uint32& v ) const                 { return hiBit(u) | loBits(v); }
    ///Mersenne Twister internal function
    uint32 twist(const uint32& m, const uint32& s0, const uint32& s1) const  { return m ^ (mixBits(s0,s1)>>1) ^ (loBit(s1) ? MAGIC : 0U); }

};


///RNG is a global instantiation of a Random Number Generator
extern Random RNG;


#endif  //RANDOM_H

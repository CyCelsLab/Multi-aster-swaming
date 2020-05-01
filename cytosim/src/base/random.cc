//RCS: $Id: random.cc,v 2.16 2005/02/14 19:43:34 nedelec Exp $
//==========================================================================

#include "random.h"
#include "assert_macro.h"

//allocate the default instantiation of the class
//RNG = Random Number Generator
Random RNG;

//------------------------creators-----------------------
Random::Random( const uint32  oneSeed ) 
{ 
    assert( sizeof( uint32 ) >= 4 );
    seed(oneSeed); 
}
//-------------------------------------------------------
Random::Random( uint32 *const bigSeed ) 
{ 
    assert( sizeof( uint32 ) >= 4 );
    seed(bigSeed); 
}
//-------------------------------------------------------
Random::Random() 
{
    assert( sizeof( uint32 ) >= 4 );
    seed(); 
}
//-------------------------------------------------------
void Random::reload()
{
  /** Generate N new values in state
      Made clearer and faster by Matthew Bellew (matthew.bellew@home.com) */
	register uint32 *p = state;
	register int i;
	for( i = N - M; i--; ++p )
		*p = twist( p[M], p[0], p[1] );
	for( i = M; --i; ++p )
		*p = twist( p[M-N], p[0], p[1] );
	*p = twist( p[M-N], p[0], state[0] );

	left = N, pNext = state;
}
//-------------------------------------------------------
Random::uint32 Random::hash( time_t t, clock_t c )
{
  /** Get a uint32 from t and c
      Better than uint32(x) in case x is floating point in [0,1]
      Based on code by Lawrence Kirby (fred@genesis.demon.co.uk) */

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
//-------------------------------------------------------
real Random::preal( const real& n )       // real number in [0,n)
{ 
  return preal() * n; 
}
//-------------------------------------------------------
float Random::pfloat()
{
  //This assumes IEEE Standard 754 Floating point numbers
  //32 bits: 1 for sign, 8 for exponents, 23 for fraction
  register uint32 result = pint();
  register uint32 E = 126;
  while (( result < BIT31 ) && ( E > 94 )) {
    result <<= 1; --E;
  }
  result = ( result & 0x7FFFFFFFU ) >> 8 | ( E << 23 );
  return *(float*)&result;
}
//-------------------------------------------------------
float Random::sfloat()
{
  //This assumes IEEE Standard 754 Floating point numbers
  //32 bits: 1 for sign, 8 for exponents, 23 for fraction
  register uint32 result = pint();
  register uint32 S = result & BIT31;
  result <<= 1;
  uint32 E = 126;
  while (( result < BIT31 ) && ( E > 94 )) {
    result <<= 1; --E;
  }
  result = S | (( result & 0x7FFFFFFFU ) >> 8 ) | ( E << 23 );
  return *(float*)&result;
}
//-------------------------------------------------------
double Random::pdouble()
{
  //This assumes IEEE Standard 754 Floating point numbers
  //64 bits: 1 for sign, 11 for exponents, 52 for Fraction
  //set all the 64 bits to random, using two random uint32:
  //register uint64 result = ( uint64(pint()) << 32 ) | pint();
  register uint64 result = pint();   *(uint32*)&result = pint();
  register uint32 E = 1022;
  while (( result < BIT63 ) && ( E > 959 )) {
    result <<= 1; --E;
  }
  result >>= 11;
  *(uint32*)&result = ( *(uint32*)&result & 0x000FFFFFU ) | ( E << 20 );
  return *(double*)&result;
}
//-------------------------------------------------------
double Random::sdouble()
{
  //This assumes IEEE Standard 754 Floating point numbers
  //64 bits: 1 for sign, 11 for exponents, 52 for Fraction
  //set all the 64 bits to random, using two random uint32:
  register uint32 E = pint();
  register uint32 S = E & BIT31;         //sign-bit
  register uint64 result = E << 1;       //stack other bits in the fraction
  *(uint32*)&result = pint();            //fill the upper part
  E = 1022;
  while (( result < BIT63 ) && ( E > 959 )) {
    result <<= 1; --E;
  }
  result >>= 11;
  *(uint32*)&result = S | ( *(uint32*)&result & 0x000FFFFFU ) | ( E << 20 );
  return *(double*)&result;
}
//-------------------------------------------------------
real Random::preal_exc()                   // real number in )0,1)
{
  register real x = preal();
  while( x == 0 ) x = preal();
  return x;
}
//-------------------------------------------------------
real Random::preal_exc( const real& n )  // real number in )0,n)
{ 
  return preal_exc() * n; 
}
//-------------------------------------------------------
real Random::real_uniform_range(const real low, const real high)
{
  return low + preal() * ( high - low );
}
//-------------------------------------------------------
real Random::gauss()        
{
  /** signed real number, following a normal law N(0,1)
      using a trick from Numerical Recipe */

  static real bufferValue = 0;
  static bool bufferFlag  = false;

  if ( bufferFlag ) { 
    bufferFlag = false;
    return bufferValue;
  } else {
    register real x, y, normsq, fac;
    do {
       x = sreal();
       y = sreal();
       normsq = x * x + y * y;
    } while (( normsq >= 1.0 ) || ( normsq == 0 ));
    fac = sqrt( -2 * log( normsq ) / normsq );
    bufferValue = fac * x;
    bufferFlag  = true;
    return fac * y;
  }
}
//-------------------------------------------------------
real Random::gauss2()        // signed real number, following a normal law N(0,1)
{
  //this version uses cos() and sin() and might be a bit slower.
  //const real PI = 3.14159265358979323846264338327950288;

  static real bufferValue = 0;
  static bool bufferFlag  = false;

  if ( bufferFlag ) {
    bufferFlag = false;
    return bufferValue;
  } else {
    real angle = real( pint() ) * 1.46291807926715968105133780430979e-9;
    //the constant is 2*pi/2^32
    real norm  = sqrt( -2 * log( preal_exc() ));
    bufferValue = norm * cos( angle );
    bufferFlag  = true;
    return norm * sin( angle );
  }
}
//-------------------------------------------------------
Random::uint32 Random::pint_inc_true( const uint32 n )    // integer in [0,n] for n < 2^32
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
//-------------------------------------------------------
Random::uint32 Random::pint_exc_true( const uint32 n )       // integer in [0,n) for n < 2^32
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
//-------------------------------------------------------
Random::uint32 Random::pint_ratio(const uint32 n, const int ratio[]) 
     //returns an integer in [0 n], with the ratios given in the array of ints
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
//-------------------------------------------------------
real Random::choose(const real x, const real y)
{
  if ( flip() )
      return x; 
  else
      return y; 
}
//--------------------------------------------------------
real Random::choose(const real x, const real y, const real z)
{
  switch( pint_exc(3) ) {
  case 0:  return x;
  case 1:  return y;
  default: return z;
  }
}
//--------------------------------------------------------
real Random::choose(const real x, const real y, const real z, 
                    const real t)
{
  switch( pint_exc(4) ) {
  case 0:  return x;
  case 1:  return y;
  case 2:  return z;
  default: return t;
  }
}
//--------------------------------------------------------
real Random::choose(const real x, const real y, const real z, 
                    const real t, const real u)
{
  switch( pint_exc(5) ) {
  case 0:  return x;
  case 1:  return y;
  case 2:  return z;
  case 3:  return t;
  default: return u;
  }
}
//--------------------------------------------------------
real Random::choose(const real x, const real y, const real z, 
                    const real t, const real u, const real v)
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
//--------------------------------------------------------
real Random::choose(const real x, const real y, const real z, 
                    const real t, const real u, const real v,
                    const real w)
{
  switch( pint_exc(7) ) {
    case 0:  return x;
    case 1:  return y;
    case 2:  return z;
    case 3:  return t;
    case 4:  return u;
    case 5:  return v;
    default: return w;
  }
}
//--------------------------------------------------------
real Random::choose(const real x, const real y, const real z, 
                    const real t, const real u, const real v,
                    const real w, const real a)
{
  switch( pint_exc(8) ) {
    case 0:  return x;
    case 1:  return y;
    case 2:  return z;
    case 3:  return t;
    case 4:  return u;
    case 5:  return v;
    case 6:  return w;
    default: return a;
  }
}
//--------------------------------------------------------
real Random::choose(const real x, const real y, const real z, 
                    const real t, const real u, const real v,
                    const real w, const real a, const real b)
{
  switch( pint_exc(9) ) {
    case 0:  return x;
    case 1:  return y;
    case 2:  return z;
    case 3:  return t;
    case 4:  return u;
    case 5:  return v;
    case 6:  return w;
    case 7:  return a;
    default: return b;
  }
}
//--------------------------------------------------------
void Random::seed( uint32 oneSeed )
{
  // last attempt to get a correct seeding:
  if ( oneSeed == 0 ) oneSeed = 1;
  
  // Seed the generator with a simple uint32
  register uint32 *s;
  register int i;
  for( i = N, s = state;
       i--;
       *s = oneSeed & 0xffff0000,
	 *s++ |= ((oneSeed *= 69069U)++ & 0xffff0000 ) >> 16,
	 (oneSeed *= 69069U)++ ) {}  // hard to read, but fast
  reload();
}
//--------------------------------------------------------
void Random::seed( uint32 *const bigSeed )
{
  /** Seed the generator with an array of 624 uint32's
      There are 2^19937-1 possible initial states.  This function allows
      any one of those to be chosen by providing 19937 bits.  The lower
      31 bits of the first element, bigSeed[0], are discarded.  Any bits
      above the lower 32 in each element are also discarded. Theoretically,
      the rest of the array can contain any values except all zeroes.
      Just call seed() if you want to get array from /dev/urandom
  */  
  register uint32 *s = state, *b = bigSeed;
  register int i = N;
  for( ; i--; *s++ = *b++ & 0xffffffff ) {}
  reload();
}
//--------------------------------------------------------
void Random::seed()
{
  /** Seed the generator with an array from /dev/urandom if available
      Otherwise use a hash of time() and clock() values
  */
  
  // First try getting an array from /dev/urandom
  FILE* urandom = fopen( "/dev/urandom", "rb" );
  if ( urandom ) {
    register uint32 *s = state;
    register int i = N;
    register bool success = true;
    while( success && i-- ) {
      success = fread( s, sizeof(uint32), 1, urandom );
      *s++ &= 0xffffffff;  // filter in case uint32 > 32 bits
    }
    fclose(urandom);
    if ( success ) {
      // There is a 1 in 2^19937 chance that a working urandom gave
      // 19937 consecutive zeroes and will make the generator fail
      // Ignore that case and continue with initialization
      reload();
      return;
    }
  }
  
  // Was not successful, so use time() and clock() instead
  seed( hash( time(NULL), clock() ));
}
//--------------------------------------------------------
Random::uint32  Random::seedTimer()
{
  uint32 s = hash( time(NULL), clock() );
  seed( s );
  return s ; 
}
//--------------------------------------------------------
void Random::save( uint32* saveArray ) const  // to array of size SAVE
{
  register uint32 *sa = saveArray;
  register const uint32 *s = state;
  register int i = N;
  for( ; i--; *sa++ = *s++ ) {}
  *sa = left;
}
//--------------------------------------------------------
void Random::load( uint32 *const loadArray )  // from such array
{
  register uint32 *s = state;
  register uint32 *la = loadArray;
  register int i = N;
  for( ; i--; *s++ = *la++ ) {}
  left = *la;
  pNext = &state[N-left];
}




//RCS: $Id: testrandom.cc,v 2.7 2005/02/18 18:25:18 nedelec Exp $
//==========================================================================
#include "random.h"
#include <cstring>

typedef Random::uint32 uint32;
typedef Random::int32   int32;
typedef Random::uint64 uint64;

const int CNT = 100;


float Convert1(uint32 x)
{
  return ( float((int32)x) * 4.6566128730773926e-10 );
}

float Convert2(uint32 x)
{
  return ( float(x) * 2.3283064365386963e-10 );
}


void test_int()
{
  int j;

  for(j=0; j<21; j++)	
    printf(" %10u%s", RNG.pint(), (j%7)==6 ? "\n" : "");

  printf("\n");

  for(j=0; j<90; j++)	
    printf(" %2u%s", RNG.pint_inc(99), (j%30)==29 ? "\n" : "");

  printf("\n");

  for(j=0; j<100; j++)	
    printf(" %3i%s", RNG.sint_inc(99), (j%20)==19 ? "\n" : "");

  printf("\n");

  for(j=0; j<42; j++)
    printf(" %10.7f%s", RNG.sreal(), (j%7)==6 ? "\n" : "");

  printf("\n");

  for(j=0; j<42; j++)
    printf(" %8f%s", RNG.preal(), (j%7)==6 ? "\n" : "");

  printf("\n");
}


float convertFix(uint32 x)
{
  //This assumes IEEE Standard 754 Floating point numbers
  //32 bits: 1 for sign, 8 for exponents, 23 for fraction
  const uint32 BIT31    = 0x80000000U;
  const uint32 EXPON    = 127 << 23;
  uint32 result = EXPON | ( x >> 9 );
  return *(float*)&result - 1.0;
}



void testbits()
{
  const int SCALE=2;
  float x;
  for( int ii=0; ii <= SCALE; ++ii) {
    x = ii / float(SCALE);
    printf(" %f :", x);
    printBits(&x, 4);
   // x = -ii / float(SCALE);
   // printf("%f :", x);
   // printBits(&x, 4);
  }
  
  double y;
  for( int ii=0; ii <= 20; ++ii ) {
    y = convertFix( RNG.pint() );
    printf(" %f :", y);
    printBits(&y,8);
  }
}

void test_test( real prob )
{
  const int MAX = 50000000;
  int cnt = 0;
  for(int jj=0; jj < MAX; ++jj)
    cnt += RNG.test( prob );
  printf("prob = %f measured = %f cnt = %i\n", prob, cnt / double(MAX), cnt );
}  

void test_real()
{
  double x;
  for(int kk=0; kk < 10; ++kk) {
    x = RNG.preal();
    printf(" %+f", x);
  }
  printf("\n");
  for(int kk=0; kk < 10; ++kk) {
    x = RNG.sreal();
    printf(" %+f", x);
  }
  printf("\n");  
}


int main(int argc, char * argv[])
{
  double xx = 2.3283064365386962890625e-10;
  printBits(&xx,8);
  xx = 4.656612873077392578125e-10;
  printBits(&xx,8);

  if ( argc == 1 ) {
    for(int kk=0; kk < 22; ++kk)
      test_test( kk / 20.0 );
  } else {
    for(int kk=0; kk < 22; ++kk)
      test_real();
  }
      
  
  printf("done\n");
  return EXIT_SUCCESS;
}


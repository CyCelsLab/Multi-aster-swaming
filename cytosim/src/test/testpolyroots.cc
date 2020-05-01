//
// C++ Implementation: testpolyroots
//
// Description: 
//
//
// Author: Dietrich Foethke <foethke@embl.de>, (C) 2005
//
// Copyright: See COPYING file that comes with this distribution
//
// CVS: $Id: testpolyroots.cc,v 1.3 2005/04/08 10:42:02 foethke Exp $


#include "types.h"
#include "poly_roots.h"
#include <stdio.h>


int main(int argc, char *argv[])
{
  using namespace PolyRoots;
  
  double   sol[4];
  complex  zSol[4];
  int num, zNum;
  double   a, b, c, d, e;
    
/*  num  =  solveCubic( 2, -4, -22, 24, sol );
  zNum = zSolveCubic( 2, -4, -22, 24, zSol, zSol+1, zSol+2 );
  printf("Solutions: real: %d, complex: %d\n", num, zNum);
  if( num == 3 )
    printf( "sol1: %f, sol2: %f, sol3: %f\n", sol[0], sol[1], sol[2] );
  printf( "zSol1: %f+%f i, zSol2: %f+%f i, zSol3: %f+%f i\n", zSol[0].dat[0], zSol[0].dat[1], zSol[1].dat[0], zSol[1].dat[1], zSol[2].dat[0], zSol[2].dat[1] );*/
  

/*  a = 1;
  b = 0;
  c = 0;
  d = 0;
  e = 1;*/
/*  a = -20;
  b = 5;
  c = 17;
  d = -29;
  e = 87;*/
/*  a = 3;
  b = 6;
  c = -123;
  d = -126;
  e = 1080;*/
  a = 1.;
  b = 37.5;
  c = 500.;
  d = 2812.5;
  e = 5615.23;
    
  zNum = zSolveQuartic( a, b, c, d, e, zSol );
  printf("number of complex solutions: %d\n", zNum);
  for( int ii = 0; ii < zNum; ii++ ) {
    printf( "sol%d: %e+%ei\n", ii, zSol[ii].dat[0], zSol[ii].dat[1] );
  }
  
  printf("\n\n");
    
  num = solveQuartic( a, b, c, d, e, sol );
  printf("number of real solutions: %d\n", num);
  for( int ii = 0; ii < num; ii++ ) {
    printf( "sol%d: %e\n", ii, sol[ii]);
  }
}

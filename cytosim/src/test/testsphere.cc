//RCS: $Id: testsphere.cc,v 2.0 2004/08/16 16:44:59 nedelec Exp $
//----------------------------- testsphere.cc ---------------------------
// Test class PointsOnSphere

#include "pointsonsphere.h"

const int mx=6;

int main() {

  PointsOnSphere S( 1 );

  int sum_iterations = 0;
  
  for( int nbp = 1; nbp < 50; ++nbp ) {
    
    fprintf(stderr, "%4i pts :", nbp);
    for(int m=0; m < mx; ++m) {
    
      S.distributePoints( nbp );
      fprintf(stderr, " %6.4f", S.minimumDistance());
      sum_iterations += S.nbIterationsToConverge();
    
    }
    fprintf(stderr, " energy %7.2f", S.finalEnergy());
    fprintf(stderr, " iter. %7i\n", sum_iterations);
  }
    
  return 0;
}

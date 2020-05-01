//RCS: $Id: testmap.cc,v 2.3 2005/01/10 17:27:47 foethke Exp $
//-----------------------------tesmap.cc-------------------------
#include "map.h"
#include "random.h"

#include <cstdio>
#include "smath.h"

void trivialTest()
{    
    real min[] = {0,0,0};
    real max[] = {3,3,3};
    int size[] = {3,3,3};
    
    Map<1, real> map1;
    map1.setDimensions( min, max, size );
    map1.allocateEdgeTypeArray();
    map1.allocateRegionArray();
    map1.dump();
    
    printf("\n");
    
    Map<2, real> map2;
    map2.setDimensions( min, max, size );
    map2.allocateEdgeTypeArray();
    map2.allocateRegionArray();
    map2.dump();
    
    printf("\n");

    map2.deleteAll();
    map2.setDimensions( min, max, size );
    map2.setPeriodic();
    map2.allocateEdgeTypeArray();
    map2.allocateRegionArray();
    map2.dump();
    
    printf("\n");
    
    Map<3, real> map3;
    map3.setDimensions( min, max, size );
    map3.allocateCellArray();
    map3.allocateEdgeTypeArray();
    map3.allocateRegionArray();
    map3.dump();
}

void realtest()
{
    
    printf("Real test...");

    real min[] = {-10,-10};
    real max[] = { 10, 10};
    int size[] = { 10, 10};
    
    Map<2, float> map;
    map.setDimensions(min, max, size);
    map.setPeriodic();
    map.allocateCellArray();
    
    real w[2];

    map.clear();

    for(int cc=0; cc<120*120; ++cc) {
        w[0] = 12*RNG.sreal();
        w[1] = 12*RNG.sreal();
        ++map( w );
    }

    FILE * test = fopen("testmap.out","w");
    map.write(test, 0);
    fclose(test);
    printf("wrote file testmap.out\n");
}


int main()
{
  trivialTest();
  realtest();
  return 0;
}

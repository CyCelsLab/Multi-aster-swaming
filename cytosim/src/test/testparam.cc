//RCS: $Id: testparam.cc,v 2.1 2005/03/13 22:42:34 nedelec Exp $
//-------------------------testparam.cc------------------------

#include "testparam.h"
#include <cstdlib>

ParamTest test, copy;

int main()
{
  test.printDescriptions();
  test.parseFile("data.good");
  /*
  test.parseFile("data.bad");
  
  test.parseLine("along=2");
  test.parseLine("along*=2.0");
  test.parseLine("areal*=2.0");
  test.parseLine("array[10]=2.0");
  test.parseLine("array[3=2.0");
  test.parseLine("array[3]*=1.5");
  test.parseLine("array[1]*=2:3");
  test.parseLine("array[1]*=2.1:3.2:3");
  test.parseLine("array[1]*=2:3");
*/
  test.printFile(stdout);
  copy = test;
  copy.printFile(stdout);
  //printf("\n"); test.printFile(stdout, 1);
  return EXIT_SUCCESS;
}

//RCS: $Id: testparam.h,v 2.1 2004/12/06 19:46:02 nedelec Exp $
#include "types.h"
#include "parameter_list.h"

class ParamTest : public ParameterList
{
 public:

  char          string[256];
  char          string2[256];
  int           aint;
  int           array[6];
  long          along;
  unsigned long aulong;

  float         afloat;
  double        adouble;
  real          areal;

  ParamTest() {
    link( "aint",        &aint,        1,  "1");
    link( "string",      string,       256, "tata");
    link( "string2",     string2,      256, "toto");
    link( "array",       array,        6,  "2");
    link_deprecated( "old",         array);
    link( "along",       &along,       1,  "1");
    link( "aulong",      &aulong,      1,  "1");
    link( "afloat",      &afloat,      1,  "1");
    link( "adouble",     &adouble,     1,  "1");
    link( "areal",       &areal,       1,  "1");
    }

    virtual ~ParamTest() { }
};

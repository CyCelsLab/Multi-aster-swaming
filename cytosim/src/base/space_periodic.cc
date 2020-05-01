//RCS $Id: space_periodic.cc,v 2.0 2004/08/16 16:38:31 nedelec Exp $
#include "space_periodic.h"

#if (DIM == 1)

real SpacePeriodic::volume() const
{
  return 2 * mSize[0];
}

#endif


//--------------------------------------------------------------------

#if (DIM == 2)

real SpacePeriodic::volume() const
{
  return 4 * mSize[0] * mSize[1];
}

#endif

//--------------------------------------------------------------------

#if (DIM == 3)

real SpacePeriodic::volume() const
{
  return 8 * mSize[0] * mSize[1] * mSize[2];
}


#endif

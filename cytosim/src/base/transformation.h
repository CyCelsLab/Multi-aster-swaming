//RCS: $Id: transformation.h,v 2.0 2004/08/16 16:39:21 nedelec Exp $
//a Transformation is just a matrix of dimension DIM x DIM
//which can be used to store rotations / symmetries, etc.

#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "main.h"

#ifndef DIM
  #error DIM is not defined
#endif


#if (DIM==1)
   #include "matrix1.h"
   typedef Matrix11 Transformation;
#endif

#if (DIM==2)
   #include "matrix2.h"
   typedef Matrix22 Transformation;
#endif

#if (DIM==3)
   #include "matrix3.h"
   typedef Matrix33 Transformation;
#endif


#endif

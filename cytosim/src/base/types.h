//RCS: $Id: types.h,v 2.14 2005/04/18 10:23:58 nedelec Exp $
//SYNOPSIS: we define and use "real" to be able to easily change
//the floating point Precision, which should normally be double
//EPSILON is a lower limit of the precision achieved


#ifndef TYPES_H
#define TYPES_H

#include <limits>

//To use float, the keyword SINGLE should be defined
//it is STRONGLY ADVISED to use double precision: disable the line below
//#define SINGLE

#ifdef SINGLE
    typedef float real;           ///<real is used as a synonym of double or float
#else
    typedef double real;          ///<real is used as a synonym of double or float
#endif

///<EPSILON is used because of round-off error
const real EPSILON = 10. * std::numeric_limits<real>::epsilon();


///functions returning an error code should return NO_ERROR when things are normal.
/**NO_ERROR should be ZERO, error codes should be non-zero integers.*/
const int NO_ERROR = 0;

///some string constant: default names for files
const char   DATA_OUT[] = "data.out";

///this is the string that identifies the beginning of a frame in a [result.out] file
const char   FRAME_START_TAG[] = "#frm ";

///this is the string that identifies the end of a frame in a [result.out] file
const char   FRAME_END_TAG[] = "#end ";

/// positive long integer type for the names of the objects
typedef unsigned long Name;

/// bindingKey is a type used to prevent internal links to microtubules
/** keys are compared to see if a binding of a hand to a MT rod is allowed or not */
typedef unsigned long BindingKey;

/// type used to count iterations of the simulation engine
typedef unsigned long Step;

///the standard size for allocation of character arrays (should be at least 64)
const unsigned int STRING_SIZE = 512;
  
#endif

//RCS: $Id: assert_macro.h,v 2.5 2005/01/09 00:11:18 nedelec Exp $
// We use the standard assert macro for debugging.
// assert( x ) stops the program if x is false

// You can disable/enable assertions everywhere by choosing
// one of the two lines below:

#ifndef ASSERT_MACRO_H
#define ASSERT_MACRO_H

// disable the line below the get the assertions
// enable DUMMY_ASSERT for faster executable / to compile on the cluster:
#define DUMMY_ASSERT

#ifdef DUMMY_ASSERT
    #define assert(ignore) ((void) 0)
// to disable assertions for the intel compiler, version <= 8.0,
// define NDEBUG. Redefining the assert macro doesn't work for 
// older versions of icc. In version 8.1 it seems to be working.
    #define NDEBUG
#else
    #include <cassert>
    #define ASSERT
#endif


#endif

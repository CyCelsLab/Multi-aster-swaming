//RCS: $Id: main.h,v 2.18 2005/04/18 10:22:29 nedelec Exp $
//This contains some global compilation options for Cytosim

#ifndef MAIN_H
#define MAIN_H

#include "types.h"
#include "assert_macro.h"
#include "smath.h"

//========================CHANGE BELOW ONLY=============================
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv


//---------------------------Number of dimensions in space: 1, 2 or 
#define DIM 2

//---------------------------Use pseudo yz coordinates in 1D. Used only for
//---------------------------display purposes for now. Only makes sense in 1D
#if (DIM==1)
  #define PSEUDO_YZ
#endif

//---------------------------Hydrodynamic correction for microtubule's mobility:
// log( hydrodynamic cut-off length / diameter of the Microtubules )

const real HYDRO = log( 2.0 / 0.025 );

//---------------------------Approximate free motors as uniform (OFF)
//#define UNIFORM_FREE_COMPLEX

//---------------------------Hand's attachements depend on the distance (OFF)
//#define NEAR_ATTACH_FASTER

//---------------------------Slow algorithms for finding close rods (OFF)
//#define SLOW_SAFE_ATTACH

//---------------------------Include a test for the Hand's attachements
//#define TEST_ATTACH

//---------------------------The cutting of Microtub is done according to curvature
//#define CUT_WITH_CURVATURE

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//=========================CHANGE ABOVE ONLY============================
//----------------------------------------------------------------------
//------------------------- GLOBAL VARIABLES ---------------------------
//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

// set finish_now to 1 to save and exit gracefully from sim (ignored in play)
extern bool finish_now;

// For the creation of initial state, we might get stuck in a loop, for
// example if we try to fit a 10 um Microtub inside a 5um box. Normally,
// we exit from such loops after MAX_TRIALS_BEFORE_STUCK unsuccesful trials.
const unsigned long MAX_TRIALS_BEFORE_STUCK = 100000;    

// Five debug switches can be interactively changed in play from the keyboard
// this can be used during development phases, to test alternative implementations
extern int user_control_keys[5];

// Current format version id. for writing [result.out] file
const int fileFormatVersion = 20;


#endif //MAIN_H

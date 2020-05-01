//RCS: $Id: microtub.h,v 2.34 2005/04/22 18:55:26 nedelec Exp $
// ----------------------------------microtub.h---------------------------
// ----------------------------FLEXIBLE MICROTUBULES----------------------

#ifndef MICROTUB_H
#define MICROTUB_H

#include "array.h"
#include "vecteur.h"
#include "object_list.h"
#include "space.h"
#include "fiber.h"
#include "sim_param.h"
#include "lattice.h"



///Enum to specify the MT end state
enum MTDynamicState {
  MT_GROWING    = 0,   ///<  Growing state
  MT_SHRINKING  = 1,   ///<  Shrinking state
  MT_PAUSE      = 2    ///<  Fix length
};

///conversion factor from dimers to Length: 8 nm per dimer, 13 protofilaments
const real tubulinLength = 8.0e-3 / 13.0;

//--------------------- some incomplete declarations:
class Hand;
class Grafted;
class Aster;
class PointMicrotub;
class Matrix;
class MicrotubOrganizer;

///Microtub are flexible Fiber with dynamic growth/shrink properties
class Microtub : public Fiber
{
  friend class   PointInterpolated;

 public:

  ///concentration of monomers, normalized in [0,1] by dividing by MP.mtmonlength,
  static real    tubulin;

  ///space in which microtubules live
  static Space * space;

  ///color of filament, used for display only
  int            mtColor;

 private:
  ///the MicrotubOrganizer to which the MT belongs: Aster / Nucleus
  MicrotubOrganizer * mtOrganizer;

  ///dynamic state of ends
  MTDynamicState mtState[3];

  ///record interation number at the time of last dynamic state change
  Step           mtStateChange[3];

  ///record position of last dynamic state change
  Vecteur        mtStateChangeWhere[3];

  ///list of attached hands
  NodeList       mtHands;

  ///a grafted used to implement boxglue: high friction on the box wall
  Grafted *      mtboxglue_grafted;

  //--------------------- methods:

  ///resets the pointers to arrays, should be called only by the constructors
  void        microtubConstructorBase();

  ///resets the values of member variables
  void        reset();

 public:

  //---------------------Constructor/Destructor

  ///Constructor - base
  Microtub()    { microtubConstructorBase(); }

  ///Constructor - creates a microtubule according to MP.mtinit
  Microtub(int dummy);

  ///Constructor - creates a new microtubule with a center of gravity at position w
  Microtub(const Vecteur w);

  ///Constructor - creates a microtubule with given specifications (position w, ...)
  Microtub(const Vecteur w, const Vecteur dir, const real len, const MTEnd which = MT_MINUS_END);

  ///Destructor
  virtual ~Microtub();

  //---------------------

  ///we override the next() derived from Node to fix the type
  Microtub *         next()    const          { return static_cast<Microtub*>( son ); }
  ///we override the prev() derived from Node to fix the type
  Microtub *         prev()    const          { return static_cast<Microtub*>( dad ); }


  //---------------------

  ///get microtubule Organizer, i.e object to which MT belongs
  const MicrotubOrganizer * getOrganizer()             const  { return mtOrganizer; }

  ///set microtubule organizer
  void        setOrganizer(MicrotubOrganizer * a)             { mtOrganizer = a; }

  //---------------------

  ///get MT_GROWING/MT_SHRINKING state of one end
  MTDynamicState  getDynamicState(const MTEnd which)   const  { return mtState[which]; }

  ///set state of one end
  void        setDynamicState(const MTEnd which, const MTDynamicState new_state);


  ///Returns the first Hand attached to the microtubule
  Hand   *    firstHand()                              const  { return (Hand*)(mtHands.getFront()); }

  ///Link the Hand to this microtubule: it is now bound to it
  void        linkHand( Hand * ha );

  ///update the PointInterpolated for all attached Hand
  void        updateHands();

  //---------------------
  ///calculates the growth rate, in units of um/s of added microtubule length = speed.
  real        growthRate(real givenGrowthRate, real forceOnEnd, real currentMTLength) const;

  ///calculate the transition rate of a MT end, as a function of the rate in data.in, the growth and the position
  real        transitionRate(const MTDynamicState currentState, const real & givenTransitionRate, const real growthRate, const Vecteur & positionOfTip) const;

  ///change microtubule plus-end for Dynamic instability: grow/shrink and transitions
  void        stepPlusEndForDynamic();

  ///change microtubule minus-end for Dynamic instability: grow/shrink and transitions
  void        stepMinusEndForDynamic();

  //---------------------

  ///monte-carlo step
  void        step();

  //---------------------
  ///true if confinement may apply
  bool        isConfined()                            const   { return MP.mtconfine && space->isConfined(); }

  ///true if point pp is inside the space
  bool        insidePoint(const int pp)               const   { return space->isInside( whereP(pp) ); }

  ///project the point pp on the nearest space boundary
  void        projectPoint(const int pp, real r[DIM]) const   { space->project( whereP(pp), r ); }

  ///number of points inside the space
  int         insidePosition(const Space * s = space) const   { return PointSet::insidePosition(s); }

  ///a pointer to the space where the mt is living
  const Space* getSpace()                             const   { return space; }

  /// set the box glue
  void        setBoxGlue();
  //---------------------

  ///Direction for MT when setting-up initial configuration
  static Vecteur    initDirection(int mode = digit(MP.mtinit, 2));

  ///Length for MT when setting-up the initial configuration
  static real       initLength(int mode = digit(MP.mtinit, 3));

  //lattice used to say if a site is occupied
  Lattice           lattice;

  //int lattice site from the real value in ab ((int) ab/stepsize)
  int               getLatticeSite(real ab) const;

  ///Tests if lattice site at ab/MP.mtdimersize is free
  bool        siteFree(real ab) const;

  ///Tests if lattice site at nSite is free
  bool        siteFree(int nSite) const { return lattice.isFree( nSite ); }

  ///Returns density of hands at end, averaging over len
  real        getHandDensityAtEnd(MTEnd end, real len) const;

  ///returns the iteration counter of the last mt state change
  Step        getLastStateChange( const MTEnd which )    const { return mtStateChange[ which ]; }

  ///returns the position of the last mt state change
  Vecteur     getPosLastStateChange( const MTEnd which ) const { return mtStateChangeWhere[ which ]; }

  ///Sets a lattice site to free
  void        setSiteFree(real ab);

  ///Sets a lattice site to free
  void        setSiteOccupied(real ab);

  //---------------------

#ifdef PSEUDO_YZ
private:
  ///sets the fake y coordinate for 1D simulations
  ///called only in microtubConstructorBase()
  void        setPseudoYZ();
  ///array of 0/1 to specify if a site is occupied or not
  static Array<int> pseudoOccupancyArray;

public:
  ///pseudo coordinates for use in 1D
  real        pseudoY;
  real        pseudoZ;
  /// index in pseudoOccupancyArray
  int         pseudoIndex;
#endif

  //---------------------I/O

  ///write 'this' microtubule to IO
  void        write();

  ///read 'this' microtub from IO
  void        read();

  ///basic consistency check (debug)
  int        looksWrong() const;
};


#endif

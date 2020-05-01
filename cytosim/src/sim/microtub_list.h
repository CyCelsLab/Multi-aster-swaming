//RCS: $Id: microtub_list.h,v 2.7 2005/04/19 18:18:15 nedelec Exp $
//---------------------------------microtublist.h-----------------------------

#ifndef MICROTUB_LIST_H
#define MICROTUB_LIST_H


#include "object_list.h"
#include "microtub.h"
#include "map.h"

///a list of Microtub
class MicrotubList {

 protected:

  ///NameList where names of Microtub are stored
  NameList       mtNameList;

  ///list of Microtub
  ObjectList     microtubs;

 public:

    ///creator
    MicrotubList() {
    microtubs.setNameList( & mtNameList );
  }

  ///destructor
  virtual ~MicrotubList() {
    MSG(13, "MicrotubList::destructor %p\n", this );
    eraseMicrotub();
    deleteMTGrid();
  }

  ///total number of Microtub
  int  nbMicrotub() const {
    return microtubs.size();
  }

  ///first Microtub in the list
  Microtub * firstMicrotub() const {
      return static_cast<Microtub*>(microtubs.getFront());
  }

  ///return a random Microtub in the list (slow)
  Microtub * randomMicrotub() const;

  ///mix the list
  void mixMicrotub() {
    microtubs.mixWell();
  }



  ///number of Microtub that do not belong to an Aster/Nucleus
  int  nbFreeMicrotub() const;

  ///total length of Microtub, also calculates Microtub::tubulin
  real totalMTLength();

  ///the total length divided by the number of Microtub
  real meanMTLength();


  //-------------------------------------------------------------

  // Functions for performing one simulation step:

  ///nucleate new Microtub is needed
  void nucleate();
  ///force the dynamic state of some Microtub
  void forceDynamicState();
  /// Monte-Carlo step for every Microtub
  void stepMicrotub();

  ///register a new Microtub in the list
  Microtub *  linkMicrotub( Microtub * );
  ///find a Microtub by its name
  Microtub *  findMicrotub( Name, int createIfNotFound = 0 );
  ///read a Name, find a corresponding Microtub and return it
  Microtub *  readMicrotubName();

  ///call modulo for all Microtub
  void moduloPositionMicrotub();

  ///rename Microtub in a contiguous manner
  void  renameMicrotub() {
    mtNameList.renameAll();
  }

  ///delete all Microtub
  void eraseMicrotub() {
    microtubs.erase();
    mtNameList.forgetAll();
  }

  ///delete Microtub which have not been updated by the last IO file read
  void cleanUpReadMicrotub(const int frameToKeep ) {
    microtubs.deleteIfFlagDifferentThan( frameToKeep );
  }

  ///write Microtub to IO file
  void writeMicrotub();


  //-------------------------------------------------------------
  // Functions for setting the grid for divide&conquer attachment

  ///step at which rod dispath was done
  Step   rod_dispatch_step;

  ///grid for divide-and-conquer strategies:
  static Map<DIM, Array<PointMicrotub*> > mtGrid;

  ///the current rod being painted:
  static PointMicrotub * paint_rod;

  ///the painting function:
  static void paintGrid( const int x, const int y, const int z );

  ///set the grid for the divide-and-conquer attachement, may return error code
  int  setMTGrid();

  ///delete the grid
  void deleteMTGrid();

  ///return true if the grid is ready
  bool isMTGridReady();

  // Attachments of motors, normal algorithm:
  ///distribute the MT rods on the grid
  void dispatchRods();

  ///given a Hand at a given position, find close rods and test attachement
  bool tryToAttach( const Vecteur &, Hand& );

  //------CHAITANYA
  // Record lengths, positions and identity of each MT.
  void countMTlenghtdyn();
  //CHAITANYA: for evaluation using statistics.out file
  //float maxLenP = 0, maxLenN=0, totLenP=0, totLenN=0, xMax=0;
  void getMTtotlen();

  // Functions assuming uniform free motors concentration:

#ifdef UNIFORM_FREE_COMPLEX
  ///attachement function for option UNIFORM_FREE_COMPLEX
  Microtub * uniformMicrotubSite(const real, real&) const;
  int attachToEnds(const int, const int, const real);
  int attachToMicrotubEverywhere(const int, const int, const real);
  int attachToMicrotubEnds(const int, const int, const real);
  void tryToAttachUniform(const int, const int, const real, const real, int& );
  void tryToAttachUniform();
#endif

  ///return the closest point to a given position
  PointMicrotub * closestPoint(const Vecteur &);

  ///test function for attachment scheme, at a particular position
  int  testAttach(Vecteur place, real dist);
  ///test function for attachement scheme
  void testAttach();

  ///CHAITANYA:Countthe capture status (if dist-to-nucleus-center <= MP.nuradius)
  int countMTcapt( int captTyp );

};


#endif  //MICROTUB_LIST_H

//RCS: $Id: grafted_list.h,v 2.4 2005/04/10 21:44:15 nedelec Exp $

#include "object_list.h"
#include "grafted.h"

///contain two list of Grafted: free and bound
class GraftedList {

  friend class Grafted;
  
 protected:

  ///NameList holding the names of the grafteds
  NameList       ghNameList;
  
  ///List holding the non-attached grafteds
  ObjectList     freeGrafteds;
  
  ///List holding the attached grafteds
  ObjectList     boundGrafteds;
  
 public:  
 
  ///creator
  GraftedList() {
      freeGrafteds.setNameList( &ghNameList );
      boundGrafteds.setNameList( &ghNameList );
    }

  ///destructor
  virtual ~GraftedList() {
    MSG(13, "GraftedList::destructor %p\n", this );
    eraseGrafted(); 
  }

  ///total number of Grafted
  int  nbGrafted() const {
    return freeGrafteds.size() + boundGrafteds.size(); 
  }

  ///number of Grafted of a given type
  int  nbGraftedOfType(int) const;

  ///perform Monte-Carlo step for all Grafted
  void stepGrafted();
  
  ///register a new Grafted in the appropriate list
  Grafted *  linkGrafted( Grafted * );
    
  ///find a Grafted from the name
  Grafted *  findGrafted( Name, int createIfNotFound = 0 );

  ///returns the first free Grafted
  Grafted *  firstFreeGrafted() const {
      return static_cast<Grafted*>( freeGrafteds.getFront() );
  }
  
  ///returns the first bound Grafted
  Grafted *  firstBoundGrafted() const {
      return static_cast<Grafted*>( boundGrafteds.getFront() );
  }
  
  ///mix the lists
  void mixGrafted() {
    freeGrafteds.mixWell();
    boundGrafteds.mixWell();
  }
  
  ///delete all Grafted
  void eraseGrafted() {
    boundGrafteds.erase();
    freeGrafteds.erase();
    ghNameList.forgetAll();
  }

  ///remove Grafted which have not been updated in the last read
  void cleanUpReadGrafted(const int frameToKeep ) {
    freeGrafteds.deleteIfFlagDifferentThan( frameToKeep );
    boundGrafteds.deleteIfFlagDifferentThan( frameToKeep );
  }
  
  ///calls moduloPosition for all Grafteds
  void moduloPositionGrafted();

  ///write Grafted on file IO
  void writeGrafted();

  ///debug function
  void checkGrafted();
};

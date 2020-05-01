//RCS: $Id: aster_list.h,v 2.4 2005/04/10 21:43:02 nedelec Exp $

#ifndef ASTER_LIST_H
#define ASTER_LIST_H

#include "object_list.h"
#include "aster.h"

///a list of Aster
class AsterList {
  
 protected:

  ///the NameList where names of Aster are stored
  NameList      asNameList;
  
  ///List of Asters
  ObjectList    asters;

 public:  
 
  ///creator
  AsterList() {
    asters.setNameList( &asNameList );
  }

  ///destructor
  virtual ~AsterList() {
    MSG(13, "AsterList::destructor %p\n", this );
    eraseAster(); 
  }

  ///total number of Aster
  int  nbAster() const {
    return asters.size(); 
  }

  ///Monte-Carlo step for the Asters
  void stepAster();

  ///register a new Aster in the list
  Aster *  linkAster( Aster * );
  
  ///find Aster in the list according to Name
  Aster *  findAster( Name, int createIfNotFound = 0 );
  
  ///first Aster in the list
  Aster *  firstAster() {
      return static_cast<Aster*>( asters.getFront() ); 
  }

  ///mix the list
  void mixAster() {
    asters.mixWell();
  }
  
  ///rename Aster to make Name contiguous
  void renameAster() {
    asNameList.renameAll();
  }

  ///delete every Aster in the list
  void eraseAster() {
    asters.erase();
    asNameList.forgetAll();
  }
  
  ///delete Aster which have not been updated during last IO file read
  void cleanUpReadAster(const int frameToKeep) {
    asters.deleteIfFlagDifferentThan( frameToKeep );  
  }
  
  ///read a Name, find a corresponding Aster, and return it (obsolete)
  Name readAsterName();           //obsolete
  
  ///calls moduloPosition for all asters
  void moduloPositionAster();
  
  ///write Aster list to IO file
  void writeAster();

  ///debug function
  void checkAster();
};


#endif //ASTER_LIST_H

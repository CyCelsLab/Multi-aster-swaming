//RCS: $Id: solid_list.h,v 2.4 2005/04/10 21:48:33 nedelec Exp $

#ifndef SOLID_LIST_H
#define SOLID_LIST_H

#include "object_list.h"
#include "solid.h"

///a list of Solid
class SolidList {

 protected:

  ///the NameList where names of Solid are stored
  NameList        soNameList;
  
  ///the list of Solid
  ObjectList      solids;

 public:  
 
  ///creator
  SolidList() {
    solids.setNameList( &soNameList ); 
  }

  ///destructor
  virtual ~SolidList() {
    MSG(13, "SolidList::destructor %p\n", this );
    eraseSolid(); 
  }

  ///total number of Solid
  int  nbSolid() const {
    return solids.size(); 
  }
  
  ///first Solid
  Solid * firstSolid() const {
      return static_cast<Solid*>( solids.getFront() );
  }

  ///mix the list
  void mixSolid() {
    solids.mixWell();
  }
  
  ///Monte-Carlo step for every Solid
  void stepSolid();
 
  ///register a Solid into the list
  Solid *  linkSolid( Solid * );
  
  ///find a Solid from its name
  Solid *  findSolid( Name, int createIfNotFound = 0 );
  
  ///read a name from IO file, find the corresponding Solid in the list, and return it
  Solid *  readSolidName();

  ///call modulo for all Microtub
  void moduloPositionSolid();

  ///rename all the Solid to make names contiguous
  void  renameSolid() {
    soNameList.renameAll();
  }

  ///delete every Solid in the list
  void eraseSolid() {
    solids.erase(); 
    soNameList.forgetAll();
  }
  
  ///delete Solid which have not been updated during last IO file read
  void cleanUpReadSolid(const int frameToKeep) {
    solids.deleteIfFlagDifferentThan( frameToKeep ); 
  }

  ///write every Solid to IO file
  void writeSolid();
};


#endif  //SOLID_LIST_H

//RCS: $Id: nucleus_list.h,v 2.4 2005/04/10 21:46:45 nedelec Exp $

#ifndef NUCLEUS_LIST_H
#define NUCLEUS_LIST_H

#include "node.h"
#include "name_list.h"
#include "object_list.h"
#include "nucleus.h"


///a list of nuclei
class NucleusList {

 protected:
  
  ///the NameList where names of nuclei are stored
  NameList    nuNameList;
  
  ///list of nuclei
  ObjectList  nuclei;

 public:

  ///creator
  NucleusList() { nuclei.setNameList( &nuNameList ); }

  ///destructor
  virtual ~NucleusList()
  {
    MSG(13, "NucleusList::destructor %p\n", this );
    eraseNucleus();
  }
  
  ///total number of nuclei
  int  nbNucleus() const { return nuclei.size(); }
  
  ///Monte-Carlo step for the nuclei
  void stepNucleus();
  
  ///register a new nucleus in the list  
  Nucleus * linkNucleus( Nucleus * );
  
  ///find nucleus in the list according to Name
  Nucleus * findNucleus( Name, int createIfNotFound = 0 );
  //  Nucleus * readNucleusName( FILE* );

  ///first nucleus in the list
  Nucleus * firstNucleus()
  {
    return static_cast<Nucleus*>( nuclei.getFront() );
  }
  
  ///delete every nucleus in the list
  void eraseNucleus()                  
  {
    nuclei.erase(); 
    nuNameList.forgetAll();
  }
    
  ///calls moduloPosition for all nuclei
  void moduloPositionNucleus();
  
  ///write nucleus list to IO file
  void writeNucleus();
  
  ///delete nuclei which have not been updated during last IO file read
  void cleanUpReadNucleus(const int frameToKeep) {
    nuclei.deleteIfFlagDifferentThan( frameToKeep );  
  }
  
};

#endif  //NUCLEUS_LIST_H

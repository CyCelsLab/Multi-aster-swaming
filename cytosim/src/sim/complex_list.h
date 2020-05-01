//RCS: $Id: complex_list.h,v 2.7 2005/04/20 12:19:56 clausen Exp $

#ifndef COMPLEX_LIST_H
#define COMPLEX_LIST_H

#include "object_list.h"
#include "complex.h"

///Complex are stored in 3 lists: free, bound (1 side bound), bridge (2 side bound)
class ComplexList {

  friend class Complex;
  
 protected:

  ///NameList holding the names of the complex
  NameList      cxNameList;
  
  ///list of complex which are not bound to any microtub
  ObjectList    freeComplex;
  
  ///list of complex bound to one and only one microtub
  ObjectList    boundComplex;
  
  ///list of complex bound to two microtubs
  ObjectList    bridgeComplex;
  
 public:  
 
  ///Creator
  ComplexList() {
     freeComplex.setNameList( &cxNameList );
    boundComplex.setNameList( &cxNameList );
   bridgeComplex.setNameList( &cxNameList );
  }

  ///destructor
  virtual ~ComplexList() {
    MSG(13, "ComplexList::destructor %p\n", this );
    eraseComplex(); 
  }

  ///total number of Complex
  int  nbComplex() const { 
    return freeComplex.size() + boundComplex.size() + bridgeComplex.size(); 
  }

  ///total number of free complexes
  int  nbFreeComplex() const   { return freeComplex.size();   }

  ///total number of bound complexes
  int  nbBoundComplex() const  { return boundComplex.size();  }

  ///total number of bridged complexes
  int  nbBridgeComplex() const { return bridgeComplex.size(); }

  ///number of complex of a given type
  int  nbComplexOfType(int) const;

  ///Monte-Carlo step for every Complex
  void stepComplex();

  ///register a Complex in the appropriate list
  Complex *  linkComplex( Complex * );
  
  ///find a Complex from its name
  Complex *  findComplex( Name, int createIfNotFound = 0 );
    
  ///first free Complex
  Complex *  firstFreeComplex() const {
    return static_cast<Complex*>( freeComplex.getFront() );
  }

  ///first bound Complex (one Hand attached)
  Complex *  firstBoundComplex() const {
    return static_cast<Complex*>( boundComplex.getFront() );
  }

  ///first bridge Complex (both Hand attached)
  Complex *  firstBridgeComplex() const {
    return static_cast<Complex*>( bridgeComplex.getFront() );
  }
  
  ///return a free complex of a specified type, or create one
  Complex *  freeComplexOfType(const int type);
  
  ///mix the lists
  void mixComplex() {
    freeComplex.mixWell();
    boundComplex.mixWell();
    bridgeComplex.mixWell();    
  }
  
  ///delete all Complex
  void eraseComplex() {
    bridgeComplex.erase(); 
    boundComplex.erase();
    freeComplex.erase();
    cxNameList.forgetAll();
  }

  ///remove Complex which have not been updated during last file read
  void cleanUpReadComplex(const int frameToKeep ) {
   bridgeComplex.deleteIfFlagDifferentThan( frameToKeep );
    boundComplex.deleteIfFlagDifferentThan( frameToKeep );
     freeComplex.deleteIfFlagDifferentThan( frameToKeep );
  }
  
  ///call modulo for all complex
  void moduloPositionComplex();

  ///write all Complex on IO file
  void writeComplex();

  ///debug function
  void checkComplex();
};


#endif //COMPLEX_LIST_H

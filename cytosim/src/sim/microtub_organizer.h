//RCS: $Id: microtub_organizer.h,v 1.11 2005/04/22 12:31:18 nedelec Exp $

#ifndef MICROTUB_ORGANIZER_H
#define MICROTUB_ORGANIZER_H

#include "object.h"
#include "microtub.h"
#include "array.h"

///class MicrotubOrganizer is an array of pointers to Microtub
/** MicrotubOrganizer is an automatic Array of Microtub pointers,
Zero values correspond to nucleation sites which have no Microtubs attached to them.
Functions allow for Microtub to be destroyed or nucleated during the simulation.
The mechanical links are handled in the derived classes: Aster & Nucleus
*/

class MicrotubOrganizer
{
private:
  
  ///array of pointers to the Microtub attached to this organizer
  Array<Microtub *>  moMT;  ///< list of pointers to Microtub

public:
    
  ///default constructor
  MicrotubOrganizer() {}
    
  ///destructor
  virtual   ~MicrotubOrganizer();

  ///the number of nucleation site on the aster, i.e. max number of Microtub it can contain
  int        maxNbMicrotub()             const  { return moMT.size(); }
  
  ///actual number of Microtub attached to the Organizer
  int        nbMicrotub()                const  { return moMT.nbNonZeroValues(); }
  
  ///return true if there is a Microtub attached in spot indx
  bool       hasMicrotub(const int indx) const  { if (indx < moMT.size()) return ( moMT[indx] != 0 ); return false; }
  
  ///return the Microtub attached to this aster in spot indx
  Microtub * getMicrotub(const int indx) const  { assert(indx < moMT.size()); return moMT[indx]; }
  
  ///find and return the index where the Microtub is stored, or -1 if the Microtub was not found
  int        findIndex(Microtub * mt)    const  { return moMT.find( mt ); }
    
  ///add Microtub at the end of the list, skipping any free spot
  int        registerMicrotub(Microtub * mt);
  
  ///add Microtub in the given spot indx, return the microtubule that was in this spot previously
  Microtub * registerMicrotub(Microtub * mt, const int indx);
  
  ///true if the Microtub can be deleted (if too short), otherwise the Microtub will be rescued
  virtual bool isDeletable(Microtub * mt) const { return true; }
  
  ///remove Microtub, returning the position where it was attached, or -1 if not found in this Organizer
  /** this can be overriden in a derived class, for example to allow replacement of the tube */
  virtual int deregisterMicrotub(Microtub * mt);
  
  ///remove all registered Microtubs
  void       deregisterAllMicrotubs();
  
  ///link the all Microtub referenced by the Organizer in the global sim. if needed
  void       linkMicrotubsIfNeeded();
  
  ///unlink all the Microtub referenced by the Organizer from the global sim
  void       unlinkMicrotubs();
  
  
  ///return the center of gravity from all MT ends specified
  Vecteur    getPosition(const MTEnd) const;
  
  ///return the center of gravity of all MTs
  virtual Vecteur    getPosition() const;
  
  ///move the aster with its associated solid and microtubs
  virtual void       translatePosition( const Vecteur & T );
  
  ///transform the aster with its associated solid and microtubs
  virtual void       transformPosition( const Transformation & T );
  
  ///inside function for space: returns how many MT points are inside
  virtual int        insidePosition( const Space * s = Microtub::space ) const;
  
  ///modulo function for periodic space
  virtual void       moduloPosition();
  
  
  
  ///debug function
  int        looksWrong()  const;
};



#endif //MICROTUB_ORGANIZER_H

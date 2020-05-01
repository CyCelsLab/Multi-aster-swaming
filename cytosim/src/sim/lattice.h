//RCS: $Id: lattice.h,v 2.8 2005/04/12 20:46:50 nedelec Exp $
// -------------------------------lattice.h--------------------------------

#ifndef LATTICE_H
#define LATTICE_H

#include "assert_macro.h"

enum LAHandType {  
  laNoHand = 0,        ///< no hand is bound at site: lasites[] = '0'
  laHandGeneric = 1    ///< any hand is attached: lasites[] = '1'
};

class Lattice
{

  ///the allocation of new points in the arrays is done by chuncks
  static const int LACHUNCK = 100;

 private:

  ///allocation size of array pssites[]
  int             laallocated;
  
  ///number of points in the set
  int             lamaxsites;

 protected:

  /// lasites[] is of size laallocated, and contain lamaxsites valid sites
  char     *      lasites;                 //usable range = [ 0 to lamaxsites [

 public:


  //--------------------------------------------------------------

  ///Constructor
  Lattice();
  
  ///Destructor
  virtual ~Lattice()                     { deallocate(); }

  ///allocate() is called by addSites for ex.,  returns size if memory was allocated
  virtual int   allocate(int size);     
  
  ///deallocate() free all memory allocated by allocate()
  virtual void  deallocate();
  
  ///No more sites
  void          clearSites()             { lamaxsites = 0; }  
  
  /// set all points to zero (debug/test)
  void          resetSites();

  ///Sets the site nSite to occupied
  void          setOccupied(int nSite, LAHandType handType);

  ///Sets the site nSite to occupied
  void          setOccupied(int nSite);
  
  ///Sets the site nSite to free
  void          setFree(int nSite);

  ///Tests if nSite if free
  enum LAHandType haType(int nSite) const;

  ///Tests if nSite if free
  bool           isFree(int nSite) const;

  ///Tests if nSite if occupied
  bool           isOccupied(int nSite) const;

  ///returns number of points
  int           nbSites()     const      { return lamaxsites; }
  
  ///Returns the address of the array of sites
  char    *     sites()                     { return lasites; }
  
  ///Returns the address of the array of sites
  char const *  getSites()            const { return lasites; }

  ///the address where site ii starts
  char    *     siteAddr(int ii)         { assert((ii>=0) && (ii<lamaxsites)); return lasites + ii; }
  
  ///write to file
  void          write();
  
  ///reads file
  void          read();
};

#endif

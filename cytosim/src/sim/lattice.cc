//RCS: $Id: lattice.cc,v 2.9 2005/04/27 13:28:56 clausen Exp $
//-------------------------------lattice.cc----------------------------------

#include "lattice.h"

#include "exceptions.h"
#include "iowrapper.h"
#include "iomessages.h"
#include "sim_param.h"

//----------------------------------------------------------------------------
Lattice::Lattice()
{
  laallocated = 0;
  lasites     = 0;
  lamaxsites  = 0;
}


//----------------------------------------------------------------------------
/** allocate(size) ensures that the set can hold <size> sites
 it returns the size if new memory was allocated 
 Sets newly allocated sites to unoccupied
*/
int Lattice::allocate(int size)
{
  if ( size < laallocated ) {
    lamaxsites = size;
    return 0;
  }

  int laallocated_new = size + LACHUNCK;

  char *lasites_new = new char [ laallocated_new ];

  if ( lasites_new == 0 )
    MSG.error("Lattice::allocate()", "memory allocation failed");

  if ( lasites ) {               
    // copy the current 
    for (int p = 0; p < lamaxsites; p++ )
      lasites_new[ p ] = lasites[ p ];
    delete [] lasites;
    // set new sites to unoccupied
    for (int p = lamaxsites; p < laallocated_new; p++ )
      lasites_new[ p ] = 0;
  } else {
    for (int p = 0; p < laallocated_new; p++ )
      lasites_new[ p ] = 0;
  }
  lasites      = lasites_new;
  laallocated  = laallocated_new;
  lamaxsites   = size;
  
  MSG(12, "Lattice::allocate %lx ask %d done %d\n", this, size, laallocated_new);
  return laallocated_new;
}


//---------------------------------------------------------------------------
void Lattice::deallocate()
{
  if ( lasites )
    delete[] lasites;   
  lasites = 0;
  laallocated = 0;
}

//----------------------------------------------------------------------------
void Lattice::setFree( int nSite )
{
  //If nSite < 0 we are beyind the MT end - assume we are at end
  if ( nSite < 0 ) {
    printf("nSite < 0. ERROR\n");
    nSite = 0;
  }

  if ( nSite < lamaxsites ) {
    if ( lasites[ nSite ] == 0 ) {
      printf("Freeing already free site: %i\n",nSite);
      exit(0);
    }
    lasites[ nSite ] = 0;
  }
}

//----------------------------------------------------------------------------
void Lattice::setOccupied( int nSite )
{
  //If nSite < 0 we are beyind the MT end - assume we are at end
  if ( nSite < 0 ) {
    printf("nSite < 0. ERROR\n");
    nSite = 0;
  }

  if ( !isFree ( nSite ) ) {
    printf("Setting already occoupied site: %i\n",nSite);
    exit(0);
  }

  if ( nSite + 1 > lamaxsites ) {
    allocate ( nSite + 1 );
  }

  // set to laHandGeneric
  lasites[ nSite ] = 1;
}

//----------------------------------------------------------------------------
void Lattice::setOccupied( int nSite, LAHandType haType)
{
  //If nSite < 0 we are beyind the MT end - assume we are at end
  if ( nSite < 0 ) {
    printf("nSite < 0. ERROR\n");
    nSite = 0;
  }

  if ( nSite + 1 > lamaxsites ) {
    allocate ( nSite + 1);
  }

  if ( haType == laHandGeneric ) {
    lasites[ nSite ] = 1;
  }
}

//----------------------------------------------------------------------------
LAHandType Lattice::haType( int nSite ) const
{
  //cant bind outside of lattice (not on MT!)
  if ( nSite<0 || nSite>=lamaxsites ) {
    return laNoHand;
  }

  if ( lasites[ nSite ] == 1 ) {
    return laHandGeneric;
  }

  return laNoHand;
}

//----------------------------------------------------------------------------
bool Lattice::isFree( int nSite ) const
{
  //If we try to bind outside of defined lattice it means we have not set that
  // site to occupied. Hence we assume it is free.
  if ( nSite >= lamaxsites || nSite < 0) {
    return 1;
  }

  return lasites[ nSite ] == 0;
}

//----------------------------------------------------------------------------

bool Lattice::isOccupied( int nSite ) const
{
  return !isFree( nSite );
}

//----------------------------------------------------------------------------
void Lattice::resetSites()
{
  assert( lasites );
  for( int q = 0; q < laallocated; ++q )
    lasites[ q ] = 0;
}

//****************************************************************************
//----------------------------------------------------------------------------
void Lattice::write()
{
  IO.writeUInt16( lamaxsites );
  for( int p = 0; p < lamaxsites ; ++p )
    IO.writeUInt8( *(lasites + p) );
}

//----------------------------------------------------------------------------
void Lattice::read()
{
  lamaxsites = 0;
  try {

    lamaxsites = IO.readUInt16();
    MSG(80, "Lattice::read( s%lx ) maxsites %i\n", this, lamaxsites );
    allocate( lamaxsites );
    resetSites();        //not absolutely necessary, but nicer for debugging

    for( int p = 0; p < lamaxsites ; ++p )
      lasites[p] = IO.readUInt8();

    //we do not call modulo, as this is dependent on the objects
    
  } catch( IOException e ) {
    //MSG("S%lx %i ", name, r);
    e.addBeforeMessage("Lattice::read : ");
    lamaxsites = 0;
    throw e;
  }
}

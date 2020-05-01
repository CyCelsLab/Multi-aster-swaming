//RCS: $Id: microtub_organizer.cc,v 1.12 2005/04/10 21:46:35 nedelec Exp $
//==============================aster.cc=============================

#include "microtub_organizer.h"
#include "sim.h"
#include "exceptions.h"
#include "iowrapper.h"

//-------------------------------------------------------------------
Vecteur MicrotubOrganizer::getPosition() const
{
  Vecteur result = VZERO;
  int nbPts = 0;
  //we take the center of gravity of all the Microtub-points
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) {
      for( int pt = 0; pt < mt->nbPoints(); ++pt )
        result += mt -> whereP(pt);
      nbPts += mt -> nbPoints();
    }
  }
  if ( nbPts ) result /= nbPts;
  return result;
}


//-------------------------------------------------------------------
Vecteur MicrotubOrganizer::getPosition(const MTEnd end) const
{
  Vecteur result = VZERO;
  int nbMTs = 0;
  //we take the center of gravity of all Microtub's ends
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) {
      result += mt -> whereEnd(end);
      ++nbMTs;
    }
  }
  if ( nbMTs ) result /= nbMTs;
  return result;
}

//-------------------------------------------------------------------
void MicrotubOrganizer::translatePosition( const Vecteur & T )
{  
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) mt -> translatePosition( T );
  }
}

//-------------------------------------------------------------------
void MicrotubOrganizer::transformPosition( const Transformation & T )
{
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) mt -> transformPosition( T );
  }
}

//-------------------------------------------------------------------
int MicrotubOrganizer::insidePosition( const Space * s ) const
{
  int result = 0;
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) result += mt -> insidePosition( s );
  }
  return result;
}

//-------------------------------------------------------------------
void MicrotubOrganizer::moduloPosition()
{
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) mt -> moduloPosition();
  }
}

//===================================================================
int MicrotubOrganizer::registerMicrotub(Microtub * mt)
{
  assert( mt != 0 );
  mt -> setOrganizer( this );
  return moMT.pushFront( mt );
}

//-------------------------------------------------------------------
Microtub * MicrotubOrganizer::registerMicrotub(Microtub * mt, const int indx)
{
  assert( indx >= 0 );
  
  //allocate if needed
  if ( indx >= moMT.size() )
    moMT.setSizeAndReset( indx+1 );
  
  //get the microtubule at the specified position
  Microtub * old_mt = moMT[ indx ];
  
  //if there was a different MT there, de-register it
  if ( old_mt ) {
    if ( old_mt == mt ) {
      if ( mt ) mt -> setOrganizer( this );
      return 0;
    }
    //MSG.warning("MicrotubOrganizer", "registerMicrotub() called for an non-empty spot containing a different MT");
    //No warning: this may happen in play, when reading an aster in which some MT were exchanged
    old_mt -> setOrganizer( 0 );
  }
  
  //link the given microtub
  moMT[ indx ] = mt;
  if ( mt ) mt -> setOrganizer( this );
  return old_mt;
}

//-------------------------------------------------------------------
int MicrotubOrganizer::deregisterMicrotub(Microtub * mt)
{
  assert( mt != 0 );

  int indx = moMT.find( mt );
  if ( indx >= 0 ) {
    assert( moMT[ indx ] == mt );
    moMT[ indx ] = 0;
    mt -> setOrganizer( 0 );
    return indx;
  }
  return -1;
}

//-------------------------------------------------------------------
void MicrotubOrganizer::deregisterAllMicrotubs()
{
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    if ( moMT[ mti ] ) {
      moMT[ mti ] -> setOrganizer( 0 );
      moMT[ mti ] = 0;
    }
  }
}

//-------------------------------------------------------------------
void MicrotubOrganizer::linkMicrotubsIfNeeded()
{
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) {
      
      //we link the Microtub in the sim::MicrotubList
      if ( ! mt -> isLinked() )
        sim.linkMicrotub( mt );
      
      //we also fix the Organizer if needed
      if ( mt -> getOrganizer() == 0 )
        mt -> setOrganizer( this );
    }
  }
}

//-------------------------------------------------------------------
void MicrotubOrganizer::unlinkMicrotubs()
{
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) {
      
      //we unlink the Microtub in the sim::MicrotubList
      if ( mt -> isLinked() )
        mt -> pop();
      
    }
  }
}


//-------------------------------------------------------------------
MicrotubOrganizer::~MicrotubOrganizer()
{
  MSG(13, "MicrotubOrganizer::destructor  %p\n", this);

  //we de-register the Microtub, without deleting them
  deregisterAllMicrotubs();
}

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
int MicrotubOrganizer::looksWrong() const
{
  for( int mti = 0; mti < moMT.size(); ++mti ) {
    Microtub * mt = moMT[ mti ];
    if ( mt ) {
      assert( mt -> getOrganizer() == this );
      assert( mt -> getName() > 0 );
    }
  }
  return NO_ERROR;
}


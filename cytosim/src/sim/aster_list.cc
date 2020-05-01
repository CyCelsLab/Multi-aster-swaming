//RCS: $Id: aster_list.cc,v 2.5 2005/01/11 16:58:32 foethke Exp $

#include "sim.h"
#include "iowrapper.h"

//-------------------------------------------------------------------
void AsterList::stepAster()
{
  //in the current model, nothing need to be done for this object
  //for(Aster * as=firstAster(); as ; as=as->next() )
  //  as -> step();
}


//-------------------------------------------------------------------
Aster * AsterList::linkAster( Aster * as )
{
  asters.pushFront( as );
  as->linkMicrotubsIfNeeded();
  if ( as->getSolid() ) {
    if (! as->getSolid() -> isLinked() )
      sim.linkSolid( as->getSolid() );
    assert( as->getName() == as->getSolid()->getName() );
  }
  
  MSG(70, "AsterList::link     as%lx ( total %i )\n", as->getName(), nbAster());
  return as;
}


//-------------------------------------------------------------------
Aster * AsterList::findAster( Name n, int createIfNotFound )
{
  MSG(90,"findAster( %i )\n", n);
  Aster * as = static_cast<Aster*>( asNameList.nodeWithThisName( n ) );
  if ( ( as == 0 ) && n )
    if ( createIfNotFound ) {
      as = new Aster();
      as -> setName( n );
      linkAster( as );
    }
  //else MSG.warning("AsterList::findAster", "cannot find %lx", n);
  return as;
}


//-------------------------------------------------------------------
Name AsterList::readAsterName()
  //read the name, and arrange the 'translation' to pointer
  //obsolete: this is not used if IO.getFileFormat() >= 18
{
  assert( IO.getFileFormat() < 18 );
  return IO.readUInt16();
}


//-------------------------------------------------------------------
void AsterList::moduloPositionAster()
{
  for(Aster * as=firstAster(); as ; as=as->next() )
    as->moduloPosition();
}


//-------------------------------------------------------------------
void AsterList::writeAster()
{
  IO.writeRecordTag("AS");
  IO.writeIntAscii( nbAster(), asNameList.nextName());

  for(Aster * as=firstAster(); as ; as=as->next() )
    as->write();
}


//-------------------------------------------------------------------
void AsterList::checkAster()
{
    for(Aster * as=firstAster(); as ; as=as->next() ) {
      assert( as->looksWrong() == NO_ERROR );
    }
}

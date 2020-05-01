//RCS: $Id: solid_list.cc,v 2.5 2005/02/02 19:30:38 nedelec Exp $
//============================solidlist.cc=======================

#include "sim.h"
#include "iowrapper.h"

//-------------------------------------------------------------------
void SolidList::stepSolid()
{  
  //we have to do nothing for this kind of object.
}


//-------------------------------------------------------------------
Solid * SolidList::linkSolid(Solid * so)
{
  //We link the solid's grafted if they are not linked already:
  for( int ii = 0; ii < so->nbGrafteds(); ++ii )
    if ( ! so->getGrafted( ii )-> isLinked() )
      sim.link( so->getGrafted( ii ) );

  solids.pushFront( so ); 
  MSG(70, "SolidList::link     so%lx ( total %i )\n", so->getName(), nbSolid());

  return so;
}


//-------------------------------------------------------------------
Solid * SolidList::findSolid( Name n, int createIfNotFound )
  //local storage of name list in a buffer for performance
{
  //MSG(90, "findSolid( %i )\n", n );
  Solid * so = static_cast<Solid*>( soNameList.nodeWithThisName( n ) );
  if ( ( so == 0 ) && n )
    if ( createIfNotFound ) {
      so = new Solid();
      so -> setName( n );
      linkSolid( so );
    }
  //else MSG.warning("SolidList::findSolid", "cannot find k%lx", n);
  return so;
}


//-------------------------------------------------------------------
Solid * SolidList::readSolidName()
  //read the name, and arrange the 'translation' to microtubule pointer
{
  Name name = IO.readUInt16();
  if ( name == 0 ) return 0;
  return findSolid( name, 0 );
}

//-------------------------------------------------------------------
void SolidList::moduloPositionSolid()
{
  for(Solid * soi=firstSolid(); soi ; soi=soi->next())
    soi->moduloPosition();
}

//-------------------------------------------------------------------
void SolidList::writeSolid()
{
  IO.writeRecordTag("SO");
  IO.writeIntAscii( nbSolid(), soNameList.nextName());

  for(Solid * soi=firstSolid(); soi ; soi=soi->next())
    soi->write();
}



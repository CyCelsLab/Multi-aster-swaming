//RCS: $Id: nucleus_list.cc,v 2.5 2005/01/08 10:28:09 nedelec Exp $

#include "nucleus_list.h"
#include "grafted.h"
#include "iowrapper.h"
#include "sim.h"


void NucleusList::stepNucleus()
{
  //there's nothing to be done for this kind of object.

}

//----------------------------------------------------------------------------

Nucleus * NucleusList::linkNucleus(Nucleus * nu)
{
  //link the nucleus' Microtubules
  MSG(70, "NucleusList::link     a%lx ( total %i )\n", nu->getName(), nbNucleus());
  nu->linkMicrotubsIfNeeded();

  // link the nucleus' grafted if they are not linked already:
  for( int ii = 0; ii < nu->nbGrafteds(); ++ii )
    if ( ! nu->getGrafted( ii )-> isLinked() )
      sim.link( nu->getGrafted( ii ) );

  nuclei.pushFront( nu );

  return nu;
}

//----------------------------------------------------------------------------

Nucleus * NucleusList::findNucleus( Name n, int createIfNotFound )
  //local storage of name list in a buffer for performance
{
  MSG(90, "findNucleus( %i )\n", n );
  Nucleus * nu = static_cast<Nucleus*>( nuNameList.nodeWithThisName( n ) );
  if ( ( nu == 0 ) && n )
    if ( createIfNotFound ) {
      nu = new Nucleus();
      nu -> setName( n );
      linkNucleus( nu );
    }
  //else MSG.warning("NucleusList::findNucleus", "cannot find k%lx", n);
  return nu;
}

//----------------------------------------------------------------------------

void NucleusList::moduloPositionNucleus()
{
  for(Nucleus * nui=firstNucleus(); nui; nui=nui->next())
    nui->moduloPosition();
}

//----------------------------------------------------------------------------

void NucleusList::writeNucleus()
{
  IO.writeRecordTag("NU");
  IO.writeIntAscii( nbNucleus(), nuNameList.nextName());

  for(Nucleus * nui=firstNucleus(); nui; nui=nui->next())
    nui->write();
}

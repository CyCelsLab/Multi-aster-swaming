//RCS: $Id: grafted_list.cc,v 2.3 2005/03/16 20:20:55 nedelec Exp $
//============================graftedlist.cc=======================

#include "sim.h"
#include "iowrapper.h"

//------------------------------------------------------------------
void GraftedList::stepGrafted()
{
  Grafted * ghi, * gh;

  //call step() for free grafted
  ghi = firstFreeGrafted();
  gh  = ghi;
  while( gh ) {
    ghi = ghi->next();
    gh->stepFree();
    gh = ghi; 
  }
  
  // call step() for all attached hands
  ghi = firstBoundGrafted();
  gh  = ghi;
  while( gh ) {
    ghi = ghi->next();
    gh->stepBound();
    gh = ghi; 
     
  }
}

//------------------------------------------------------------------
int GraftedList::nbGraftedOfType(int whichtype) const
{
  int result = 0;
  Grafted * ghi;
  for( ghi = firstFreeGrafted(); ghi ; ghi=ghi->next() )
    if ( ghi->getType()==whichtype) result++;
  for( ghi = firstBoundGrafted();  ghi ; ghi=ghi->next() )
    if ( ghi->getType()==whichtype) result++;
  return result;
}



//------------------------------------------------------------------
Grafted * GraftedList::linkGrafted(Grafted * gh)
{
  gh->setState();

  if ( gh->ghHand.isAttached() )
    boundGrafteds.pushFront(gh); 
  else
    freeGrafteds.pushFront(gh); 

  MSG(70, "GraftedList::link g%lx ( total %i )\n", gh->getName(), nbGrafted());
  return gh;
}


//------------------------------------------------------------------
Grafted * GraftedList::findGrafted( Name n, int createIfNotFound )
{
  MSG(90,"findGrafted( name %i, create %i )\n", n, createIfNotFound );
  Grafted * gh = static_cast<Grafted*>( ghNameList.nodeWithThisName( n ) );
  assert( ( gh == 0 ) || ( gh->getName() == n ));
  if ( ( gh == 0 ) && n )
    if ( createIfNotFound ) {
      gh = new Grafted();
      gh -> setName( n );
      linkGrafted( gh );
    }
  //else MSG.warning("GraftedList::findGrafted", "cannot find g%lx", n);
  return gh;
}


//------------------------------------------------------------------
void GraftedList::moduloPositionGrafted()
{
  Grafted * ghi;
  for( ghi = firstFreeGrafted(); ghi ; ghi=ghi->next() )
    ghi->moduloPosition();
  
  for( ghi = firstBoundGrafted();  ghi ; ghi=ghi->next() )
    ghi->moduloPosition();
}


//------------------------------------------------------------------
void GraftedList::writeGrafted()
{
  IO.writeRecordTag("GH");
  IO.writeIntAscii(nbGrafted(), ghNameList.nextName());
  
  Grafted * ghi;
  for( ghi = firstFreeGrafted(); ghi ; ghi=ghi->next() ) {
    ghi->write();
    assert( ghi->isFree() );
  }
 
  for( ghi = firstBoundGrafted();  ghi ; ghi=ghi->next() ) {
    ghi->write();
    assert( ghi->isAttached() );
  }
}



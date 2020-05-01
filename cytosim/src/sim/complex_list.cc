//RCS: $Id: complex_list.cc,v 2.5 2005/04/22 11:44:06 nedelec Exp $
//==================================complexlist.cc===========================

#include "sim.h"
#include "iowrapper.h"

//-------------------------------------------------------------------------
//-------------------------- STEPS FOR COMPLEX ----------------------------
//-------------------------------------------------------------------------


//elementary steps for all the complex in the box
void ComplexList::stepComplex()
{
  Complex * cxi, * cx;

  /*
  MSG(9,"ComplexList::stepComplex entry :  free %5i, bound %5i, bridge %5i\n", 
	sim.freeComplex.size(),
	sim.boundComplex.size(),
	sim.bridgeComplex.size());
  */

  //--------------------------BIND------------------------------------
  //binding can send complexes from free to bound, and bound to linking
  //To avoid processing the same complex twice, we process the list
  //in the inverse order : bound, then free
  
  cxi = firstBoundComplex();
  cx  = cxi;
  while( cx ) {
    cxi = cxi->next();
    cx->tryToAttachBound();
    cx = cxi; 
  }

#ifdef UNIFORM_FREE_COMPLEX

  sim.tryToAttachUniform();

#else

  cxi = firstFreeComplex();
  cx  = cxi;
  while( cx ) {
    cxi = cxi->next();
    cx->tryToAttachFree();
    cx = cxi; 
  }

#endif

  //---------------------------------STEP------------------------------  
  //below, unbinding can occur, which send complexes from bound to free
  //and linking to bound. To avoid testing for unbinding twice, we process 
  //bound before linking. We process free last, to allow motor which just
  //unbound to diffuse, this reduces recapture events

  cxi = firstBoundComplex();
  cx  = cxi;
  while( cx ) {
    cxi = cxi->next();
    cx->stepBound();
    cx = cxi; 
  }

  cxi = firstBridgeComplex();
  cx  = cxi;
  while( cx ) {
    cxi = cxi->next();
    cx->stepBridge();
    cx = cxi; 
  }
  
  cxi = firstFreeComplex();
  cx  = cxi;
  while( cx ) {
    cxi = cxi->next();
    cx->stepFree();
    cx = cxi;
  }
}


//============================================================================
//============================================================================
//============================================================================

//------------------------------------------------------------------------
int ComplexList::nbComplexOfType(int whichtype) const
{
  int result = 0;
  Complex * cxi;

  for(cxi=firstFreeComplex(); cxi ; cxi = cxi->next() )
    if ( cxi->getType() == whichtype ) result++;

  for(cxi=firstBoundComplex(); cxi ; cxi = cxi->next() )
    if ( cxi->getType() == whichtype ) result++;

  for(cxi=firstBridgeComplex(); cxi ; cxi = cxi->next() )
    if ( cxi->getType() == whichtype ) result++;

  return result;
}


//------------------------------------------------------------------------
Complex * ComplexList::freeComplexOfType(const int type)
{
  Complex * cx = firstFreeComplex();
  while ( cx && ( cx->getType() != type ))
    cx = cx->next();
  
  //if this failed, we create a complex of the requested type:
  if ( cx == 0 ) 
    cx = linkComplex( new Complex( type ));
  
  return cx;
}


//------------------------------------------------------------------------
Complex * ComplexList::linkComplex(Complex * cx)
{
  cx -> setState();
          
  switch( cx->nbBoundHands() ) {
    case 0:   freeComplex.pushFront(cx);  break;
    case 1:  boundComplex.pushFront(cx);  break;
    case 2: bridgeComplex.pushFront(cx);  break;
  }

  //MSG(70, "Complexlist::linkComplex cx%lx ( total %i )\n", cx->getName(), nbComplex());
  return cx;
}

//------------------------------------------------------------------------
Complex * ComplexList::findComplex( Name n, int createIfNotFound )
{
  //MSG(90,"findComplex( %i ) in %lx\n", n, all );
  Complex * cx = static_cast<Complex*>( cxNameList.nodeWithThisName( n ) );  
  assert( ( cx == 0 ) || ( cx->getName() == n ));
  if ( ( cx == 0 ) && n )
    if ( createIfNotFound ) {
      //MSG("create c%lx\n", n);
      cx = new Complex();
      cx -> setName( n );
      linkComplex( cx );
    }
      //else MSG.warning("ComplexList::findComplex", "cannot find m%lx", n);
  return cx;
}


//------------------------------------------------------------------------
void ComplexList::moduloPositionComplex()
{
  Complex * cxi;
  for(cxi=firstFreeComplex(); cxi ; cxi = cxi->next() )
    cxi->moduloPosition();
  
  for(cxi=firstBoundComplex(); cxi ; cxi = cxi->next() )
    cxi->moduloPosition();
  
  for(cxi=firstBridgeComplex(); cxi ; cxi = cxi->next() )
    cxi->moduloPosition();
}

//------------------------------------------------------------------------
void ComplexList::writeComplex()
{
  IO.writeRecordTag("CX");
  IO.writeIntAscii( nbComplex(), cxNameList.nextName());

  Complex * cxi;
  for(cxi=firstFreeComplex(); cxi ; cxi = cxi->next() )
      cxi->write();
  
  for(cxi=firstBoundComplex(); cxi ; cxi = cxi->next() )
      cxi->write();
  
  for(cxi=firstBridgeComplex(); cxi ; cxi = cxi->next() )
      cxi->write();
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



void ComplexList::checkComplex()
{
  Complex * cxi; 
  assert( freeComplex.looksWrong() == NO_ERROR );
  for(cxi=firstFreeComplex(); cxi ; cxi = cxi->next() ) {
    if ( ! cxi->isFree() )
      MSG("wrong free cx %lx\n", cxi->getName() );
  }

  assert( boundComplex.looksWrong() == NO_ERROR );
  for(cxi=firstBoundComplex(); cxi ; cxi = cxi->next() ) {
    if ( ! cxi->isAttached() )
      MSG("wrong bound cx %lx\n", cxi->getName() );
  }

  assert( bridgeComplex.looksWrong() == NO_ERROR );
  for(cxi=firstBridgeComplex(); cxi ; cxi = cxi->next() ) {
    if ( ! cxi->isBridge() )
      MSG("wrong link cx %lx\n", cxi->getName() );
  }
}

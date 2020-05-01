//RCS: $Id: complex.cc,v 2.19 2005/04/22 18:54:02 nedelec Exp $
//================================Complex.cc==============================
//------------------------------------------------------------------------
//========================================================================


// Neha Khetan, 10 December 2016 - modified the write file as earlier error i nreading the complexes while running the play command. Now solved.
#include "complex.h"
#include "sim.h"
#include "exceptions.h"
#include "iowrapper.h"

Space * Complex::space = 0;

//------------------------------------------------------------------------
void Complex::complexConstructorBase(const int type, const int hatype1, const int hatype2, const Vecteur w)
{
  assert( type >= 0 );
  assert( type < MP.MAX );
  cxType     = type;
  cxPos      = w;
  cxHand1.handConstructorBase(HA_BASE_CX_SIDE1, this, hatype1);
  cxHand2.handConstructorBase(HA_BASE_CX_SIDE2, this, hatype2);
  cxStateUpdateTime = sim.iterationCnt();
}


//------------------------------------------------------------------------
Complex::Complex(int type, int hatype1, int hatype2, Vecteur w )
{
  complexConstructorBase( type, hatype1, hatype2, w );
}


//------------------------------------------------------------------------
Complex::Complex(int type, Vecteur w )
{
  complexConstructorBase( type, MP.cxhatype1[type], MP.cxhatype2[type], w );
}

//------------------------------------------------------------------------
Complex & Complex::operator = (const Complex & c)
{
  name              = c.name;
  cxType            = c.cxType;
  cxPos             = c.cxPos;
  cxHand1           = c.cxHand1;
  cxHand1.haBase    = this;
  cxHand2           = c.cxHand2;
  cxHand2.haBase    = this;
  cxState           = getState();
  cxStateUpdateTime = c.cxStateUpdateTime;

  assert( cxHand1.haBaseType == HA_BASE_CX_SIDE1 );
  assert( cxHand2.haBaseType == HA_BASE_CX_SIDE2 );
  return *this;
}


//------------------------------------------------------------------------
Complex::Complex(int type)
{
  unsigned long ouf=0;
  Vecteur place = VZERO;
  assert( ( type>=0 ) && ( type < MP.MAX ) );

  do {

    place = initPosition( space, MP.cxinit[type] );
    if (++ouf>MAX_TRIALS_BEFORE_STUCK)
        throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in createComplex");

  } while ( ! space->isInside(place) );

  *this = Complex(type, place);
}


//------------------------------------------------------------------------
Complex::~Complex()
{
  MSG(13, "Complex::destructor %lx\n", name);
  if (cxHand1.isAttached()) cxHand1.detach( DETACH_DESTRUCT );
  if (cxHand2.isAttached()) cxHand2.detach( DETACH_DESTRUCT );
}




//========================================================================
//------------------------------------------------------------------------
void Complex::stepBound()
{
  assert( nbBoundHands() == 1 );
  if ( cxHand1.isAttached() )
    cxHand1.step();
  else
    cxHand2.step();
}

//------------------------------------------------------------------------
void Complex::stepBridge()
{
  assert( isBridge() );
  static Vecteur V;
  
  V = cxHand2.where() - cxHand1.where();
  Space::modulo( V );

  //correct the force if cxlength > 0
  if ( MP.cxlength[ cxType ] > 0 ) {
    real d = V.norm();
    if ( d > EPSILON ) 
      V *= ( 1.0 - MP.cxlength[ cxType ] / d );
  }
  
  cxHand1.stepUnderForce( V );
  V.oppose();
  cxHand2.stepUnderForce( V );
  
}

//------------------------------------------------------------------------
void Complex::stepFree()
{
  assert( isFree() );

  cxPos.addRandom( MP.cxdiff_dt[ cxType ] );
  
  if ( ! insidePosition() ) {
    static Vecteur V;
    space->project( cxPos, V );
    cxPos = 2 * V - cxPos;
  }
}

//------------------------------------------------------------------------
void Complex::tryToAttachFree()
{
  assert( isFree() );
  
  if ( RNG.flip() )                //give equal chance to both sides
    sim.tryToAttach(cxPos, cxHand1 );
  else
    sim.tryToAttach(cxPos, cxHand2 );
}

//------------------------------------------------------------------------
void Complex::tryToAttachBound()
{
  assert( nbBoundHands() == 1 );
  
  if ( cxHand1.isAttached() )
    sim.tryToAttach(cxHand1.where(), cxHand2 );
  else
    sim.tryToAttach(cxHand2.where(), cxHand1 );
}



//------------------------------------------------------------------------
//exchanges the two hands...without changing anything physically
void Complex::swapHands()
{
  if ( cxHand1.isLinked() )  cxHand1.pop();
  if ( cxHand2.isLinked() )  cxHand2.pop();

  static Hand tmp;       //static for efficiency

  tmp     = cxHand2;
  cxHand2 = cxHand1;
  cxHand1 = tmp;

  cxHand1.haBaseType = HA_BASE_CX_SIDE1;
  cxHand2.haBaseType = HA_BASE_CX_SIDE2;

  if ( cxHand1.isAttached() ) cxHand1.getMT() -> linkHand( & cxHand1 );
  if ( cxHand2.isAttached() ) cxHand2.getMT() -> linkHand( & cxHand2 );

  tmp.clear();       //to escape the destructor
}

//------------------------------------------------------------------------
void Complex::swapHandTypes()
{
  int type = cxHand1.haType;
  cxHand1.haType = cxHand2.haType;
  cxHand2.haType = type;
}

//------------------------------------------------------------------------
Vecteur Complex::calculatePosition()
{
  if ( cxHand2.isAttached() ) cxPos = cxHand2.where();
  else
    if ( cxHand1.isAttached() ) cxPos = cxHand1.where();
  return cxPos;
}


//------------------------------------------------------------------------
Vecteur Complex::getStretch() const
//returns distance between first and second hand
{
  Vecteur tmp = cxHand2.where() - cxHand1.where();
  Space::modulo( tmp );
  return ( tmp );
}


//------------------------------------------------------------------------
void Complex::setState()
{
  cxHand1.setState();
  cxHand2.setState();
  cxState           = getState();
  cxStateUpdateTime = sim.iterationCnt();
}

//------------------------------------------------------------------------
void Complex::updateState( int why ) 
{
  assert( cxHand1.haBaseType == HA_BASE_CX_SIDE1 );
  assert( cxHand2.haBaseType == HA_BASE_CX_SIDE2 );
  assert( cxHand1.haBase == this );
  assert( cxHand2.haBase == this );

  int newstate = getState();
  if ( newstate != cxState )
    {
      cxState           = newstate;
      cxStateUpdateTime = sim.iterationCnt();

      if ( isLinked() )
        switch ( nbBoundHands() )  {
          case 0: sim.freeComplex.popAndPush( this );  break;
          case 1: sim.boundComplex.popAndPush( this ); break;
          case 2: sim.bridgeComplex.popAndPush( this ); break;
            //FIX EXPERIMENTAL: complex flipping for Gohta simulations mars 2005. F. Nedelec
            //The motor flips when it binds a second time:
            //if (( MP.magic == MAGIC_TYPE ) && ( MP.initcode == 40 )) Complex::swapHandTypes();
            //break;
          default: MSG.error("Complex::updateState","invalid nbBoundHands()");
        }
    }
}

//------------------------------------------------------------------------
bool  Complex::isStateNew() const
{
  return cxStateUpdateTime == sim.iterationCnt(); 
}


//------------------------------------------------------------------------
//-----------------------------------FILE---------------------------------
//------------------------------------------------------------------------

void Complex::write()
//prints 'this' on IO, using global coordinates: mt., rod nb, abscices.
{
  //we use a different 'twoLetters' ID: cX instead of cx 
  //to distinguish standard (2 bytes) names from long names
  assert( name > 0 );

  if ( name < 65536 ) {
    IO.writeRecordTag( "cx" );
    IO.writeUInt32( name );
  } else {
    IO.writeRecordTag( "cX" );
    IO.writeUInt32( name );  
  }
      
  IO.writeUInt16( cxType );
  cxHand1.write();
  cxHand2.write();
  if ( isFree() ) 
    IO.writeReal32Vect(DIM, cxPos.getXYZ());
}

//------------------------------------------------------------------------
void Complex::read()
//reads a complex in global coordinates.
{
  setFlag( sim.frameInBuffer() );

  MSG(80, "Complex::read( c%lx )\n", name);
  try {

    cxType = IO.readUInt32();
    cxHand1.read();
    cxHand2.read();
    
    if ( isFree() ) {
      IO.readReal32Vect(cxPos);																											// Neha, modified Dec. 2016.
    } else {
      calculatePosition();
    }
    
  } catch( IOException e ) {
    e.addBeforeMessage("Complex::read : ");
    throw e;
  }
}



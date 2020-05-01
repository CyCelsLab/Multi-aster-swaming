//RCS: $Id: grafted.cc,v 2.21 2005/04/22 18:54:29 nedelec Exp $
//==============================grafted.cc=============================
#include "grafted.h"
#include "sim.h"
#include "exceptions.h"
#include "iowrapper.h"

Space * Grafted::space = 0;

//========================================================================

void Grafted::graftedConstructorBase( const int type, const Vecteur w, const GraftedType kind )
{
  ghBase.clear();
  ghBaseType        = kind;
  if ( kind == GH_CORTEX )
    space->project( w, ghPos );
  else
    ghPos = w;
  ghHand.handConstructorBase(HA_BASE_GRAFTED, this, type);
}

//------------------------------------------------------------------------
Grafted::Grafted( const int type, const Vecteur w, const GraftedType kind )
{
  assert(( kind == GH_BASE ) || ( kind == GH_CORTEX ) || ( kind == GH_DIFFUSE ));
  graftedConstructorBase( type, w, kind );
}


//------------------------------------------------------------------------
Grafted::Grafted( const int type, const Solid * so, const int point_index )
{
  graftedConstructorBase( type );
  ghBase.setTo( so, point_index );
  ghBaseType = GH_SOLID;
}

//------------------------------------------------------------------------
Grafted::Grafted( const int type, const Nucleus * nu, const int point_index )
{
  graftedConstructorBase( type );
  ghBase.setTo( nu, point_index );
  ghBaseType = GH_NUCLEUS;
}


//------------------------------------------------------------------------
//copy assignment, copies the members of ghHand as well.
Grafted & Grafted::operator =(const Grafted & g)
{
  name             = g.name;
  ghPos            = g.ghPos;
  ghHand           = g.ghHand;
  ghHand.haBase    = this;
  ghBase           = g.ghBase;
  ghBaseType       = g.ghBaseType;


  assert( ghHand.haBaseType == HA_BASE_GRAFTED );
  assert( ghHand.haBase == this );
  return *this;
}


//------------------------------------------------------------------------
Grafted::Grafted(int type)
{
  unsigned long ouf=0;
  Vecteur w;

  do {

    //type = RNG.pint( MP.MAX, MP.ghmax );
    w = initPosition( space, MP.ghinit[type], type );      

    if (++ouf>MAX_TRIALS_BEFORE_STUCK)
      throw StuckException("reached MAX_TRIALS_BEFORE_STUCK in createGrafted");
    
  } while ( ! space->isInside( w ) );
  
  graftedConstructorBase( type, w ); 
  
  assert( ghHand.haBaseType == HA_BASE_GRAFTED );
  assert( ghHand.haBase     == this );
}


//------------------------------------------------------------------------
Grafted::~Grafted()
{
  MSG(13, "Grafted::destructor %p\n", this );
  if ( ghHand.isAttached() )
    ghHand.detach( DETACH_DESTRUCT );
}


//========================================================================
//========================================================================
//========================================================================

void Grafted::stepFree()
{
  assert( isFree() );      
  
  //Cortical Grafted diffuse in the membrane (experimental)
  switch( ghBaseType ) {
    case GH_DIFFUSE:
      ghPos.addRandom(MP.cxdiff_dt[getType()]);
      if ( space->isOutside(ghPos.getXYZ()) ) {
        static Vecteur V;
        space->project( ghPos, V );
        ghPos = 2 * V - ghPos;
      } 
      break;
    case GH_CORTEX:
      ghPos.addRandom(MP.cxdiff_dt[getType()]);
      space->project(ghPos);        
      break;
    default:
      break;
  }
    
  //attachment trial:
  sim.tryToAttach(whereGrafted(), ghHand );
}

//------------------------------------------------------------------------
void Grafted::stepBound()
{
  static Vecteur stretch;
  assert( isAttached() );
  
  if ( isDiffusing() ) {
    ghHand.step();
    return;
  }
  
  stretch = whereGrafted() - ghHand.where();
  Space::modulo( stretch );
  
  ghHand.stepUnderForce( stretch );
  
  /*
  //Attached cortical Grafted are easily dragged by their Hand (experimental)
  if ( isCortical() ) {
    ghPos -= stretch;
    space->project(ghPos);
    //printf("Dragging of cortical based on stretch: %2.6f \n", stretch.norm() );
  }
  */

  /* Neha Khetan, 19 October 2018 
  Implementing the motor slip based on the stretch of the motor
  Slip causes the displacement of the Anchor point i.e. base of the Grafted gets displaced.
  Also, as in Grover et al. Motor 
  */
  if( MP.haSlip[ ghHand.getType()] )
    {
    //printf("Gh %d \t haSlip: %i \n", ghHand.getType() , MP.haSlip[ ghHand.getType()] );
    ghPos -= stretch;
    space->project(ghPos);
    }
}

//------------------------------------------------------------------------
Vecteur Grafted::whereGrafted()
{ 
  //attached to a solid/Microtubule
  if ( isLink() ) 
    ghPos = ghBase.where();
      
  return ghPos; 
}

//------------------------------------------------------------------------
Vecteur Grafted::getStretch()
{
  Vecteur tmp = whereGrafted() - ghHand.where();
  Space::modulo( tmp );
  return tmp;
}

//------------------------------------------------------------------------
void Grafted::updateState( const int why )
{
  assert( ghHand.haBaseType == HA_BASE_GRAFTED );
  assert( ghHand.haBase == this );

  if ( isStateNew() ) {
    //printf("M%-5lx update %i time %6i\n", name, isAttached(), sim.iterationCnt() );
    if ( isLinked() ) {
      if ( isAttached() ) 
        sim.boundGrafteds.popAndPush( this ); 
      else
        sim.freeGrafteds.popAndPush( this );
    }
  }
}

//========================================================================
void Grafted::write()
{
  assert( name > 0 );
  
  //write a name, or a long name, using a different record-tag gH or gh
  if ( name < 2<<32 ) {
    IO.writeRecordTag( "gh" );
    IO.writeUInt32( name );        //Neha, MARCH 2016
  } else {
    IO.writeRecordTag( "gH" );
    IO.writeUInt32( name );  
  }

  //write the hand
  ghHand.write();
  
  //write the base type in clear ASCII
  IO.writeCharAscii( ghBaseType );

  if ( isLink() ) {
  
    //for a grafted on a solid/nucleus: write the object's Name and point
    //IO.writeUInt16( ghBase.getPS()->getName() );
    //IO.writeUInt16( ghBase.getPoint() );
	IO.writeUInt32( ghBase.getPS()->getName() );//CA: 2012dec
	IO.writeUInt32( ghBase.getPoint() );//CA: 2012dec
	  
  } else {
    
    //for an immobile grafted, write the coordinates
    IO.writeReal32Vect( DIM, ghPos.getXYZ() );
    
  }
}

//------------------------------------------------------------------------
void Grafted::read()
{
  setFlag( sim.frameInBuffer() );

  try {

    MSG(80, "Grafted::read( g%lx )\n", name );
    ghHand.read();
    ghPos.clear();
    ghBase.clear();

    Name name = 0;
    SimObject * gps = 0;
    
    if (IO.getFileFormat() > 18 ) {
      
      ghBaseType = GraftedType(IO.readCharAscii());
      
      switch ( ghBaseType ) {
        
        case GH_BASE:
         break;
                    
        case GH_CORTEX:
          break;
          
        case GH_DIFFUSE:
          break;
          
        case GH_SOLID:
		  if ( name < 65536 ) {
		    name = IO.readUInt16();
		  } else {
		    name = IO.readUInt32();
		  }
		  //name = IO.readUInt16();//CA: to overcome unit16 limit 2012dec
          gps  = sim.findSolid( name, 0 );
          if ( gps == 0 )
            throw IOException("Solid not found");
          break;
          
        case GH_NUCLEUS:
		  if ( name < 65536 ) {
		    name = IO.readUInt16();
		  } else {
		    name = IO.readUInt32();
		  }
		  //name = IO.readUInt16();
          gps  = sim.findNucleus( name, 0 );
          if ( gps == 0 )
            throw IOException("Nucleus not found");
          break;
        
         default:
          throw IOException("invalid GraftedType");

      }
      
    } else {
    
      if ( IO.getFileFormat() > 16 ) {
		if ( name < 65536 ) {
		  name = IO.readUInt16();
		} else {
		  name = IO.readUInt32();
		}
		  //name = IO.readUInt16();
        if ( name ) {
          gps = sim.findSolid( name, 0 );
          if ( gps == 0 )
            throw IOException("Solid not found");
        }
      }
      
    }
    
    if ( gps == 0 ) {
      IO.readReal32Vect( ghPos );
    } else {
      ghBase.setTo( gps, IO.readUInt16() );
    }
    
  } catch( IOException e ) {
    e.addBeforeMessage("Grafted::read() : ");
    throw e;
  }
}


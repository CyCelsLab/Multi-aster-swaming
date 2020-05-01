//RCS: $Id: sim.cc,v 2.18 2005/04/27 13:28:56 clausen Exp $
//------------------------------------sim.cc-----------------------------------

#include "sim.h"
#include "exceptions.h"

#include "sim_solve.cc"
#include "sim_test.cc"
#include "sim_initial.cc"
#include "sim_file.cc"
#include "sim_stats.cc"
#include "scalar_field.cc"
#include "set_field.cc"

//--------------------------  globals variables -------------------------------

SIM              sim;
int              user_control_keys[5]  = { 0, 0, 0, 0, 0 };
bool             finish_now            = false;

//=============================================================================
//=============================================================================
//=============================================================================


void SIM::stepAll()
//elementary step for all the objects in 'this' sim
{
  ++stepCnt;

  //all lists are mixed, to ensure equal chances to all, and valid results...

  mixSolid();
  mixMicrotub();
  mixAster();
  mixComplex();
  mixGrafted();

  //specialEvent();

  stepAster();
  stepMicrotub();

  //STEP THE FIELD
  //CHAITANYA:
  
  //setstabilField();//here we set initial values for the fields of stabiliz (rescue and catastrophe frequencies)

  //calculate the motion of the Microtub & other SimObject
  if ( MP.mtmobile )
		solve();
		 
	

  //TODO: if the microtubule never move, we could do the dispatch only once:

  //prepare for Hand binding
  if ( nbFreeHands() > 0 )
    dispatchRods();

#ifdef TEST_ATTACH
  //for debug: continuously test Attachment algorithm (slow)
  if ( nbFreeHands() > 0 )
    testAttach();
#endif

  stepSolid();
  stepNucleus();
  stepGrafted();
  stepComplex();

 
  // QUick & dirty way of doing it,....shortcut for Gliding Assay ( for KJ ) , to account for stepping rate
  // of dynein in the discrete nodel. haCnt is global counter for time - it holds on stepGrafted function()
  // TO IMPROVE, for right logic and implementation - a counter as f( time , hand ) instead of just f( t ).
  
  // dthastep as implemnet as default is for Athale 2014 dynein ( hard-wired ). Has been made a
  // this loop works when
  
  /*
  if ( haCnt < 9 ){
	haCnt+= 1;
  }else{
	stepGrafted();
	//printf("Time: %.3f \t move motor \n", simTime() );
        haCnt = 0;
       }
  */	
  //hadtstepValue = sim.updatehadtstepValue();
  //printf(" real time %.3f  \t havalue %f \n",  simTime()  ,  hadtstepValue );
  //printf(" real time %.3f  \t havalue %.3f \n",  simTime()  ,  hadtstepValue );
  //hadtstepValue = sim.updatehadtstepValue();
  /*
  if ( simTime() == 0 ){
		hadtstepValue = 0.00;

		printf(" real time %i  \t havalue %.3f \n",  simTime()  ,  hadtstepValue );
		hadtstepValue =  sim.updatehadtstepValue();
		printf(" real time %i  \t havalue %.3f \n",  simTime()  ,  hadtstepValue );
   };
  */
  //hadtstepValue = sim.updatehadtstepValue();
  //if( simTime() ==  hadtstepValue ){
  //if( simTime() ==  ( sim.iterationCnt() -1 )*MP.dthastep ){
  //			printf(" Motor steps! \n");			
  //			//stepGrafted();
  //			hadtstepValue = sim.updatehadtstepValue();
 //			printf(" Motor steps! real time %.3f  \t havalue %.3f \n", simTime()  ,  sim.iterationCnt()*MP.dthastep );
 // }
		
 
  /*
  //printf("cnt: %f the updated step %f \n", MP.dt*sim.iterationCnt(), hadtstepValue );


  hadtstepValue = hadtstepValue +  MP.dthastep;
//if( MP.dt*sim.iterationCnt() == hadtstepValue ){
		
  if( sim.iterationCnt() == 1 ){
	stepGrafted();
	hadtstepValue = MP.dthastep;	
		 
  
  if( MP.dt*sim.iterationCnt() == hadtstepValue ){
			stepGrafted();
			//hadtstepValue = hadtstepValue +  MP.dthastep;
			printf(" Motor steps: iternum %i \t factor %f \t havalue %f \n", sim.iterationCnt()  , MP.dt*sim.iterationCnt() , hadtstepValue );
  }else
      printf("Wait \n");
	
  */	

  //CHAITANYA: This is where the diffision step is taken for diffusive fields.
  //setstabilField();// preferably use a method that steps and updates it.
  //STABILIZ FIELD STEP--meth still to be implemented
 // printf("ficellsize= %f\t boxsize= %f, boxshape= %d, mtdynregulation= %d\n", MP.ficellsize, MP.boxsize[0], MP.boxshape, MP.mtdynregulation);
  //stepStabField();


}




//=============================================================================
//=============================================================================


void SIM::specialEvent()
{
  /*
   //this if for iva, to make the spindle longer
   if ( MP.initcode == 50 ) {
     if ( MP.cxlength[0] < 20 )
       MP.cxlength[0] = 2 + simTime() / 10;
     }
   */
  /*
  if ((simTime() > MP.boxtime[2]) && (simTime() < MP.boxtime[3])) {
  MSG(3, "Cutting spindle...\n");
  Aster * asi = firstAster();
  Aster * as = asi;
  while( as ) {
      asi = asi->next();

      if ( as->getName() > 2 ) delete( as );

      as = asi;
  }
  */
}


//=============================================================================
//=============================================================================


void SIM::resetCounters()
{
  finish_now         = false;
  stepCnt            = 0;
  starting_time      = 0;
  frame_in_buffer    = -1;  // NO frame to begin with.
  nb_frames_to_skip  = 0;
  frame_info[0]      = '\0';
  haCnt      = 0;
}

void SIM::clearSpace()
{
  if ( Microtub::space ) delete Microtub::space; Microtub::space = 0;
  if ( Solid::space )    delete Solid::space;    Solid::space    = 0;
  if ( Nucleus::space )  delete Nucleus::space;  Nucleus::space  = 0;
  if ( Grafted::space )  delete Grafted::space;  Grafted::space  = 0;
  if ( Complex::space )  delete Complex::space;  Complex::space  = 0;
}

void SIM::setSpace()
{
  clearSpace();
  Microtub::space=newSpace(MP.boxshape, MP.boxsize, MP.boxinflate);
     Solid::space=newSpace(MP.boxshape, MP.boxsize, MP.boxinflate);
   Nucleus::setSpace();
   Grafted::space=newSpace(MP.boxshape, MP.boxsize, MP.boxinflate);
   Complex::space=newSpace(MP.boxshape, MP.boxsize, MP.boxinflate);
}

int SIM::getReady()
{
  if (0 == MP.nbValuesRead())
    return 1;
  MP.computeDerivedValues();
  resetCounters();
  setSpace();
  return setMTGrid();
}

bool SIM::isReady()
{
  //we test the Microtub's space, and the attachment grid
  return Microtub::space && MicrotubList::isMTGridReady();
}


int SIM::nbObjects() const
{
    return nbMicrotub()
    + nbSolid()
    + nbGrafted()
    + nbComplex()
    + nbAster()
    + nbNucleus();
}

int SIM::nbFreeHands() const
{
    return freeGrafteds.size()
    + 2*freeComplex.size()
    + boundComplex.size();
}


void SIM::eraseState()
{
  eraseNucleus();
  eraseAster();
  eraseGrafted();
  eraseComplex();
  eraseMicrotub();
  eraseSolid();
}

void SIM::moduloPosition()
{
  if ( Microtub::space == 0 )
    MSG.error("SIM::moduloPosition()", "Microtub::space is not set");

  if ( Microtub::space -> isPeriodic() ) {
    //moduloPositionGrafted();  //that is probably not necessary
    moduloPositionComplex();
    //moduloPositionAster();   //this is redundant with Microtub / Solid
    moduloPositionMicrotub();
    moduloPositionSolid();
    moduloPositionNucleus();
  }
}

//link() is overloaded here, but that was not possible in the ObjectLists

Microtub *  SIM::link(Microtub * x) { return linkMicrotub(x); }
Complex  *  SIM::link(Complex * x)  { return linkComplex(x); }
Grafted  *  SIM::link(Grafted * x)  { return linkGrafted(x); }
Aster    *  SIM::link(Aster * x)    { return linkAster(x); }
Solid    *  SIM::link(Solid * x)    { return linkSolid(x); }
Nucleus  *  SIM::link(Nucleus * x)  { return linkNucleus(x); }
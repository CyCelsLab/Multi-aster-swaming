//RCS: $Id: microtub_attach3.cc,v 2.3 2005/02/18 18:25:57 nedelec Exp $

//======================================================================
//===========================DISPATCH===================================
//======================================================================

// A very simple algorithm for finding attachments of Hands to Tube
// Brutal force, but simple and most probably correct!
// Use this to compare with other more astute methods

void MicrotubList::setMTGrid() {}
void MicrotubList::dispatchRods() {}
void MicrotubList::deleteMTGrid() {}
bool MicrotubList::isMTGridReady() { return true; }

//======================================================================
//===========================ATTACH=====================================
//======================================================================

//return 1 if it attached the motor, 0 if it did not
bool MicrotubList::tryToAttach( const Vecteur & place, Hand& ha )
{
  //int er, excluded;
  int targetCnt = -1;

  int p;
  Microtub * mt;
  PointMicrotub * trod;
  
  microtubs.mixWell();
  
  real attach_dist = MP.haattachdist[ ha.getType() ];

  //we test all rods:
  for( mt=firstMicrotub(); mt ; mt=mt->next() )
    for( p = 0; p < mt->nbSegments(); ++p ) {
      trod = mt->getRod( p );

      //prevent illegal attachments to be made, base on key match
      //if ( trod->keyMatch( ha.bindingKey() )) continue;

      //we now compute the distance from the hand to the rod,
      //and compare it to the maximum attach distance of the hand.
      real dist = trod->distanceToRod(place);
      
      if ( dist > attach_dist ) continue;     
      
#ifdef NEAR_ATTACH_FASTER
      
      //NEW: probability increase for closer targets, following a 1/R law
      //the random test for binding is eventually done by calling RNG
      if ( dist * RNG.preal() > MP.haattachrate_dt[ ha.getType() ] ) {
        // we can now attach the hand to the target rod:
        if ( ha.attach( *trod ))
          return true;
        else
          ha.clear();
      }
      
#else    
      
      //the random test for binding is eventually done by calling RNG
      if ( RNG.test( MP.haattachrate_dt[ ha.getType() ] )) {
        // we can now attach the hand to the target rod:
        if ( ha.attach( *trod )) 
          return true;
        else
          ha.clear();
      }
      
#endif
      
      //we limit the number of trial according to the parameter MP.hamaxtargets.
      //With the default value MP.hamaxtargets=1, we give two chances to an hand from
      //a free complex, and one to the remaining hand from an already bound complex.
      //Usually, there is only one MT rod close enough to bind to, in which
      //case the hand has only this one chance irrespective of the complex state.
      
      //initialize targetCnt according to whether the hand is in a free complex or not:
      if ( targetCnt < 0 )
        targetCnt = ha.otherBound();
      //test if the number of target considered is above the specified limit:
      if ( ++targetCnt > MP.hamaxtargets )  {             //one more target
                                                          //if above the limit, we return without binding:
        return false;
      }
    }
  return false;
}


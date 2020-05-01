//RCS: $Id: microtub_attach2.cc,v 2.2 2005/04/19 18:17:52 nedelec Exp $

//======================================================================
//=========================== UNIFORM VERSION ==========================
//======================================================================
// A set to function to implement a Monte-Carlo approach to binding the 
// motors complexes, under the assumption that free motors have an uniform
// distribution:
//
// tryToAttachUniform() keeps track of the number of free complexes.
//
// tryToAttachUniform(int cxty, int side, real mtLength, real vol, int nbFree)
// estimates how many motors should bind from the number of available motors,
// and the relevant parameters.
//
// attachUniform( int cxty, int side, real totalLength )
// finds a microtubule, and attaches hand side of a complex of type cxty



//--------------------------------------------------------------------
Microtub * MicrotubList::uniformMicrotubSite(const real totalLength, real & ab) const
{
  assert( totalLength > 0 );
  
  ab = totalLength * RNG.preal();
  
  //find the microtubule corresponding to abscisse 'ab'
  Microtub * mt = firstMicrotub();
  do {
    
    if ( ab < mt->length() ) 
      break;
    else
      ab -= mt->length();
    
    mt  = mt->next();
    if ( mt == 0 ) mt = firstMicrotub();  //that should not happen, but safer

  } while ( ab > 0 );

  return mt;
}
  
//--------------------------------------------------------------------
int MicrotubList::attachToMicrotubEverywhere( const int cx_type, const int cx_side, const real totalLength )
{
  real ab;
  Microtub * mt = uniformMicrotubSite( totalLength, ab );
  
  //try to find a free complex of the requested type in the free list:
  Complex * cx = sim.freeComplexOfType(cx_type);
  assert( cx );

  //attach the complex to the microtubule:
  switch( cx_side ) {
    case 1: cx->attachTo1( mt, ab, MT_MINUS_END ); break;
    case 2: cx->attachTo2( mt, ab, MT_MINUS_END ); break;
  }
  return 1;
}



//--------------------------------------------------------------------
int MicrotubList::attachToMicrotubEnds( const int cx_type, const int cx_side, const real endBindingLength )
{
  //find a random microtubule:
  Microtub * mt = randomMicrotub();
  
  //try to find a free complex of the requested type in the free list:
  Complex * cx = sim.freeComplexOfType(cx_type);
  assert( cx );
  
  real len = endBindingLength * RNG.preal();
  if ( len > mt->length() ) len = mt->length();
  //attach the complex to the microtubule:
  switch( cx_side ) {
    case 1: cx->attachTo1( mt, len, MT_PLUS_END ); break;
    case 2: cx->attachTo2( mt, len, MT_PLUS_END ); break;
  }
  
  return 1;
}



//--------------------------------------------------------------------

//--------------------------------------------------------------------
void MicrotubList::tryToAttachUniform(const int cx_type, const int cx_side, const real totalLength, 
                                      const real totalVolume, int & nbFree )
{
  //printf("tryToAttachUniform( %.2f cx%i ha%i)\n", esp, cxty, side );
  int ha_type = 0;
  switch( cx_side ) {
    case 1: ha_type = MP.cxhatype1[cx_type];
    case 2: ha_type = MP.cxhatype1[cx_type];
  }
  
  //estimate the capture volume:
  real captureVolume = 1; //default case for DIM==1
  if ( DIM == 3 ) captureVolume = PI * sqr(MP.haattachdist[ha_type]);
  if ( DIM == 2 ) captureVolume = 2 * MP.haattachdist[ha_type];
  
  //multiply by the length of MT on which the motor can bind:
  if ( MP.haendattachdist[ha_type] > 0 )
    captureVolume *= MP.haendattachdist[ha_type] * nbMicrotub();
  else
    captureVolume *= totalLength;
  
  //estimate the number of hands in this region that would be tested:
  real bind = nbFree * captureVolume / totalVolume;
  
  //we would like to test <bind> times for MP.haattachrate_dt[haty], but
  //we correct the probability for the integer error.
  int  nb  = (int) ceil( bind );
  if ( nb == 0 ) return;
  
  //estimate the probability for each test:
  real prob = MP.haattachrate_dt[ha_type] * bind / real( nb );
  
  //TODO: we could use a Poisson derivative from random.h
  //we will test <nb> times the probability <prob>
  while ( --nb >= 0 ) if ( RNG.test( prob ) ) {
    if ( MP.haendattachdist[ha_type] > 0 )
      nbFree -= attachToMicrotubEnds( cx_type, cx_side, MP.haendattachdist[ha_type] );
    else
      nbFree -= attachToMicrotubEverywhere( cx_type, cx_side, totalLength );
  }
}



//--------------------------------------------------------------------
void MicrotubList::tryToAttachUniform()
{
  static int virgin = 1;
  static int nbFree[ ParamSim::MAX ];
  static int nbComplex[ ParamSim::MAX ];
  static real volume;

  if ( nbMicrotub() == 0 ) return;

  int cxty;

  if ( virgin ) { 
    //we initialize the number of free hands:
    volume = Complex::space->volume();   //this is slow, so we do it only once
    for( cxty = 0; cxty < MP.MAX; ++cxty ) {
      nbFree[ cxty ]    = MP.cxmax[ cxty ];
      nbComplex[ cxty ] = MP.cxmax[ cxty ];
    }

    //we update the number of free complex:
    Complex * cxi;
    for(cxi=sim.firstBoundComplex(); cxi ; cxi=cxi->next() )
      --nbFree[ cxi->getType() ];

    for(cxi=sim.firstBridgeComplex(); cxi ; cxi=cxi->next() )
      --nbFree[ cxi->getType() ];

    //for( cxty = 0; cxty < MP.MAX; ++cxty )
      //printf("attachUniform: nbFree[%i]=%i\n", ty, sim.nbFreeComplex[ty]);
  }

  real mtLength = totalMTLength();

  for( cxty = 0; cxty < MP.MAX; ++cxty ) {
    if ( nbFree[ cxty ] > 0 )
      tryToAttachUniform( cxty, 1, mtLength, volume, nbFree[cxty] );
    if ( nbFree[ cxty ] > 0 )
      tryToAttachUniform( cxty, 2, mtLength, volume, nbFree[cxty] );
  }
}


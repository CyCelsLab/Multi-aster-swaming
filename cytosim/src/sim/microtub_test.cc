//RCS: $Id: microtub_test.cc,v 2.3 2005/01/09 00:20:04 nedelec Exp $

//============================================================================
//===                        TEST  ATTACHMENT                             ====
//============================================================================

//test that attachement has equal probability to all targets.
int MicrotubList::testAttach(Vecteur place, real dist)
{
  int somethingIsWrong = 0;

#ifdef TEST_ATTACH
  const int NB = 128;
  int p, targets = 0;
  real dis;
  Microtub * mt;
  PointMicrotub * rod = 0;  

  //go through all the rods to find which are close enough from place:
  
  for( Microtub * mt=firstMicrotub(); mt; mt=mt->next() ) 
    for( p = 0; p < mt->nbSegments(); ++p ) {
      rod = mt->getRod( p );
      dis = mt->distanceToRod( rod, place );
      if ( dis < dist ) {
        ++targets; 
        rod->attachCnt = 0;
        MSG("testAttach target %i: M%-3lx rod %2i at %.3f um\n", targets, mt->getName(), p, dis);
      }
      else
        rod->attachCnt = -1;          //this rod is too far for attach
    }
      
  MSG("testAttach: spot has %i targets\n", targets);
  

  //create a test motor of dummy type:
  int hatype = MAGIC_TYPE;
  Hand ha( HA_BASE_NONE, 0, hatype );
  MP.haattachdist[ hatype ] = dist;
  MP.haattachrate_dt[ hatype ] = 1;
  
  
  //call tryTyAttach NB times to check to which rods the motor binds:
  for( int cnt = 0; cnt < NB; ++cnt )
    if ( tryToAttach( place, ha ) ) {
      mt  = ha.getMT();
      p   = ha.getPoint1();

      if ( mt->getRod( p ) -> attachCnt == -1 ) {
      
        MSG("testAttach false positive M%lx %2i : %3i %%", mt->getName(), p, (100*rod->attachCnt)/NB );
        sim.link( new Grafted( 2, place ) );
        somethingIsWrong = 1;
      
      } else {
        
        mt->getRod( p ) -> attachCnt++;
        
        if ( false ) {
          MSG("%i tryToAttach -> ", cnt);
          MSG(" bound to M%lx at %.2f\n", mt->getName(), ha.abscissa() );
        }
      }
      
      ha.detach(DETACH_DESTRUCT);                   //in testAttach()
    } 
        
  //----------------------check the results for all rods:
  
  real expected = NB/float(targets);
  real standard = expected + 2*sqrt(NB/float(targets)) ;
  
  for( Microtub * mt=firstMicrotub(); mt; mt=mt->next() ) 
    for( int p = 0; p < mt->nbSegments(); ++p ) {
      int cnt = mt->getRod( p )->attachCnt;
      if ( cnt >= 0 ) {
        
        MSG("testAttach: M%-3lx rod %2i : %3i %3i %%", mt->getName(), p, cnt, (100*cnt)/NB );
        if ( cnt > standard ) MSG(" : unfair ?");
        if ( cnt == 0 ) {
          MSG(" : missed");
          sim.link( new Grafted( 3, place ) );
          somethingIsWrong = 2;
        }
        MSG("\n");
      }
    }
#endif
    return somethingIsWrong;
}

//------------------------------------------------------------------------------
//test that attachement has equal probability to all targets.
void MicrotubList::testAttach() 
{
  int verboseLevel = MSG.getVerboseLevel();
  Vecteur where;
  
  for( int cnt=0; cnt < 100; ++cnt ) {
    where = Microtub::space->randomPlaceInVolume();
    MSG.setVerboseLevel(0);
    if ( testAttach( where, MP.MHGL ) ) {
      MSG.setVerboseLevel(1);
      testAttach( where, MP.MHGL );
    }
  }
  
  MSG.setVerboseLevel( verboseLevel );
}

//RCS: $Id: meanlength.cc,v 2.6 2005/04/22 18:56:36 nedelec Exp $
//-----------------------------meanlength.cc----------------------------
//this program simulates microtubule dynamic, without any border effect.
//use it to estimate the mt mean-length before starting a simulation

#include <sys/time.h>
#include <sys/types.h>
#include <sys/times.h>
#include <climits>
#include <cstring>
#include "sim.h"

#include "map.h"

//===========================================================================
const int LMAX=500;
int deleted=0;
Step stepCnt;

//the function MesureMTLength is copied from <analyse1.cc>

void MesureMTLength()
{
  Microtub * mt;
  static real ml_old=0, time_old=0, mminus_old=0;
  int state[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  int max=0, min=0;
  real mtop, top[LMAX], mtip, tip[LMAX];
  real Vm, Vp, l, mminus, ml=0, mlsqr=0, var=0;
  int t, w;
  int MAX=minT(sim.nbMicrotub()/10, LMAX);
  

  for(t=0; t<LMAX; t++) {
    top[t] = 0; 
    tip[t] = MP.mtmaxlength; 
  }

  for (mt=sim.firstMicrotub(); mt ; mt=mt->next() ) {
    
      ++state[ mt->getDynamicState(MT_PLUS_END)    ];
      ++state[ mt->getDynamicState(MT_MINUS_END)+2 ];

      l = mt->length();
      mminus += mt->abscissa( MT_MINUS_END );
      ml += l;
      mlsqr += sqr(l);
      if ( l < MP.mtminlength+MP.mtrodlength ) ++min;
      if ( l > MP.mtmaxlength-MP.mtrodlength ) ++max;

      for(t=0; t<MAX; t++)
        if ( l > top[t] ) {
          for(w=MAX-1; w>t; --w) top[w]=top[w-1];
          top[t]=l;
          break;
        };

      for(t=0; t<MAX; t++)
        if ( l < tip[t] ) {
          for(w=MAX-1; w>t; --w) tip[w]=tip[w-1];
          tip[t]=l;
          break;
        }
  }
      
  mtop=0;
  mtip=0;

  if ( MAX > 0 ) {
    for(t=0; t<MAX; t++) mtop+=top[t];
    mtop/=MAX;
    for(t=0; t<MAX; t++) mtip+=tip[t];
    mtip/=MAX;
  }

  if ( sim.nbMicrotub() ) {
    ml     /= real(sim.nbMicrotub());
    mminus /= real(sim.nbMicrotub());
    mlsqr  /= real(sim.nbMicrotub());
    var     = mlsqr - sqr(ml);
    if ( var > 0 ) var = sqrt( var ); else var = 0;
  }
  
  int nbmts = sim.nbMicrotub();

  MSG("F %3i %8.2fs ", sim.frameInBuffer(), MP.dt*stepCnt );
  MSG("%4i mt ", sim.nbMicrotub());
  MSG("( %4u %4u )+G/S ", (100*state[0])/nbmts, (100*state[1])/nbmts);
  MSG("( %4u %4u )-G/S ", (100*state[2])/nbmts, (100*state[3])/nbmts);

  MSG("%5.2f +/- %.3f um", ml, var);
  MSG(" (%.2f %.2f)B/T%i ", mtip, mtop, MAX);
  MSG(" (%4i %4i)@m/M", min, max);

  if (sim.simTime() != time_old) {
    Vm = ( mminus - mminus_old ) / ( sim.simTime() - time_old );
    Vp = ( ml - ml_old ) / ( sim.simTime() - time_old ) - Vm;
    MSG("  %6.2f %6.2f  nm/s\n", 1000*Vm, 1000*Vp);
  }
  else
    MSG("\n");

  ml_old = ml;
  mminus_old = mminus;
  time_old = sim.simTime();
}


//===========================================================================




void simLength()
{
  int  frame = 0;
  real steps_per_frame = real( MP.nbiter ) / MP.record;
  Step next_record = 0; //when we should record the next frame
  
  MSG("total real time %.3f s, dt=%.2e s, %i intervals of %.3f s\n",
           MP.nbiter*MP.dt, MP.dt, MP.record, MP.dt * steps_per_frame );
    
  //we mimick cytosim, but only performing the Microtub's operations
  stepCnt = 0;
  
  do {
      
      if ( stepCnt >= next_record ) {
          
          if ( stepCnt >= (Step)MP.nbiter ) break;
          
          MesureMTLength();
          sim.incFrameInBuffer();
          
          next_record = Step( ++frame * steps_per_frame );
          if ( next_record > (Step)MP.nbiter )
              next_record = MP.nbiter;     
      }
      
      
      sim.stepMicrotub();
      ++stepCnt;

  } while ( finish_now == false );
}




//===========================================================================
//===========================================================================
//===========================================================================

void microtubTheory(int nbmt)
{
  real total = MP.nbiter*MP.dt;
  MSG("simulated real time = %.2fs\n", total);
  MSG("mtdynamic = %i\n", MP.mtdynamic[0]);
}

//==========================================================================
//==========================================================================


int main(int argc, char* argv[])
{
  MSG.shutUp();
  if (NO_ERROR != MP.parseFile(MP.datafile)) 
    MP.parseFile(DATA_OUT);
  MP.massageParameters();

  MSG("============================INITIAL STATE============================\n");

  sim.getReady();
  sim.populateState();

  MSG.setVerboseLevel(1);

  microtubTheory( sim.nbMicrotub());
  MSG("===============================================================\n");
  if ( MP.mtdynamic[0] )    simLength();

  //fclose(mlout);
  return EXIT_SUCCESS;
}




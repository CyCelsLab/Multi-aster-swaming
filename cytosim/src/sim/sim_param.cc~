//RCS: $Id: sim_param.cc,v 1.22 2005/03/31 13:41:45 foethke Exp $

#include "sim_param.h"
#include "main.h"
#include "smath.h"
#include "iomessages.h"
#include "random.h"

//Declaration of the global (unique) parameter object:
ParamSim MP;   // short for MyParam

const int ParamSim::MAX;

//An event with a rate K is modelled discretely with a time step dt.
// at every time step, a random number uniforsm in [0,1] is compared
// to K * dt. This works well if the probability of having two event
// in the same time step is small, i.e. if K * dt is small compared to 1.
// Below is a precision threshold: a warning is issued if K * dt > ACCEPTABLE_RATE
const real ACCEPTABLE_RATE = 0.25;


//calculate derived parameters from the primary ones, set by the user
//TODO: ParamSim::computeDerivedValues() should throw exceptions instead of calling MSG.error()
void ParamSim::computeDerivedValues()
{
  int ty;
  
  //compatibilityOperations() contains the compatibility codes for older versions:
  compatibilityOperations();
  
  //if ( dimension != DIM )
  //  MSG.warning(lastFileParsed(),"ignoring entry : dimension = %i", dimension);
  
  dimension = DIM;
  
  //-----------------------seed the Random Number Generator
  
  if ( randseed )
    RNG.seed( randseed ); 
  else {
    randseed = RNG.seedTimer();
    MSG(4, "Time generated random seed 0x%lx\n", randseed );
  }
  
  //--------------------------box----------------------------------
    
  boxkmratio = boxkm / km;

  if ( dt   <= 0 )  MSG.error( lastFileParsed(), "dt <= 0" );
  if ( km   <= 0 )  MSG.error( lastFileParsed(), "km <= 0" );
  if ( kT   <  0 )  MSG.error( lastFileParsed(), "kT < 0" );
  if ( visc <= 0 )  MSG.error( lastFileParsed(), "visc <= 0" );
  
  for(ty = 0; ty < 6; ++ty) {
    if ( boxtime[ ty ] < 0 ) MSG.error( lastFileParsed(), "boxtime[%i] < 0", ty);
  }
  
  //-----------------------microtubules-----------------------------
  
  if ( mtrodlength <= 0 ) 
    MSG.error( lastFileParsed(), "mtrodlength <= 0" );
  
  if ( mtdynforce < 0 )
    MSG.error( lastFileParsed(), "mtdynforce < 0" );
  
  mtdynstretch = mtdynforce / km;
  
  for(ty=0; ty<4; ++ty) {
    if ( mtdyntrans[ty] < 0 ) 
      MSG.error(lastFileParsed(), "mtdyntrans %i is negative", ty);
  }
  
  //---------------------------hands------------------------------
  
  assert( MAGIC_TYPE < MAX );
  //MAGIC_TYPE never detach:
  hamodel[ MAGIC_TYPE ]         = 0;
  haspeed[ MAGIC_TYPE ]         = 0;
  haforce[ MAGIC_TYPE ]         = km;
  hadetachmodel[ MAGIC_TYPE ]   = 0;
  hadetachrate[ MAGIC_TYPE ]    = 0;
  haenddetachrate[ MAGIC_TYPE ] = 0;
  haattachrate[ MAGIC_TYPE ]    = 0.2/dt;
  
  //Max Hand attachment Length
  MHGL = 0;
  
  for (ty=0; ty<MAX; ty++) {
    
    if ( haattachdist[ty] < 0 )
      MSG.error(lastFileParsed(),"haattachdist < 0 for hand type %i",ty);
        
    if (( ty != MAGIC_TYPE ) && ( haattachdist[ty] > MHGL ))
      MHGL = haattachdist[ty];
    
    hamaxstretch[ty]        = haforce[ty] / km;
    hamaxstretch_sq[ty]     = sqr( hamaxstretch[ty] );
    habreakstretch[ty]      = sqrt( DIM ) * hamaxstretch[ty];
    habreakstretch_sq[ty]   = sqr( habreakstretch[ty] );
    haspeed_dt[ty]          = dt * haspeed[ty];
    
    if ( haforce[ty] )
      havarspeed_dt[ty]     = fabs( dt * haspeed[ty] * km / haforce[ty] );
    
    haattachrate_dt[ty]     = integratedProbability( haattachrate[ty], dt );
    hadetachrate_dt[ty]     = integratedProbability( hadetachrate[ty], dt );
    haenddetachrate_dt[ty]  = integratedProbability( haenddetachrate[ty], dt );
    
    //if haenddetach[] has been set, set haenddetachrate[] to a very high value
    if ( haenddetach[ty] ) {
      haenddetachrate[ty]    = 10.0 / dt;
      haenddetachrate_dt[ty] = 1.0;
    }
  }  
  
  //MAGIC_TYPE attached as fast as possible:
  haattachdist[ MAGIC_TYPE ]  = 1;
  
  //if MHGL==0, the attachement will not work for mouse-grab in play:
  if (MHGL == 0) MHGL = 0.1;

  //------------------------grafted motors----------------------

  for (ty=0; ty<MAX; ty++) {
  }
  
  //------------------------motor complexes----------------------
  
  for (ty=0; ty<MAX; ty++)
    cxdiff_dt[ty]      = sqrt( 6.0 * cxdiff[ty] * dt );
  
  //we ensure that a static complex is made of two static hands:
  cxhatype1[ MAGIC_TYPE ] = MAGIC_TYPE;
  cxhatype2[ MAGIC_TYPE ] = MAGIC_TYPE;
  
  //----------------------------solids---------------------------
  
  if ( sosize == 0 )
    MSG.error(lastFileParsed(),"sosize == 0");
  
  //----------------------------asters---------------------------
  
  askmratio = askm / km;
  
  //----------------------------nuclei---------------------------
  
  for ( int ii=0; ii<3; ii++)
    nukmratio[ii] = nukm[ii] / km;
}



//===========================================================================
//===========================================================================
//===========================================================================


//TODO: all test on validity of parameter values should be automatized!
void ParamSim::issueWarnings()
{
  int ty;
  real ef;
  
  //------------------------------box-------------------------------

  if ( boxshape == 0 )
    MSG.warning( lastFileParsed(), "boxshape == %i (SHAPE_UNSET)", boxshape );
    
 //------------------------------scalar_field-----------------------
 
 if (ficellsize < 0)
     MSG.warning( lastFileParsed(),"ficellsize= %d < 0", ficellsize);
 
 
    
  //--------------------------microtubules--------------------------
  
  if (mtminlength > mtmaxlength)
    MSG.error(lastFileParsed(),"mtminlength > mtmaxlength");
  
  for( int ii = 0; ii < MAX; ++ii ) {
    
    if (mtinitlength[ii] < mtminlength)
      MSG.warning(lastFileParsed(),"mtinitlength[%i] < mtminlength", ii);
  
    if (mtinitlength[ii] > mtmaxlength) 
      MSG.warning(lastFileParsed(),"mtinitlength[%i] > mtmaxlength", ii);
  }
  
  for(ty=0; ty<4; ++ty) {

    if ( mtdyntrans[ty] * dt > ACCEPTABLE_RATE )
      MSG.error(lastFileParsed(), "mtdyntrans*dt %i is too high : reduce dt", ty);
  
    if ( mtdynspeed[ty] * dt > 0.1 )
      MSG.error(lastFileParsed(), "mtdynspeed*dt %i is too high : reduce dt", ty);
  }
  
  //the warning below are just conventions, cytosim does not care of the order
  if ( mtdynspeed[0] < 0 )
    MSG.warning(lastFileParsed(), "grow speed of minus end = mtdynspeed[0] < 0");
  
  if ( mtdynspeed[1] > 0 )
    MSG.warning(lastFileParsed(), "shrink speed of minus end = mtdynspeed[1] > 0");

  if ( mtdynspeed[2] < 0 )
    MSG.warning(lastFileParsed(), "grow speed of plus end = mtdynspeed[2] < 0");
  
  if ( mtdynspeed[3] > 0 )
    MSG.warning(lastFileParsed(), "shrink speed of plus end = mtdynspeed[3] > 0");
  
          
  //--------------------------hands---------------------------------
  
  for (ty=0; ty < MAX; ty++) 
    if ( ty != MAGIC_TYPE ) {
    
    if ( hamodel[ty] )
      MSG.warning(lastFileParsed(),"hamodel[%i]=%i: you should use hamodel=0 in most cases", ty, hamodel[ty]);
    
    //TODO: we have to see if the hand is inside a complex with cxlength >0
    if ( haattachdist[ty] > 1.01 * hamaxstretch[ty] ) {
      MSG.warning(lastFileParsed(),"km * haattachdist / haforce = %.2f > 1 for hand %i",
                haattachdist[ty]/hamaxstretch[ty],ty);
    }
    
    if ( hamaxstretch[ty] == 0 ) 
      MSG.warning(lastFileParsed(),"ha <%i> has maxstretch == 0",ty);
        
    if ( haattachrate[ty] && ( haattachrate_dt[ty] * (1+hamaxtargets) > ACCEPTABLE_RATE )) 
      MSG.warning(lastFileParsed(),"haattachrate[%i]*dt*(1+hamaxtargets) = %.2f is high: reduce dt", 
                  ty, haattachrate_dt[ty]*(1+hamaxtargets));
    
    if ( hadetachrate[ty] && ( hadetachrate_dt[ty] > ACCEPTABLE_RATE ))
      MSG.warning(lastFileParsed(),"dt * hadetachrate[%i] = %.2f is high: reduce dt", 
                  ty, hadetachrate_dt[ty]);
    
    if (( haenddetach[ty] == 0 ) && ( haenddetachrate_dt[ty] > ACCEPTABLE_RATE ))
      MSG.warning(lastFileParsed(),"dt * haenddetachrate[%i] = %.2f is high: reduce dt", 
                  ty, haenddetachrate_dt[ty] );
    
    // we compare the stretch in a link due to the equipartition theorem
    // to the maximum stretch that the motor can sustain before detaching.
    if ( sqrt( DIM * kT / km ) > hamaxstretch[ty] )
      MSG.warning(lastFileParsed(),"Thermal stretch of hand %i= %.4f > maxStretch %.4f nm : increase haforce or decrease km",
                  ty, 1000 * sqrt( DIM*kT/km ), 1000 * hamaxstretch[ty] );

    if (hadetachrate[ty]) {
      //we compare the distance traveled during the time 1/detachrate,
      //and how much force this would induce compared to haforce.
      //this is an effective efficiency of the motor. Normally pretty high
      ef = fabs( km * haspeed[ty] / ( hadetachrate[ty] * haforce[ty] ));
      if (( haspeed[ty] != 0 ) && ( ef < 2. ))
        MSG.warning(lastFileParsed(),"Low efficacy for ha %i : km * speed / ( force * offrate ) = %.2f", ty, ef);
    }
    
    //we compare the force reached after one time step at max speed, 
    //with the maximum force of the motors.
    if ( fabs( haspeed[ty] * dt * km ) > 0.5 * haforce[ty] )
      MSG.warning(lastFileParsed(),"Force achieved in one step = %.1f > maxforce[%i]/2 = %.1f : reduce dt",
                  fabs(haspeed[ty] * dt * km), ty, haforce[ty] );
    }
  //--------------------------complex--------------------------------
  for (ty=0; ty<MAX; ty++) {
    
    if ((cxhatype1[ty] < 0) || (cxhatype1[ty] >= MAX))
      MSG.error(lastFileParsed(),"wrong mohatype1 for complex <%i>",ty);
    
    if ((cxhatype2[ty] < 0) || (cxhatype2[ty] >= MAX))
      MSG.error(lastFileParsed(),"wrong mohatype2 for complex <%i>",ty);
    
    if (cxlength[ty] < 0)
      MSG.error(lastFileParsed(), "cxlength[%i] < 0", ty);
    
    if (cxlength[ty] > 0) {
      if ( haattachdist[ cxhatype1[ ty ] ] < cxlength[ ty ] )
        MSG.warning(lastFileParsed(),"cxlength[%i] > haattachdist[ hand 1 ]", ty );
      if ( haattachdist[ cxhatype2[ ty ] ] < cxlength[ ty ] )
        MSG.warning(lastFileParsed(),"cxlength[%i] > haattachdist[ hand 2 ]", ty );
      
      if ( hamaxstretch[ cxhatype1[ ty ] ] < cxlength[ ty ] )
        MSG.warning(lastFileParsed(),"cxlength[%i] > maxstretch[ hand 1 ]", ty );
      if ( hamaxstretch[ cxhatype2[ ty ] ] < cxlength[ ty ] )
        MSG.warning(lastFileParsed(),"cxlength[%i] > maxstretch[ hand 2 ]", ty );
    }
  }
  //---------------------------solids-------------------------------
  
  if ( soptsmax < 1 )
    MSG.error(lastFileParsed(), "soptsmax < 1");
  
  //---------------------------asters-------------------------------

  if (( asmax > 0 ) && ( assize[0] > mtminlength ))
    MSG.error(lastFileParsed(), "assize[0] = %.2f > mtminlength = %.2f", assize[0], mtminlength);
  
  if ( asmtmax[0] < 1 )
    MSG.error(lastFileParsed(), "asmtmax[0] < 1");
  if ( asghmax[0] < 0 )
    MSG.error(lastFileParsed(), "asghmax[0] < 0");
}

//===========================================================================
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//===========================================================================
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//===========================================================================

void ParamSim::compatibilityOperations()
//to keep compatibility with older versions
{
  //--------------- Compatibility with public Cytosim 1.0 -------------------
  // if only two values where read for mtdynspeed and mtdyntrans,
  // then we assume that these values where for the plus end, 
  // and that the minus-end is non-dynamic:
  if (( 2 == nbValuesRead("mtdynspeed")) && ( 2 == nbValuesRead("mtdyntrans"))) {
    MSG.warning(lastFileParsed(), "the values for MT dynamic instability are used for the plus-end");
    for( int ii = 1; ii >= 0; --ii ) {
      mtdynspeed[ii+2] = mtdynspeed[ii];
      mtdynspeed[ii] = 0;
      mtdyntrans[ii+2] = mtdyntrans[ii];
      mtdyntrans[ii] = 0;
    }
    showThisManyValues("mtdynspeed", 4);
    showThisManyValues("mtdyntrans", 4);
  }
}


void ParamSim::printSome()
{
  int ty;
  
  MSG(4, "\n");
  
  MSG(5, "MHGL = %.2f nm\n", 1000 * MHGL );
  
  real ptsmobdt  = HYDRO * dt / ( 4 * PI * visc * mtrodlength );
  //real ptsdiffdt = sqrt( 6.0 * kT * ptsmobdt );
  
  //this is the calibration of the link's angular rigidity:
  real mtrigidkm_typical = mtrigid / cub( mtrodlength );
  MSG(4, "Hardness of matrix = links %7.3f   MT rigidity %7.3f\n", km * ptsmobdt, mtrigidkm_typical * ptsmobdt);
  MSG(4, "Force to fold a MT at right angle = %.2f pN\n", mtrodlength * mtrigidkm_typical/2.0);
  MSG(4, "Sampling error for a MT curved at (5 um radius) = %.2f nm\n", 1000*sqr( mtrodlength ) / ( 8 * 5 ));
  MSG(4, "Stretch of a link due to Brownian energy: sqrt(%i*kT/km) = %.2f nm", DIM, 1000*sqrt( DIM * kT / km ));
  MSG(4, "  ( %.2f pN )\n", sqrt( DIM * kT * km ));
  MSG(4, "\n");
  
  //---------------------------filaments------------------------------
  
  MSG(4, "MT -grow %7.2f -shrk %7.2f   +grow %7.2f +shrink %7.2f nm\n",
      1000*mtdynspeed[0] * dt,1000*mtdynspeed[1] * dt,
      1000*mtdynspeed[2] * dt,1000*mtdynspeed[3] * dt);
  
  MSG(4, "MT -cat. %7.5f -res. %7.5f   +cat. %7.5f +res.   %7.5f\n",
      mtdyntrans[0] * dt, mtdyntrans[1] * dt, 
      mtdyntrans[2] * dt, mtdyntrans[3] * dt);
  
  //------------------------------hands-------------------------------
  MSG(4,"MAXTARGET = %i\n", hamaxtargets);
  
  
  MSG(4,"hands dx[]={");
  for(ty=0; ty<MAX; ty++) MSG(4, " %5.2f", 1000*haspeed_dt[ty]);
  MSG(4,"} nm\n");
  
  
  MSG(4,"     pON[]={");
  for(ty=0; ty<MAX; ty++) MSG(4," %5.3f", haattachrate_dt[ty]);
  MSG(4,"} / step\n");
  
  
  MSG(4,"    pOFF[]={");
  for(ty=0; ty<MAX; ty++) MSG(4," %5.3f", hadetachrate_dt[ty]);
  MSG(4,"} / step\n");
  
  
  MSG(4," pOFFend[]={");
  for(ty=0; ty<MAX; ty++) MSG(4," %5.3f", haenddetachrate_dt[ty]);
  MSG(4,"} / step\n");
  
  //------------------------motor complexes----------------------
  
  MSG(4,"Cx diffus.={");
  for (ty=0; ty<MAX; ty++) MSG(4," %5.1f", 1000*cxdiff_dt[ty]);
  MSG(4,"} nm / step\n");
}

//===========================================================================
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//===========================================================================



void ParamSim::massageParameters()
{
  if ( nbValuesRead() == 0 )
    MSG.error(lastFileParsed(),"no values set or read from parameter file");
  
  computeDerivedValues();
  issueWarnings();
}

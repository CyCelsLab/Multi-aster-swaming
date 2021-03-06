//RCS: $Id: sim_initial.cc,v 2.20 2005/04/27 13:28:56 clausen Exp $
//---------------------------siminitial.cc-----------------------------//

#include "smath.h"
#include "iowrapper.h"
#include "sim_param.h"


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ============================================================================================
//------ Neha Khetan, 2014 September
//------ populateMouse, for initialization of MTOCs as modelled for  --------
//       the centripetal movement of the asters in meiotic oocytes   --------
// ============================================================================================

// 1. MTOCs : some % can be initialized at the nuclear envelope
// 2. Immobilized motors can be initialized in a gradient
// 3. MTs initialized in a gradient
// 4. MTs initialized within the nucleus radius





//========================================================================
//============================POPULATE====================================
//========================================================================


void SIM::populateState()
{
  //--------  PART 1: input from file, controlled by MP.initfile
  switch( MP.initfile ) {

  case 0:
    break;

  case 1:
    readState(START_IN);
    frame_in_buffer = 0;
    return;   //we will not add any object to the state read from file

  case 2:
    readState(START_IN);
    frame_in_buffer = 0;
    break;

  case 3:       //this can be used to finish an interupted simulation
    printf( "system( cp sim.out sim1.out )\n");
    system( "cp sim.out sim1.out" );
    //load the last frame from result.out:
    IO.openInputFile(RESULT_OUT);
    while ( NO_ERROR == readState() );
    IO.closeInputFile();
    //reopen the file with "append" mode
    IO.openOutputFile(1);
    skipNextFrame();
    ++frame_in_buffer;
    MSG(4, "appending simulation from %.3f s", starting_time);
    return;    //we will not add any object to the state read from file

  default:
    MSG.warning("SIM::populate", "ignored parameter initfile=%i", MP.initfile);
  }


  //--------  PART 2: add random objects with the constructors

  //the shape should be set. If the assertion fails, it might be
  //that getReady() was not called.
  assert( Microtub::space->isSet() );

  try {
    switch ( digit( MP.initcode, 3 )) {

      //normal initial conditions defined below in sim_initial.cc
      case 0:
        switch ( digit( MP.initcode, 2 )) {
          case 0: populateRandom();  break;
          
          case 1: populateMouse(); break;
          // nk on 19 July: as populateXenopus is missing I'm modifying this as populateMouse
          //case 1: populateXenopus(); break;
          case 2: populatePombe();   break;
          case 3: populateElegans(); break;
          case 4: populateSpindle(); break;
          case 5: populatePombe2();  break;
          case 6: populateXenbeads(); break;
          case 7: populateChromArray(); break;
          case 8: populateChrompatt(); break;
          case 9: populateDNApatt(); break;
          default:
            MSG.warning("SIM::populate","wrong 2d digit of MP.initcode");
            break;
        }
        break;

      //test initial conditions defined in sim_test.cc
      case 1: populateTest1();   break;
      case 2: populateTest2();   break;
      case 3: populateTest3();   break;
      case 4: populateTest4();   break;
      case 5: populateTest5();   break;
      default:
        MSG.warning("SIM::populate","wrong 3d digit of MP.initcode");
        break;
    }

  } catch( StuckException e ) {
    MSG("SIM::Populate %s\n", e.getMessage() );
  }

  //we should rename all objects, to have a clean start
  renameMicrotub();
}

//========================================================================
//========================================================================
//========================================================================

// Can we ing used for Gliding assay
void SIM::populateRandom()
{
  int i, type;
  printf("populate random \n");
  // normal asters: minus-end at center
  for (i=0; i<MP.asmax; ++i)
    link( new Aster( MT_MINUS_END ));

  //inverted aster structures:
  for (i=0; i<MP.chmax; ++i)
    link( new Aster( MT_PLUS_END ));

  //nucleation is delayed if boxtime[] is set properly:

  // For MT lengths of length set by initlength
  if ( MP.boxtime[0] + MP.boxtime[1] == 0 )
    for(i=0; i<MP.mtmax; ++i)
      link( new Microtub(0) );
 /*  
 
	// Neha Khetan, 1 August 2016
	//  MTs exponentially distributed for Gliding assays 
	
	Vecteur place , dir;                     // position of the MT 
 	Microtub * mt;
	
	double length;
	double lambda;
	if ( MP.boxtime[0] + MP.boxtime[1] == 0 ) {
		for( i=0; i < MP.mtmax; ++i){
			//printf("hi! Neha, I am here!!\n");
			place = Microtub::initPosition( Microtub::space, 0 );			
			dir   = Microtub::initDirection( 0 );
			lambda  = ( 1/MP.mtlengthexpmean ) ;
			//length = -log( RNG.real_uniform_range(0, 1) )/ lambda  ; //
			//length = -log( RNG.preal())/ lambda  ; //  
			length  =  RNG.exponential(lambda);		 
			printf("%3.6lf\n", length);
  			link( new Microtub( place, dir, length ));
  			}
		}

*/






  //add grafted of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.ghmax[type]; ++i)
      link( new Grafted(type) );

#ifndef UNIFORM_FREE_COMPLEX
  //add motors of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.cxmax[type]; ++i)
      link( new Complex(type) );
#endif

  //add solids
  for(i=0; i < MP.somax; ++i)
    link( new Solid(0) );

  //add nucleus
  for(i=0; i < MP.numax; ++i)
    link( new Nucleus(0) );

  //CHAITANYA: Initialize the field here - ie set the size and call ranGTPField.readFromFile()

  setstabilField(MP.superMod);//here we set initial values for the fields of stabiliz (rescue and catastrophe frequencies)

  //printf("field param: %i\n", MP.fieldTyp[0]);
  //  FILE * cfile = fopen("mtcat.out", "w");
  //  for (int kk=0; kk < sim.catMap.nbCells(); kk++){
  //	  fprintf(cfile, "%i, %f\n", kk, sim.catMap.cell( kk ) );}
  //  fclose( cfile);
  //  printf("@microtub_list wrote file\n");

}
//========================================================================
//------------------------------------------------------------------------
//========================================================================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ============================================================================================
//------ Neha Khetan, 2014 September
//------ populateMouse, for initialization of MTOCs as modelled for  --------
//       the centripetal movement of the asters in meiotic oocytes   --------
// ============================================================================================

// 1. MTOCs : some % can be initialized at the nuclear envelope
// 2. Immobilized motors can be initialized in a gradient
// 3. MTs initialized in a gradient
// 4. MTs initialized within the nucleus radius



void SIM::populateMouse()
{
  
  int i, type;
  printf("populate Mouse!!! :) \n");
 
  // ----------------------------------------------------------------------------
  //                              MICROTUBULES
  // ----------------------------------------------------------------------------
  
 /*	
  //  Jan 2015: MT initialized in a sigmoid gradient (function same as used in case of dynein localization) 

	  Vecteur pos , dir;
	  Microtub * mt;
	  
	  FILE * mtGrad = fopen("shortMicro_Grad.txt", "w");

	  if ( MP.boxtime[0] + MP.boxtime[1] == 0 ) {
		for(i=0; i<MP.mtmax; ++i){
			if( MP.mtGradient == 1){ 		  	
	
				Vecteur posmt;
				real pRmt, pNmt; 					
		
				int boolPosmt = 0;
	
				while (boolPosmt == 0) {
				
					posmt = sim.mtNucl();			 // Retrurns the coordinates random in volume
					
	   				pRmt  = RNG.real_uniform_range(0, 1);
					pNmt  =	sim.mtMap( posmt[0],posmt[1]);		// this places  gradfted motors as a func. of nuclear distance,//			Calculating the prob. of positioning a grafte
					
					
				
					if (pRmt < pNmt ) {
							boolPosmt = 1;
					}
				}
					dir   = Microtub::initDirection( 9 );		
					link( new Microtub( posmt, dir, MP.mtinitlength[0], MT_MINUS_END));
					fprintf(mtGrad,"%2.2f\t%2.2f\t%2.2f\t%2.2f\n",posmt[0],posmt[1],(sqrt(pow(posmt[0],2) + pow(posmt[1],2))), pNmt);
			}else {
				link( new Microtub(0) );
			}
		}

	fclose(mtGrad);	  
	}



	//  MTs inside the nucleus
	
	Vecteur place , dir;                     // position of the MT 
 	Microtub * mt;
	unsigned long ouf = 0;

	if ( MP.boxtime[0] + MP.boxtime[1] == 0 ) {
		for( i=0; i < MP.mtmax; ++i){
			place = Vecteur::randSphere(MP.nuradius);
			dir   = Microtub::initDirection( 0 );	
  			link( new Microtub( place, dir, MP.mtinitlength[0], MT_MINUS_END));
  			}
		}
  */








  // ----------------------------------------------------------------------------
  //                              MOTOR COMPLEXES
  // ----------------------------------------------------------------------------


#ifndef UNIFORM_FREE_COMPLEX

  //add motors of all types
  for (type=0; type < MP.MAX; ++type)
 	 for(i=0; i<MP.cxmax[type]; ++i)
      link( new Complex(type) );

#endif

  //add solids
  for(i=0; i < MP.somax; ++i)
    link( new Solid(0) );
    
    
  //add nucleus
  Nucleus* nuc = 0;
  real nuctrans[] = { MP.initsize[2], MP.initsize[3], MP.initsize[4] };
  for(i=0; i < MP.numax; ++i) {
    nuc = link( new Nucleus(0) );
    nuctrans[0] = nuctrans[0] + i * 6.0;
    nuc->translatePosition( nuctrans );
  }  
  
 // ----------------------------------------------------------------------------------
 //                              IMMOBILIZED MOTORS ( Modified:  18 October 2014, Neha Khetan )
 // ----------------------------------------------------------------------------------
 //add grafted motors:  As a Gradient of motors as a function of space set up ghMap  

for (type=0; type < MP.MAX; ++type){

	switch ( MP.ghGradient[type] ){

		// Grafted is located uniformly in the cell		
		case 0:
		{
			for(i=0; i<MP.ghmax[type]; ++i){
				link( new Grafted(type)  );
			}
		} break;
		
		
		// Initialized in a gradient originating from the center of the cell in a shape as defined by the user
		// Implemented for Mouse multi-aster centering work
		case 1:
			{
			
			FILE * motGrad = fopen("motorGrad.txt", "w");
			for(i=0; i<MP.ghmax[type]; ++i){
				Vecteur place;			
				int boolPosgh = 0;
		
				do{

					real pRgh, pNgh; 					//random and gh-positioning probabilities					
					place = sim.ghNucl();			    // Retrurns the coordinates random in volume// 
					pRgh = RNG.real_uniform_range(0, 1);
					pNgh =	sim.ghMap( place[0] , place[1]  );		//Calculating the prob. of positioning a grafted  //pNgh =	sim.fieldModel( posGh[0],posGh[1]);
					
					if (pRgh < pNgh )
						 {
							boolPosgh = 1;
							fprintf(motGrad,"%2.2f\t%2.2f\t%2.2f\t%2.2f\n", place[0],place[1],(sqrt(pow(place[0],2) + pow(place[1],2))), pNgh);
						 }
					}while ( boolPosgh == 0);
										
					link( new Grafted( type , place ) );
					

							
				};
				fclose(motGrad);
			} break;
		
    // Motors at the cortex
		case 2:
		   {
		   printf("Case 2: I am Anchor (grafted) at the cortex for %d \n", MP.ghmax[type] );
			for(i=0; i<MP.ghmax[type]; ++i)
				{			
					Vecteur place;			
					place = Vecteur::randSphere(MP.boxsize[0]);
					link( new Grafted( type, place , GH_CORTEX));			
			    }
			} break;

	// Walking Motors at the cortex , 13 October 2017
		case 3:
			{
		  printf("Case3 : I am Motor (grafted)  at the cortex for %d \n", MP.ghmax[type] );
			for(i=0; i<MP.ghmax[type]; ++i)
				{			
					Vecteur place;			
					place = ( MP.boxsize[0] + MP.boxedge * RNG.preal() ) * Vecteur::randNormed();
					link( new Grafted( type, place ));			
			    }
			}break;
	   
}
}
	  


  
  
// Modified on 16 Sept- nk
setstabilField(MP.superMod);    //here we set initial values for the fields of stabiliz (rescue and catastrophe frequencies)
    
// ---------------------------------------------------------------------------------------
//              MTOCs or ASTERS ( Normal asters: minus end at the center) (22 July 2014)
// --------------------------------------------------------------------------------------
// ------------ perAster is the percentage of asters initialized at the NE if asinit ==25 else normal nucleation

  if (MP.asinit == 25 ){  		
    for (i=0; i<MP.asmax; ++i)
   		{
   		if ( i < (MP.asmax * MP.perAster))
        {
  				MP.asinit = 23;
  				link( new Aster( MT_MINUS_END ) );     
      				}
      			else {
      				MP.asinit = 24;
      				link( new Aster( MT_MINUS_END ) );     
  					}
    		}
    
  } else {
  for (i=0; i<MP.asmax; ++i)
      link( new Aster( MT_MINUS_END ));
  	}
  }



     


  
//------------------- Neha Khetan, signing off here....bye bye bye!! --------------

//------------------------------------------------------------------------
void SIM::populateXenopus()
{
}


//------------------------------------------------------------------------
void SIM::populatePombe()
{
  int i, type;

  //nucleation is delayed if boxtime[] is set properly:
  if ( MP.boxtime[0] + MP.boxtime[1] == 0 )
    for(i=0; i<MP.mtmax; ++i)
      link( new Microtub(0) );

  //add grafted of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.ghmax[type]; ++i)
      link( new Grafted(type) );

#ifndef UNIFORM_FREE_COMPLEX
  //add motors of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.cxmax[type]; ++i)
      link( new Complex(type) );
#endif

  //add nucleus
  Nucleus* nuc = 0;
  real nuctrans[] = { MP.initsize[2], MP.initsize[3], MP.initsize[4] };
  for(i=0; i < MP.numax; ++i) {
    nuc = link( new Nucleus(0) );
    nuctrans[0] = nuctrans[0] + i * 6.0;
    nuc->translatePosition( nuctrans );
  }
}
// to add the numnber of asters
//------------------------------------------------------------------------
void SIM::populatePombe2()
{
  Vecteur dir, pos;
  dir.setRandomSphere();

  pos = Microtub::space->randomPlaceInVolume();

  Microtub * mt1 = link( new Microtub( pos,  dir, MP.mtinitlength[0] ));
  Microtub * mt2 = link( new Microtub( pos,  dir, MP.mtinitlength[0] ));
  Microtub * mt3 = link( new Microtub( pos, -dir, MP.mtinitlength[0] ));
  Microtub * mt4 = link( new Microtub( pos, -dir, MP.mtinitlength[0] ));

  Complex * co;
  real mid = ( MP.initsize[0] + MP.initsize[1] )/2.0;
  co = link( new Complex(0) );
  co->attachTo1( mt1, mid, MT_MINUS_END );
  co->attachTo2( mt2, mid, MT_MINUS_END );
  co = link( new Complex(0) );
  co->attachTo1( mt2, mid, MT_MINUS_END );
  co->attachTo2( mt1, mid, MT_MINUS_END );

  co  = linkComplex( new Complex(MAGIC_TYPE, VZERO) );
  co->attachTo1( mt1, MP.initsize[0], MT_MINUS_END);
  co->attachTo2( mt3, MP.initsize[1], MT_MINUS_END);
  co  = linkComplex( new Complex(MAGIC_TYPE, VZERO) );
  co->attachTo1( mt3, MP.initsize[0], MT_MINUS_END);
  co->attachTo2( mt1, MP.initsize[1], MT_MINUS_END);

  co  = linkComplex( new Complex(MAGIC_TYPE, VZERO) );
  co->attachTo1( mt2, MP.initsize[0], MT_MINUS_END);
  co->attachTo2( mt4, MP.initsize[1], MT_MINUS_END);
  co  = linkComplex( new Complex(MAGIC_TYPE, VZERO) );
  co->attachTo1( mt4, MP.initsize[0], MT_MINUS_END);
  co->attachTo2( mt2, MP.initsize[1], MT_MINUS_END);
  
  
  
  
  
  


}


//------------------------------------------------------------------------
void SIM::populateElegans()
{
  MSG("populateElegans()\n");
  //Cleo: you can modify here
}

//CHAITANYA (14-06-2006)
//--- based on the idea of clusters of beads distributed in some fixed pattern
//--- Pattern is a Y-axis oriented column of beads
//----Overrides the settings of MP.numax == max. no. of nuclei
//------------------------------------------------------------------------
void SIM::populateXenbeads()
{
  int i, type;
  //int maxN;
  printf("populate xenbeads \n");
  // normal asters: minus-end at center
  for (i=0; i<MP.asmax; ++i)
    link( new Aster( MT_MINUS_END ));

  //inverted aster structures:
  for (i=0; i<MP.chmax; ++i)
    link( new Aster( MT_PLUS_END ));

  //nucleation is delayed if boxtime[] is set properly:
  if ( MP.boxtime[0] + MP.boxtime[1] == 0 )
    for(i=0; i<MP.mtmax; ++i)
      link( new Microtub(0) );

  //add grafted of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.ghmax[type]; ++i)
      link( new Grafted(type) );

#ifndef UNIFORM_FREE_COMPLEX
  //add motors of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.cxmax[type]; ++i)
      link( new Complex(type) );
#endif

  //add solids
  for(i=0; i < MP.somax; ++i)
    link( new Solid(0) );

///----
  //add nucleus
  Nucleus* nuc = 0;
  maxN = (int) floor( MP.boxsize[1]/MP.nuradius );//(2*boxsi/diameter)
  printf("@sim_init maxNucl: %i\n", maxN);
  Vecteur nuctrans((real)0.0,  -(MP.boxsize[1]+MP.nuradius), (real)0.0);
  //real nucPos[][];
  for(i=0; i < maxN; ++i) {
    nuc = link( new Nucleus(0) );
    nuctrans[1] = nuctrans[1] + (2* MP.nuradius);
    nuc->translatePosition( nuctrans );
  }
  printf("@sim_initial Max.nuclei: %i\n", nbNucleus());

  //CHAITANYA: Initialize the field here - ie set the size and call ranGTPField.readFromFile()

  setstabilField(MP.superMod);//here we set initial values for the fields of stabiliz (rescue and catastrophe frequencies)

}

///-----------------------------+++++++++++++++++
//CHAITANYA (22-06-2006)
//--- based on the idea of clusters of beads distributed in some fixed pattern
//--- the bead patterns and optimal distance are tested
//------------------------------------------------------------------------
void SIM::populateChromArray()
{
  int i, type;
  //int maxN;
  printf("populate chromArray \n");
  // normal asters: minus-end at center
  for (i=0; i<MP.asmax; ++i)
    link( new Aster( MT_MINUS_END ));

  //inverted aster structures:
  for (i=0; i<MP.chmax; ++i)
    link( new Aster( MT_PLUS_END ));

  //nucleation is delayed if boxtime[] is set properly:
  if ( MP.boxtime[0] + MP.boxtime[1] == 0 )
    for(i=0; i<MP.mtmax; ++i)
      link( new Microtub(0) );

  //add grafted of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.ghmax[type]; ++i)
      link( new Grafted(type) );

#ifndef UNIFORM_FREE_COMPLEX
  //add motors of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.cxmax[type]; ++i)
      link( new Complex(type) );
#endif

  //add solids
  for(i=0; i < MP.somax; ++i)
    link( new Solid(0) );

///----
  //add nucleus
  Nucleus* nuc = 0;
  maxN = (int) floor( 2* MP.boxsize[1]/MP.nuradius );//(2*boxsi/radius)
  printf("@sim_init maxNucl: %i\n", maxN);
  //one column of beads
  Vecteur nuctransA((real)(0.0- MP.arrayDist/2),  -(MP.boxsize[1]+MP.nuradius), (real)0.0);
  //real nucPos[][];
  for(i=0; i < floor(maxN/2); ++i) {
    nuc = link( new Nucleus(0) );
    nuctransA[1] = nuctransA[1] + (2* MP.nuradius);
    nuc->translatePosition( nuctransA );
  }
  Vecteur nuctransB((nuctransA[0] + MP.arrayDist),  -(MP.boxsize[1]+MP.nuradius), (real)0.0);
  for(i=0; i < floor(maxN/2); ++i) {
      nuc = link( new Nucleus(0) );
      nuctransB[1] = nuctransB[1] + (2* MP.nuradius);
      nuc->translatePosition( nuctransB );
  }

  printf("@sim_initial Max.nuclei: %i\n", nbNucleus());

  //CHAITANYA: Initialize the field here - ie set the size and call ranGTPField.readFromFile()

  setstabilField(MP.superMod);//here we set initial values for the fields of stabiliz (rescue and catastrophe frequencies)

}
///-----------------------------+++++++++++++++++
///-----------------------------+++++++++++++++++
//CHAITANYA (14-06-2006)
//--- based on the idea of clusters of beads distributed in some fixed pattern
//--- the bead patterns and optimal distance are tested
//------------------------------------------------------------------------
void SIM::populateChrompatt()
{
  int i, type;
  //int maxN;
  printf("populate chromArray \n");
  // normal asters: minus-end at center
  for (i=0; i<MP.asmax; ++i)
    link( new Aster( MT_MINUS_END ));

  //inverted aster structures:
  for (i=0; i<MP.chmax; ++i)
    link( new Aster( MT_PLUS_END ));

  //nucleation is delayed if boxtime[] is set properly:
  if ( MP.boxtime[0] + MP.boxtime[1] == 0 )
    for(i=0; i<MP.mtmax; ++i)
      link( new Microtub(0) );

  //add grafted of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.ghmax[type]; ++i)
      link( new Grafted(type) );

#ifndef UNIFORM_FREE_COMPLEX
  //add motors of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.cxmax[type]; ++i)
      link( new Complex(type) );
#endif

  //add solids
  for(i=0; i < MP.somax; ++i)
    link( new Solid(0) );

///----
  //add nucleus
  Nucleus* nuc = 0;
  maxN = MP.numax;
  printf("@sim_init maxNucl: %i\n", maxN);
  //one column of beads
  Vecteur nuctrans;
  if( MP.pattorient == 0 && maxN >1){
	  nuctrans[0]=(real)MP.nuradius*maxN;
	  }
  else if(MP.pattorient == 1 && maxN >1){
	  nuctrans[1]=(real)MP.nuradius*maxN;
	  }
  else
	  printf("@sim_initial.cc, Array orientation not set\n");
  //real nucPos[][];
  for(int ii=0; ii < maxN; ++ii) {
	  //for(int  x=0; x < floor(maxN/4); ++x){
		  nuc = link( new Nucleus(0) );
		  if( MP.pattorient == 0 && maxN >1){
		    nuctrans[0] = nuctrans[0] - (2* MP.nuradius);
		    }
		  else if(MP.pattorient == 1 && maxN >1){
		    nuctrans[1] = nuctrans[1] - (2* MP.nuradius);
		    }
		  nuc->translatePosition( nuctrans );
	  //}
  }

  printf("@sim_initial Max.nuclei: %i\n", nbNucleus());

  //CHAITANYA: Initialize the field here - ie set the size and call ranGTPField.readFromFile()

  setstabilField(MP.superMod);//here we set initial values for the fields of stabiliz (rescue and catastrophe frequencies)

}

///-----------------------------+++++++++++++++++
//CHAITANYA (04-02-2008)
//--- even distribution of beads clusters in a grid pattern
//--- nbNucleus decides if space filled
//
//-- NEXT STEP- read in pattern from array-design files and take
//-- coordinates automatically
//
//------------------------------------------------------------------------
void SIM::populateDNApatt()
{
  int i, type;
  //int maxN;
  printf("populate DNApatt in regular grid \n");
  // normal asters: minus-end at center
  for (i=0; i<MP.asmax; ++i)
    link( new Aster( MT_MINUS_END ));

  //inverted aster structures:
  for (i=0; i<MP.chmax; ++i)
    link( new Aster( MT_PLUS_END ));

  //nucleation is delayed if boxtime[] is set properly:
  if ( MP.boxtime[0] + MP.boxtime[1] == 0 )
    for(i=0; i<MP.mtmax; ++i)
      link( new Microtub(0) );

  //add grafted of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.ghmax[type]; ++i)
      link( new Grafted(type) );

#ifndef UNIFORM_FREE_COMPLEX
  //add motors of all types
  for (type=0; type < MP.MAX; ++type)
    for(i=0; i<MP.cxmax[type]; ++i)
      link( new Complex(type) );
#endif

  //add solids
  for(i=0; i < MP.somax; ++i)
    link( new Solid(0) );

///----
  //add nucleus
  Nucleus* nuc = 0;
  maxN = (int) pow( floor( 2*MP.boxsize[1]/( 2*MP.nuradius + MP.arrayDist ) ),2);//(box/(r+dist))^2

  printf("@sim_init maxNucl: %i\n", maxN);
  //along grid pattern-- initial value in regular square grid pattern
  int nL = (int) floor(pow( maxN, 0.5 )); //length of nucleus box
  real offsetS =  (nL/2) * (2*MP.nuradius + MP.arrayDist);

  Vecteur nuctransA; //vector to store nuclear translation
  //real nucPos[][];
  for(int i=0; i < nL; i++) {
	  for(int j=0; j < nL; j++) {
	    nuc = link( new Nucleus(0) ); //create a new nucleus object
    	nuctransA[1] = (real)offsetS - j*( 2*MP.nuradius + MP.arrayDist );
    	nuc->translatePosition( nuctransA );
	}
	nuctransA[0] = (real)(-offsetS) + i*( 2*MP.nuradius + MP.arrayDist );
  }

  printf("@sim_initial Max.nuclei: %i\n", nbNucleus());

  //CHAITANYA: Initialize the field here - ie set the size and call ranGTPField.readFromFile()

  setstabilField(MP.superMod);//here we set initial values for the fields of stabiliz (rescue and catastrophe frequencies)

}

//------------------------------------------------------------------------
void SIM::populateSpindle()
{
  Vecteur pos, dir = -VY;
  Microtub * mt;
  Grafted * gh;

  //the width of the spindle is about 5 micro-meters wide
  real width = 2.5 / (MP.mtmax-1);

  //set-up the Kinetochore fibers
  for( int ii=0; ii<MP.mtmax; ++ii) {
    int pp = (2*(ii%2)-1)*ii+ii%2; //this is: center, right, left, right, left, etc.
    pos.set( width * pp, 0, 0);

    //setup a microtubule pointing up and down:
    for( int d = -1; d < MP.asmax; d += 2 ) {
      dir = d * VY;
      mt = link( new Microtub( pos, dir, MP.mtinitlength[0], MT_PLUS_END ));
      mt->setDynamicState( MT_PLUS_END, MT_PAUSE );
      mt->setType(MAGIC_TYPE);  //mark them with a special type
      pos = mt->whereEnd( MT_PLUS_END );
      gh = link( new Grafted( MAGIC_TYPE, pos));
      gh->attachTo( mt, 0, MT_PLUS_END );
      gh = link( new Grafted( MAGIC_TYPE, pos+dir));
      gh->attachTo( mt, 1, MT_PLUS_END );
    }
  }

#ifndef UNIFORM_FREE_COMPLEX
//add motors of all types
for (int type=0; type < MP.MAX; ++type)
  for(int ii=0; ii < MP.cxmax[type]; ++ii)
    link( new Complex(type) );
#endif

// normal asters: minus-end at center
for (int ii=0; ii<MP.asmax; ++ii)
  link( new Aster( MT_MINUS_END ));

}

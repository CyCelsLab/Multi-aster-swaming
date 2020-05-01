//RCS: $Id: microtub_list.cc,v 2.12 2005/04/19 18:18:00 nedelec Exp $
//---------------------------microtublist.cc--------------------------------
#include "sim.h"
#include "iowrapper.h"

#include "smath.h"


#ifdef SLOW_SAFE_ATTACH
  #include "microtub_attach3.cc"
#else
  #include "microtub_attach1.cc"
#endif

#ifdef UNIFORM_FREE_COMPLEX
  #include "microtub_attach2.cc"
#endif

#include "microtub_test.cc"

//--------------------------------------------------------------------
Microtub * MicrotubList::randomMicrotub() const
{
  //find a random Microtub:
  Microtub * mt = firstMicrotub();
  for(int ii = RNG.pint_exc( nbMicrotub() ); ii > 0; --ii )
    mt = mt->next();
  assert( mt );
  return mt;
}



//-------------------------------------------------------------------------
//--------------------------------MT STEPS --------------------------------
//-------------------------------------------------------------------------

//you can set the Dynamic state of microbutules according to their
//position in space, or other factors. (CHAITANYA)
void MicrotubList::forceDynamicState()
{
  switch( MP.mtdynregulation ) {
  case 0:
    break;

  // microtubule grow inside a central disc
  case 1:
    for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() ) {
      real pen = mti-> whereEnd( MT_PLUS_END ).norm();
      if ( pen > MP.boxsize[5] )
        mti -> setDynamicState( MT_PLUS_END, MT_SHRINKING );
    }
    break;

  //microtubule grow inside and shrinks outside
  //CHAITANYA: the state in which the self-assembly of chromatin mediated spindles happens
  case 2:
    for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
      if ( mti -> getDynamicState( MT_PLUS_END ) != MT_PAUSE ) {
        if ( Microtub::space->isInside( mti -> whereEnd( MT_PLUS_END )) )
          mti -> setDynamicState( MT_PLUS_END, MT_GROWING );
        else
          mti -> setDynamicState( MT_PLUS_END, MT_SHRINKING );
      }
    break;

  //microtubule grow outside (surface) and shrinks inside
  case 3:
    for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
      if ( mti -> getDynamicState( MT_PLUS_END ) != MT_PAUSE ) {
        if ( Microtub::space->isInside( mti -> whereEnd( MT_PLUS_END )) )
          mti -> setDynamicState( MT_PLUS_END, MT_SHRINKING );
        else
          mti -> setDynamicState( MT_PLUS_END, MT_GROWING );
      }
    break;

  }
}


//-------------------------------------------------------------------------
void MicrotubList::nucleate()
{
  //-----------------nucleation to keep a constant number of Tubes:

  // boxtime[0] is the initial time delay, before any nucleation starts
  if ( sim.simTime() < MP.boxtime[0] )
    return;

  int desiredNbOfTubes = MP.mtmax;

  // between boxtime[0] and boxtime[0]+boxtime[1], there is a linear increase:
  if ( sim.simTime() < MP.boxtime[0] + MP.boxtime[1] ) {
    desiredNbOfTubes = int( desiredNbOfTubes * ( sim.simTime() - MP.boxtime[0] ) / MP.boxtime[1]);
  }

  // calculate the number of Microtub that need to be nucleated:
  desiredNbOfTubes -= nbFreeMicrotub();

  // if no nucleation is needed, we exit
  if ( desiredNbOfTubes <= 0 ) return;

  // report if this is going to be the first time we nucleate
  static int virgin=1;
  if ( virgin ) {
    virgin = 0;
    MSG(11, "MT nucleation starts...\n");
  }

  MSG(11, "MT nucleation: total = %i, to be nucleated %i\n", nbMicrotub(), desiredNbOfTubes );

  // nucleation of new Microtub is done according to MP.mtnucleation
  for (; desiredNbOfTubes > 0; --desiredNbOfTubes )
    switch( MP.mtnucleation ) {

      case 1:       //-----------------nucleation according to MP.mtinit:
        linkMicrotub( new Microtub(1) );
        break;

      case 2: {     //-----------------nucleation on the complex:
        //since the list are mixed continuously, firstComplex should be a random one
        Complex * cx = sim.firstFreeComplex();
        if ( cx )
          linkMicrotub( new Microtub( cx->calculatePosition() ));

      } break;

      case 3: {     //-----------------nucleation on the solids:
        //since the list are mixed continuously, firstSolid should be a random one
        Solid * so = sim.firstSolid();
        if ( so ) {
          Microtub * mt = new Microtub( so->calculatePosition(),
                                        Vecteur::randNormed(),
                                        Microtub::initLength(), MT_MINUS_END );
          linkMicrotub( mt );
        }
      } break;

      case 4: {     //-----------------nucleation on the free points of solids:

        Solid * so = sim.firstSolid();
        if ( so && so->nbGrafteds() ) {
          Grafted * gh = so->getGrafted( RNG.pint_exc( so->nbGrafteds() ));
          if ( gh && gh->isFree() ) {
            Microtub * mt = new Microtub( so->calculatePosition(),
                                          Vecteur::randNormed(),
                                          Microtub::initLength(), MT_MINUS_END );
            linkMicrotub( mt );
            gh->attachTo( mt, 0.1, MT_MINUS_END );
          }
        }
      } break;

      //CHAITANYA
      case 5: {     //-----------------nucleation on the free points of solids:
        Vecteur vect;
      	Nucleus * nu = sim.firstNucleus();
      	for( Nucleus* nu = sim.firstNucleus(); nu; nu=nu->next() ) {
			//the position of the nucleus
			//for( int n = 1; n < floor(sim.maxN / (1+MP.mtmax)); n++){
				vect = nu->getPosition();
				Microtub * mt = new Microtub( vect, Vecteur::randNormed(),
							Microtub::initLength(), MT_MINUS_END );
				linkMicrotub( mt );
			//}
			nu++;
		 }

	   } break;
	   

	   // Neha Khetan, 11 March 2015, 
	   // Nucleation of MTs in a MT gradient ( The sigmoid gradient as used in the MTOC motility for gh initialization)
		case 6: {
	   			Vecteur pos , dir;
  				
	   			if( MP.mtGradient == 1){ 		  	
					
					Vecteur posmt;
					real pRmt, pNmt; 					//random and gh-positioning probabilities
		
				    int boolPosmt = 0;
					while (boolPosmt == 0) {
				
						posmt = sim.mtNucl();			 // Retrurns the coordinates random in volume
		   			    dir = Microtub::initDirection( 9  );		
						pRmt = RNG.real_uniform_range(0, 1);
						pNmt =	sim.mtMap( posmt[0],posmt[1]);		// this places  gradfted motors as a func. of nuclear distance,//Calculating the prob. of positioning a grafte
				
						if (pRmt < pNmt ) {
							boolPosmt = 1;
						}
					}
					    
                   Microtub * mt = new Microtub( posmt, dir, MP.mtinitlength[0], MT_MINUS_END );
		           linkMicrotub( mt );
                 }
                 	
		    			
				}break; 
				
				
        // Neha Khetan, 4 Jan 2016, 
        // Nucleation of MTs within the nucleus            
		case 7: {
				Vecteur place , dir;                     // position of the MT , direction of the MT
						
				place = Vecteur::randSphere( MP.nuradius );
				dir   = Microtub::initDirection( 0 );	
				Microtub * mt = new Microtub( place, dir, MP.mtinitlength[0], MT_MINUS_END );
		        linkMicrotub( mt );
			  	
			  	
			  	}break;
				
      default:      
      case 0:
        break;
    }
}

//-------------------------------------------------------------------------

//elementary steps for all the microtubules in 'this' box
void MicrotubList::stepMicrotub()
{

  //CHAITANYA: termination criterion.
  countMTcapt( MP.mtcaptMod ); //CHAITANYA: ADD THIS parameter to the data.in and param.h files. This is only a temp fix
  //nucleate new MTs if needed
  nucleate();

  //impose the dynamic state (Growing/Shrinking) if needed
  if ( MP.mtdynregulation )
    forceDynamicState();

  //calculate the free tubulin monomer concentration
  totalMTLength();

  if ( Microtub::tubulin < 0 )  //you may want that tubulin concentration is positive
    MSG.warning("MicrotubList::stepMicrotub","negative tubulin concentration\n");

  //call step() for every Microtub
  Microtub * mti = firstMicrotub();
  Microtub * mt  = mti;
  while( mt )	{
    mti = mti->next();
    mt -> step();
    mt = mti;
  }

  //CHAITANYA- evaluates the simple assymetry right and left mass
  //countMTlenghtdyn();
  //getMTtotlen();

  ////
/*  if(sim.simTime() == 5*MP.dt){
  	    //printf("res.size %i, cat.size %i \n", sim.resMap.nbCells(), sim.catMap.nbCells());

  	    FILE * map1 = fopen("resmap.out", "w");
  	    sim.resMap.write(map1, 0);
  	    fclose( map1 );

  	    FILE * map2 = fopen("catmap.out", "w");
  	    sim.catMap.write(map2, 0);
  	    fclose( map2 );
	    printf( "@microtub_list: wrote cat & res maps!\n");
     }
     */
  ////



}



//=========================================================================
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
//=========================================================================


int MicrotubList::nbFreeMicrotub() const
{
  int result = 0;
  for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
    result += ( mti->getOrganizer() == 0 );
  return result;
}

//-------------------------------------------------------------------------
real MicrotubList::totalMTLength()
{
  real result = 0;

  for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
    result += mti->length();

  //we update the (normalized) monomer concentration:
  Microtub::tubulin = 1.0 - result / MP.mtmonlength;

  return result;
}

//-------------------------------------------------------------------------
real MicrotubList::meanMTLength()
{
  int nbmt = microtubs.size();

  if (nbmt)
    return totalMTLength() / (real) nbmt;
  else
    return 0;
}



//CHAITANYA: monitor steric interaction of microtubules with nucleus
//ver. 2008-05-13
// 09/02/2006 NOTE: Need to modify this to involve actual capture in terms of steric interaction
// Right now this function only looks for the position of the +tips of the MTs. IT ALSO ASSUMES all MTs are aster nucleated!!!
int MicrotubList::countMTcapt( int captTyp )
{
	//Param captTyp
	// captTyp=0 NO CAPTURE
	// captTyp = 1, eucleidean dist betw MT plus-tip and chromatin-edge, SIMULATION ENDS if captured
	// captTyp = 2, (13/05/2007) euclidean dist, SHRINKS- sets vg,vs,fcat,fres==0.
	// captTyp = 3, (12/05/2008) euclidean dist, PAUSES- sets vg,vs,fcat,fres==0.

	// TO DO: add a probabilistic capture and release rng>Pcapt, PAUSE, else Catastrophe (collision based shrinkage)
	int nbmtcapt=0;
	    switch( captTyp ) {

		  case 0:
		    break;

		  case 1:{
		  //based on euclidean distance between EDGE of nucleus and location of MT-tip (+ end)
		  real minDist = MP.nuradius;
		    for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
				{
					real distNuc = 0.0;
					real ptip[ DIM ];

				    if (DIM == 1){
					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
					   distNuc = pow( ( pow((MP.nuVect[0] - ptip[0]), 2)) , 0.5);
				    }

					else if ( DIM == 2){
					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
					   ptip[1] = mti->whereEnd(MT_PLUS_END).YY;
					   distNuc = pow( (pow( (MP.nuVect[0] - ptip[0]), 2) + pow((MP.nuVect[1] - ptip[1]),2) ), 0.5);
				    }
				    else if( DIM == 3 ){
					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
					   ptip[1] = mti->whereEnd(MT_PLUS_END).YY;
					   ptip[2] = mti->whereEnd(MT_PLUS_END).ZZ;
					   distNuc = pow( (pow((MP.nuVect[0] - ptip[0]),2) + pow((MP.nuVect[1] - ptip[1]),2) + pow((MP.nuVect[2] - ptip[2]),2) ), 0.5);

				   }
					//finding if the microtubule is within the radius of the nucleus/bead
					//distance between +tip and center of bead


					distNuc = pow( ( pow( (MP.nuVect[0] - ptip[0]), 2) + pow( (MP.nuVect[1] - ptip[1]), 2) + pow( (MP.nuVect[2] - ptip[2]), 2) ), 0.5);
					//printf("distNuc = %f, nucRadius = %f \n", distNuc, MP.nuradius );

					if ( distNuc <= MP.nuradius ){
					   nbmtcapt = nbmtcapt + 1;
					   if (minDist > distNuc)
					   		minDist = distNuc;
					//   printf("distNuc = %f, nucRadius = %f, Nr.Capt MTs= %i \n", distNuc, MP.nuradius, nbmtcapt );
				    }
				}
			if(nbmtcapt >= 1){
			  finish_now = true;//CHAI: Once the chromosomes are "reached", capture is automatic (for now)
			  printf("Time %f, Nr.Capt MTs= %i, minDist= %f \n", sim.simTime(), nbmtcapt, minDist );
	  		}
		  }//case1
		  break;

		  case 2:{
			  //SHRINKS ON CONTACT
		  		  //based on euclidean distance between EDGE of nucleus and location of MT-tip (+ end)
		  		  real minDist = MP.nuradius;
		  		    for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
		  				{
		  					real distNuc = 0.0;
		  					real ptip[ DIM ];

		  				    if (DIM == 1){
		  					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
		  					   distNuc = pow( ( pow((MP.nuVect[0] - ptip[0]), 2)) , 0.5);
		  				    }

		  					else if ( DIM == 2){
		  					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
		  					   ptip[1] = mti->whereEnd(MT_PLUS_END).YY;
		  					   distNuc = pow( (pow( (MP.nuVect[0] - ptip[0]), 2) + pow((MP.nuVect[1] - ptip[1]),2) ), 0.5);
		  				    }
		  				    else if( DIM == 3 ){
		  					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
		  					   ptip[1] = mti->whereEnd(MT_PLUS_END).YY;
		  					   ptip[2] = mti->whereEnd(MT_PLUS_END).ZZ;
		  					   distNuc = pow( (pow((MP.nuVect[0] - ptip[0]),2) + pow((MP.nuVect[1] - ptip[1]),2) + pow((MP.nuVect[2] - ptip[2]),2) ), 0.5);

		  				   }
		  				    //finding if the microtubule is within the radius of the nucleus/bead
		  				   	//distance between +tip and center of bead
   	  					    //distNuc = pow( ( pow( (MP.nuVect[0] - ptip[0]), 2) + pow( (MP.nuVect[1] - ptip[1]), 2) + pow( (MP.nuVect[2] - ptip[2]), 2) ), 0.5);
		  					//printf("distNuc = %f, nucRadius = %f \n", distNuc, MP.nuradius );

		  					if ( distNuc <= MP.nuradius ){
		  					   nbmtcapt = nbmtcapt + 1;
		  					   //printf("distNuc = %f, nucRadius = %f, Nr.Capt MTs= %i \n", distNuc, MP.nuradius, nbmtcapt );
							   //-- 15-11-2007 CHAITANYA-- here we set the MT state to SHRINKING!!
							   mti->setDynamicState(MT_PLUS_END, MT_SHRINKING);
							   //12-jun-2008 Chaitanya
		  				    }
		  				}
		  			//END of ALL MICROTUBULE LIST
		  			//if(nbmtcapt >= 1){
		  			//  finish_now = true;//CHAI: Once the chromosomes are "reached", capture is automatic (for now)
		  			//printf("Time %f, Nr.Capt MTs= %i\n", sim.simTime(), nbmtcapt );
		  	  		//}
		  }//case2
		  break;

		  case 3:{
			  //PAUSES ON CONTACT
			 //based on euclidean distance between EDGE of nucleus and location of MT-tip (+ end)
	 		 real minDist = MP.nuradius;
	 	  	 for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() ){
				 real distNuc = 0.0;
		 		 real ptip[ DIM ];

			    if (DIM == 1){
  					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
  					   distNuc = pow( ( pow((MP.nuVect[0] - ptip[0]), 2)) , 0.5);
  				    }
				else if ( DIM == 2){
  					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
  					   ptip[1] = mti->whereEnd(MT_PLUS_END).YY;
  					   distNuc = pow( (pow( (MP.nuVect[0] - ptip[0]), 2) + pow((MP.nuVect[1] - ptip[1]),2) ), 0.5);
  				    }
			    else if( DIM == 3 ){
  					   ptip[0] = mti->whereEnd(MT_PLUS_END).XX;
  					   ptip[1] = mti->whereEnd(MT_PLUS_END).YY;
  					   ptip[2] = mti->whereEnd(MT_PLUS_END).ZZ;
  					   distNuc = pow( (pow((MP.nuVect[0] - ptip[0]),2) + pow((MP.nuVect[1] - ptip[1]),2) + pow((MP.nuVect[2] - ptip[2]),2) ), 0.5);
  				   }
		 		  				    //finding if the microtubule is within the radius of the nucleus/bead
		 		  				   	//distance between +tip and center of bead
		    	  					    //distNuc = pow( ( pow( (MP.nuVect[0] - ptip[0]), 2) + pow( (MP.nuVect[1] - ptip[1]), 2) + pow( (MP.nuVect[2] - ptip[2]), 2) ), 0.5);
		 		  					//printf("distNuc = %f, nucRadius = %f \n", distNuc, MP.nuradius );

				if ( distNuc <= MP.nuradius ){
  					   nbmtcapt = nbmtcapt + 1;
  					   //printf("distNuc = %f, nucRadius = %f, Nr.Capt MTs= %i \n", distNuc, MP.nuradius, nbmtcapt );
					   //--CHAITANYA-- here we set the MT state to PAUSE
					   //12-jun-2008 Chaitanya Athale
					   mti->setDynamicState(MT_PLUS_END, MT_PAUSE);

  				    }
  				}
  			//END of ALL MICROTUBULE LIST
  			//if(nbmtcapt >= 1){
  			//  finish_now = true;//CHAI: Once the chromosomes are "reached", capture is automatic (for now)
  			//printf("Time %f, Nr.Capt MTs= %i\n", sim.simTime(), nbmtcapt );
  	  		//}
		  }//case3
		  break;

		}////switch
	return nbmtcapt;
}



//-------------------------------------------------------------------------
Microtub * MicrotubList::linkMicrotub( Microtub * mt )
{
  MSG(70, "Microtub::link   mt%lx ( total %i )\n", mt->getName(), nbMicrotub());
  microtubs.pushFront( mt );
  return mt;
}


//-------------------------------------------------------------------------
Microtub * MicrotubList::findMicrotub( Name n, int createIfNotFound )
  //local storage of name list in a buffer for performance
{
  MSG(90,"findMicrotub( M%lx  create %i )\n", n, createIfNotFound );
  Microtub * mt = static_cast<Microtub*>(mtNameList.nodeWithThisName( n ));
  if (( mt == 0 ) && n )
    if ( createIfNotFound ) {
      mt = new Microtub();
      if ( mt == 0 ) {
        fprintf(stderr, "MicrotubList::findMicrotub() memory allocation failed\n");
        exit(1);
      }
      mt -> setName( n );
      linkMicrotub( mt );
    }
  //else MSG.warning("MicrotubList::findMicrotub", "cannot find M%lx", n);
  return mt;
}


//-------------------------------------------------------------------------
//read the name, and arrange the 'translation' to microtubule pointer
Microtub * MicrotubList::readMicrotubName() {

  Name name = IO.readUInt16();
  if ( name == 0 ) return 0;
  return findMicrotub( name, 0 );

}

//CHAITANYA: Print the length of each microtubule and its +end location at every time point. This data will be used for evaluation
void MicrotubList::countMTlenghtdyn()
{
	//Initialize the max and min values for a left-right gradient along X-axis
	//printf("aster.X= %f \n", MP.asVect[0]);
	float maxLenP = 0, maxLenN=0, totLenP=0, totLenN=0, xMax=0;
    //int numMTcapt =countMTcapt();
    //printf("numAtChrom= %i\n",MP.numAtChr);
	//printf("Time, %f\n", sim.simTime() );
	 for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() ){
		 real x = mti->whereEnd(MT_PLUS_END).XX;
		 //printf("%lx, %f, %f\n",  mti->getName(),  mti->whereEnd(MT_PLUS_END).XX, mti->length());
		 if (mti->whereEnd(MT_PLUS_END).XX <= MP.asVect[0]){ // negative X-axis
		     totLenN= totLenN + mti->length();
			 if(mti->length() > maxLenN)
			    maxLenN = mti->length();
		 }
	     else if (mti->whereEnd(MT_PLUS_END).XX > MP.asVect[0]){
			 totLenP= totLenP + mti->length();
			 if(mti->length() > maxLenP)
			    maxLenP = mti->length();
		 }
		 if(x > xMax)
		    xMax = x;

	     //printf("negMaxLen= %f, posMaxLen= %f\n",  maxLenN,  maxLenP);
	     }
	 printf("%f,  %f,  %f,  %f,  %f\n", sim.simTime(), totLenP, totLenN, maxLenP, maxLenN);
	 // Chromosomes are at a dist MP.xdist2Chrom from the clamped aster at (0,0) in the positive direction along X-axis
	 // terminates simulation if X-coord of MTs in positive direction== dist2Chrom
	 //if (MP.numAtChr== 1 && xMax >= MP.xdist2Chrom){
	//	printf("Chromos. reached at time[s]= %f, xMax= %f, maxLenP= %f, maxLenN= %f\n", sim.simTime(),  xMax, maxLenP, maxLenN);
	 //   finish_now = true;//CHAI: Once the chromosomes are "reached", capture is automatic (for now)
	 //}
	 //else
	 /*
	 if(numMTcapt >=MP.numAtChr)
	 {
		 printf("Chromos. reached at time[s]= %f, xMax= %f, maxLenP= %f, maxLenN= %f, numMTcapt= %i\n", sim.simTime(),  xMax, maxLenP, maxLenN, numMTcapt);
		 finish_now = true;//CHAI: Once the chromosomes are "reached", capture is automatic (for now)
	 }
	 */
}

////////////-------------------------------
//CHAITANYA: Print the length of each microtubule and its +end location at every time point. This data will be used for evaluation
void MicrotubList::getMTtotlen()
{
	float maxLenP = 0, maxLenN=0, totLenP=0, totLenN=0, xMax=0;

    for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() ){
		 real x = mti->whereEnd(MT_PLUS_END).XX;
		 //printf("%lx, %f, %f\n",  mti->getName(),  mti->whereEnd(MT_PLUS_END).XX, mti->length());
		 if (mti->whereEnd(MT_PLUS_END).XX <= MP.asVect[0]){ // negative X-axis
		     totLenN= totLenN + mti->length();
			 if(mti->length() > maxLenN)
			    maxLenN = mti->length();
		 }
	     else if (mti->whereEnd(MT_PLUS_END).XX > MP.asVect[0]){
			 totLenP= totLenP + mti->length();
			 if(mti->length() > maxLenP)
			    maxLenP = mti->length();
		 }
		 if(x > xMax)
		    xMax = x;
		 }
	//printf("Time: %f, %f, %f, %f, %f\n", sim.simTime());
}



///////------------------------------------


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//----------------------------   WRITE  -------------------------------------
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
void MicrotubList::moduloPositionMicrotub()
{
  for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
    mti->moduloPosition();
}

void MicrotubList::writeMicrotub()
{
  IO.writeRecordTag("MT");
  IO.writeIntAscii(nbMicrotub(), mtNameList.nextName() );
  for( Microtub * mti=firstMicrotub(); mti; mti=mti->next() )
    mti->write();
}



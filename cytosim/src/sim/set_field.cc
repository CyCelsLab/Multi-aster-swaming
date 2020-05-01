//RCS: $Id: sim_initial.cc,v 2.20 2005/04/27 13:28:56 clausen Exp $
//---------------------------set_field.cc-----------------------------//
// Neha, modified on 07 - 03 - 2016 :to incorporate the gradient as in RWD model
// -->  last modified on 13 dec 2014: to incorporate motor gradient without any normalization
// 		case 2 of ghModelType is for sigmoid shape of motor density w/o normalization

//---------------------------------------------- Neha, 18 october 2014...map for grafted motors

#include "smath.h"
#include "iowrapper.h"

//========================================================================
//============================set_field.cc====================================
//========================================================================

/**
This set of functions sets the field of choice based on the init_code that decides
how many nuclei and where they are to be placed
*/







//====================================================================================================================
//---------------------------------------------- Neha, 18 october 2014...map for grafted motors
// To initialze the dynein motors 
Vecteur SIM::ghNucl(){
	Vecteur posGh;
	posGh = Grafted::space->randomPlaceInVolume();
	return posGh;
}     

//=================================================================================================================================
// To create a motor distribution function

real SIM:: ghMap(real xx, real yy){
	
	real rVal, rValn, pVal;						// coordinates of dynein random in volume in cell
	rVal = sqrt(pow(xx,2) + pow(yy,2));			// now the rVal will range from edge of chrom. to the periphery 

	
	real ymax, ymin, gradVal;					//variables for max, min (normalization values), gradient Value


	switch(MP.ghModelType){
    
    	case 0:{
		// Gradient extends from cell upto the chromatin edge, Max at the chromatin edge, the max. values continues inside the nucleus.

		    if (rVal < (MP.nuradius)  ){
				pVal  = 1;
				}
    		else{
        		rValn = rVal - MP.nuradius;
        		
        		
        		real pmax  = 1/( 1 + exp(( 0    -  MP.ghGradInfl)/ MP.ghGradSlope ));
				real pmin  = 1/( 1 + exp((( MP.boxsize[0]-rValn)- MP.ghGradInfl)/MP.ghGradSlope));
				real p     = ( 1/( 1 + exp(( rValn - MP.ghGradInfl)/MP.ghGradSlope) ) );
    			
    			pVal= (p- pmin)/(pmax-pmin);
        		
        		}
        	return pVal;
        		

    			}break;

   	  	// min - max normalized wherein the gradient range from chromatin periphery to the cell periphery
   	  	// used in BJ 2015 (Feb) and PLoS comp bio (Nov. 2015) submission
   	  	case 1:{
    			real pmax = 1/(1+ exp((0-MP.ghGradInfl)/MP.ghGradSlope));
				real pmin = 1/(1+ exp((MP.boxsize[0]-MP.ghGradInfl)/MP.ghGradSlope));
				real p = 1/(1+ exp((rVal- MP.ghGradInfl)/MP.ghGradSlope));
    			
    			pVal= (p- pmin)/(pmax-pmin);
    			
				return pVal;
				}
	   			 break;
	    
		// modifying on 13 dec 2014: modfying the case 1 such that no normalization is done.
	    case 2:{	
    			
				real pVal = 1/(1+ exp((rVal- MP.ghGradInfl)/MP.ghGradSlope));
				return pVal;
				}
	    break; 
	    
	    
	    	
	    
	    default:{
	    	printf("Hi! Something is wrong!! Invalid motor gradient model\n");  
	    	}
	    break;
	    }

}



//===========================================     End of Neha, 19 october 2014 for motor Gradient





//====================================================================================================================
//---------------------------------------------- Neha, 18 october 2014...map for grafted motors



// To initialze the MT filaments

Vecteur SIM::mtNucl(){
	Vecteur posmt;
	posmt = Microtub::space->randomPlaceInVolume();
	return posmt;
}     

//=================================================================================================================================
// To create a motor distribution function

real SIM:: mtMap(real xx, real yy){
	
	real rVal, rValn, pVal;						// coordinates of dynein random in volume in cell
	rVal = sqrt(pow(xx,2) + pow(yy,2));			// now the rVal will range from edge of chrom. to the periphery 

	
	real ymax, ymin, gradVal;					//variables for max, min (normalization values), gradient Value


	switch(MP.mtModelType){
    
    
   // ------------------------- TO DO!!!  
    	///NEED to modify as in cse of motor and MT gradient on 7 - march 2016
    	case 0:{
		// Gradient extends from cell upto the chromatin edge, Max at the chromatin edge, the max. values continues inside the nucleus.
		//if the location is within the chromatin region, set pNgh == 1
		    if (rVal < (MP.nuradius)  ){
				pVal  = 1;
				}
    		else{
        		rValn = rVal - MP.nuradius;
    			pVal = (1/( 1 + exp(( rValn - MP.mtGradInfl)/MP.mtGradSlope) ) );//gradient 0-1 (peak at centre)
    			}
		return pVal;
		}
	    break;

   	  	
   	  	case 1:{	
    			real pmax = 1/(1+ exp((0-MP.mtGradInfl)/MP.mtGradSlope));
				real pmin = 1/(1+ exp((MP.boxsize[0]-MP.mtGradInfl)/MP.mtGradSlope));
				real p = 1/(1+ exp((rVal- MP.mtGradInfl)/MP.mtGradSlope));
    			
    			pVal= (p- pmin)/(pmax-pmin);
    			
				return pVal;
				}
	   			 break;
	    
		// modifying on 13 dec 2014: modfying the case 1 such that no normalization is done.
	    case 2:{	
    			
				real pVal = 1/(1+ exp((rVal- MP.mtGradInfl)/MP.mtGradSlope));
				return pVal;
				}
	    break; 
	    
	    default:{
	    	printf("Hi! Something is wrong!! Invalid motor gradient model\n");  
	    	}
	    break;
	    }

}































//------------------------------------------------------------------------
 //CHAITANYA: populates the maps for rescue and catastrophe gradients
void SIM::setstabilField(int superMode){
  //TESTING ANGLE CALCULATION- used in sim_stats.cc #7
  //float angl = atan2( 5.0, 0.0 ) *180/M_PI;
  //printf("ANGLE x0,y5= %f\n", angl);
  Vecteur bounds = Microtub::space->getBoundingRect();
  //  MSG("set3stabilField()\n");
  //printf("bounds[0] %f \n", bounds[0] );
  //--- for single nucleus
  //Vecteur nucV;
  //nucV = Vecteur(MP.nuVect);
  //printf("nucl pos vector nucV[0]: %f, bound-X: %f, boxsiz.x: %f\n", nucV[1], bounds[1], MP.boxsize[1]);
  //this function simply tests if map.h is being correctly called
  //catMap.test();
  Vecteur bd_l, bd_r, cellC;
  int num_cells[DIM];

  for (int ii=0; ii<DIM; ii++){
	  bd_l[ ii ] = -bounds[0];
	  bd_r[ ii ] =  bounds[0];
	  num_cells[ ii ] = int( round( MP.boxsize[ii]/MP.ficellsize ) );
	  //printf("%f\n", bd_l[ii]);
  }

  distMap.setDimensions( bd_l, bd_r, num_cells );
  distMap.setConfined();
  distMap.allocateCellArray();
  distMap.clear();

  catMap.setDimensions( bd_l, bd_r, num_cells );
  catMap.setConfined();
  catMap.allocateCellArray();
  catMap.clear();

  resMap.setDimensions( bd_l, bd_r, num_cells );
  resMap.setConfined();
  resMap.allocateCellArray();
  resMap.clear();

  int maxcells=catMap.nbCells();

  real xc, yc, zc, nuDist, xDist, epos[DIM];
  Vecteur v, vect;
  int nn=0;

  //FILE * catfil = fopen("catmap.out", "w");
  //FILE * resfil = fopen("resmap.out", "w");
  //FILE * vfile = fopen("dynVals.out", "w");

  //-----------
  //for Multiple NUCLEI

  //printf("max grad exp: %f, exp-ext: %f\n ",fieldModel( 1, 100, 100 ),fieldModel( 4, 100,100 ) );
  real ndist; //minimal distance to the nuclei
  int nnu=0; //nucleus number
  for( Nucleus* nu = firstNucleus(); nu; nu=nu->next() ) {
	  vect = nu->getPosition(); //position of the nucleus
	  //printf("@set_field Nuc: %i, X: %f, Y: %f\n", nnu, vect[0], vect[1]);
  	  //-----------------------------------------
  	  //------------SCAN EACH FIELD CELL AND SET VALUES of DISTANCE-MAP ------------
  	  //-----------------------------------------

  	  for( int ii=0; ii < maxcells; ii++){
		    ndist=0;
  			v = distMap.positionFromIndex( ii, 0 );

  			//Correcting for the offset in the grid
  			//v[0] = v[0] - MP.ficellsize/2;
  			//v[1] = v[1] + MP.ficellsize/2;
			for( int d=0; d<DIM; d++){
				ndist += pow( (v[d]-vect[d]), 2) ;
			}
  		    ndist = pow( ( ndist ), 0.5);
  		    //if( v[0] == -25 && nnu == 0 ){
			//	printf("grid_x: %f, y: %f, d: %f\n", v[0], v[1], ndist );
  			//}

  		    if( superMode == 0){
			//minimal distance method
				if( distMap[ ii ] == 0 || distMap[ ii ] > ndist ){
					distMap[ ii ] = ndist;
				}else{
					distMap[ ii ] = distMap[ ii ];
				}

			}else if( superMode == 1){
				//summation of distance method
				distMap[ ii ] = distMap[ ii ] + ndist;
			}
			else if( superMode == 2){
				//average of distance method
				if(distMap[ ii ]==0){distMap[ ii ] = ndist;}
				else{distMap[ ii ] = (distMap[ ii ]+ndist)/2;}
			}
			//if( ii == 0 )
	        //    printf("dist at ii=0: %f\n", distMap[ ii ] );
  		}
  	    nnu++;
	}//all nuclei scanned

  real maxV = distMap.getMaxval();
  real minV = distMap.getMinval();
  //printf("abs_max: %f, abs_min Dist: %f\n", maxV, minV);
  for( int m=0; m< maxcells; m++){


	 if(maxV > 0 ){
		 //08-aug-2007 scale simply with max distance
		 distMap[ m ] = distMap[ m ]/maxV;
		 //distMap[ m ] = (distMap[ m ]-minV)/(maxV-minV);
	 }
	 else{
		distMap[m] = distMap[m]/1;
	 }

	 //08-aug-2007
	 //scale to boxsize
	 //distMap[ m ] = distMap[ m ]/(2*MP.boxsize[0]);

  }
  real maxVnorm = distMap.getMaxval();
  real minVnorm = distMap.getMinval();
  //printf("new_max: %f, min: %f\n", maxVnorm, minVnorm );

  //printf("maxdist: %f, nNucl: %i, maxCells: %i\n", maxV, nnu, maxcells);
  int catBin=0;

  //
  //---- Using DISTANCE MAP, SET VALUES FOR CATASTROPHE AND RESCUE MAPS -----
  //

  for( int q=0; q< maxcells; q++){
	  ////For catMap we need to invert the gradient!!! - by inverting min-max values
	  if( MP.fieldTyp[0] == 0){
			MP.mtdynamic[0]   =     1;
			//printf("No gradient\n");
	  }
	  else if( MP.fieldTyp[0] == 1 ||  MP.fieldTyp[0] == 2 ||MP.fieldTyp[0] == 3 ||  MP.fieldTyp[0] == 4 ||  MP.fieldTyp[0] == 6   ||  MP.fieldTyp[0] == 7 ||  MP.fieldTyp[0] == 8 ||  MP.fieldTyp[0] == 9 ||  MP.fieldTyp[0] == 10 ||  MP.fieldTyp[0] == 11 || MP.fieldTyp[0] == 12){
		  //1: Exponential, 2: Linear gradient, 3: exponential-with-power (0.1); 6: Verhulst growth curve; 7: sigmoid; 8:optimized step (12mu) + exp(-x/10); 9: mitotic inside DNA, stabilized outside
		  // 10: Hill reaction (approx of PDE solution)
		  //catMap[ q ] = fieldModel( MP.fieldTyp[0], distMap[ q ], maxV );
	   	  //resMap[ q ] = fieldModel( MP.fieldTyp[0], distMap[ q ], maxV );
		  catMap[ q ] = (( MP.mtcatagradient[1] - MP.mtcatagradient[0] ) * fieldModel( MP.fieldTyp[0], distMap[ q ] , maxV)) + MP.mtcatagradient[0];
		  resMap[ q ] = (( MP.mtresgradient[1] - MP.mtresgradient[0] ) * fieldModel( MP.fieldTyp[0], distMap[ q ] , maxV)) + MP.mtresgradient[0];


		  MP.mtdynamic[0]    =    5;
	  }
	  else if( MP.fieldTyp[0] == 5 ){
		  //RADIAL Step gradient
		  catMap[ q ] = fieldModel( MP.fieldTyp[0], distMap[ q ] , maxV) + MP.mtcatagradient[1];
		  resMap[ q ] = MP.mtresgradient[0];
		  MP.mtdynamic[0]    =    5;
	  }

   	  else{
		throw IOException("MP.fieldTyp not implemented yet in SIM::sim_initial");
	  }
   }
   //printf("catmax: %f, catmin: %f\n", catMap.getMaxval(), catMap.getMinval() );

  //catMap=distMap;
  //resMap=distMap;

  //fclose( vfile);
  //printf("@microtub_list wrote vals near nucleus\n");

  //fclose( catfil );
  //fclose( resfil );
  //printf("res.size %i, cat.size %i \n", resMap.nbCells(),catMap.nbCells());
if( MP.fieldTyp[1] == 1){
	FILE * map1 = fopen("resmap.out", "w");
	resMap.write(map1, 0);
	fclose( map1 );

	FILE * map2 = fopen("catmap.out", "w");
	catMap.write(map2, 0);
	fclose( map2 );
	printf( "wrote cat & res maps!\n");
	}



}


//------------------------------------------------------------------------
 //CHAITANYA (22-10-2008): provides the different models of gradient-- based on distance from nucleus and model type- RETURN transition rate
 //the value is used thus setup the MAP.CATASTROPHE & MAP.RESCUE
 //RETURN VALUE is the 0-1 normalized value in the gradient.
 //THE FUNCTION THAT CALLs it, e.g. stabilizField will in turn scale it to the appropriate values of rescue and cat.
real SIM::fieldModel( int model, real dNucl, real maxDist ){
	//@@
	//@@ dNucl is a scaled (0 to 1) distance from the edge
	//@@ model is the case to be chosen of gradient model
	//@@ maxDist is the threshold distance in the case of the step function
	//@@
	//@@ grad is the value returned from each of these models
	//@@
	//step: at critical distance val_1 else val_2
	//linear: between certain values val_1 and val_2
	//exponential: using Maiwen's approximation of RanGTP-ImpBeta distribution

	real grad, dr;
				//% Analytical approximation of (RanGTP)-(Importin-beta) complex by Caudron et.al. (2005) based on
				//% Brown and Kholodenko (1999).
	real relNRad = (MP.stepDist + MP.nuradius) / maxDist;//relative distance of nucl edge + step


	switch(model){
	  case 0:{
		  grad=0;
		  return grad;
	  }
	    break;

   	  case 1:{
	    	//Caudron et al. 2005 Science. Exponential distribution, scaled to max and min rescue and catastrophe frequencies
			//% Analytical approximation of (RanGTP)-(Importin-beta) complex by Caudron et.al. (2005) based on
			//% Brown and Kholodenko (1999).

			//-- parameters
			real D2 = 13,             // um^2/s
				koff2 = 4.5 / pow(10, 4), // s^(-1)
				kon3 = 0.3,               // uM^(-1) s^(-1)
				RB_0 = 0.35,              // uM  //CHAITANYA: how many molecules does that come down to in a small vol? Variability? MC-Gillespie approach reqd.?
				RanBP1=1.2,               // uM Total Ran-Bind Prot 1 conc.
				L = MP.boxsize[0],                  // um Limit of cytoplasm
				lamda, A1, A2, RB_r,
				ll = MP.nuradius;     //um length from center of chromatin to edge

			//--- CORRECTION CHAITANYA: CHANGES: 25-02-2008
		    //dr = dNucl  - (MP.stepDist + MP.nuradius)/maxDist;
		    //printf("@fieldModel, dist-rescaled: %f, total-dist: %f\n", dNucl, (2 * MP.boxsize[0]));
			if( dNucl <= relNRad  ){
				grad = 1;
				}
			else{
				dr = ll + (dNucl * maxDist ) - (MP.stepDist + MP.nuradius ) + MP.ficellsize;

   				lamda = sqrt( ( koff2 + ( kon3 * RanBP1 ) )/ D2 );// %variable for calculation
		    	A1 = ( RB_0 * ll * exp( lamda * ll ) )/( 1 + ( ( (lamda*L + 1) * exp(2 * lamda * (ll-L)) )/( (lamda*L)  -1 ) ) );
		    	A2 = (( lamda*L + 1 )* RB_0 * ll * exp( lamda * ( ll - 2*L ) ))/( ( lamda*L - 1) + ((lamda*L + 1 ) * exp( 2*lamda*(ll-L))) );
		    	//the gradient
		    	grad = (( A1 * exp(-lamda * dr)) + (A2 * exp( lamda * dr)))/dr; // %at dist (r)
				grad = grad/(RB_0);//scaled betw 0 and 1
			}
			/*
				dr = dNucl*maxDist - (MP.nuradius + MP.stepDist);
				grad = exp( (1 - dr)/1 );//TESTING a simple exponential gradient with max/2 at dNucl=7.5
			}*/
			return grad;
			} break;

	    /*
	    case 2:{
		//linear gradient scaled linearly to distance from edge of nucleus
			//if ( dNucl <= (MP.stepDist + MP.nuradius)/maxDist )
			//	dr =0;
			//else
			//dr = dNucl  - (MP.stepDist + MP.nuradius) / maxDist;
		    //grad = dr;
		    //--- CORRECTION CHAITANYA: 2007-08-01
		    dr = dNucl  - (MP.stepDist + MP.nuradius)/maxDist;
		    if ( dr < 0 ){
				grad = 1;//max. value if within the nucleus
			   }
		    else
		    	grad = 1-dr;

			return grad;
			} break;
         
	    case 3:{
  		  //STEP gradient: returns 0 or 1: used to get the MP.catgradient[] value in function setStabilField
  		  //dNucl = (dNucl*MP.boxsize[0]);// - MP.stepDist - MP.nuradius;

		  if ( dNucl < (MP.nuradius + MP.stepDist)/maxDist ){
		  	grad = 1;//max. value if within the nucleus
			}
		  else{
		    grad = 0;
		    }

 		  return grad;
		} break;
		*/
   	    case 4:{
	    	//Modified gradient of Kholodenko by power factor.

	    	//-- parameters
			real D2 = 13,             // um^2/s
				koff2 = 4.5 / pow(10, 4), // s^(-1)
				kon3 = 0.3,               // uM^(-1) s^(-1)
				RB_0 = 0.35,              // uM  //CHAITANYA: how many molecules does that come down to in a small vol? Variability? MC-Gillespie approach reqd.?
				RanBP1=1.2,               // uM Total Ran-Bind Prot 1 conc.
				L = MP.boxsize[0],                  // um Limit of cytoplasm
				lamda, A1, A2, RB_r, dr,
				ll =MP.nuradius;     //um length from center of chromatin to edge

			dNucl = (dNucl *  2 * MP.boxsize[0]);// - MP.stepDist - MP.nuradius;
			//printf("@fieldModel, dist-rescaled: %f, total-dist: %f\n", dNucl, (2 * MP.boxsize[0]));
			if (dNucl <= (MP.nuradius + MP.stepDist) ){
			   dr = MP.nuradius + MP.stepDist;   }
			else if(dNucl > (MP.nuradius + MP.stepDist)){
			   dr = dNucl - MP.stepDist; }

   			lamda = sqrt( ( koff2 + ( kon3 * RanBP1 ) )/ D2 );// %variable for calculation
		    A1 = ( RB_0 * ll * exp( lamda * ll ) )/( 1 + ( ( (lamda*L + 1) * exp(2 * lamda * (ll-L)) )/( (lamda*L)  -1 ) ) );
		    A2 = (( lamda*L + 1 )* RB_0 * ll * exp( lamda * ( ll - 2*L ) ))/( ( lamda*L - 1) + ((lamda*L + 1 ) * exp( 2*lamda*(ll-L))) );
		    //the gradient
		    grad = (( A1 * exp(-lamda * dr)) + (A2 * exp( lamda * dr)))/dr; // %at dist (r)
			grad = pow( (grad/(RB_0)), 0.1 );//scaled betw 0 and 1 AND power of 0.1
			return grad;
			} break;

   	    case 5:{
	    	//Data based on experimental fitting of RanImp lifetime to mean aster lengths via RCC1 (Caudron et al 2005).
			//dNucl = (dNucl*2*MP.boxsize[0]);// - MP.stepDist - MP.nuradius;
			//printf("@fieldModel, dist-rescaled: %f, total-dist: %f\n", dNucl, (2 * MP.boxsize[0]));
			//if (dNucl <= (MP.nuradius + MP.stepDist) ){ //CHANGES: 25-02-2008
			if( dNucl <= relNRad ){
			   dr = MP.nuradius + MP.stepDist;   }
			else{
			   dr = maxDist * (dNucl - relNRad); }

			grad = (0.0964*pow(dr,2.6782) )/(pow(132.1616, 2.6782) + pow(dr, 2.6782) );
			return grad;
		} break;

		case 6:{
		  	//Gradient based on - INVERSE- Verhulst's Growth with carrying capacity limit. It takes the form P(t)=K*Po*exp(r*t)/(K+P0*(exp(rt)-1))
			//dNucl = (dNucl*MP.boxsize[0]);// - MP.stepDist - MP.nuradius;
			//Parameters of the Verhulst or Logistic Growth Eqn.
			real K=MP.sigmaCap, rr=MP.sigmaGro, Po=1;//K=carrying capacity, rr=growth rate, Po=intial population size
			//printf("@fieldModel, dist-rescaled: %f, total-dist: %f\n", dNucl, (2 * MP.boxsize[0]));

			if (dNucl*maxDist <= (MP.nuradius + MP.stepDist )  ){
			   dr = 1;
			   grad = 1;//max. value if within the nucleus
			   }
			else{// if(dNucl > (MP.nuradius + MP.stepDist)/maxDist ){
			   dr = (dNucl*maxDist) - MP.nuradius - MP.stepDist;
			   //scale to 0:100
			   //grad = 1 - (Po * exp( rr * dNucl ))/( K + Po * (exp( rr * dNucl ) - 1 ) );// Po <= grad <= K and
			   grad = 1 - (Po * exp( rr * dr ))/( K + Po * (exp( rr * dr ) - 1 ) );// Po <= grad <= K and

			   // 1 - (term) -- leads to decay not growth
			   }
			return grad;
		} break;


		case 7:{
				  	//Sigmoid gradient- 1-Parameter=inflection point sigmoid=1/( 1+exp(-x) )
				  	//Modified sigmoid gradient: 1/( 1 + exp(-( x - inflPt )) ).

					if (dNucl*maxDist <= (MP.nuradius + MP.stepDist )  ){
					   dr = 1;
					   grad = 1;//max. value if within the nucleus
					   }
					else{// if(dNucl > (MP.nuradius + MP.stepDist)/maxDist ){
					   dr = (dNucl*maxDist) - MP.nuradius - MP.stepDist;
					   //scale to 0:100
					   //grad = 1 - (Po * exp( rr * dNucl ))/( K + Po * (exp( rr * dNucl ) - 1 ) );// Po <= grad <= K and
					   grad = 1 - ( 1/( 1 + exp(-dr)));// Po <= grad <= K and

					   // 1 - (term) -- leads to decay not growth
					   }
					return grad;
		} break;

		//--11-aug-2008-----CHAITANYA---- (modified over Nov-2007- including exponential fraction)
		case 8:{
		  	//Gradient based on ultrasensitive gradient- step + exponential -- APPROX. HILL REAC OUTPUT
			real cDist = MP.stepDist; //critical distance fitted to the ultrasenstive gradient output
   			real dNucl_r = dNucl * maxDist;//absolute value of dist from nucleus

   			if ( dNucl < (MP.nuradius + cDist)/maxDist ){
					grad = 1;//max. value if within the nucleus & critical-dist
 			   }
			else{
				dr = dNucl_r - cDist - MP.nuradius + MP.ficellsize;//ABSOL DIST scaled for separation from surface
				grad = exp(- dr/MP.expFrac );//exponential decay outside
			   }

			return grad;

		} break;

		case 9:{
			real dNucl_r = dNucl * maxDist;//absolute value of dist from nucleus
			//Simple exponential gradient and WITHIN CHROMOSOMES MITOTIC STABILIZATION VALUES
			if ( dNucl <= MP.nuradius/maxDist  ){
				grad = 0;
			}
			else{
				dr = dNucl_r - MP.nuradius + MP.ficellsize;//dist from edge + grid-size (to prevent zero-vals)
				grad = exp(- dr/10 );//exponential decay outside
			}

			return grad;

		} break;
		case 10:{
			//CHAITANYA: 2008-Apr-15- HILL REACTION- 2 parameters, theta and nHill
			//Gradient based on ultrasensitive gradient- step + exponential -- APPROX. HILL REAC OUTPUT
			/* MATLAB code:
			theta = constVals(1);
			nHill = constVals(2);
			y = (xdata.^nHill)./( theta^nHill + xdata.^nHill );
			F= 1-y;
			*/

			real theta = MP.thetaHill; //critical distance fitted to the ultrasenstive gradient output
			real nHill = MP.nHill; //critical distance fitted to the ultrasenstive gradient output

			real dNucl_r = dNucl * maxDist;//absolute value of dist from nucleus
			if ( dNucl <= ( MP.nuradius/maxDist )){
				grad = 1;//max. value if within the nucleus & critical-dist
			    }
			else{
				dr = dNucl_r - MP.nuradius;
				grad = 1 - ( pow( dr,nHill )/( pow(theta,nHill) + pow(dr,nHill) ));//Hill reaction decay (inverted) outside

			}
			return grad;
		} break;

		/* Polynomial produces wierd curves! Remove!!
		case 11:{
			// CHAITANYA: 2008-Aug-08- polynomial of arbitrary degree
			// Gradient based on fit using gradientApprox and MATLAB function polyval
			// here we have been using a polynomial of order 6
			// Create Sim.param variable polyCoeff[] of size 7, with coefficients for
			// x^0 to x^6
			//MP.polyCoeff; //Array of 7 elements for X-vals at powers 6-0
            real dNucl_r = dNucl * maxDist;//absolute value of dist from nucleus

            //printf("polyval[0] %f\n", MP.polyCoeff[0]);

            grad = 0;



			if ( dNucl <= ( MP.nuradius/maxDist )){
				grad = 1;//max. value if within the nucleus & critical-dist
			    }
			else{
				dr = dNucl_r - MP.nuradius;

				for( int nn=0; nn < 7; nn++){
					grad =+ MP.polyCoeff[nn] * pow( dr, nn );//Hill reaction decay (inverted) outside
				}



					//The polynomial gradient
			}
			return grad;
		} break;

		case 12:{
					// CHAITANYA: 2008-Aug-08- polynomial of arbitrary degree
					// Gradient based on fit using gradientApprox and MATLAB function polyval
					// here we have been using a polynomial of order 6
					// Create Sim.param variable polyCoeff[] of size 7, with coefficients for
					// x^0 to x^6
					//MP.polyCoeff; //Array of 7 elements for X-vals at powers 6-0
		            real dNucl_r = dNucl * maxDist;//absolute value of dist from nucleus

		            //printf("polyval[0] %f\n", MP.polyCoeff[0]);

		            grad = 0;



					if ( dNucl <= ( MP.nuradius/maxDist )){
						grad = 1;//max. value if within the nucleus & critical-dist
					    }
					else{
						dr = dNucl_r - MP.nuradius;

						for( int nn=0; nn < 7; nn++){
							grad =+ MP.polyCoeff[nn] * pow( dr, nn );//Hill reaction decay (inverted) outside
						}


							//The polynomial gradient
					}
					return grad;
		} break;
		*/
		case 11:{
				  	//MODIFIED SIGMOID GRADIENT ( crit-dist, base-value, slope )
				    // 1/( 1/slope + exp(-( x - inflPt )) ).

				    real inflPt = MP.sigVals[0];
					real baseVal = MP.sigVals[1];
					real slopeVal = MP.sigVals[2];

					if (dNucl*maxDist <= (MP.nuradius + MP.stepDist )  ){
					   dr = 1;
					   grad = 1;//max. value if within the nucleus
					   }
					else{// if(dNucl > (MP.nuradius + MP.stepDist)/maxDist ){
					   dr = (dNucl*maxDist) - MP.nuradius - MP.stepDist;
					   //scale to 0:100
					   // SIMPLE SIGMOID: grad = 1 - ( 1/( 1 + exp(-dr)));// Po <= grad <= K and
					   grad = 1 - (1./( 1/baseVal + exp(-( dr - inflPt )/slopeVal) ) );

							   // 1 - (term) -- leads to decay not growth
							   }
							return grad;
		} break;

		case 12:{
				  	//4-parameter SIGMOID GRADIENT ( crit-dist, base-value, slope, peak-val )
				    // 1/( 1/slope + exp(-( x - inflPt )) ).

				    real inflPt = MP.sigVals[0];
					real baseVal = MP.sigVals[1];
					real slopeVal = MP.sigVals[2];
					real peakVal = MP.sigVals[3];

					if (dNucl*maxDist <= (MP.nuradius)  ){
					   dr = 0;//1
					   grad = peakVal - (1./( 1/baseVal + exp(-( dr - inflPt )/slopeVal) ) );
					   //peakVal;//1;//max. value if within the nucleus
					   }
					else{// if(dNucl > (MP.nuradius + MP.stepDist)/maxDist ){
					   dr = (dNucl*maxDist) - MP.nuradius;
					   //scale to 0:100
					   // SIMPLE SIGMOID: grad = 1 - ( 1/( 1 + exp(-dr)));// Po <= grad <= K and
					   grad = peakVal - (1/( 1 + exp(-( dr - inflPt )/slopeVal) ) );

							   // 1 - (term) -- leads to decay not growth
							   }
				return grad;
		} break;

        case 2:{
			  	//CA/NK: 2014-10-31: logistic function inv. to go from nucl centre outwards, normalized to min-max (0-1) 
                //All parameters are derived from ghMAP parameters!! this OVERRIDES the mtgrad field type parameters... simply calling
                //fieldTyp   2 is enough
                
                // commenting out on : 03- 07 - 2016
                /*
  			    real rVal;
  			    rVal =   dNucl * maxDist;//abs dist from the centre
				real pmax = 1/(1+ exp((0-MP.ghGradInfl)/MP.ghGradSlope));
				real pmin = 1/(1+ exp((MP.boxsize[0]-MP.ghGradInfl)/MP.ghGradSlope));
				//neha's model: also used up to localize graftedmotors (ghmap)
				real p = 1/(1+ exp((rVal- MP.ghGradInfl)/MP.ghGradSlope));
    			
    			grad= (p- pmin)/(pmax-pmin);
    			*/
    			
    			real rVal;
    			rVal =   dNucl * maxDist;  //abs dist from the centre
    			if (rVal <= MP.nuradius ){
					grad  = 1 ;
					}
	    		else{

        		real rValn;
        		rValn = rVal - MP.nuradius;
        		
        		real pmax = 1/(1+ exp((0 - MP.ghGradInfl)/MP.ghGradSlope));
				real pmin = 1/(1+ exp(( ( MP.boxsize[0]- MP.nuradius) -MP.ghGradInfl)/MP.ghGradSlope));
				real p = 1/(1+ exp((rValn - MP.ghGradInfl)/MP.ghGradSlope));
    			
    			grad= (p- pmin)/(pmax-pmin);
    			};
    			
							   
			   return grad;
		} break;
		
		
		
		case 3:{
			  	
  			    // modifying on 13 dec 2014: modfying the case 2 such that no normalization is done.
  			    
  			    real rVal;
  			    rVal =   dNucl * maxDist;//abs dist from the centre
				real grad = 1/(1+ exp((rVal- MP.ghGradInfl)/MP.ghGradSlope));			   
			   return grad;
		} break;
		
		
		
	  default:
        MSG.warning("SIM::populate", "ignored parameter initfile=%i", MP.initfile);
      }
    return grad;

	}

//



//

/*
void SIM::scaleGradient(real nuDist){
		//------SET THE GRADIENT
		///** Nuclear distance is distance to center of nucleus! scale with the radius!
		if( MP.fieldTyp == 0){
			MP.mtdynamic[0]   =     1;
		    //printf("No gradient\n");
			}
		else if( MP.fieldTyp == 1 ||  MP.fieldTyp == 2 ||  MP.fieldTyp == 10 ||  MP.fieldTyp == 11 ){
			//RADIAL, Exponential or Linear gradient
			catMap[ ii ] = ( ( MP.mtcatagradient[1] - MP.mtcatagradient[0] ) * fieldModel( MP.fieldTyp, nuDist ) )  + MP.mtcatagradient[0];
			resMap[ ii ] = ( ( MP.mtresgradient[1] - MP.mtresgradient[0] ) * fieldModel( MP.fieldTyp, nuDist ) ) + MP.mtresgradient[0];
			MP.mtdynamic[0]    =    5;
			}
		else if( MP.fieldTyp == 3 ){
			//RADIAL Step gradient
			catMap[ ii ] = MP.mtcatagradient[ (int) fieldModel( 3, nuDist ) ];
			resMap[ ii ] = MP.mtresgradient[ (int) fieldModel( 3, nuDist ) ];
			MP.mtdynamic[0]    =    5;
			}
		else if( MP.fieldTyp == 4 ){
			//BILATERAL exponential gradient
			catMap[ ii ] = ( ( MP.mtcatagradient[1] - MP.mtcatagradient[0] ) * fieldModel( 10, xDist ) )  + MP.mtcatagradient[0];
			resMap[ ii ] = ( ( MP.mtresgradient[1] - MP.mtresgradient[0] ) * fieldModel( 10, xDist ) ) + MP.mtresgradient[0];
			MP.mtdynamic[0]    =    5;
			}
		else if( MP.fieldTyp == 5 ){
			//BILATERAL linear gradient-- WITH STEP LENGTH TO SCALE
			catMap[ ii ] = ( ( MP.mtcatagradient[1] - MP.mtcatagradient[0] ) * fieldModel( 2, xDist ) )  + MP.mtcatagradient[0];
			resMap[ ii ] = ( ( MP.mtresgradient[1] - MP.mtresgradient[0] ) * fieldModel( 2, xDist ) ) + MP.mtresgradient[0];
			MP.mtdynamic[0]    =    5;
			}
		else if( MP.fieldTyp == 6){
			//BILATERAL step gradient
			catMap[ ii ] = ( ( MP.mtcatagradient[1] - MP.mtcatagradient[0] ) * fieldModel( 3, xDist ) )  + MP.mtcatagradient[0];
			resMap[ ii ] = ( ( MP.mtresgradient[1] - MP.mtresgradient[0] ) * fieldModel( 3, xDist ) ) + MP.mtresgradient[0];
			MP.mtdynamic[0]    =    5;
			}
		else{
			throw IOException("MP.fieldTyp not implemented yet in SIM::sim_initial");
			}
}
*/


//-------------
// in order to find the minium, maximum and mean: returns a Vecteur (to avoid using STANDARD TEMPLATE LIBRARY FUNCTIONS)
// {MIN, MAX, MEAN} from a 1-D array
Vecteur SIM::findMinmax( real arr[], int aSize ){
	Vecteur minmax;
	real meanVal;
	minmax[0] = arr[0];
	minmax[1] = arr[0];
	for (int m=0; m < aSize; m++){
		//minimum
		if (minmax[0] > arr[ m ]){
			minmax[0] = arr[ m ];
		}
		//maximum
		if (minmax[1] < arr[ m ]){
			minmax[1] = arr[ m ];
		}
		meanVal +=arr[m];
	}
	meanVal= meanVal/aSize;
	minmax[2] = meanVal;
	return minmax;
}
//--------------------


//write out the selected field
void SIM::writeMap()
{
  //write all 3 gradients

  IO.writeRecordTag("mp");
  //printf("@IO.writeRecordTag('mp'); \n");
  resMap.write();
  catMap.write();
  distMap.write();



}


//read in the written map information, using same methods as
//setstabilField(int superMode)

void SIM::readMap(){

	Vecteur bounds = Microtub::space->getBoundingRect();//bounds
	Vecteur b_l, b_r, cellP;//left, right and cell-pos vectors
  	int n_cells[DIM], max_cells[3], indx[3];//maximum no. of cells, index counter

	// read as is written in result.out
	for (int ii=0; ii<DIM; ii++){
		b_l[ ii ] = -bounds[0];
		b_r[ ii ] =  bounds[0];
		n_cells[ ii ] = int( round( MP.boxsize[ii]/MP.ficellsize ) );
		//printf("%f\n", b_l[ii]);
	  }

	  distMap.setDimensions( b_l, b_r, n_cells );
	  distMap.setConfined();
	  distMap.allocateCellArray();
	  distMap.clear();

	  catMap.setDimensions( b_l, b_r, n_cells );
	  catMap.setConfined();
	  catMap.allocateCellArray();
	  catMap.clear();

	  resMap.setDimensions( b_l, b_r, n_cells );
	  resMap.setConfined();
	  resMap.allocateCellArray();
	  resMap.clear();

	  max_cells[0] = IO.readUInt32();

	  //printf("@set_field.cc, maxIndx: %i: \n", max_cells );

	  for (int ii = 0; ii < max_cells[0] ; ++ii ) {

		indx[0]  = IO.readUInt32(); //index
		//printf("@set_field::readMap: nn = %i", indx);

		resMap[ indx[0] ] =  IO.readReal32() ;
	    //printf("@set_field::readMap: cellVal: %5.5f: \n", distMap[ indx ]);
	    }

	  max_cells[1] = IO.readUInt32();

	  for (int ii = 0; ii < max_cells[1] ; ++ii ) {

		indx[1]  = IO.readUInt32(); //index
		//printf("@set_field::readMap: nn = %i", indx);

		catMap[ indx[1] ] =  IO.readReal32() ;
		  //printf("@set_field::readMap: cellVal: %5.5f: \n", distMap[ indx ]);
	    }

	  max_cells[2] = IO.readUInt32();

	  for (int ii = 0; ii < max_cells[2] ; ++ii ) {

		indx[2]  = IO.readUInt32(); //index
		//printf("@set_field::readMap: nn = %i", indx);

		distMap[ indx[2] ] =  IO.readReal32() ;
		  //printf("@set_field::readMap: cellVal: %5.5f: \n", distMap[ indx ]);
	    }

		//printf("MAXval distMap= %5.5f , catMap= %5.5f, resMap= %5.5f\n", distMap.getMaxval(),catMap.getMaxval(),resMap.getMaxval());

	}




